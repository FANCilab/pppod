"""Export Suite2p calcium outputs to ALF/ONE multiphoton datasets."""

from __future__ import annotations

import csv
import re
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Iterable

import numpy as np


REQUIRED_PLANE_FILES = {
    "F.npy",
    "Fneu.npy",
    "spks.npy",
    "stat.npy",
    "ops.npy",
}


@dataclass(frozen=True)
class SessionInfo:
    subject: str
    date: str
    session: str
    local_raw_path: Path
    timing_path: Path
    time_zero: float
    frame_start: int
    raw_di3_count: int
    export_frame_count: int

    @property
    def frame_stop(self) -> int:
        return self.frame_start + self.export_frame_count

    @property
    def frame_slice(self) -> slice:
        return slice(self.frame_start, self.frame_stop)


@dataclass(frozen=True)
class RawSessionInfo:
    subject: str
    date: str
    session: str
    local_raw_path: Path
    timing_path: Path
    time_zero: float
    di3_count: int
    complete_cycle_count: int


@dataclass(frozen=True)
class PlaneInfo:
    path: Path
    plane_index: int
    fov_index: int

    @property
    def fov_name(self) -> str:
        return f"FOV_{self.fov_index:02d}"


def single_session_calcium_export(s2prunpath, dataoutroot, *, stim_type, fix_interleave=True):
    """Export one Suite2p run folder to ALF/ONE calcium imaging datasets.

    Parameters
    ----------
    s2prunpath : str or pathlib.Path
        Suite2p run folder containing run-level ``db.npy`` and ``plane*``
        folders. Network paths embedded in ``db.npy`` are parsed only as
        strings and remapped to the local ``2P`` data mirror.
    dataoutroot : str or pathlib.Path
        Root output folder. Files are written below
        ``subject/date/session/alf/FOV_NN``.
    stim_type : {"sparse_noise", "gratings"}
        Stimulus type used to select the target raw session when Suite2p
        concatenated multiple sessions. Only sparse-noise export is currently
        implemented.
    fix_interleave : bool, default=True
        If true, treat Suite2p plane folders as mROI FOVs and correct each
        FOV's timestamps from DI3/None ScanImage cycle intervals.

    Returns
    -------
    list[pathlib.Path]
        Paths to all written ALF datasets.
    """

    import sparse  # Required by ONE for .sparse_npz mask datasets.

    s2p_run_path = Path(s2prunpath)
    data_out_root = Path(dataoutroot)
    stim_type = _normalize_stim_type(stim_type)

    run_db = _load_dict_npy(s2p_run_path / "db.npy", "run-level Suite2p db.npy")
    raw_sessions = _raw_sessions_from_db(run_db, s2p_run_path)
    session = _select_session(raw_sessions, stim_type, fix_interleave=fix_interleave)
    planes = _find_planes(s2p_run_path)

    if not planes:
        raise FileNotFoundError(f"No Suite2p plane folders found in {s2p_run_path}")

    _validate_fov_frame_counts_match_sessions(planes, raw_sessions)
    iscell_by_plane = _load_iscell_by_plane(s2p_run_path, planes)

    if fix_interleave:
        frame_times_by_fov, timing_bad_frames_by_fov = _corrected_fov_times(
            session.timing_path,
            planes,
            n_frames=session.export_frame_count,
            time_zero=session.time_zero,
        )
    else:
        frame_times = _load_session_di3_times(
            session.timing_path,
            session.export_frame_count,
            time_zero=session.time_zero,
        )
        timing_bad = np.zeros(session.export_frame_count, dtype=bool)
        frame_times_by_fov = {plane.fov_name: frame_times for plane in planes}
        timing_bad_frames_by_fov = {plane.fov_name: timing_bad for plane in planes}

    written_files: list[Path] = []
    for plane in planes:
        written_files.extend(
            _export_plane(
                plane=plane,
                iscell=iscell_by_plane[plane.path],
                session=session,
                data_out_root=data_out_root,
                sparse_module=sparse,
                frame_times=frame_times_by_fov[plane.fov_name],
                timing_bad_frames=timing_bad_frames_by_fov[plane.fov_name],
            )
        )

    return written_files


def _export_plane(
    plane: PlaneInfo,
    iscell: np.ndarray,
    session: SessionInfo,
    data_out_root: Path,
    sparse_module,
    *,
    frame_times: np.ndarray,
    timing_bad_frames: np.ndarray | None = None,
):
    plane_path = plane.path
    missing = sorted(name for name in REQUIRED_PLANE_FILES if not (plane_path / name).is_file())
    if missing:
        raise FileNotFoundError(
            f"Suite2p plane folder is missing required files: {plane_path}; "
            f"missing={missing}"
        )

    F = np.load(plane_path / "F.npy", allow_pickle=False)
    Fneu = np.load(plane_path / "Fneu.npy", allow_pickle=False)
    spks = np.load(plane_path / "spks.npy", allow_pickle=False)
    stat = np.load(plane_path / "stat.npy", allow_pickle=True)
    ops = _load_dict_npy(plane_path / "ops.npy", "Suite2p ops.npy")

    _validate_trace_shapes(plane_path, F, Fneu, spks, iscell, stat)

    full_n_frames = F.shape[1]
    F = F[:, session.frame_slice]
    Fneu = Fneu[:, session.frame_slice]
    spks = spks[:, session.frame_slice]

    n_rois, n_frames = F.shape

    frame_times = np.asarray(frame_times, dtype=np.float64).reshape(-1)
    if frame_times.size != n_frames:
        raise ValueError(
            f"{plane.fov_name} corrected timing length does not match sliced trace frames: "
            f"{frame_times.size} vs {n_frames}"
        )

    bad_frames = _bad_frames_from_ops(ops, full_n_frames)[session.frame_slice].astype(bool, copy=False)

    if timing_bad_frames is not None:
        timing_bad_frames = np.asarray(timing_bad_frames, dtype=bool).reshape(-1)
        if timing_bad_frames.size != n_frames:
            raise ValueError(
                f"{plane.fov_name} timing_bad_frames length does not match sliced trace frames: "
                f"{timing_bad_frames.size} vs {n_frames}"
            )
        bad_frames = bad_frames | timing_bad_frames

    cell_classifier = _cell_classifier_from_iscell(iscell)
    roi_types = iscell[:, 0].astype(np.int16, copy=False)
    stack_pos = _stack_pos_from_stat(stat, plane.fov_index)
    mean_image = _mean_image_from_ops(ops)

    sparse_masks = _sparse_masks_from_stat(stat, mean_image.shape[:2], sparse_module)

    alf_folder = (
        data_out_root
        / session.subject
        / session.date
        / session.session
        / "alf"
        / plane.fov_name
    )
    alf_folder.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    written.append(_save_npy(alf_folder / "mpci.times.npy", frame_times.astype(np.float64)))
    written.append(_save_npy(alf_folder / "mpci.badFrames.npy", bad_frames.astype(bool)))
    written.append(_save_npy(alf_folder / "mpci.ROIActivityF.npy", F.T.astype(np.float32, copy=False)))
    written.append(
        _save_npy(
            alf_folder / "mpci.ROINeuropilActivityF.npy",
            Fneu.T.astype(np.float32, copy=False),
        )
    )
    written.append(
        _save_npy(
            alf_folder / "mpci.ROIActivityDeconvolved.npy",
            spks.T.astype(np.float32, copy=False),
        )
    )
    written.append(_save_npy(alf_folder / "mpciROIs.cellClassifier.npy", cell_classifier))
    written.append(_save_npy(alf_folder / "mpciROIs.mpciROITypes.npy", roi_types))
    written.append(_save_tsv(alf_folder / "mpciROITypes.names.tsv", _roi_type_rows()))
    written.append(_save_npy(alf_folder / "mpciROIs.stackPos.npy", stack_pos.astype(np.int64)))
    written.append(_save_sparse_mask(alf_folder / "mpciROIs.masks.sparse_npz", sparse_masks, sparse_module))
    written.append(_save_npy(alf_folder / "mpciMeanImage.images.npy", mean_image))

    return written


def _load_dict_npy(path: Path, description: str) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {description}: {path}")

    value = np.load(path, allow_pickle=True)
    if isinstance(value, np.ndarray) and value.shape == ():
        value = value.item()

    if not isinstance(value, dict):
        raise ValueError(f"{description} did not contain a dictionary: {path}")

    return value


def _normalize_stim_type(stim_type: str) -> str:
    normalized = str(stim_type).strip().lower()
    if normalized == "sparse_noise":
        return normalized
    if normalized == "gratings":
        raise NotImplementedError(
            'Calcium export for stim_type="gratings" is not implemented yet.'
        )
    raise ValueError('stim_type must be "sparse_noise" or "gratings"')


def _raw_sessions_from_db(db: dict, s2p_run_path: Path) -> tuple[RawSessionInfo, ...]:
    data_paths = _data_paths_with_known_edge_cases(_as_list(db.get("data_path")), s2p_run_path)
    if not data_paths:
        raise ValueError('Run-level db.npy must contain non-empty "data_path"')

    run_date_folder = None
    try:
        parsed = [_parse_subject_date_session(str(path_value)) for path_value in data_paths]
    except ValueError:
        run_date_folder, run_date = _date_from_run_path(s2p_run_path)
        parsed = [
            _parse_subject_date_session(str(path_value), fallback_date=run_date)
            for path_value in data_paths
        ]
    subjects = {item[0] for item in parsed}
    dates = {item[1] for item in parsed}

    if len(subjects) != 1 or len(dates) != 1:
        raise ValueError(
            "Suite2p run must correspond to exactly one subject/date; "
            f"subjects={sorted(subjects)}, dates={sorted(dates)}"
        )

    local_2p_root = _infer_local_2p_root(s2p_run_path)
    raw_sessions: list[RawSessionInfo] = []
    for path_value, (subject, date, session) in zip(data_paths, parsed):
        local_raw_path = _local_raw_path(
            path_value,
            local_2p_root,
            fallback_date_folder=run_date_folder,
        )
        timing_path = _find_timing_csv(local_raw_path)
        raw_sessions.append(
            RawSessionInfo(
                subject=subject,
                date=date,
                session=session,
                local_raw_path=local_raw_path,
                timing_path=timing_path,
                time_zero=_session_time_zero(local_raw_path),
                di3_count=_count_di3_times(timing_path),
                complete_cycle_count=_count_complete_di3_none_cycles(timing_path),
            )
        )

    return tuple(raw_sessions)


def _data_paths_with_known_edge_cases(data_paths: list, s2p_run_path: Path) -> list:
    """Return Suite2p db data paths, correcting known mishandled runs.

    NM036/20260507 was processed as a fictional raw session "3", but the TIFFs
    came from raw sessions 1 and 2. Keep this hardcoded and exact so all normal
    exporter code paths keep using the db.npy metadata unchanged.

    NM037/20260609/Processed contains TIFFs from raw sessions 1, 2, and 3, but
    its Suite2p db.npy lists only sessions 1 and 2. Session 1 is the sparse-noise
    session; the exporter still selects only that session after validating the
    full concatenated frame count against all three raw sessions.
    """

    if _is_nm036_20260507_fictional_session3_run(data_paths, s2p_run_path):
        return [
            "Z:/Data/2P/NM036/20260507/1",
            "Z:/Data/2P/NM036/20260507/2",
        ]

    if _is_nm037_20260609_missing_session3_run(data_paths, s2p_run_path):
        return [
            "Z:/Data/2P/NM037/20260609/1",
            "Z:/Data/2P/NM037/20260609/2",
            "Z:/Data/2P/NM037/20260609/3",
        ]

    return data_paths


def _is_nm036_20260507_fictional_session3_run(data_paths: list, s2p_run_path: Path) -> bool:
    if [str(value).replace("\\", "/").rstrip("/") for value in data_paths] != ["Z:/Data/2P/NM036/3"]:
        return False

    parts = [part.lower() for part in _split_path_parts(str(s2p_run_path))]
    required_tail = ["2p", "nm036", "20260507", "processed", "1"]
    return len(parts) >= len(required_tail) and parts[-len(required_tail):] == required_tail


def _is_nm037_20260609_missing_session3_run(data_paths: list, s2p_run_path: Path) -> bool:
    normalized_data_paths = [str(value).replace("\\", "/").rstrip("/") for value in data_paths]
    if normalized_data_paths != [
        "Z:/Data/2P/NM037/20260609/1",
        "Z:/Data/2P/NM037/20260609/2",
    ]:
        return False

    parts = [part.lower() for part in _split_path_parts(str(s2p_run_path))]
    required_tail = ["2p", "nm037", "20260609", "processed"]
    return len(parts) >= len(required_tail) and parts[-len(required_tail):] == required_tail


def _select_session(
    raw_sessions: tuple[RawSessionInfo, ...],
    stim_type: str,
    *,
    fix_interleave: bool,
) -> SessionInfo:
    if stim_type != "sparse_noise":
        raise NotImplementedError(f'Calcium export for stim_type="{stim_type}" is not implemented.')

    matches = [
        raw_session
        for raw_session in raw_sessions
        if any(path.is_file() for path in _raw_file_candidates(raw_session.local_raw_path, "SparseNoise_Log.bin"))
    ]
    if not matches:
        raise FileNotFoundError(
            "Could not identify sparse-noise raw session: none of the local raw "
            "sessions contains SparseNoise_Log.bin"
        )
    if len(matches) > 1:
        raise ValueError(
            "Could not identify sparse-noise raw session unambiguously: multiple "
            f"sessions contain SparseNoise_Log.bin: {[item.session for item in matches]}"
        )

    target = matches[0]
    frame_start = 0
    for raw_session in raw_sessions:
        if raw_session is target:
            break
        frame_start += raw_session.di3_count

    # Never silently drop imaging frames.  For mROI timing, a missing DI3/None
    # pair is detected later and raises instead of shortening the export.
    export_frame_count = target.di3_count
    if export_frame_count <= 0:
        raise ValueError(
            f"Target session {target.local_raw_path} has no exportable imaging frames"
        )

    return SessionInfo(
        subject=target.subject,
        date=target.date,
        session=target.session,
        local_raw_path=target.local_raw_path,
        timing_path=target.timing_path,
        time_zero=target.time_zero,
        frame_start=frame_start,
        raw_di3_count=target.di3_count,
        export_frame_count=export_frame_count,
    )


def _parse_subject_date_session(raw_path: str, fallback_date: str | None = None) -> tuple[str, str, str]:
    parts = _split_path_parts(raw_path)
    two_p_indices = [i for i, part in enumerate(parts) if part.lower() == "2p"]
    if not two_p_indices:
        raise ValueError(f'Could not find "2P" folder in db data_path: {raw_path}')

    i_two_p = two_p_indices[-1]
    if i_two_p + 3 >= len(parts):
        if fallback_date is None or i_two_p + 2 >= len(parts):
            raise ValueError(
                "db data_path must contain .../2P/<subject>/<date>/<session> "
                "or, with a Suite2p run-path date fallback, .../2P/<subject>/<session>: "
                f"{raw_path}"
            )

        subject = parts[i_two_p + 1]
        date = fallback_date
        session = parts[i_two_p + 2]
    else:
        subject = parts[i_two_p + 1]
        date = _format_date(parts[i_two_p + 2])
        session = parts[i_two_p + 3]


    if not subject or not session:
        raise ValueError(f"Empty subject or session parsed from db data_path: {raw_path}")
    if not re.fullmatch(r"\d+", session):
        raise ValueError(f"Session folder must be numeric for ONE output: {session}")

    return subject, date, session


def _format_date(token: str) -> str:
    match = re.fullmatch(r"((?:19|20)\d{2})-?(\d{2})-?(\d{2})", token)
    if not match:
        raise ValueError(f"Date folder must be YYYYMMDD or YYYY-MM-DD, got: {token}")
    return f"{match.group(1)}-{match.group(2)}-{match.group(3)}"


def _date_from_run_path(s2p_run_path: Path) -> tuple[str, str]:
    matches: dict[str, str] = {}
    for part in _split_path_parts(str(s2p_run_path)):
        try:
            formatted = _format_date(part)
        except ValueError:
            continue
        matches.setdefault(formatted, part)

    if len(matches) != 1:
        raise ValueError(
            "Suite2p run path must contain exactly one YYYYMMDD or YYYY-MM-DD "
            f"date folder for db data_path fallback; found {len(matches)}: {sorted(matches)}"
        )
    formatted, folder = next(iter(matches.items()))
    return folder, formatted


def _split_path_parts(path_string: str) -> list[str]:
    normalized = path_string.strip().replace("\\", "/")
    return [part for part in normalized.split("/") if part and part not in {".", ".."}]


def _infer_local_2p_root(s2p_run_path: Path) -> Path:
    resolved_parts = s2p_run_path.resolve().parts
    for idx in range(len(resolved_parts) - 1, -1, -1):
        if resolved_parts[idx].lower() == "2p":
            return Path(*resolved_parts[: idx + 1])
    raise ValueError(
        f'Could not infer local 2P root from Suite2p run path "{s2p_run_path}". '
        'The run path must live under a folder named "2P" so remote db paths can '
        "be remapped without hardcoded machine-specific roots."
    )


def _local_raw_path(raw_data_path, local_2p_root: Path, fallback_date_folder: str | None = None) -> Path:
    parts = _split_path_parts(str(raw_data_path))
    two_p_indices = [i for i, part in enumerate(parts) if part.lower() == "2p"]
    if not two_p_indices:
        raise ValueError(f'Could not remap raw path because it has no "2P" component: {raw_data_path}')

    suffix = parts[two_p_indices[-1] + 1 :]
    if len(suffix) == 2 and fallback_date_folder is not None:
        suffix = [suffix[0], fallback_date_folder, suffix[1]]
    return local_2p_root.joinpath(*suffix)


def _find_timing_csv(local_raw_path: Path) -> Path:
    candidates = [
        local_raw_path / "timestamped_digital.csv",
        local_raw_path.parent / "timestamped_digital.csv",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate

    raise FileNotFoundError(
        "Could not find local timestamped_digital.csv for raw path. "
        f"Tried: {[str(path) for path in candidates]}"
    )


SESSION_TIME_ZERO_FILES = (
    "photodiode.csv",
    "timestamped_digital.csv",
    "encoder_log.csv",
    "times_grid.csv",
)


def _session_time_zero(raw_path: Path) -> float:
    """Return the shared session time zero for a raw session folder.

    If any raw files use an absolute acquisition clock, use the earliest absolute
    first sample and ignore already-relative files.  If all files are relative,
    use the earliest first sample.
    """

    first_values: list[float] = []
    for filename in SESSION_TIME_ZERO_FILES:
        for candidate in _raw_file_candidates(raw_path, filename):
            if not candidate.is_file():
                continue
            first = _first_numeric_csv_value(candidate)
            if first is not None:
                first_values.append(float(first))

    if not first_values:
        raise FileNotFoundError(
            f"Could not determine session time zero for raw path {raw_path}; "
            f"tried {list(SESSION_TIME_ZERO_FILES)} in the raw folder and parent folder"
        )

    absolute_values = [value for value in first_values if abs(value) > 1e5]
    if absolute_values:
        return float(min(absolute_values))
    return float(min(first_values))


def _raw_file_candidates(raw_path: Path, filename: str) -> tuple[Path, Path]:
    return (raw_path / filename, raw_path.parent / filename)


def _first_numeric_csv_value(path: Path) -> float | None:
    with open(path, newline="") as file_obj:
        reader = csv.reader(file_obj)
        for row in reader:
            if not row:
                continue
            try:
                value = float(row[0])
            except ValueError:
                continue
            if np.isfinite(value):
                return value
    return None


def _to_session_relative_times(
    times: np.ndarray,
    time_zero: float,
    *,
    description: str,
) -> np.ndarray:
    """Convert absolute timestamps to session-relative time without dropping rows.

    If the input is already relative, it is returned unchanged.  The decision is
    based on which representation is closer to a session-start timebase.
    """

    values = np.asarray(times, dtype=np.float64).reshape(-1)
    if values.size == 0:
        raise ValueError(f"{description} must be non-empty")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{description} contains non-finite values")
    if not np.isfinite(time_zero):
        raise ValueError(f"Session time zero for {description} is non-finite: {time_zero}")

    shifted = values - float(time_zero)
    raw_scale = float(np.nanmedian(np.abs(values)))
    shifted_scale = float(np.nanmedian(np.abs(shifted)))

    # Absolute acquisition clocks are usually much farther from zero than the
    # corresponding session-relative values.  Relative inputs remain unchanged.
    if shifted_scale < raw_scale:
        values = shifted

    if values.size > 1 and np.any(np.diff(values) <= 0):
        raise ValueError(f"{description} must be strictly increasing")
    return values.astype(np.float64, copy=False)


def _count_di3_times(timing_path: Path) -> int:
    return int(_load_all_di3_times(timing_path).size)


def _count_complete_di3_none_cycles(timing_path: Path) -> int:
    di3_times, _ = _load_di3_none_pairs(timing_path)
    return int(di3_times.size)


def _validate_fov_frame_counts_match_sessions(
    planes: list[PlaneInfo],
    raw_sessions: tuple[RawSessionInfo, ...],
) -> None:
    expected_frames = sum(raw_session.di3_count for raw_session in raw_sessions)
    for plane in planes:
        F = np.load(plane.path / "F.npy", allow_pickle=False)
        if F.ndim != 2:
            raise ValueError(f"F.npy must be 2D [nROIs, nFrames]: {plane.path / 'F.npy'}")
        if F.shape[1] != expected_frames:
            counts = [(item.session, item.di3_count) for item in raw_sessions]
            raise ValueError(
                f"{plane.path / 'F.npy'} has {F.shape[1]} frames but raw sessions contain "
                f"{expected_frames} DI3 onsets in db['data_path'] order: {counts}"
            )


def _load_session_di3_times(
    timing_path: Path,
    n_frames: int,
    *,
    time_zero: float,
) -> np.ndarray:
    times = _load_all_di3_times(timing_path)
    if times.size != n_frames:
        raise ValueError(
            f"DI3 frame-time count in {timing_path} does not match Suite2p frames; "
            f"found {times.size}, expected {n_frames}. Refusing to trim data."
        )
    return _to_session_relative_times(
        times.astype(np.float64, copy=False),
        time_zero,
        description=f"DI3 frame times from {timing_path}",
    )

def _corrected_fov_times(
    timing_path: Path,
    planes: list[PlaneInfo],
    *,
    n_frames: int,
    time_zero: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Return corrected mROI frame times and timing-derived bad-frame masks.

    Normal case:
      n DI3 onsets == n None edges == n Suite2p frames.

    Accepted terminal edge case:
      n DI3 onsets == n Suite2p frames
      n complete DI3/None cycles == n Suite2p frames - 1

    In the terminal edge case, the final DI3 frame is preserved. The missing
    final cycle duration is estimated from recent complete cycles, and the
    final frame is marked bad for every FOV. This preserves all Suite2p frames
    without silently trimming data.
    """

    all_di3_times = _load_all_di3_times(timing_path)
    paired_di3_times, none_times = _load_di3_none_pairs(timing_path)

    if all_di3_times.size == 0:
        raise ValueError(f"No DI3 imaging onsets found in {timing_path}")

    if all_di3_times.size != n_frames:
        raise ValueError(
            f"DI3 onset count in {timing_path} does not match Suite2p frames; "
            f"found {all_di3_times.size}, expected {n_frames}. Refusing to trim data."
        )

    timing_bad = np.zeros(n_frames, dtype=bool)

    if paired_di3_times.size == n_frames:
        di3_times = paired_di3_times

    elif paired_di3_times.size == n_frames - 1:
        if paired_di3_times.size == 0:
            raise ValueError(
                f"{timing_path} has one DI3 onset but no complete DI3/None cycle, "
                "so the missing cycle duration cannot be estimated."
            )

        if not np.allclose(
            paired_di3_times,
            all_di3_times[: paired_di3_times.size],
            rtol=0.0,
            atol=1e-9,
        ):
            raise ValueError(
                f"{timing_path} has one missing DI3/None cycle, but it is not a clean "
                "terminal missing-None edge case. Refusing to guess timing."
            )

        di3_times = all_di3_times
        timing_bad[-1] = True

    else:
        raise ValueError(
            f"Complete DI3/None cycle count in {timing_path} does not match Suite2p frames; "
            f"found {paired_di3_times.size}, expected {n_frames}. "
            "Only the terminal case with exactly one missing final None edge is allowed."
        )

    cycle_dt_complete = none_times - paired_di3_times
    if cycle_dt_complete.size == 0:
        raise ValueError(f"No complete DI3/None cycles available in {timing_path}")

    if np.any(cycle_dt_complete <= 0):
        bad = np.flatnonzero(cycle_dt_complete <= 0)[:5]
        raise ValueError(
            f"DI3/None cycles must have positive duration in {timing_path}; "
            f"bad cycle indices include {bad.tolist()}"
        )

    if paired_di3_times.size == n_frames:
        cycle_dt = cycle_dt_complete
    else:
        recent = cycle_dt_complete[-min(100, cycle_dt_complete.size):]
        estimated_final_dt = float(np.median(recent))

        if not np.isfinite(estimated_final_dt) or estimated_final_dt <= 0:
            raise ValueError(
                f"Could not estimate missing final DI3/None duration in {timing_path}: "
                f"estimated value was {estimated_final_dt}"
            )

        cycle_dt = np.concatenate(
            [
                cycle_dt_complete.astype(np.float64, copy=False),
                np.asarray([estimated_final_dt], dtype=np.float64),
            ]
        )

        warnings.warn(
            f"{timing_path} has {n_frames} DI3 onsets but only {paired_di3_times.size} "
            f"complete DI3/None cycles. Preserving the final Suite2p frame, estimating "
            f"its missing cycle duration as {estimated_final_dt:.9g} s from recent cycles, "
            "and marking the final frame bad in mpci.badFrames.",
            RuntimeWarning,
            stacklevel=2,
        )

    if di3_times.size != n_frames or cycle_dt.size != n_frames:
        raise ValueError(
            f"Internal mROI timing length mismatch for {timing_path}: "
            f"di3_times={di3_times.size}, cycle_dt={cycle_dt.size}, n_frames={n_frames}"
        )

    fov_heights = _fov_heights_from_ops(planes)
    height_values = np.asarray([fov_heights[plane.fov_name] for plane in planes], dtype=np.float64)
    total_height = float(np.sum(height_values))
    if not np.isfinite(total_height) or total_height <= 0:
        raise ValueError(f"Invalid total FOV Ly from Suite2p ops.npy files: {total_height}")

    corrected_relative: dict[str, np.ndarray] = {}
    timing_bad_by_fov: dict[str, np.ndarray] = {}

    cumulative_height = 0.0
    for plane, height in zip(planes, height_values):
        offset = cycle_dt * (cumulative_height / total_height)
        corrected_absolute = di3_times + offset

        corrected_relative[plane.fov_name] = _to_session_relative_times(
            corrected_absolute.astype(np.float64, copy=False),
            time_zero,
            description=f"{plane.fov_name} corrected frame times from {timing_path}",
        )

        timing_bad_by_fov[plane.fov_name] = timing_bad.copy()
        cumulative_height += float(height)

    return corrected_relative, timing_bad_by_fov

def _fov_heights_from_ops(planes: list[PlaneInfo]) -> dict[str, int]:
    heights: dict[str, int] = {}
    for plane in planes:
        ops = _load_dict_npy(plane.path / "ops.npy", "Suite2p ops.npy")
        if "Ly" not in ops:
            raise ValueError(f'{plane.path / "ops.npy"} must contain "Ly" for mROI timing')
        ly = int(np.asarray(ops["Ly"]).reshape(-1)[0])
        if ly <= 0:
            raise ValueError(f'Invalid "Ly" in {plane.path / "ops.npy"}: {ly}')
        heights[plane.fov_name] = ly
    return heights


def _load_di3_none_pairs(timing_path: Path) -> tuple[np.ndarray, np.ndarray]:
    di3_times: list[float] = []
    none_times: list[float] = []
    pending_di3: float | None = None
    di3_label_mode = _infer_di3_label_mode(timing_path)
    with open(timing_path, newline="") as file_obj:
        reader = csv.reader(file_obj)
        for row in reader:
            label = _timestamped_digital_label(row, di3_label_mode)
            if label not in {"DI3", "None"}:
                continue
            try:
                timestamp = float(row[0])
            except ValueError as exc:
                raise ValueError(f"Invalid {label} timestamp in {timing_path}: {row}") from exc
            if not np.isfinite(timestamp):
                raise ValueError(f"Non-finite {label} timestamp in {timing_path}: {row}")
            if label == "DI3":
                pending_di3 = timestamp
            elif label == "None" and pending_di3 is not None:
                di3_times.append(pending_di3)
                none_times.append(timestamp)
                pending_di3 = None

    di3_values = np.asarray(di3_times, dtype=np.float64)
    none_values = np.asarray(none_times, dtype=np.float64)
    if di3_values.size and np.any(np.diff(di3_values) <= 0):
        raise ValueError(f"DI3/None pair DI3 times must be strictly increasing in {timing_path}")
    return di3_values, none_values



def _find_planes(s2p_run_path: Path) -> list[PlaneInfo]:
    plane_pattern = re.compile(r"^plane(\d+)$", re.IGNORECASE)
    candidates: list[tuple[int, Path]] = []
    for child in s2p_run_path.iterdir():
        if not child.is_dir():
            continue
        match = plane_pattern.match(child.name)
        if match:
            candidates.append((int(match.group(1)), child))

    candidates.sort(key=lambda item: item[0])
    planes: list[PlaneInfo] = []
    fallback_fov_index = 0
    for plane_index, plane_path in candidates:
        plane_db_path = plane_path / "db.npy"
        if plane_db_path.is_file():
            plane_db = _load_dict_npy(plane_db_path, "plane-level Suite2p db.npy")
            fov_index = int(plane_db.get("iroi", fallback_fov_index))
        else:
            fov_index = fallback_fov_index
        planes.append(PlaneInfo(path=plane_path, plane_index=plane_index, fov_index=fov_index))
        fallback_fov_index += 1

    _assert_unique_fov_names(planes)
    return planes


def _assert_unique_fov_names(planes: Iterable[PlaneInfo]) -> None:
    by_name: dict[str, list[Path]] = defaultdict(list)
    for plane in planes:
        by_name[plane.fov_name].append(plane.path)

    duplicates = {name: paths for name, paths in by_name.items() if len(paths) > 1}
    if duplicates:
        raise ValueError(
            "Multiple plane folders map to the same fov_name; this exporter expects one "
            f"Suite2p plane folder per FOV. duplicates={duplicates}"
        )


def _infer_timing_stream_count(run_db: dict, planes: list[PlaneInfo], local_raw_path: Path) -> int:
    candidates: list[tuple[str, int]] = []
    candidates.extend(_nplanes_candidates_from_mapping(run_db, "run db.npy"))
    candidates.extend(_nplanes_candidates_from_scanimage_tiff(local_raw_path))

    for plane in planes:
        plane_db_path = plane.path / "db.npy"
        if plane_db_path.is_file():
            plane_db = _load_dict_npy(plane_db_path, "plane-level Suite2p db.npy")
            candidates.extend(
                _nplanes_candidates_from_mapping(plane_db, f"{plane.path.name}/db.npy")
            )

        ops_path = plane.path / "ops.npy"
        if ops_path.is_file():
            ops = _load_dict_npy(ops_path, "Suite2p ops.npy")
            candidates.extend(
                _nplanes_candidates_from_mapping(ops, f"{plane.path.name}/ops.npy")
            )

    unique_by_source = {}
    for source, value in candidates:
        unique_by_source.setdefault((source, value), None)
    candidates = list(unique_by_source)
    values = {value for _, value in candidates}

    if len(values) > 1:
        details = ", ".join(f"{source}={value}" for source, value in candidates)
        raise ValueError(
            "Conflicting Suite2p/ScanImage nplanes metadata for timing stream count: "
            f"{details}"
        )

    if values:
        n_streams = values.pop()
        source_names = sorted({source for source, _ in candidates})
        source_families = {_nplanes_source_family(source) for source in source_names}
        if len(source_families) == 1:
            warnings.warn(
                f"Using timing stream count nplanes={n_streams} from "
                f"{', '.join(source_names)}; "
                "no independent Suite2p/ScanImage nplanes source was available to verify it.",
                RuntimeWarning,
                stacklevel=2,
            )
        return n_streams

    warnings.warn(
        "Could not find Suite2p/ScanImage nplanes metadata "
        "(Suite2p nplanes or SI.hStackManager.numSlices/zs); falling back to one "
        "timing stream.",
        RuntimeWarning,
        stacklevel=2,
    )
    return 1


def _nplanes_source_family(source: str) -> str:
    lower = source.lower()
    if lower.endswith("si.hstackmanager.numslices"):
        return "scanimage_num_slices"
    if lower.endswith("si.hstackmanager.zs"):
        return "scanimage_zs"
    if lower.endswith("nplanes"):
        return "suite2p_nplanes"
    return source


def _nplanes_candidates_from_mapping(mapping: dict, source: str) -> list[tuple[str, int]]:
    candidates: list[tuple[str, int]] = []
    for key_path, value in _iter_metadata_values(mapping):
        key_path_lower = key_path.lower()
        if key_path_lower.endswith("nplanes"):
            candidates.extend(_candidate_ints(value, f"{source}:{key_path}"))
        elif key_path_lower.endswith("si.hstackmanager.numslices"):
            candidates.extend(_candidate_ints(value, f"{source}:{key_path}"))
        elif key_path_lower.endswith("si.hstackmanager.zs"):
            candidates.extend(_candidate_count_from_zs(value, f"{source}:{key_path}"))
        elif isinstance(value, str):
            candidates.extend(_nplanes_candidates_from_header_text(value, source))
    return candidates


def _iter_metadata_values(value, prefix: str = ""):
    if isinstance(value, dict):
        for key, child in value.items():
            key_text = str(key)
            child_prefix = f"{prefix}.{key_text}" if prefix else key_text
            yield from _iter_metadata_values(child, child_prefix)
        return

    if isinstance(value, (list, tuple)):
        for idx, child in enumerate(value):
            child_prefix = f"{prefix}[{idx}]" if prefix else f"[{idx}]"
            yield from _iter_metadata_values(child, child_prefix)
        return

    if isinstance(value, np.ndarray) and value.dtype == object and value.size <= 100:
        for idx, child in enumerate(value.reshape(-1)):
            child_prefix = f"{prefix}[{idx}]" if prefix else f"[{idx}]"
            yield from _iter_metadata_values(child, child_prefix)
        return

    yield prefix, value


def _nplanes_candidates_from_scanimage_tiff(local_raw_path: Path) -> list[tuple[str, int]]:
    tiff_path = _find_first_tiff(local_raw_path)
    if tiff_path is None:
        return []

    try:
        import tifffile
    except ImportError:
        warnings.warn(
            "Cannot verify ScanImage SI.hStackManager metadata from TIFF header because "
            "the optional tifffile package is unavailable.",
            RuntimeWarning,
            stacklevel=2,
        )
        return []

    try:
        with tifffile.TiffFile(tiff_path) as tif:
            software_tag = tif.pages[0].tags.get("Software")
            software = software_tag.value if software_tag is not None else ""
    except Exception as exc:
        warnings.warn(
            f"Could not read ScanImage metadata from TIFF header {tiff_path}: {exc}",
            RuntimeWarning,
            stacklevel=2,
        )
        return []

    if not isinstance(software, str):
        software = str(software)

    return _nplanes_candidates_from_header_text(software, f"{tiff_path.name}:Software")


def _find_first_tiff(local_raw_path: Path) -> Path | None:
    for pattern in ("*.tif", "*.tiff", "*.TIF", "*.TIFF"):
        matches = sorted(local_raw_path.glob(pattern))
        if matches:
            return matches[0]
    return None


def _nplanes_candidates_from_header_text(text: str, source: str) -> list[tuple[str, int]]:
    candidates: list[tuple[str, int]] = []
    fields = {
        "SI.hStackManager.numSlices": _candidate_ints,
        "SI.hStackManager.zs": _candidate_count_from_zs,
    }
    for field, parser in fields.items():
        match = re.search(rf"(?:^|\n){re.escape(field)}\s*=\s*([^\r\n]+)", text)
        if match:
            candidates.extend(parser(match.group(1).strip(), f"{source}:{field}"))
    return candidates


def _candidate_ints(value, source: str) -> list[tuple[str, int]]:
    array = np.asarray(value)
    if array.size != 1:
        return []
    try:
        number = float(array.reshape(-1)[0])
    except (TypeError, ValueError):
        try:
            number = float(str(value).strip().strip("'\""))
        except ValueError:
            return []

    if not np.isfinite(number) or number < 1 or not number.is_integer():
        return []
    return [(source, int(number))]


def _candidate_count_from_zs(value, source: str) -> list[tuple[str, int]]:
    if isinstance(value, str):
        zs_text = value.strip().strip("'\"")
        numbers = re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", zs_text)
        if numbers:
            return [(source, len(numbers))]

    array = np.asarray(value)
    if array.size:
        return [(source, int(array.size))]

    return []

def _load_iscell_by_plane(s2p_run_path: Path, planes: list[PlaneInfo]) -> dict[Path, np.ndarray]:
    """Load Suite2p iscell labels for each exported FOV.

    The Suite2p combined GUI is the authoritative curation source when it
    exists. Single-FOV Suite2p runs may not have a combined folder; in that
    specific case, use the curated plane0/iscell.npy labels directly.

    Mapping rule:
      combined/iscell.npy row i belongs with combined/stat.npy row i.
      combined/stat.npy row i has an "iplane" field indicating the source plane.
      Within each iplane group, rows are kept in combined/stat row order and
      assigned to that plane's ROI order.

    For multi-plane runs, planeX/iscell.npy may be used only to verify row
    order via the unchanged classifier score column, not as a label source.
    """

    combined_dir = s2p_run_path / "combined"
    combined_iscell_path = combined_dir / "iscell.npy"
    combined_stat_path = combined_dir / "stat.npy"

    if not combined_dir.exists() and len(planes) == 1:
        plane = planes[0]
        return {plane.path: _load_single_plane_iscell(plane)}

    if not combined_iscell_path.is_file():
        raise FileNotFoundError(
            "Missing combined/iscell.npy. Cell labels must come from the visually "
            "curated Suite2p combined view for multi-plane runs. Single-plane runs "
            f"without a combined folder use plane-level labels. Missing: {combined_iscell_path}"
        )
    if not combined_stat_path.is_file():
        raise FileNotFoundError(
            "Missing combined/stat.npy. It is required to map combined/iscell.npy "
            f"back to FOVs: {combined_stat_path}"
        )

    combined_iscell = np.load(combined_iscell_path, allow_pickle=False)
    combined_stat = np.load(combined_stat_path, allow_pickle=True)

    if combined_iscell.ndim != 2:
        raise ValueError(f"combined/iscell.npy must be 2D, got {combined_iscell.shape}")
    if combined_stat.ndim != 1:
        raise ValueError(f"combined/stat.npy must be a 1D object array, got {combined_stat.shape}")
    if combined_iscell.shape[0] != combined_stat.shape[0]:
        raise ValueError(
            "combined/iscell.npy and combined/stat.npy row counts differ: "
            f"{combined_iscell.shape[0]} vs {combined_stat.shape[0]}"
        )

    plane_by_index = {plane.plane_index: plane for plane in planes}
    combined_plane_indices = np.asarray(
        [
            _combined_stat_plane_index(roi_stat, idx, combined_stat_path)
            for idx, roi_stat in enumerate(combined_stat)
        ],
        dtype=np.int64,
    )

    unknown_planes = sorted(set(combined_plane_indices.tolist()) - set(plane_by_index))
    if unknown_planes:
        raise ValueError(
            "combined/stat.npy references plane indices not exported as FOVs: "
            f"{unknown_planes}. Exported plane indices: {sorted(plane_by_index)}"
        )

    iscell_by_plane: dict[Path, np.ndarray] = {}

    for plane in planes:
        F = np.load(plane.path / "F.npy", allow_pickle=False)
        if F.ndim != 2:
            raise ValueError(f"F.npy must be 2D [nROIs, nFrames]: {plane.path / 'F.npy'}")

        n_plane_rois = F.shape[0]
        combined_indices = np.flatnonzero(combined_plane_indices == plane.plane_index)

        if combined_indices.size != n_plane_rois:
            raise ValueError(
                f"combined/stat.npy has {combined_indices.size} ROI(s) for plane{plane.plane_index}, "
                f"but {plane.path / 'F.npy'} has {n_plane_rois} ROI(s). Refusing to guess mapping."
            )

        labels = combined_iscell[combined_indices, :].copy()

        plane_iscell_path = plane.path / "iscell.npy"
        if plane_iscell_path.is_file():
            plane_iscell = np.load(plane_iscell_path, allow_pickle=False)

            if plane_iscell.ndim != 2:
                raise ValueError(f"{plane_iscell_path} must be 2D, got {plane_iscell.shape}")
            if plane_iscell.shape[0] != n_plane_rois:
                raise ValueError(
                    f"{plane_iscell_path} has {plane_iscell.shape[0]} rows but "
                    f"{plane.path / 'F.npy'} has {n_plane_rois} ROIs"
                )

            if combined_iscell.shape[1] >= 2 and plane_iscell.shape[1] >= 2:
                score_diff = np.abs(
                    labels[:, 1].astype(np.float64, copy=False)
                    - plane_iscell[:, 1].astype(np.float64, copy=False)
                )
                if not np.all(np.isfinite(score_diff)):
                    raise ValueError(
                        f"Non-finite iscell classifier score differences while verifying "
                        f"combined-to-plane mapping for plane{plane.plane_index}"
                    )
                if np.max(score_diff) > 1e-12:
                    raise ValueError(
                        f"combined/iscell.npy rows for plane{plane.plane_index} do not appear "
                        f"to be in the same ROI order as {plane_iscell_path}: "
                        f"max classifier-score difference is {np.max(score_diff):.6g}. "
                        "Refusing to guess mapping."
                    )
        else:
            if combined_indices.size > 1 and not np.all(np.diff(combined_indices) == 1):
                raise ValueError(
                    f"{plane_iscell_path} is missing and combined/stat.npy rows for "
                    f"plane{plane.plane_index} are not contiguous. Cannot verify row order."
                )
            warnings.warn(
                f"{plane_iscell_path} is missing. Using combined/iscell.npy labels grouped by "
                f"combined/stat.npy iplane only; no plane-level score diagnostic was possible.",
                RuntimeWarning,
                stacklevel=2,
            )

        iscell_by_plane[plane.path] = labels

    return iscell_by_plane


def _load_single_plane_iscell(plane: PlaneInfo) -> np.ndarray:
    plane_iscell_path = plane.path / "iscell.npy"
    if not plane_iscell_path.is_file():
        raise FileNotFoundError(
            "Missing plane-level iscell.npy for single-plane Suite2p run without "
            f"a combined folder: {plane_iscell_path}"
        )

    F = np.load(plane.path / "F.npy", allow_pickle=False)
    if F.ndim != 2:
        raise ValueError(f"F.npy must be 2D [nROIs, nFrames]: {plane.path / 'F.npy'}")

    iscell = np.load(plane_iscell_path, allow_pickle=False)
    if iscell.ndim != 2:
        raise ValueError(f"{plane_iscell_path} must be 2D, got {iscell.shape}")
    if iscell.shape[0] != F.shape[0]:
        raise ValueError(
            f"{plane_iscell_path} has {iscell.shape[0]} rows but "
            f"{plane.path / 'F.npy'} has {F.shape[0]} ROIs"
        )

    warnings.warn(
        f"{plane.path.parent / 'combined'} is missing for a single-plane Suite2p run; "
        f"using curated labels from {plane_iscell_path}.",
        RuntimeWarning,
        stacklevel=2,
    )
    return iscell.copy()


def _combined_stat_plane_index(roi_stat, roi_index: int, stat_path: Path) -> int:
    try:
        return int(np.asarray(_stat_value(roi_stat, "iplane")).reshape(-1)[0])
    except Exception as exc:
        raise KeyError(
            f'{stat_path} entry {roi_index} is missing required "iplane" field'
        ) from exc


def _stat_signature_to_unique_index(stat: np.ndarray, stat_path: Path) -> dict[tuple, int]:
    signatures: dict[tuple, int] = {}
    for roi_index, roi_stat in enumerate(stat):
        signature = _stat_signature(roi_stat)
        if signature in signatures:
            raise ValueError(
                f"Duplicate ROI mask signature in {stat_path}: ROI {signatures[signature]} "
                f"and ROI {roi_index}. Cannot verify combined iscell mapping."
            )
        signatures[signature] = roi_index
    return signatures


def _stat_signature(roi_stat) -> tuple:
    ypix = np.asarray(_stat_value(roi_stat, "ypix"), dtype=np.int64).reshape(-1)
    xpix = np.asarray(_stat_value(roi_stat, "xpix"), dtype=np.int64).reshape(-1)
    lam = np.asarray(_stat_value(roi_stat, "lam"), dtype=np.float64).reshape(-1)
    if ypix.size != xpix.size or ypix.size != lam.size:
        raise ValueError(
            "Suite2p stat ypix/xpix/lam lengths do not match while verifying combined iscell mapping"
        )
    if ypix.size == 0:
        raise ValueError("Suite2p stat contains an empty ROI mask")
    order = np.lexsort((xpix, ypix))
    return (
        tuple(ypix[order].tolist()),
        tuple(xpix[order].tolist()),
        tuple(np.round(lam[order], 8).tolist()),
    )


def _load_all_di3_times(timing_path: Path) -> np.ndarray:
    times: list[float] = []
    di3_label_mode = _infer_di3_label_mode(timing_path)
    with open(timing_path, newline="") as file_obj:
        reader = csv.reader(file_obj)
        for row in reader:
            if _timestamped_digital_label(row, di3_label_mode) == "DI3":
                try:
                    times.append(float(row[0]))
                except ValueError as exc:
                    raise ValueError(f"Invalid DI3 timestamp in {timing_path}: {row}") from exc

    if not times:
        raise ValueError(f"No DI3 frame times found in {timing_path}")

    values = np.asarray(times, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"Non-finite DI3 frame times found in {timing_path}")
    if values.size > 1 and np.any(np.diff(values) <= 0):
        raise ValueError(f"DI3 frame times must be strictly increasing in {timing_path}")
    return values


def _infer_di3_label_mode(timing_path: Path) -> str:
    """Infer where timestamped_digital.csv stores DI3 labels.

    Most sessions use column 2 directly (timestamp, DI3/None).  Some exports
    instead put the port name in column 2 and the DI channel in column 3, e.g.
    timestamp, DIPort2, DI3.  Detect that file format from the rows instead of
    naming individual sessions.
    """

    has_primary_di3 = False
    has_secondary_di3 = False

    with open(timing_path, newline="") as file_obj:
        reader = csv.reader(file_obj)
        for row in reader:
            if len(row) >= 2 and row[1].strip() == "DI3":
                has_primary_di3 = True
            if len(row) >= 3 and row[1].strip() == "DIPort2" and row[2].strip() == "DI3":
                has_secondary_di3 = True
            if has_primary_di3:
                return "primary"

    if has_secondary_di3:
        return "dipport2_secondary"
    return "primary"


def _timestamped_digital_label(row: list[str], di3_label_mode: str) -> str | None:
    if len(row) < 2:
        return None

    primary = row[1].strip()
    if primary in {"DI3", "None"}:
        return primary

    if di3_label_mode == "dipport2_secondary" and primary == "DIPort2":
        secondary = row[2].strip() if len(row) >= 3 else ""
        if secondary == "DI3":
            return "DI3"
        if secondary in {"", "None"}:
            return "None"

    return None


def _load_frame_times(
    timing_path: Path,
    n_frames: int,
    *,
    stream_index: int = 0,
    n_streams: int = 1,
    time_zero: float | None = None,
) -> np.ndarray:
    if n_streams < 1:
        raise ValueError(f"n_streams must be positive, got {n_streams}")
    if stream_index < 0 or stream_index >= n_streams:
        raise ValueError(f"stream_index {stream_index} is outside n_streams={n_streams}")

    all_times = _load_all_di3_times(timing_path)
    selected = all_times[stream_index::n_streams]

    if selected.size != n_frames:
        raise ValueError(
            f"DI3 timing count in {timing_path} for stream {stream_index + 1}/"
            f"{n_streams} is {selected.size}, expected {n_frames}; refusing to trim data"
        )

    if time_zero is None:
        time_zero = float(all_times[0])

    return _to_session_relative_times(
        selected.astype(np.float64, copy=False),
        time_zero,
        description=f"DI3 frame times from {timing_path}",
    )


def _validate_trace_shapes(plane_path: Path, F, Fneu, spks, iscell, stat) -> None:
    if F.ndim != 2:
        raise ValueError(f"F.npy must be 2D [nROIs, nFrames]: {plane_path / 'F.npy'}")
    if Fneu.shape != F.shape:
        raise ValueError(f"Fneu.npy shape {Fneu.shape} does not match F.npy shape {F.shape}")
    if spks.shape != F.shape:
        raise ValueError(f"spks.npy shape {spks.shape} does not match F.npy shape {F.shape}")
    if iscell.ndim != 2 or iscell.shape[0] != F.shape[0]:
        raise ValueError(
            f"iscell.npy shape {iscell.shape} is incompatible with {F.shape[0]} ROIs"
        )
    if stat.ndim != 1 or stat.shape[0] != F.shape[0]:
        raise ValueError(f"stat.npy shape {stat.shape} is incompatible with {F.shape[0]} ROIs")


def _bad_frames_from_ops(ops: dict, n_frames: int) -> np.ndarray:
    bad_frames = ops.get("badframes")
    if bad_frames is None:
        return np.zeros(n_frames, dtype=bool)

    bad_frames = np.asarray(bad_frames, dtype=bool).reshape(-1)
    if bad_frames.size != n_frames:
        raise ValueError(
            f'ops["badframes"] has {bad_frames.size} entries, expected exactly {n_frames}; '
            "refusing to trim or pad bad-frame metadata"
        )
    return bad_frames


def _cell_classifier_from_iscell(iscell: np.ndarray) -> np.ndarray:
    if iscell.shape[1] >= 2:
        return iscell[:, 1].astype(np.float64, copy=False)
    warnings.warn(
        "iscell.npy has only one column; writing mpciROIs.cellClassifier as 0/1 labels.",
        RuntimeWarning,
        stacklevel=2,
    )
    return iscell[:, 0].astype(np.float64, copy=False)


def _stack_pos_from_stat(stat: np.ndarray, plane_index: int) -> np.ndarray:
    stack_pos = np.empty((stat.shape[0], 3), dtype=np.int64)
    for idx, roi_stat in enumerate(stat):
        med = np.asarray(_stat_value(roi_stat, "med"))
        if med.size < 2:
            raise ValueError(f'stat[{idx}] is missing centroid field "med"')
        y, x = med[:2]
        stack_pos[idx] = [int(round(float(y))), int(round(float(x))), plane_index]
    return stack_pos


def _mean_image_from_ops(ops: dict) -> np.ndarray:
    if "meanImg" not in ops:
        raise ValueError('ops.npy must contain "meanImg" for mpciMeanImage.images.npy')

    mean_image = np.asarray(ops["meanImg"], dtype=np.float64)
    if mean_image.ndim != 2:
        raise ValueError(f'ops["meanImg"] must be 2D, got shape {mean_image.shape}')

    return mean_image[:, :, np.newaxis, np.newaxis]


def _sparse_masks_from_stat(stat: np.ndarray, image_shape: tuple[int, int], sparse_module):
    height, width = image_shape
    coords_chunks: list[np.ndarray] = []
    data_chunks: list[np.ndarray] = []

    for roi_index, roi_stat in enumerate(stat):
        ypix = np.asarray(_stat_value(roi_stat, "ypix"), dtype=np.int64).reshape(-1)
        xpix = np.asarray(_stat_value(roi_stat, "xpix"), dtype=np.int64).reshape(-1)
        lam = np.asarray(_stat_value(roi_stat, "lam"), dtype=np.float32).reshape(-1)

        if ypix.size != xpix.size or ypix.size != lam.size:
            raise ValueError(
                f"stat[{roi_index}] ypix/xpix/lam lengths do not match: "
                f"{ypix.size}, {xpix.size}, {lam.size}"
            )
        if np.any(ypix < 0) or np.any(ypix >= height) or np.any(xpix < 0) or np.any(xpix >= width):
            raise ValueError(f"stat[{roi_index}] mask pixels exceed mean image shape {image_shape}")

        roi_coords = np.vstack(
            [
                np.full(ypix.size, roi_index, dtype=np.int64),
                ypix,
                xpix,
            ]
        )
        coords_chunks.append(roi_coords)
        data_chunks.append(lam)

    if coords_chunks:
        coords = np.concatenate(coords_chunks, axis=1)
        data = np.concatenate(data_chunks).astype(np.float32, copy=False)
    else:
        coords = np.empty((3, 0), dtype=np.int64)
        data = np.empty((0,), dtype=np.float32)

    return sparse_module.COO(coords, data, shape=(stat.shape[0], height, width))


def _save_sparse_mask(path: Path, mask, sparse_module) -> Path:
    with open(path, "wb") as file_obj:
        sparse_module.save_npz(file_obj, mask)
    return path


def _stat_value(roi_stat, key: str):
    if not isinstance(roi_stat, dict) and hasattr(roi_stat, "item"):
        try:
            roi_stat = roi_stat.item()
        except ValueError:
            pass
    if isinstance(roi_stat, dict):
        return roi_stat.get(key)
    if hasattr(roi_stat, key):
        return getattr(roi_stat, key)
    try:
        return roi_stat[key]
    except Exception as exc:
        raise ValueError(f'Suite2p stat entry is missing "{key}"') from exc


def _save_npy(path: Path, value: np.ndarray) -> Path:
    np.save(path, value)
    return path


def _save_tsv(path: Path, rows: list[tuple[int, str]]) -> Path:
    with open(path, "w", newline="") as file_obj:
        writer = csv.writer(file_obj, delimiter="\t", lineterminator="\n")
        writer.writerow(["roi_values", "roi_names"])
        writer.writerows(rows)
    return path


def _roi_type_rows() -> list[tuple[int, str]]:
    return [(0, "nonCell"), (1, "neuron")]


def _as_list(value) -> list:
    if value is None:
        return []
    if isinstance(value, (str, bytes, Path, PureWindowsPath, PurePosixPath)):
        return [value]
    try:
        return list(value)
    except TypeError:
        return [value]
