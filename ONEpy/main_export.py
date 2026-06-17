"""Find valid Suite2p runs and export session data to ALF."""

from __future__ import annotations

import argparse
import traceback
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from .calcium_export import (
        _load_dict_npy,
        _raw_sessions_from_db,
        _session_time_zero,
        _to_session_relative_times,
        single_session_calcium_export,
    )
    from .sparse_noise_export import (
        _align_movie_and_times,
        _load_frame_times as _load_sparse_frame_times,
        _load_sparse_noise_movie,
        _save_sparse_noise_raw_bin,
        _sparse_noise_events,
        single_session_sparse_noise_export,
    )
    from .s2p_path_utils import s2p_filter_run_paths, s2p_run_paths
    from .wheel_export import (
        _find_raw_file,
        _load_unwrapped_encoder_counts,
        _velocity,
        single_session_wheel_export,
    )
except ImportError:
    from calcium_export import (
        _load_dict_npy,
        _raw_sessions_from_db,
        _session_time_zero,
        _to_session_relative_times,
        single_session_calcium_export,
    )
    from sparse_noise_export import (
        _align_movie_and_times,
        _load_frame_times as _load_sparse_frame_times,
        _load_sparse_noise_movie,
        _save_sparse_noise_raw_bin,
        _sparse_noise_events,
        single_session_sparse_noise_export,
    )
    from s2p_path_utils import s2p_filter_run_paths, s2p_run_paths
    from wheel_export import (
        _find_raw_file,
        _load_unwrapped_encoder_counts,
        _velocity,
        single_session_wheel_export,
    )


def export_valid_s2p_calcium(
    root: Path,
    out: Path,
    *,
    stim_type: str,
    encoder_counts_per_revolution=65536,
    sparse_noise_edges_deg=(-135.0, 135.0, -40.0, 40.0),
) -> int:
    """Find filtered Suite2p runs under ``root`` and export them to ``out``."""

    root = Path(root)
    out = Path(out)

    print(f"Scanning for Suite2p runs under: {root}")
    ready_paths = s2p_run_paths(root)
    print(f"Found {len(ready_paths)} Suite2p-ready run path(s).")

    valid_paths, session_info, data_path_lists = s2p_filter_run_paths(ready_paths)
    print(f"Filtered to {len(valid_paths)} subject/date-compatible run path(s).")

    grouped_runs: dict[tuple[str, str, str], list[tuple[Path, list[str]]]] = defaultdict(list)
    for run_path, info, data_paths in zip(valid_paths, session_info, data_path_lists):
        subject, date, session = info
        grouped_runs[(subject, date, session)].append((Path(run_path), list(data_paths)))

    n_success = 0
    failures: list[tuple[Path, BaseException]] = []

    for session_key in sorted(grouped_runs):
        subject, date, session = session_key
        run_candidates = sorted(grouped_runs[session_key], key=lambda item: _run_preference_key(item[0]))
        if len(run_candidates) > 1:
            print(f"\n[DUPLICATES] {subject} {date} session {session}")
            for candidate, _ in run_candidates:
                print(f"  candidate: {candidate}")

        session_exported = False
        candidate_failures: list[tuple[Path, BaseException]] = []
        for run_path, data_paths in run_candidates:
            print(f"\n[EXPORT] {subject} {date} session {session}: {run_path}")
            try:
                written = []
                written.extend(single_session_calcium_export(run_path, out, stim_type=stim_type))
                if _is_multi_session_data_path_list(data_paths):
                    written.extend(
                        _export_wheel_and_sparse_for_raw_sessions(
                            run_path,
                            out,
                            stim_type=stim_type,
                            encoder_counts_per_revolution=encoder_counts_per_revolution,
                            sparse_noise_edges_deg=sparse_noise_edges_deg,
                        )
                    )
                else:
                    written.extend(
                        single_session_wheel_export(
                            run_path,
                            out,
                            encoder_counts_per_revolution=encoder_counts_per_revolution,
                        )
                    )
                    if stim_type == "sparse_noise":
                        written.extend(
                            single_session_sparse_noise_export(
                                run_path,
                                out,
                                stimulus_edges_deg=sparse_noise_edges_deg,
                            )
                        )
            except Exception as exc:
                candidate_failures.append((Path(run_path), exc))
                print(f"[FAILED] {run_path}")
                print(f"  {type(exc).__name__}: {exc}")
                traceback.print_exc()
                continue

            for alf_folder in _alf_folders_from_written(written):
                _validate_alf_timebase(alf_folder)

            n_success += 1
            session_exported = True
            print(f"[OK] Wrote {len(written)} file(s) to {out / subject / date / session / 'alf'}")
            break

        if not session_exported:
            failures.extend(candidate_failures)
            print(f"[FAILED SESSION] {subject} {date} session {session}")

    print("\nSummary")
    print(f"  Ready runs: {len(ready_paths)}")
    print(f"  Valid subject/date-compatible runs: {len(valid_paths)}")
    print(f"  Successful exports: {n_success}")
    print(f"  Failed exports: {len(failures)}")

    if failures:
        print("\nFailures")
        for run_path, exc in failures:
            print(f"  {run_path}: {type(exc).__name__}: {exc}")
        return 1

    return 0


def _is_multi_session_data_path_list(data_paths: list[str]) -> bool:
    sessions = []
    for value in data_paths:
        parts = [part for part in str(value).replace("\\", "/").split("/") if part]
        if parts:
            sessions.append(parts[-1])
    return len(set(sessions)) > 1


def _export_wheel_and_sparse_for_raw_sessions(
    run_path: Path,
    out: Path,
    *,
    stim_type: str,
    encoder_counts_per_revolution: float,
    sparse_noise_edges_deg,
) -> list[Path]:
    run_db = _load_dict_npy(run_path / "db.npy", "run-level Suite2p db.npy")
    raw_sessions = _raw_sessions_from_db(run_db, run_path)
    written: list[Path] = []

    for raw_session in raw_sessions:
        print(
            f"[EXPORT wheel] {raw_session.subject} {raw_session.date} "
            f"session {raw_session.session}"
        )
        written.extend(
            _export_wheel_raw_session(
                raw_session,
                out,
                encoder_counts_per_revolution=encoder_counts_per_revolution,
            )
        )

        if stim_type == "sparse_noise":
            try:
                written.extend(
                    _export_sparse_noise_raw_session(
                        raw_session,
                        out,
                        stimulus_edges_deg=sparse_noise_edges_deg,
                    )
                )
            except FileNotFoundError as sparse_exc:
                print(f"[SKIP sparse_noise session {raw_session.session}] {sparse_exc}")

    return written


def _export_wheel_raw_session(
    raw_session,
    out: Path,
    *,
    encoder_counts_per_revolution: float,
) -> list[Path]:
    encoder_path = _find_raw_file(raw_session.local_raw_path, "encoder_log.csv")
    timestamps, position_counts = _load_unwrapped_encoder_counts(encoder_path)
    timestamps = _to_session_relative_times(
        timestamps,
        _session_time_zero(raw_session.local_raw_path),
        description=f"wheel timestamps from {encoder_path}",
    )

    if timestamps.size == 0:
        raise ValueError(f"No wheel samples found in {encoder_path}")
    if not np.all(np.isfinite(position_counts)):
        raise ValueError(f"Wheel positions contain non-finite values: {encoder_path}")
    if timestamps.size != position_counts.size:
        raise ValueError(
            f"Wheel timestamps and positions differ in length for {encoder_path}: "
            f"{timestamps.size} vs {position_counts.size}"
        )
    if np.any(np.diff(timestamps) <= 0):
        raise ValueError(
            f"Wheel timestamps are not strictly increasing; refusing to drop duplicates: {encoder_path}"
        )

    scale = 2 * np.pi / float(encoder_counts_per_revolution)
    position = position_counts.astype(np.float64) * scale
    velocity = _velocity(position, timestamps)

    alf_folder = out / raw_session.subject / raw_session.date / raw_session.session / "alf"
    alf_folder.mkdir(parents=True, exist_ok=True)
    return [
        _save_npy(alf_folder / "_ibl_wheel.position.npy", position),
        _save_npy(alf_folder / "_ibl_wheel.timestamps.npy", timestamps.astype(np.float64)),
        _save_npy(alf_folder / "_ibl_wheel.velocity.npy", velocity),
    ]


def _export_sparse_noise_raw_session(
    raw_session,
    out: Path,
    *,
    stimulus_edges_deg,
) -> list[Path]:
    sparse_noise_path = _find_raw_file(raw_session.local_raw_path, "SparseNoise_Log.bin")
    times_grid_path = _find_raw_file(raw_session.local_raw_path, "times_grid.csv")

    stim_map = _load_sparse_noise_movie(sparse_noise_path)
    frame_times = _load_sparse_frame_times(times_grid_path, raw_session.local_raw_path)
    stim_map, frame_times = _align_movie_and_times(stim_map, frame_times)
    sparse_times, sparse_positions = _sparse_noise_events(
        stim_map,
        frame_times,
        stimulus_edges_deg,
    )

    alf_folder = out / raw_session.subject / raw_session.date / raw_session.session / "alf"
    raw_passive_folder = out / raw_session.subject / raw_session.date / raw_session.session / "raw_passive_data"

    alf_folder.mkdir(parents=True, exist_ok=True)
    raw_passive_folder.mkdir(parents=True, exist_ok=True)

    return [
        _save_npy(alf_folder / "_ibl_passiveRFM.times.npy", frame_times.astype(np.float64)),
        _save_sparse_noise_raw_bin(
            raw_passive_folder / "_iblrig_RFMapStim.raw.bin",
            sparse_noise_path,
            expected_shape=stim_map.shape,
        ),
        _save_npy(alf_folder / "_ibl_sparseNoise.times.npy", sparse_times.astype(np.float64)),
        _save_npy(alf_folder / "_ibl_sparseNoise.xy.npy", sparse_positions.astype(np.float64)),
    ]


def _save_npy(path: Path, value: np.ndarray) -> Path:
    np.save(path, value)
    return path


def _alf_folders_from_written(written: list[Path]) -> list[Path]:
    folders: set[Path] = set()
    for path in written:
        parts = Path(path).parts
        for idx, part in enumerate(parts):
            if part == "alf":
                folders.add(Path(*parts[: idx + 1]))
                break
    return sorted(folders)


def _validate_alf_timebase(alf_folder: Path) -> None:
    stim_path = alf_folder / "_ibl_passiveRFM.times.npy"
    wheel_path = alf_folder / "_ibl_wheel.timestamps.npy"
    fov_time_paths = sorted(alf_folder.glob("FOV_*/mpci.times.npy"))

    stim_times = _load_time_vector(stim_path) if stim_path.is_file() else None
    wheel_times = _load_time_vector(wheel_path) if wheel_path.is_file() else None
    fov_times = [(path.parent.name, _load_time_vector(path)) for path in fov_time_paths]

    if stim_times is not None:
        for fov_name, mpci_times in fov_times:
            if not _time_ranges_overlap(stim_times, mpci_times):
                raise ValueError(
                    f"No overlap between sparse-noise times and {fov_name} calcium times in {alf_folder}: "
                    f"stim=[{stim_times[0]:.6g}, {stim_times[-1]:.6g}], "
                    f"mpci=[{mpci_times[0]:.6g}, {mpci_times[-1]:.6g}]"
                )
        if wheel_times is not None and not _time_ranges_overlap(stim_times, wheel_times):
            raise ValueError(
                f"No overlap between sparse-noise times and wheel times in {alf_folder}: "
                f"stim=[{stim_times[0]:.6g}, {stim_times[-1]:.6g}], "
                f"wheel=[{wheel_times[0]:.6g}, {wheel_times[-1]:.6g}]"
            )


def _load_time_vector(path: Path) -> np.ndarray:
    values = np.asarray(np.load(path, allow_pickle=False), dtype=np.float64).reshape(-1)
    if values.size == 0:
        raise ValueError(f"Time vector is empty: {path}")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"Time vector contains non-finite values: {path}")
    if values.size > 1 and np.any(np.diff(values) <= 0):
        raise ValueError(f"Time vector is not strictly increasing: {path}")
    return values


def _time_ranges_overlap(a: np.ndarray, b: np.ndarray) -> bool:
    return max(float(a[0]), float(b[0])) <= min(float(a[-1]), float(b[-1]))


def _run_preference_key(run_path: Path) -> tuple[int, int, str]:
    parts = [part.lower() for part in run_path.parts]
    has_processed_component = any(part == "processed" for part in parts)
    return (0 if has_processed_component else 1, len(run_path.parts), str(run_path).lower())


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export valid Suite2p runs to ALF."
    )
    parser.add_argument(
        "root",
        type=Path,
        help="Root folder to scan for Suite2p runs.",
    )
    parser.add_argument(
        "out",
        type=Path,
        help="ALF output root.",
    )
    parser.add_argument(
        "--stim-type",
        required=True,
        choices=("sparse_noise", "gratings"),
        help="Stimulus type used to select the target raw session.",
    )
    parser.add_argument(
        "--encoder-counts-per-revolution",
        type=float,
        default=65536,
        help="Encoder counts per full wheel revolution for radians conversion.",
    )
    parser.add_argument(
        "--sparse-noise-edges-deg",
        type=float,
        nargs=4,
        default=(-135.0, 135.0, -40.0, 40.0),
        metavar=("LEFT", "RIGHT", "BOTTOM", "TOP"),
        help="Sparse-noise stimulus edges in visual degrees.",
    )
    args = parser.parse_args()

    return export_valid_s2p_calcium(
        args.root,
        args.out,
        stim_type=args.stim_type,
        encoder_counts_per_revolution=args.encoder_counts_per_revolution,
        sparse_noise_edges_deg=tuple(args.sparse_noise_edges_deg),
    )


if __name__ == "__main__":
    raise SystemExit(main())
