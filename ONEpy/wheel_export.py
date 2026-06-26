"""Export rotary encoder wheel data to ALF/ONE datasets."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

try:
    from .calcium_export import (
        _as_list,
        _date_from_run_path,
        _infer_local_2p_root,
        _load_dict_npy,
        _local_raw_path,
        _parse_subject_date_session,
        _session_time_zero,
        _to_session_relative_times,
    )
except ImportError:
    from calcium_export import (
        _as_list,
        _date_from_run_path,
        _infer_local_2p_root,
        _load_dict_npy,
        _local_raw_path,
        _parse_subject_date_session,
        _session_time_zero,
        _to_session_relative_times,
    )


DEFAULT_ENCODER_COUNTS_PER_REVOLUTION = 65536


def single_session_wheel_export(
    s2prunpath,
    dataoutroot,
    *,
    encoder_counts_per_revolution=DEFAULT_ENCODER_COUNTS_PER_REVOLUTION,
):
    """Export one session's rotary encoder data to ALF wheel datasets.

    The raw ``encoder_log.csv`` second column is treated as a signed 16-bit
    wrapped encoder position. It is unwrapped and converted to radians.
    Samples are never silently dropped; invalid or duplicate timestamps raise.
    """

    s2p_run_path = Path(s2prunpath)
    data_out_root = Path(dataoutroot)
    subject, date, session, raw_path = _session_from_run_db(s2p_run_path)

    encoder_path = _find_raw_file(raw_path, "encoder_log.csv")
    timestamps, position_counts = _load_unwrapped_encoder_counts(encoder_path)
    timestamps = _to_session_relative_times(
        timestamps,
        _session_time_zero(raw_path),
        description=f"wheel timestamps from {encoder_path}",
    )
    _validate_wheel_samples(timestamps, position_counts, encoder_path)

    scale = 2 * np.pi / float(encoder_counts_per_revolution)
    position = position_counts.astype(np.float64) * scale
    velocity = _velocity(position, timestamps)

    alf_folder = data_out_root / subject / date / session / "alf"
    alf_folder.mkdir(parents=True, exist_ok=True)

    written = [
        _save_npy(alf_folder / "_ibl_wheel.position.npy", position),
        _save_npy(alf_folder / "_ibl_wheel.timestamps.npy", timestamps.astype(np.float64, copy=False)),
        _save_npy(alf_folder / "_ibl_wheel.velocity.npy", velocity),
    ]
    return written


def _session_from_run_db(s2p_run_path: Path) -> tuple[str, str, str, Path]:
    run_db = _load_dict_npy(s2p_run_path / "db.npy", "run-level Suite2p db.npy")
    data_paths = _as_list(run_db.get("data_path"))
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
    sessions = {item[2] for item in parsed}
    if len(subjects) != 1 or len(dates) != 1 or len(sessions) != 1:
        raise ValueError(
            "Suite2p run must correspond to exactly one subject/date/session; "
            f"subjects={sorted(subjects)}, dates={sorted(dates)}, sessions={sorted(sessions)}"
        )

    local_2p_root = _infer_local_2p_root(s2p_run_path)
    raw_path = _local_raw_path(
        data_paths[0],
        local_2p_root,
        fallback_date_folder=run_date_folder,
    )
    subject, date, session = parsed[0]
    return subject, date, session, raw_path


def _find_raw_file(raw_path: Path, filename: str) -> Path:
    candidates = [raw_path / filename, raw_path.parent / filename]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"Could not find {filename}. Tried: {[str(path) for path in candidates]}")


def _session_start_time(raw_path: Path) -> float:
    """Backward-compatible alias for the shared session time-zero helper."""

    return _session_time_zero(raw_path)


def _first_numeric_csv_value(path: Path) -> float | None:
    with open(path, newline="") as file_obj:
        reader = csv.reader(file_obj)
        for row in reader:
            if not row:
                continue
            try:
                return float(row[0])
            except ValueError:
                continue
    return None


def _load_unwrapped_encoder_counts(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, delimiter=",", ndmin=2)
    if data.shape[1] < 2:
        raise ValueError(f"encoder_log.csv must have at least two columns: {path}")

    timestamps = data[:, 0].astype(np.float64, copy=False)
    raw_position = data[:, 1].astype(np.float64, copy=False)

    deltas = np.diff(raw_position)
    deltas[deltas > 32768] -= 65536
    deltas[deltas < -32768] += 65536
    unwrapped = np.concatenate(([raw_position[0]], raw_position[0] + np.cumsum(deltas)))

    return timestamps, unwrapped


def _validate_wheel_samples(timestamps: np.ndarray, position_counts: np.ndarray, path: Path) -> None:
    if timestamps.ndim != 1 or position_counts.ndim != 1:
        raise ValueError(f"Wheel timestamps and positions must be vectors: {path}")
    if timestamps.size != position_counts.size:
        raise ValueError(
            f"Wheel timestamps and positions differ in length for {path}: "
            f"{timestamps.size} vs {position_counts.size}"
        )
    if timestamps.size == 0:
        raise ValueError(f"No wheel samples found in {path}")
    if not np.all(np.isfinite(timestamps)) or not np.all(np.isfinite(position_counts)):
        raise ValueError(f"Wheel data contain non-finite values: {path}")
    if timestamps.size > 1 and np.any(np.diff(timestamps) <= 0):
        raise ValueError(
            f"Wheel timestamps are not strictly increasing; refusing to drop duplicates: {path}"
        )


def _velocity(position: np.ndarray, timestamps: np.ndarray) -> np.ndarray:
    velocity = np.zeros_like(position, dtype=np.float64)
    if timestamps.size < 2:
        return velocity

    dt = np.diff(timestamps)
    dp = np.diff(position)
    step_velocity = np.divide(
        dp,
        dt,
        out=np.full_like(dp, np.nan, dtype=np.float64),
        where=dt > 0,
    )
    velocity[1:] = step_velocity

    finite = np.isfinite(velocity)
    if not np.any(finite):
        return np.zeros_like(position, dtype=np.float64)

    first_good = np.flatnonzero(finite)[0]
    velocity[:first_good] = velocity[first_good]
    for idx in range(first_good + 1, velocity.size):
        if not np.isfinite(velocity[idx]):
            velocity[idx] = velocity[idx - 1]
    return velocity


def _save_npy(path: Path, value: np.ndarray) -> Path:
    np.save(path, value)
    return path
