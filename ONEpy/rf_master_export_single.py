"""Export exactly one configured Suite2p sparse-noise session to ALF.

The public ``export_rf_session_single`` function is designed for a direct
MATLAB-to-Python call. It searches only the selected subject/date subtree,
keeps support for every Suite2p plane/FOV in that run, and never writes under
the raw 2P root.
"""

from __future__ import annotations

import contextlib
import csv
import io
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping

import numpy as np

try:
    from . import calcium_export as _calcium_export_module
    from . import wheel_export as _wheel_export_module
    from .calcium_export import (
        _as_list,
        _data_paths_with_known_edge_cases,
        _load_dict_npy,
        single_session_calcium_export,
    )
    from .main_export import (
        _alf_folders_from_written,
        _export_wheel_and_sparse_for_sparse_noise_session,
        _is_multi_session_data_path_list,
        _run_preference_key,
        _sparse_noise_session_key,
        _validate_alf_timebase,
    )
    from .s2p_path_utils import s2p_run_paths
    from .sparse_noise_export import single_session_sparse_noise_export
    from .wheel_export import single_session_wheel_export
except ImportError:  # pragma: no cover - useful for direct debugging from ONEpy/
    import calcium_export as _calcium_export_module
    import wheel_export as _wheel_export_module
    from calcium_export import (
        _as_list,
        _data_paths_with_known_edge_cases,
        _load_dict_npy,
        single_session_calcium_export,
    )
    from main_export import (
        _alf_folders_from_written,
        _export_wheel_and_sparse_for_sparse_noise_session,
        _is_multi_session_data_path_list,
        _run_preference_key,
        _sparse_noise_session_key,
        _validate_alf_timebase,
    )
    from s2p_path_utils import s2p_run_paths
    from sparse_noise_export import single_session_sparse_noise_export
    from wheel_export import single_session_wheel_export


def export_rf_session_single(
    suite2p_root: str | Path,
    raw_2p_root: str | Path,
    alf_root: str | Path,
    subject: str,
    date: str,
    session: str | int,
    log_root: str | Path,
    export_options: Mapping[str, Any],
) -> dict[str, Any]:
    """Discover and export only ``subject/date/session``.

    Required ``export_options`` keys are supplied by ``masterPipeline_single``:
    ``encoder_counts_per_revolution``, ``sparse_noise_edges_deg``,
    ``force_reexport``, and ``dry_run``.
    """

    options = dict(export_options)
    required = {
        "encoder_counts_per_revolution",
        "sparse_noise_edges_deg",
        "force_reexport",
        "dry_run",
    }
    missing = sorted(required - options.keys())
    if missing:
        raise ValueError(f"Missing export option(s): {', '.join(missing)}")

    suite2p_root = Path(suite2p_root)
    raw_2p_root = Path(raw_2p_root)
    alf_root = Path(alf_root)
    log_root = Path(log_root)
    subject = str(subject).strip()
    date = normalize_date(str(date))
    session = str(session).strip()
    encoder_counts = float(options["encoder_counts_per_revolution"])
    sparse_noise_edges = tuple(float(value) for value in options["sparse_noise_edges_deg"])
    force_reexport = bool(options["force_reexport"])
    dry_run = bool(options["dry_run"])

    if not subject or not session:
        raise ValueError("subject and session must be non-empty")
    if len(sparse_noise_edges) != 4:
        raise ValueError("sparse_noise_edges_deg must contain [left, right, bottom, top]")
    if encoder_counts <= 0:
        raise ValueError("encoder_counts_per_revolution must be positive")

    assert_not_under_raw(alf_root, raw_2p_root, "ALF root")
    assert_not_under_raw(log_root, raw_2p_root, "log root")
    assert_not_under_raw(alf_root, suite2p_root, "ALF root")
    assert_not_under_raw(log_root, suite2p_root, "log root")
    log_root.mkdir(parents=True, exist_ok=True)

    alf_session = alf_root / subject / date / session / "alf"
    print("\n========== rf_master_export_single ==========")
    print(f"Started:       {datetime.now().isoformat(timespec='seconds')}")
    print(f"Selected:      {subject}/{date}/{session}")
    print(f"Raw 2P root:   {raw_2p_root}")
    print(f"Suite2p root:  {suite2p_root}")
    print(f"ALF root:      {alf_root}")
    print(f"ALF session:   {alf_session}")
    print(f"Log root:      {log_root}")
    print(f"Force export:  {force_reexport}")
    print(f"Dry run:       {dry_run}")

    candidates = discover_session_candidates(
        suite2p_root, raw_2p_root, subject, date, session
    )
    print(f"[DISCOVER] Matching ready Suite2p run(s): {len(candidates)}")
    for index, candidate in enumerate(candidates, start=1):
        print(f"  candidate {index}: {candidate}")
    if not candidates:
        raise FileNotFoundError(
            f"No ready Suite2p sparse-noise run matched {subject}/{date}/{session} "
            f"under {[str(path) for path in selected_scan_roots(suite2p_root, subject, date)]}"
        )

    if dry_run:
        return {
            "subject": subject,
            "date": date,
            "session": session,
            "alf_root": str(alf_root),
            "alf_session": str(alf_session),
            "selected_s2p_run_path": str(candidates[0]),
            "candidate_count": len(candidates),
            "export_ok": False,
            "export_status": "dry_run",
            "export_log": "",
            "message": "Dry run: discovery only, no ALF files written.",
        }

    if alf_session_is_complete(alf_session) and not force_reexport:
        message = f"Skipped export: complete ALF session already exists at {alf_session}."
        print(f"[EXPORT SKIP] {message}")
        return {
            "subject": subject,
            "date": date,
            "session": session,
            "alf_root": str(alf_root),
            "alf_session": str(alf_session),
            "selected_s2p_run_path": str(candidates[0]),
            "candidate_count": len(candidates),
            "export_ok": True,
            "export_status": "skipped_existing",
            "export_log": "",
            "message": message,
        }

    messages: list[str] = []
    for index, run_path in enumerate(candidates, start=1):
        export_log = log_root / (
            f"{safe_filename(subject)}_{date}_{safe_filename(session)}_candidate{index:02d}.log"
        )
        print(f"\n[EXPORT {index}/{len(candidates)}] {run_path}")
        ok, message = export_one_candidate(
            run_path,
            alf_root,
            export_log,
            raw_2p_root=raw_2p_root,
            encoder_counts_per_revolution=encoder_counts,
            sparse_noise_edges_deg=sparse_noise_edges,
        )
        messages.append(message)
        if ok:
            if not alf_session_is_complete(alf_session):
                raise RuntimeError(
                    f"Export returned successfully but ALF validation failed: {alf_session}"
                )
            result = {
                "subject": subject,
                "date": date,
                "session": session,
                "alf_root": str(alf_root),
                "alf_session": str(alf_session),
                "selected_s2p_run_path": str(run_path),
                "candidate_count": len(candidates),
                "export_ok": True,
                "export_status": "ok",
                "export_log": str(export_log),
                "message": message,
            }
            print(f"[DONE] Single-session ALF export validated: {alf_session}")
            return result

    raise RuntimeError(
        f"All {len(candidates)} export candidate(s) failed for "
        f"{subject}/{date}/{session}: {' || '.join(messages)}"
    )


def selected_scan_roots(suite2p_root: Path, subject: str, date: str) -> list[Path]:
    """Return supported narrow date roots; never scan the full subject tree."""

    compact_date = normalize_date(date).replace("-", "")
    subject_root = Path(suite2p_root) / subject
    candidates = (
        subject_root / compact_date,
        subject_root / normalize_date(date),
        subject_root / "Processed" / compact_date,
        subject_root / "Processed" / normalize_date(date),
    )
    roots: list[Path] = []
    seen: set[str] = set()
    for path in candidates:
        key = normalized_path_text(path)
        if key not in seen:
            roots.append(path)
            seen.add(key)
    return roots


def discover_session_candidates(
    suite2p_root: Path,
    raw_2p_root: Path,
    subject: str,
    date: str,
    session: str,
) -> list[Path]:
    """Find ready Suite2p runs only inside the configured subject/date tree."""

    matches: list[Path] = []
    existing_roots = [path for path in selected_scan_roots(suite2p_root, subject, date) if path.is_dir()]
    if not existing_roots:
        raise FileNotFoundError(
            "No selected subject/date folder exists at any supported location: "
            + ", ".join(str(path) for path in selected_scan_roots(suite2p_root, subject, date))
        )
    for scan_root in existing_roots:
        print(f"[DISCOVER] Scanning selected subject/date only: {scan_root}")
        for run_path_value in s2p_run_paths(scan_root):
            run_path = Path(run_path_value)
            try:
                with raw_2p_root_override(raw_2p_root):
                    found_subject, found_date, found_session = _sparse_noise_session_key(run_path)
            except Exception as exc:
                print(f"[EXCLUDED] {run_path} | sparse-noise session unresolved: {exc}")
                continue

            found_key = (str(found_subject), normalize_date(str(found_date)), str(found_session))
            wanted_key = (subject, normalize_date(date), str(session))
            if found_key == wanted_key:
                matches.append(run_path)
                print(f"[TARGET] {found_subject}/{found_date}/{found_session}: {run_path}")
            else:
                print(
                    f"[EXCLUDED] {run_path} | resolved {found_key[0]}/{found_key[1]}/{found_key[2]}, "
                    f"wanted {wanted_key[0]}/{wanted_key[1]}/{wanted_key[2]}"
                )

    return sorted(set(matches), key=_run_preference_key)


def export_one_candidate(
    run_path: Path,
    alf_root: Path,
    export_log: Path,
    *,
    raw_2p_root: Path,
    encoder_counts_per_revolution: float,
    sparse_noise_edges_deg: tuple[float, float, float, float],
) -> tuple[bool, str]:
    """Run one established export implementation, log it, and echo its output."""

    export_log.parent.mkdir(parents=True, exist_ok=True)
    buffer = io.StringIO()
    ok = False
    message = ""
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        print("[EXPORT LOG]")
        print(f"Started: {datetime.now().isoformat(timespec='seconds')}")
        print(f"Run path: {run_path}")
        print(f"ALF root: {alf_root}")
        try:
            written = export_sparse_noise_run(
                run_path,
                alf_root,
                raw_2p_root=raw_2p_root,
                encoder_counts_per_revolution=encoder_counts_per_revolution,
                sparse_noise_edges_deg=sparse_noise_edges_deg,
            )
            ok = True
            message = f"Wrote {len(written)} file(s)."
            print(f"[OK] {message}")
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            print(f"[FAILED] {message}")
            traceback.print_exc()
        finally:
            print(f"Finished: {datetime.now().isoformat(timespec='seconds')}")

    captured = buffer.getvalue()
    export_log.write_text(captured, encoding="utf-8")
    print(captured, end="")
    print(f"[{'OK' if ok else 'FAILED'}] {message} Log: {export_log}")
    return ok, message


def export_sparse_noise_run(
    run_path: Path,
    alf_root: Path,
    *,
    raw_2p_root: Path,
    encoder_counts_per_revolution: float,
    sparse_noise_edges_deg: tuple[float, float, float, float],
) -> list[Path]:
    """Reuse the established calcium, wheel, sparse-noise, and validation logic."""

    run_path = Path(run_path)
    alf_root = Path(alf_root)
    with raw_2p_root_override(raw_2p_root):
        run_db = _load_dict_npy(run_path / "db.npy", "run-level Suite2p db.npy")
        data_paths = [
            str(value)
            for value in _data_paths_with_known_edge_cases(_as_list(run_db.get("data_path")), run_path)
        ]

        written: list[Path] = []
        written.extend(single_session_calcium_export(run_path, alf_root, stim_type="sparse_noise"))
        if _is_multi_session_data_path_list(data_paths):
            written.extend(
                _export_wheel_and_sparse_for_sparse_noise_session(
                    run_path,
                    alf_root,
                    encoder_counts_per_revolution=encoder_counts_per_revolution,
                    sparse_noise_edges_deg=sparse_noise_edges_deg,
                )
            )
        else:
            written.extend(
                single_session_wheel_export(
                    run_path,
                    alf_root,
                    encoder_counts_per_revolution=encoder_counts_per_revolution,
                )
            )
            written.extend(
                single_session_sparse_noise_export(
                    run_path,
                    alf_root,
                    stimulus_edges_deg=sparse_noise_edges_deg,
                )
            )

        for alf_folder in _alf_folders_from_written(written):
            _validate_alf_timebase(alf_folder)
        return written


@contextlib.contextmanager
def raw_2p_root_override(raw_2p_root: str | Path):
    """Route existing exporter path remapping to the configured raw-data root.

    The established helpers infer a 2P root from the Suite2p run path. On this
    installation, processed runs and raw acquisitions live in separate trees,
    so only that inferred root is temporarily supplied explicitly.
    """

    configured_root = Path(raw_2p_root)
    original_calcium_infer = _calcium_export_module._infer_local_2p_root
    original_wheel_infer = _wheel_export_module._infer_local_2p_root
    original_load_all_di3 = _calcium_export_module._load_all_di3_times
    original_load_pairs = _calcium_export_module._load_di3_none_pairs

    def configured_infer(_s2p_run_path: Path) -> Path:
        return configured_root

    _calcium_export_module._infer_local_2p_root = configured_infer
    _wheel_export_module._infer_local_2p_root = configured_infer
    _calcium_export_module._load_all_di3_times = _load_all_di3_times_compatible
    _calcium_export_module._load_di3_none_pairs = _load_di3_none_pairs_compatible
    try:
        yield
    finally:
        _calcium_export_module._infer_local_2p_root = original_calcium_infer
        _wheel_export_module._infer_local_2p_root = original_wheel_infer
        _calcium_export_module._load_all_di3_times = original_load_all_di3
        _calcium_export_module._load_di3_none_pairs = original_load_pairs


def _timing_row_label(row: list[str]) -> str | None:
    """Accept both legacy label-column and current DIPort2/third-column CSVs."""

    if len(row) < 2:
        return None
    second = row[1].strip()
    if second in {"DI3", "None"}:
        return second
    if second == "DIPort2":
        if len(row) >= 3 and row[2].strip() == "DI3":
            return "DI3"
        return "None"
    return None


def _load_all_di3_times_compatible(timing_path: Path) -> np.ndarray:
    times: list[float] = []
    with open(timing_path, newline="") as file_obj:
        for row in csv.reader(file_obj):
            if _timing_row_label(row) != "DI3":
                continue
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


def _load_di3_none_pairs_compatible(timing_path: Path) -> tuple[np.ndarray, np.ndarray]:
    di3_times: list[float] = []
    none_times: list[float] = []
    pending_di3: float | None = None
    with open(timing_path, newline="") as file_obj:
        for row in csv.reader(file_obj):
            label = _timing_row_label(row)
            if label is None:
                continue
            try:
                timestamp = float(row[0])
            except ValueError as exc:
                raise ValueError(f"Invalid {label} timestamp in {timing_path}: {row}") from exc
            if not np.isfinite(timestamp):
                raise ValueError(f"Non-finite {label} timestamp in {timing_path}: {row}")
            if label == "DI3":
                pending_di3 = timestamp
            elif pending_di3 is not None:
                di3_times.append(pending_di3)
                none_times.append(timestamp)
                pending_di3 = None

    di3_values = np.asarray(di3_times, dtype=np.float64)
    none_values = np.asarray(none_times, dtype=np.float64)
    if di3_values.size and np.any(np.diff(di3_values) <= 0):
        raise ValueError(f"DI3/None pair DI3 times must be strictly increasing in {timing_path}")
    return di3_values, none_values


def alf_session_is_complete(alf_session: Path) -> bool:
    """Return True when one ALF session has inputs required by both mappers."""

    alf_session = Path(alf_session)
    session_dir = alf_session.parent
    required_session_files = (
        alf_session / "_ibl_passiveRFM.times.npy",
        alf_session / "_ibl_sparseNoise.times.npy",
        alf_session / "_ibl_sparseNoise.xy.npy",
        alf_session / "_ibl_wheel.position.npy",
        alf_session / "_ibl_wheel.timestamps.npy",
        alf_session / "_ibl_wheel.velocity.npy",
        session_dir / "raw_passive_data" / "_iblrig_RFMapStim.raw.bin",
    )
    if not all(path.is_file() for path in required_session_files):
        return False

    fov_dirs = sorted(path for path in alf_session.glob("FOV_*") if path.is_dir())
    if not fov_dirs:
        return False
    required_fov_files = (
        "mpci.times.npy",
        "mpci.ROIActivityF.npy",
        "mpci.ROINeuropilActivityF.npy",
        "mpci.ROIActivityDeconvolved.npy",
    )
    return all(all((fov_dir / name).is_file() for name in required_fov_files) for fov_dir in fov_dirs)


def normalize_date(value: str) -> str:
    text = str(value).strip()
    digits = text.replace("-", "")
    if len(digits) == 8 and digits.isdigit():
        return f"{digits[0:4]}-{digits[4:6]}-{digits[6:8]}"
    return text


def safe_filename(value: str) -> str:
    label = "".join(char if char.isalnum() or char in ("-", "_") else "_" for char in str(value))
    return label.strip("_") or "unknown"


def assert_not_under_raw(path: Path, raw_root: Path, label: str) -> None:
    path_text = normalized_path_text(path)
    raw_text = normalized_path_text(raw_root)
    if path_text == raw_text or path_text.startswith(raw_text.rstrip("/\\") + "\\"):
        raise ValueError(f"Refusing to use {label} under read-only raw root: {path}")


def normalized_path_text(path: Path) -> str:
    return str(path).replace("/", "\\").rstrip("\\").lower()
