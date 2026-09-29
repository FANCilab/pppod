r"""Selected RF plotting script Final for per-FOV, combined-session, and group outputs.

This script reads RF analysis results written by ALF_RF_pipeline.m:
    D:\Pipeline\Results\SUBJECT\DATE\SESSION\alf\FOV_*

The plotting groups are supplied by ALF_RF_pipeline.m from subject membership
only.  Analysis results are not expected to live inside group folders.

All figures are written to:
    D:\Pipeline\Plots\group*

Every plot is saved as both PNG and SVG. Per-neuron RF maps are split into peak/noise PASS and FAIL folders. RF-size histograms are made from valid Gaussian sigma parameters.

RF map image display note:
    MATLAB saves RF edges as [left, right, top, bottom], e.g. [-135 135 40 -40].
    Matplotlib imshow expects extent=(left, right, bottom, top).  The helper
    rf_image_extent() below performs this conversion so the RF pixels and
    Gaussian fit overlay use the same elevation convention.
"""

from __future__ import annotations

import csv
import json
import math
import re
import sys
import textwrap
import warnings
import argparse
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "text.usetex": False,
        "axes.grid": False,
    }
)
warnings.filterwarnings("ignore", message=r".*will be ignored.*", category=UserWarning)

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import LinearLocator
from matplotlib.transforms import blended_transform_factory

try:
    import h5py
except Exception:  # pragma: no cover - only used when h5py is absent on user machine
    h5py = None

try:
    from pycolorbar.bivariate import BivariateColormap
except Exception:  # pragma: no cover - pycolorbar version compatibility
    try:
        from pycolorbar.bivariate.cmap import BivariateColormap
    except Exception:  # pragma: no cover - optional plotting dependency on user machine
        BivariateColormap = None


# =============================================================================
# Session-independent configuration
# =============================================================================

PIPELINE_ROOT = Path(r"D:\Pipeline")
RESULTS_ROOT = PIPELINE_ROOT / "FinalResults"
PLOTS_ROOT = PIPELINE_ROOT / "FinalPlots"
# Plot only this group. Set to None to plot all configured groups.
SELECTED_GROUP: str | None = None
REMOTE_2P_ROOT = Path(r"\\10.233.25.135\FANCiNAS1\Data\2P")
# Set to None to process every FOV_* found under the ALF session.
FOV_NAMES: list[str] | None = None

GROUP_SUBJECTS: dict[str, set[str]] = {}
SUBJECT_GROUP: dict[str, str] = {}

MIN_EV = 0.01
MAX_P_VALUE = 0.05
MIN_PEAK_TO_NOISE = 7.7
MAX_NEURON_RF_PLOTS: int | None = 12  # None = plot all significant RF neurons
AGGREGATE_EV_X_LIMITS = (0.0, 0.3)
PLOT_SELECTED_GOOD_RF_MAPS_ONLY = True
SELECTED_GOOD_RF_MAPS_PER_GROUP = 12
SELECTED_GOOD_RF_EDGE_MARGIN_FRACTION = 0.08
SELECTED_GOOD_RF_EDGE_MARGIN_DEG = 5.0
SELECTED_GOOD_RF_COLLECTION_FOLDER = "Selected_good_RF_maps"

MAKE_GAUSSIAN_COVERAGE = True
MAKE_PER_NEURON_RF_MAPS = False
MAKE_FOV_OVERLAYS = True
MAKE_NEUTRAL_ROI_OVERLAYS = True
MAKE_SHIFT_EV_HISTOGRAM = True
MAKE_PREDICTED_VS_ACTUAL = True
MAKE_PEAK_TIME_COURSES = True
MAKE_RF_SIZE_HISTOGRAM = True
MAKE_COMBINED_SESSION_PLOTS = True
MAKE_COMBINED_FOV_OVERLAYS = True

# Trace panel downsampling for very dense plots. 1 keeps every point.
TRACE_PLOT_STRIDE = 1

# pycolorbar is used for RF-center colors in the FOV overlay plots and their
# 2D color legend. Install with: pip install pycolorbar
PYCOLORBAR_BIVARIATE_CMAP_NAME = "teuling.GRMB"
PYCOLORBAR_TEULING_DIAGONAL_TILT = 0.37
PYCOLORBAR_TEULING_OFFDIAG_TILT = 1.0
BIVARIATE_CMAP_WHEEL_RESOLUTION = 301
# FOV overlay 2D pycolorbar plots use this azimuth plotting/normalization range
# for the color map and color wheel only. Set to None to use the RF result edges.
# This does not modify the underlying RF-center data.
FOV_OVERLAY_BIVARIATE_AZIMUTH_LIMITS: tuple[float, float] | None = (-40.0, 130.0)


# =============================================================================
# Styling
# =============================================================================

PEAK_ROI_ALPHA = 1.0
ROI_MARKER_SIZE = 3
NON_PEAK_CONTOUR_LINEWIDTH = 1.1
NEUTRAL_ROI_ALPHA = 0.1
FOV_BLACK_FLOOR_PERCENTILE = 20.0

RF_ON_CMAP = LinearSegmentedColormap.from_list(
    "rf_on",
    ["#25456f", "#ffffff", "#c7342e"],
    N=256,
)
RF_OFF_CMAP = LinearSegmentedColormap.from_list(
    "rf_off",
    ["#61223b", "#ffffff", "#1f6fba"],
    N=256,
)
RF_BWR_CMAP = LinearSegmentedColormap.from_list(
    "rf_bwr",
    ["#244c9a", "#ffffff", "#c42c2c"],
    N=256,
)
RF_MAGENTA_GREEN_CMAP = LinearSegmentedColormap.from_list(
    "rf_magenta_green",
    ["#61223b", "#ffffff", "#1a9850"],
    N=256,
)
# Per-neuron RF-map colors from v4.
RF_NEGATIVE_GRAY = "#b8b8b8"
RF_ZERO_COLOR = "#ffffff"
RF_BRIGHT_GREEN = "#00ff00"
RF_PEAK_GREEN = "#126d36"
RF_PER_NEURON_CMAP = LinearSegmentedColormap.from_list(
    "rf_per_neuron_signed_lightgray_white_010_green_darkgreen",
    [
        (0.00, RF_NEGATIVE_GRAY),  # negative values: lighter gray, not black
        (0.50, RF_ZERO_COLOR),     # zero: white
        (0.83, RF_BRIGHT_GREEN),   # strong positive: exact RGB 010 green
        (1.00, RF_PEAK_GREEN),     # peak positive: same dark green as v4
    ],
    N=256,
)

COMBINED_FOV_GAP_PIXELS = 0
COMBINED_CELL_ROI_FILLED_ALPHA = 0.5
COMBINED_CELL_ROI_CONTOUR_LINEWIDTH = 0.7
ROI_CONTOUR_OVERLAY_DPI = 300


# =============================================================================
# Main orchestration
# =============================================================================


def main() -> int:
    global RESULTS_ROOT, PLOTS_ROOT, SELECTED_GROUP

    args = parse_cli_args()
    if args.results_root is not None:
        RESULTS_ROOT = Path(args.results_root)
    if args.plots_root is not None:
        PLOTS_ROOT = Path(args.plots_root)

    try:
        configure_group_subjects(parse_group_subjects_arg(args.group_subjects))
    except ValueError as exc:
        print(f"[ERROR] Invalid --group-subjects: {exc}")
        return 2

    if args.group is not None:
        group = args.group.strip()
        SELECTED_GROUP = None if group.lower() == "all" else group
    if SELECTED_GROUP is not None and SELECTED_GROUP not in GROUP_SUBJECTS:
        available = ", ".join(sorted(GROUP_SUBJECTS))
        print(f"[ERROR] Unknown --group {SELECTED_GROUP!r}. Configured groups: {available}")
        return 2

    flip_lookup = parse_flip_fov_args(args.flip_fov)
    sessions = discover_rf_sessions(RESULTS_ROOT)
    sessions = apply_flip_lookup_to_sessions(sessions, flip_lookup)

    if SELECTED_GROUP is not None:
        sessions = [info for info in sessions if info["group"] == SELECTED_GROUP]
    group_rows: dict[str, list[dict]] = {}
    group_ev_rows: dict[str, list[dict]] = {}
    selected_per_neuron_keys: set[tuple[str, str, str, str, str, int]] | None = None

    print("\n========== Selected FANCi RF plotting ==========")
    print(f"Results root:    {RESULTS_ROOT}")
    print(f"Plots root:      {PLOTS_ROOT}")
    print(f"Group filter:    {SELECTED_GROUP if SELECTED_GROUP is not None else 'all groups'}")
    print(f"Configured groups: {', '.join(f'{group}({len(subjects)})' for group, subjects in sorted(GROUP_SUBJECTS.items()))}")
    print(f"Flip FOV specs:  {sum(len(fovs) for fovs in flip_lookup.values())}")
    print(f"Skip existing:   {'no, force replot requested' if args.force_replot else 'yes'}")
    print(f"Sessions:     {len(sessions)}")
    if PLOT_SELECTED_GOOD_RF_MAPS_ONLY:
        print(
            "Per-neuron RF maps: selected good RFs only "
            f"({SELECTED_GOOD_RF_MAPS_PER_GROUP} diverse center locations per group)"
        )
        group_rows, group_ev_rows = collect_group_aggregate_rows_for_sessions(sessions)
        selected_rows_by_group = select_diverse_good_rf_rows_by_group(
            group_rows,
            n_per_group=SELECTED_GOOD_RF_MAPS_PER_GROUP,
            edge_margin_fraction=SELECTED_GOOD_RF_EDGE_MARGIN_FRACTION,
            edge_margin_deg=SELECTED_GOOD_RF_EDGE_MARGIN_DEG,
        )
        selected_per_neuron_keys = {
            per_neuron_selection_key_from_row(row)
            for rows in selected_rows_by_group.values()
            for row in rows
        }
        for group, rows in sorted(selected_rows_by_group.items()):
            output_dir = PLOTS_ROOT / group / SELECTED_GOOD_RF_COLLECTION_FOLDER
            output_dir.mkdir(parents=True, exist_ok=True)
            write_metadata(rows, output_dir / f"{group}_selected_good_RF_maps_for_per_neuron_plotting.csv")
            print(f"Selected per-neuron RF maps for {group}: {len(rows)}")
    else:
        print("Per-neuron RF maps: all significant RFs")

    for info in sessions:
        group = info["group"]
        subject = info["subject"]
        date = info["date"]
        session = info["session"]
        flip_fovs = set(info.get("flip_fovs", set()))
        results_root = Path(info.get("results_root", RESULTS_ROOT))

        summary = plot_rf_session(
            results_root=results_root,
            subject=subject,
            date=date,
            session=session,
            plots_root=PLOTS_ROOT / group,
            group=group,
            flip_fovs=flip_fovs,
            min_ev=MIN_EV,
            max_p_value=MAX_P_VALUE,
            min_peak_to_noise=MIN_PEAK_TO_NOISE,
            max_neuron_plots=MAX_NEURON_RF_PLOTS,
            selected_per_neuron_keys=selected_per_neuron_keys,
            selected_good_rf_collection_dir=(
                PLOTS_ROOT / group / SELECTED_GOOD_RF_COLLECTION_FOLDER
                if selected_per_neuron_keys is not None
                else None
            ),
            skip_fov_overlays=False,
            make_combined_session_plots=MAKE_COMBINED_SESSION_PLOTS,
            skip_existing=not args.force_replot,
        )
        if not PLOT_SELECTED_GOOD_RF_MAPS_ONLY:
            group_rows.setdefault(group, []).extend(summary["aggregate_rows"])
            group_ev_rows.setdefault(group, []).extend(summary["aggregate_ev_rows"])

    for group, rows in sorted(group_rows.items()):
        plot_group_aggregates(
            group,
            rows,
            group_ev_rows.get(group, []),
            PLOTS_ROOT / group / "Aggregate",
            min_ev=MIN_EV,
            max_p_value=MAX_P_VALUE,
            min_peak_to_noise=MIN_PEAK_TO_NOISE,
        )

    plot_subject_debug_aggregates(
        sessions,
        group_rows,
        PLOTS_ROOT,
        min_ev=MIN_EV,
        max_p_value=MAX_P_VALUE,
        min_peak_to_noise=MIN_PEAK_TO_NOISE,
    )

    return 0


def parse_cli_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot FANCi RF analysis outputs.")
    parser.add_argument("--results-root", type=str, default=None, help="Root containing subject/date/session/alf results.")
    parser.add_argument("--plots-root", type=str, default=None, help="Root for grouped plot outputs.")
    parser.add_argument("--group", type=str, default=None, help="Group to plot, e.g. groupA, groupB, groupC, groupD, or all.")
    parser.add_argument(
        "--group-subjects",
        type=str,
        required=True,
        help="Plot group definitions, formatted as groupA=SUBJECT1,SUBJECT2;groupB=SUBJECT3.",
    )
    parser.add_argument(
        "--flip-fov",
        action="append",
        default=[],
        help="FOV to plot with flipped RF azimuth, formatted as SUBJECT/DATE/SESSION/FOV_XX. Can be repeated.",
    )
    parser.add_argument("--force-replot", action="store_true", help="Rebuild session plots even when plot_summary.json already matches.")
    return parser.parse_args()


def parse_group_subjects_arg(value: str) -> dict[str, set[str]]:
    group_subjects: dict[str, set[str]] = {}
    text = str(value).strip()
    if not text:
        raise ValueError("expected at least one group definition")

    for group_spec in text.split(";"):
        group_spec = group_spec.strip()
        if not group_spec:
            continue
        if "=" not in group_spec:
            raise ValueError(f"{group_spec!r} is missing '='")
        group, subjects_text = group_spec.split("=", 1)
        group = group.strip()
        if not group:
            raise ValueError(f"{group_spec!r} has an empty group name")
        subjects = {subject.strip() for subject in subjects_text.split(",") if subject.strip()}
        if not subjects:
            raise ValueError(f"{group!r} has no subjects")
        group_subjects[group] = subjects

    if not group_subjects:
        raise ValueError("expected at least one non-empty group definition")
    return group_subjects


def configure_group_subjects(group_subjects: dict[str, set[str]]) -> None:
    global GROUP_SUBJECTS, SUBJECT_GROUP

    subject_group: dict[str, str] = {}
    for group, subjects in group_subjects.items():
        for subject in subjects:
            previous_group = subject_group.get(subject)
            if previous_group is not None:
                raise ValueError(f"subject {subject!r} appears in both {previous_group!r} and {group!r}")
            subject_group[subject] = group

    GROUP_SUBJECTS = {group: set(subjects) for group, subjects in group_subjects.items()}
    SUBJECT_GROUP = subject_group


def parse_flip_fov_args(values: list[str]) -> dict[tuple[str, str, str], set[str]]:
    lookup: dict[tuple[str, str, str], set[str]] = {}
    for raw_value in values or []:
        text = str(raw_value).strip().strip('"').strip("'")
        if not text:
            continue
        parts = [part for part in re.split(r"[\\/]+", text) if part]
        if len(parts) != 4:
            print(f"[WARN] Ignoring --flip-fov {raw_value!r}; expected SUBJECT/DATE/SESSION/FOV_XX")
            continue
        subject, date, session, fov_name = parts
        fov_name = normalize_fov_name(fov_name)
        if not fov_name.startswith("FOV_"):
            print(f"[WARN] Ignoring --flip-fov {raw_value!r}; FOV must look like FOV_XX")
            continue
        key = (subject.strip(), normalize_date_label(date), str(session).strip())
        lookup.setdefault(key, set()).add(fov_name)
    return lookup


def apply_flip_lookup_to_sessions(
    sessions: list[dict[str, Any]],
    flip_lookup: dict[tuple[str, str, str], set[str]],
) -> list[dict[str, Any]]:
    matched: set[tuple[str, str, str]] = set()
    out: list[dict[str, Any]] = []
    for info in sessions:
        key = (str(info["subject"]), normalize_date_label(info["date"]), str(info["session"]))
        updated = dict(info)
        if key in flip_lookup:
            updated["flip_fovs"] = set(updated.get("flip_fovs", set())) | set(flip_lookup[key])
            matched.add(key)
        out.append(updated)
    unmatched = set(flip_lookup) - matched
    if unmatched:
        print(f"[INFO] Flip FOV specs not used for currently discovered sessions: {len(unmatched)}")
    return out


def discover_rf_sessions(results_root: Path) -> list[dict[str, str]]:
    sessions = []
    for alf_dir in sorted(results_root.glob("*/*/*/alf")):
        fov_dirs = [
            path
            for path in alf_dir.iterdir()
            if path.is_dir()
            and path.name.upper().startswith("FOV_")
            and (path / "_FANCi_rf.maps.npy").exists()
        ]
        if not fov_dirs:
            continue
        session_dir = alf_dir.parent
        date_dir = session_dir.parent
        subject_dir = date_dir.parent
        group = group_for_subject(subject_dir.name)
        if group is None:
            print(f"[WARN] Skipping ungrouped subject under results root: {subject_dir.name}")
            continue
        sessions.append(
            {
                "group": group,
                "subject": subject_dir.name,
                "date": normalize_date_label(date_dir.name),
                "session": session_dir.name,
                "results_root": results_root,
            }
        )
    return sorted(sessions, key=lambda item: (item["group"], item["subject"], item["date"], item["session"]))


def collect_group_aggregate_rows_for_sessions(
    sessions: list[dict[str, str]],
) -> tuple[dict[str, list[dict]], dict[str, list[dict]]]:
    group_rows: dict[str, list[dict]] = {}
    group_ev_rows: dict[str, list[dict]] = {}
    for info in sessions:
        group = info["group"]
        subject = info["subject"]
        date = info["date"]
        session = info["session"]
        results_root = Path(info.get("results_root", RESULTS_ROOT))
        results_session = results_root / subject / date / session / "alf"
        if FOV_NAMES is not None:
            fov_names = FOV_NAMES
        else:
            fov_names = discover_fovs(results_session, results_session)

        records: list[dict] = []
        flip_fovs = set(info.get("flip_fovs", set()))
        for fov_name in fov_names:
            fov_results_dir = results_session / fov_name
            if not fov_results_dir.is_dir():
                continue
            try:
                records.append(
                    {
                        "fov_name": fov_name,
                        "result": load_fanci_result(fov_results_dir),
                        "flip_azimuth": fov_name in flip_fovs,
                    }
                )
            except Exception as exc:
                print(f"[WARN] Aggregate pre-scan skipped {group}/{subject}/{date}/{session}/{fov_name}: {exc}")
        if not records:
            continue

        rows = session_aggregate_rows(
            records,
            group=group,
            subject=subject,
            date=date,
            session=session,
            min_ev=MIN_EV,
            max_p_value=MAX_P_VALUE,
            min_peak_to_noise=MIN_PEAK_TO_NOISE,
        )
        group_rows.setdefault(group, []).extend(rows)
        ev_rows = session_ev_aggregate_rows(
            records,
            group=group,
            subject=subject,
            date=date,
            session=session,
            max_p_value=MAX_P_VALUE,
        )
        group_ev_rows.setdefault(group, []).extend(ev_rows)
    return group_rows, group_ev_rows


def group_for_subject(subject: str) -> str | None:
    return SUBJECT_GROUP.get(str(subject).strip())


def per_neuron_selection_key(
    group: str,
    subject: str,
    date: str,
    session: str,
    fov_name: str,
    cell_id: int,
) -> tuple[str, str, str, str, str, int]:
    return (
        str(group),
        str(subject),
        normalize_date_label(date),
        str(session),
        str(fov_name),
        int(cell_id),
    )


def per_neuron_selection_key_from_row(row: dict) -> tuple[str, str, str, str, str, int]:
    return per_neuron_selection_key(
        row["group"],
        row["subject"],
        row["date"],
        row["session"],
        row["fov"],
        int(row["rf_result_index_0based"]),
    )


def select_diverse_good_rf_rows_by_group(
    group_rows: dict[str, list[dict]],
    *,
    n_per_group: int,
    edge_margin_fraction: float,
    edge_margin_deg: float,
) -> dict[str, list[dict]]:
    selected: dict[str, list[dict]] = {}
    for group, rows in sorted(group_rows.items()):
        good_rows = [
            row
            for row in rows
            if bool(row.get("is_good", False))
            and np.isfinite(row.get("azimuth_deg", np.nan))
            and np.isfinite(row.get("elevation_deg", np.nan))
        ]
        if not good_rows:
            selected[group] = []
            continue

        interior_rows = [
            row
            for row in good_rows
            if row_is_not_too_close_to_edge(
                row,
                edge_margin_fraction=edge_margin_fraction,
                edge_margin_deg=edge_margin_deg,
            )
        ]
        edge_relaxed = len(interior_rows) < min(n_per_group, len(good_rows))
        candidates = interior_rows if interior_rows else good_rows
        picked = diverse_rf_center_sample(candidates, n_per_group)
        selected[group] = [
            annotate_selected_rf_row(row, rank=rank, edge_relaxed=edge_relaxed)
            for rank, row in enumerate(picked, start=1)
        ]
    return selected


def row_is_not_too_close_to_edge(
    row: dict,
    *,
    edge_margin_fraction: float,
    edge_margin_deg: float,
) -> bool:
    az_min = float(row["azimuth_min"])
    az_max = float(row["azimuth_max"])
    el_min = float(row["elevation_min"])
    el_max = float(row["elevation_max"])
    az_margin = max(float(edge_margin_deg), abs(az_max - az_min) * float(edge_margin_fraction))
    el_margin = max(float(edge_margin_deg), abs(el_max - el_min) * float(edge_margin_fraction))
    az = float(row["azimuth_deg"])
    el = float(row["elevation_deg"])
    return (az_min + az_margin <= az <= az_max - az_margin) and (el_min + el_margin <= el <= el_max - el_margin)


def diverse_rf_center_sample(rows: list[dict], n_select: int) -> list[dict]:
    if len(rows) <= n_select:
        return list(rows)

    az_min, az_max = float(min(row["azimuth_min"] for row in rows)), float(max(row["azimuth_max"] for row in rows))
    el_min, el_max = float(min(row["elevation_min"] for row in rows)), float(max(row["elevation_max"] for row in rows))
    az_span = max(az_max - az_min, 1e-9)
    el_span = max(el_max - el_min, 1e-9)
    points = np.asarray(
        [
            [
                (float(row["azimuth_deg"]) - az_min) / az_span,
                (float(row["elevation_deg"]) - el_min) / el_span,
            ]
            for row in rows
        ],
        dtype=float,
    )

    center = np.asarray([0.5, 0.5], dtype=float)
    peak_to_noise = np.asarray([float(row.get("peak_to_noise", np.nan)) for row in rows], dtype=float)
    peak_to_noise[~np.isfinite(peak_to_noise)] = -np.inf
    first_scores = np.linalg.norm(points - center, axis=1)
    first_idx = int(np.lexsort((-peak_to_noise, first_scores))[0])
    selected_ids = [first_idx]
    remaining = set(range(len(rows))) - {first_idx}

    while remaining and len(selected_ids) < n_select:
        chosen_points = points[selected_ids, :]
        remaining_ids = np.asarray(sorted(remaining), dtype=int)
        min_dist = np.min(
            np.linalg.norm(points[remaining_ids, None, :] - chosen_points[None, :, :], axis=2),
            axis=1,
        )
        # Primary objective: maximize distance from already selected centers.
        # Secondary objective: prefer stronger good RFs when distances tie.
        order = np.lexsort((-peak_to_noise[remaining_ids], -min_dist))
        next_idx = int(remaining_ids[order[0]])
        selected_ids.append(next_idx)
        remaining.remove(next_idx)

    return [rows[idx] for idx in selected_ids]


def annotate_selected_rf_row(row: dict, *, rank: int, edge_relaxed: bool) -> dict:
    out = dict(row)
    out["selected_rank"] = int(rank)
    out["selected_for_per_neuron_rf_plotting"] = True
    out["selection_edge_margin_relaxed"] = bool(edge_relaxed)
    return out


def plot_rf_session(
    results_root: str | Path,
    subject: str,
    date: str,
    session: str | int,
    plots_root: str | Path,
    *,
    group: str = "",
    fovs: list[str] | tuple[str, ...] | None = None,
    flip_fovs: set[str] | None = None,
    min_ev: float = MIN_EV,
    max_p_value: float = MAX_P_VALUE,
    min_peak_to_noise: float = MIN_PEAK_TO_NOISE,
    max_neuron_plots: int | None = MAX_NEURON_RF_PLOTS,
    selected_per_neuron_keys: set[tuple[str, str, str, str, str, int]] | None = None,
    selected_good_rf_collection_dir: Path | None = None,
    skip_fov_overlays: bool = False,
    make_combined_session_plots: bool = MAKE_COMBINED_SESSION_PLOTS,
    skip_existing: bool = True,
) -> dict[str, Any]:
    """Plot per-FOV and combined-session RF outputs for one subject/date/session."""

    subject = str(subject)
    date = normalize_date_label(date)
    session = str(session)
    results_root = Path(results_root)
    plots_root = Path(plots_root)

    results_session = results_root / subject / date / session / "alf"
    plots_session = plots_root / subject / date / session / "alf"
    flip_fovs = set(flip_fovs or ())

    if fovs:
        fov_names = list(fovs)
    elif FOV_NAMES is not None:
        fov_names = FOV_NAMES
    else:
        fov_names = discover_fovs(results_session, results_session)

    if not fov_names:
        raise FileNotFoundError(f"No FOV folders found under results={results_session}")

    print("\n========== FANCi RF plotting ==========")
    print(f"group/subject/date/session: {group} / {subject} / {date} / {session}")
    print(f"Results session:      {results_session}")
    print(f"Plots session:        {plots_session}")
    print(f"FOVs:                 {', '.join(fov_names)}")
    if flip_fovs:
        print(f"Flip azimuth FOVs:    {', '.join(sorted(flip_fovs))}")

    if skip_existing and session_plot_is_current(plots_session, fov_names, flip_fovs):
        print(f"[SKIP] Existing plots are current for {group}/{subject}/{date}/{session}: {plots_session}")
        combined_records = load_session_records_for_aggregates(results_session, fov_names, flip_fovs)
        aggregate_rows = session_aggregate_rows(
            combined_records,
            group=group,
            subject=subject,
            date=date,
            session=session,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
        )
        aggregate_ev_rows = session_ev_aggregate_rows(
            combined_records,
            group=group,
            subject=subject,
            date=date,
            session=session,
            max_p_value=max_p_value,
        )
        return {
            "group": group,
            "subject": subject,
            "date": date,
            "session": session,
            "plots_session": str(plots_session),
            "skipped_existing": True,
            "aggregate_rows": aggregate_rows,
            "aggregate_ev_rows": aggregate_ev_rows,
        }

    raw_suite2p = None
    if not skip_fov_overlays and (
        MAKE_FOV_OVERLAYS or MAKE_NEUTRAL_ROI_OVERLAYS or make_combined_session_plots or MAKE_COMBINED_FOV_OVERLAYS
    ):
        try:
            raw_suite2p = find_raw_suite2p_folder(subject=subject, date=date, session=session)
            print(f"Raw Suite2p folder:   {raw_suite2p}")
        except FileNotFoundError as exc:
            raw_suite2p = None
            print(f"[WARN] FOV overlay plots will be skipped: {exc}")

    summary_rows = []
    combined_records = []
    for fov_name in fov_names:
        print(f"\n=== Plotting {fov_name} ===")
        fov_results_dir = results_session / fov_name
        fov_alf_dir = results_session / fov_name
        fov_output_dir = plots_session / fov_name
        flip_azimuth = fov_name in flip_fovs
        fov_output_dir.mkdir(parents=True, exist_ok=True)

        result = load_fanci_result(fov_results_dir)
        combined_records.append(
            {
                "fov_name": fov_name,
                "result": result,
                "fov_alf_dir": fov_alf_dir,
                "fov_output_dir": fov_output_dir,
                "flip_azimuth": flip_azimuth,
            }
        )
        counts = plot_one_fov(
            result,
            fov_name=fov_name,
            fov_alf_dir=fov_alf_dir,
            output_dir=fov_output_dir,
            raw_suite2p_folder=raw_suite2p,
            group=group,
            subject=subject,
            date=date,
            session=session,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
            max_neuron_plots=max_neuron_plots,
            selected_per_neuron_keys=selected_per_neuron_keys,
            selected_good_rf_collection_dir=selected_good_rf_collection_dir,
            skip_fov_overlays=skip_fov_overlays,
            flip_azimuth=flip_azimuth,
        )
        summary_rows.append({"fov": fov_name, **counts, "output_dir": str(fov_output_dir)})
        print(f"[DONE] {fov_name}: {counts}")

    combined_counts = {}
    if make_combined_session_plots:
        try:
            combined_counts = plot_combined_session_outputs(
                combined_records,
                output_dir=plots_session / "combined",
                raw_suite2p_folder=raw_suite2p,
                subject=subject,
                date=date,
                session=session,
                group=group,
                min_ev=min_ev,
                max_p_value=max_p_value,
                min_peak_to_noise=min_peak_to_noise,
                max_neuron_plots=max_neuron_plots,
                selected_per_neuron_keys=selected_per_neuron_keys,
                skip_fov_overlays=skip_fov_overlays,
            )
            print(f"[DONE] combined: {combined_counts}")
        except Exception as exc:
            combined_counts = {"error": str(exc)}
            print(f"[WARN] Combined session plots failed: {exc}")
    elif MAKE_COMBINED_FOV_OVERLAYS:
        combined_counts = {"combined_rf_roi_overlay": 0, "errors": []}
        if skip_fov_overlays or raw_suite2p is None:
            message = "combined FOV overlay plots skipped; raw Suite2p folder unavailable or --skip-fov-overlays was set."
            combined_counts["errors"].append(message)
            print(f"[WARN] {message}")
        else:
            try:
                combined_counts["combined_rf_roi_overlay"] = plot_combined_rf_roi_overlay(
                    combined_records,
                    raw_suite2p_folder=raw_suite2p,
                    output_dir=plots_session / "combined" / "FOV_overlays",
                    subject=subject,
                    date=date,
                    session=session,
                    min_ev=min_ev,
                    max_p_value=max_p_value,
                    min_peak_to_noise=min_peak_to_noise,
                )
                print(f"[DONE] combined FOV overlays: {combined_counts}")
            except Exception as exc:
                message = f"combined RF/FOV overlay failed: {exc}"
                combined_counts["errors"].append(message)
                print(f"[WARN] {message}")

    aggregate_rows = session_aggregate_rows(
        combined_records,
        group=group,
        subject=subject,
        date=date,
        session=session,
        min_ev=min_ev,
        max_p_value=max_p_value,
        min_peak_to_noise=min_peak_to_noise,
    )
    aggregate_ev_rows = session_ev_aggregate_rows(
        combined_records,
        group=group,
        subject=subject,
        date=date,
        session=session,
        max_p_value=max_p_value,
    )

    summary_path = plots_session / "plot_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "group": group,
        "subject": subject,
        "date": date,
        "session": session,
        "results_session": str(results_session),
        "plots_session": str(plots_session),
        "raw_suite2p_folder": str(raw_suite2p) if raw_suite2p is not None else None,
        "flipped_azimuth_fovs": sorted(flip_fovs),
        "fovs": summary_rows,
        "combined": combined_counts,
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nSummary written to: {summary_path}")
    summary["aggregate_rows"] = aggregate_rows
    summary["aggregate_ev_rows"] = aggregate_ev_rows
    return summary


def session_plot_is_current(plots_session: Path, fov_names: list[str], flip_fovs: set[str]) -> bool:
    summary_path = plots_session / "plot_summary.json"
    if not summary_path.exists():
        return False
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"[WARN] Existing plot summary could not be read, replotting: {summary_path} ({exc})")
        return False

    expected_fovs = sorted(normalize_fov_name(name) for name in fov_names)
    actual_fovs = sorted(normalize_fov_name(row.get("fov", "")) for row in summary.get("fovs", []))
    expected_flips = sorted(normalize_fov_name(name) for name in flip_fovs)
    actual_flips = sorted(normalize_fov_name(name) for name in summary.get("flipped_azimuth_fovs", []))
    if actual_fovs != expected_fovs:
        return False
    if actual_flips != expected_flips:
        return False

    for fov_name in expected_fovs:
        if not (plots_session / fov_name / "plot_input_summary.json").exists():
            return False
    return True


def load_session_records_for_aggregates(
    results_session: Path,
    fov_names: list[str],
    flip_fovs: set[str],
) -> list[dict]:
    records: list[dict] = []
    for fov_name in fov_names:
        fov_results_dir = results_session / fov_name
        if not fov_results_dir.is_dir():
            continue
        records.append(
            {
                "fov_name": fov_name,
                "result": load_fanci_result(fov_results_dir),
                "fov_alf_dir": fov_results_dir,
                "fov_output_dir": None,
                "flip_azimuth": fov_name in flip_fovs,
            }
        )
    return records


def plot_one_fov(
    result: FanciResult,
    *,
    fov_name: str,
    fov_alf_dir: Path,
    output_dir: Path,
    raw_suite2p_folder: Path | None,
    group: str,
    subject: str,
    date: str,
    session: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    max_neuron_plots: int | None,
    selected_per_neuron_keys: set[tuple[str, str, str, str, str, int]] | None,
    selected_good_rf_collection_dir: Path | None,
    skip_fov_overlays: bool,
    flip_azimuth: bool,
) -> dict[str, int]:
    counts = {
        "gaussian_coverage": 0,
        "per_neuron_rf_maps": 0,
        "neutral_roi_overlays": 0,
        "rf_roi_overlays": 0,
        "shift_ev_histograms": 0,
        "predicted_vs_actual": 0,
        "peak_time_courses": 0,
        "rf_size_histograms": 0,
    }

    if MAKE_GAUSSIAN_COVERAGE:
        counts["gaussian_coverage"] = plot_gaussian_coverage(
            result,
            output_dir / "RF_gaussian_coverage",
            subject=subject,
            date=date,
            session=session,
            fov_name=fov_name,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
            flip_azimuth=flip_azimuth,
        )

    if MAKE_PER_NEURON_RF_MAPS:
        counts["per_neuron_rf_maps"] = plot_per_neuron_metadata(
            result,
            output_dir / "per_neuron_metadata",
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
            max_neuron_plots=max_neuron_plots,
            selected_per_neuron_keys=selected_per_neuron_keys,
            group=group,
            subject=subject,
            date=date,
            session=session,
            fov_name=fov_name,
            flip_azimuth=flip_azimuth,
            selected_good_rf_collection_dir=selected_good_rf_collection_dir,
        )

    if MAKE_SHIFT_EV_HISTOGRAM:
        counts["shift_ev_histograms"] = plot_shift_ev_histogram(
            result,
            output_dir / "shift_ev_diagnostics",
            subject=subject,
            date=date,
            session=session,
            fov_name=fov_name,
            max_p_value=max_p_value,
        )

    if MAKE_PREDICTED_VS_ACTUAL:
        counts["predicted_vs_actual"] = plot_predicted_vs_actual_traces(
            result,
            output_dir / "predicted_vs_actual_traces",
            subject=subject,
            date=date,
            session=session,
            fov_name=fov_name,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
        )

    if MAKE_PEAK_TIME_COURSES:
        counts["peak_time_courses"] = plot_peak_time_courses(
            result,
            output_dir / "peak_time_courses",
            subject=subject,
            date=date,
            session=session,
            fov_name=fov_name,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
        )

    if MAKE_RF_SIZE_HISTOGRAM:
        counts["rf_size_histograms"] = plot_rf_size_histogram(
            result,
            output_dir / "RF_size_histogram",
            subject=subject,
            date=date,
            session=session,
            fov_name=fov_name,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
        )

    if not skip_fov_overlays and raw_suite2p_folder is not None:
        if MAKE_NEUTRAL_ROI_OVERLAYS:
            counts["neutral_roi_overlays"] = plot_neutral_roi_overlays(
                raw_suite2p_folder,
                output_dir / "FOV_overlays",
                subject=subject,
                date=date,
                session=session,
                fov_name=fov_name,
            )
        if MAKE_FOV_OVERLAYS:
            counts["rf_roi_overlays"] = plot_rf_roi_overlay(
                result,
                raw_suite2p_folder,
                output_dir / "FOV_overlays",
                subject=subject,
                date=date,
                session=session,
                fov_name=fov_name,
                min_ev=min_ev,
                max_p_value=max_p_value,
                min_peak_to_noise=min_peak_to_noise,
                flip_azimuth=flip_azimuth,
            )

    write_fov_summary(result, output_dir / "plot_input_summary.json")
    return counts


# =============================================================================
# Result loading
# =============================================================================


class FanciResult:
    def __init__(self, results_dir: Path, arrays: dict[str, np.ndarray], mat_path: Path | None):
        self.results_dir = results_dir
        self.arrays = arrays
        self.mat_path = mat_path

    def has(self, key: str) -> bool:
        return key in self.arrays

    def get(self, key: str, default: Any = None) -> Any:
        return self.arrays.get(key, default)

    def require(self, key: str) -> np.ndarray:
        if key not in self.arrays:
            raise KeyError(f"Missing required result array: {key} in {self.results_dir}")
        return self.arrays[key]


FANCI_NPY_FILE_MAP = {
    "maps": "_FANCi_rf.maps.npy",
    "explVars": "_FANCi_rf.explainedVariance.npy",
    "explVarsFullModel": "_FANCi_rf.explainedVarianceFullModel.npy",
    "lambdas": "_FANCi_rf.lambdas.npy",
    "runLambdas": "_FANCi_rf.runLambdas.npy",
    "pValues": "_FANCi_rf.pValues.npy",
    "timestamps": "_FANCi_rf.temporalLags.npy",
    "gaussPars": "_FANCi_rf.gaussPars.npy",
    "peakToNoise": "_FANCi_rf.peakToNoise.npy",
    "runningLags": "_FANCi_rf.runningLags.npy",
    "runningKernels": "_FANCi_rf.runningKernels.npy",
    "roiIndex0": "_FANCi_rf.roiIndex.npy",
    "roiIndexMatlab": "_FANCi_rf.roiIndexMatlab.npy",
    "cellClassifier": "_FANCi_rf.cellClassifier.npy",
    "predictionsFullModel": "_FANCi_rf.predictionsFullModel.npy",
    "predictionsRFOnly": "_FANCi_rf.predictionsRFOnly.npy",
    "tracesPreprocessed": "_FANCi_rf.tracesPreprocessed.npy",
    "tracesPreprocessed_times": "_FANCi_rf.tracesPreprocessed_times.npy",
    "tracesPreprocessed_roiIndex0": "_FANCi_rf.tracesPreprocessed_roiIndex.npy",
    "tracesPreprocessed_roiIndexMatlab": "_FANCi_rf.tracesPreprocessed_roiIndexMatlab.npy",
}

# MAT fields used only as a fallback when the matching _FANCi_rf.*.npy files are absent.
MAT_FIELD_MAP = {
    "edges": "results/edges",
    "bestSubfields": "results/bestSubfields",
    "subfieldSigns": "results/subfieldSigns",
    "optimalDelays": "results/optimalDelays",
    "evRFOnly": "results/evRFOnly",
    "evShiftRFOnly": "results/evShiftRFOnly",
    "predictionsFullModel": "results/predictionsFullModel",
    "predictionsRFOnly": "results/predictionsRFOnly",
    "runningTimestamps": "results/running/timestamps",
    "minPeak_shift": "results/minPeak_shift",
}


def load_fanci_result(results_dir: Path) -> FanciResult:
    if not results_dir.is_dir():
        raise FileNotFoundError(f"Missing FOV result directory: {results_dir}")

    arrays: dict[str, np.ndarray] = {}
    for key, filename in FANCI_NPY_FILE_MAP.items():
        path = results_dir / filename
        if path.exists():
            arrays[key] = np.load(path, allow_pickle=True)

    # Prefer .npy outputs.  Only use MAT files to backfill arrays that were not
    # saved as .npy, and check both running-regressor modes.
    mat_candidates = [
        results_dir / "_FANCi_rf.withRunningRegressor.mat",
        results_dir / "_FANCi_rf.withoutRunningRegressor.mat",
    ]
    mat_path = next((path for path in mat_candidates if path.exists()), None)
    if h5py is not None:
        for candidate in mat_candidates:
            if not candidate.exists():
                continue
            with h5py.File(candidate, "r") as handle:
                for key, h5_path in MAT_FIELD_MAP.items():
                    if key not in arrays and h5_path in handle:
                        arrays[key] = read_h5_numeric(handle[h5_path])
    else:
        for candidate in mat_candidates:
            if candidate.exists():
                print(f"[WARN] h5py is not available; MAT-only arrays will not be loaded from {candidate}")

    # Backstop defaults.
    maps = arrays.get("maps")
    if maps is not None:
        n_cells = maps.shape[0]
        if "edges" not in arrays:
            arrays["edges"] = np.asarray([-135.0, 135.0, 40.0, -40.0], dtype=float)
        if "timestamps" not in arrays:
            arrays["timestamps"] = np.arange(maps.shape[3], dtype=float)
        if "roiIndex0" not in arrays:
            arrays["roiIndex0"] = np.arange(n_cells, dtype=int)
        if "roiIndexMatlab" not in arrays:
            arrays["roiIndexMatlab"] = arrays["roiIndex0"].astype(int) + 1

    for required in ("maps", "explVars", "pValues", "gaussPars", "peakToNoise"):
        if required not in arrays:
            raise FileNotFoundError(
                f"Missing required {required} file in {results_dir}. "
                "Expected latest _FANCi_rf.* outputs."
            )

    orient_loaded_arrays(arrays)
    return FanciResult(results_dir, arrays, mat_path if mat_path is not None and mat_path.exists() else None)


def read_h5_numeric(dataset) -> np.ndarray:
    arr = np.asarray(dataset)
    # MATLAB v7.3 HDF5 stores MATLAB column-major arrays in reversed dimension order.
    # Squeeze vectors, but keep matrices/tensors for later orientation against known shapes.
    arr = np.array(arr)
    if arr.dtype.kind in {"S", "U", "O"}:
        return arr
    return np.squeeze(arr)


def orient_loaded_arrays(arrays: dict[str, np.ndarray]) -> None:
    maps = np.asarray(arrays["maps"])
    if maps.ndim != 5:
        raise ValueError(f"_FANCi_rf.maps.npy must be 5D [ROI,row,col,lag,ON/OFF], got {maps.shape}")
    n_cells = maps.shape[0]
    n_lags = maps.shape[3]

    for key in (
        "explVars",
        "explVarsFullModel",
        "lambdas",
        "runLambdas",
        "pValues",
        "peakToNoise",
        "bestSubfields",
        "optimalDelays",
        "evRFOnly",
        "roiIndex0",
        "roiIndexMatlab",
        "cellClassifier",
        "tracesPreprocessed_roiIndex0",
        "tracesPreprocessed_roiIndexMatlab",
    ):
        if key in arrays:
            arrays[key] = np.ravel(np.asarray(arrays[key]))

    for key in ("timestamps", "edges", "runningLags", "runningTimestamps", "minPeak_shift"):
        if key in arrays:
            arrays[key] = np.ravel(np.asarray(arrays[key], dtype=float))

    for key in ("gaussPars",):
        if key in arrays:
            arrays[key] = orient_matrix(np.asarray(arrays[key], dtype=float), n_rows=n_cells, n_cols=None)

    if "subfieldSigns" in arrays:
        arrays["subfieldSigns"] = orient_matrix(np.asarray(arrays["subfieldSigns"], dtype=float), n_rows=n_cells, n_cols=2)

    for key in ("tracesPreprocessed", "predictionsFullModel", "predictionsRFOnly"):
        if key in arrays:
            arrays[key] = orient_time_by_cell(np.asarray(arrays[key], dtype=float), n_cells=n_cells)

    if "evShiftRFOnly" in arrays:
        arr = np.asarray(arrays["evShiftRFOnly"], dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(n_cells, -1) if arr.size % n_cells == 0 else arr.reshape(-1, 1)
        elif arr.shape[0] != n_cells and arr.shape[1] == n_cells:
            arr = arr.T
        arrays["evShiftRFOnly"] = arr

    # If MAT vectors came in as [1 x n] or [n x 1], ravel fixed them. Trim/pad happens later.
    if arrays["timestamps"].size != n_lags:
        arrays["timestamps"] = np.ravel(arrays["timestamps"])[:n_lags]


def orient_matrix(arr: np.ndarray, *, n_rows: int, n_cols: int | None) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 1:
        if n_cols is None:
            return arr.reshape(n_rows, -1)
        return arr.reshape(n_rows, n_cols)
    if arr.shape[0] == n_rows and (n_cols is None or arr.shape[1] == n_cols):
        return arr
    if arr.shape[-1] == n_rows:
        arr = arr.T
    if n_cols is not None and arr.shape[1] != n_cols and arr.shape[0] == n_cols:
        arr = arr.T
    return arr


def orient_time_by_cell(arr: np.ndarray, *, n_cells: int) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 1:
        if arr.size == n_cells:
            return arr.reshape(1, n_cells)
        return arr.reshape(-1, 1)
    if arr.shape[1] == n_cells:
        return arr
    if arr.shape[0] == n_cells:
        return arr.T
    return arr


# =============================================================================
# Shared result helpers
# =============================================================================


def result_column(result: FanciResult, key: str, n_expected: int) -> np.ndarray:
    values = np.full(n_expected, np.nan, dtype=float)
    if not result.has(key):
        return values
    incoming = np.ravel(np.asarray(result.get(key), dtype=float))
    n = min(n_expected, incoming.size)
    values[:n] = incoming[:n]
    return values


def good_rf_mask(
    result: FanciResult,
    *,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
) -> np.ndarray:
    maps = result.require("maps")
    n = maps.shape[0]
    ev = result_column(result, "explVars", n)
    p_values = result_column(result, "pValues", n)
    peak_to_noise = result_column(result, "peakToNoise", n)
    gauss = result.require("gaussPars")
    has_gaussian = np.all(np.isfinite(gauss[:, :6]), axis=1)
    return (
        np.isfinite(ev)
        & np.isfinite(p_values)
        & np.isfinite(peak_to_noise)
        & has_gaussian
        & (ev > min_ev)
        & (p_values < max_p_value)
        & (peak_to_noise > min_peak_to_noise)
    )


def significant_rf_mask(result: FanciResult, *, min_ev: float, max_p_value: float) -> np.ndarray:
    maps = result.require("maps")
    n = maps.shape[0]
    ev = result_column(result, "explVars", n)
    p_values = result_column(result, "pValues", n)
    return np.isfinite(ev) & np.isfinite(p_values) & (ev > min_ev) & (p_values < max_p_value)


def top_good_by_peak_noise(
    result: FanciResult,
    *,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    n_top: int = 5,
) -> np.ndarray:
    mask = good_rf_mask(
        result,
        min_ev=min_ev,
        max_p_value=max_p_value,
        min_peak_to_noise=min_peak_to_noise,
    )
    peak_to_noise = result_column(result, "peakToNoise", result.require("maps").shape[0])
    ids = np.flatnonzero(mask)
    if ids.size == 0:
        return ids
    order = np.argsort(peak_to_noise[ids])[::-1]
    return ids[order[:n_top]]


def ev_for_shift_plots(result: FanciResult) -> tuple[np.ndarray, str]:
    n = result.require("maps").shape[0]
    if result.has("evRFOnly"):
        return result_column(result, "evRFOnly", n), "RF-only EV from final circular-shift test"
    return result_column(result, "explVars", n), "CV RF-only explained variance (fallback; evRFOnly not found)"


def gaussian_display_pars(gauss: np.ndarray, *, flip_azimuth: bool) -> np.ndarray:
    pars = np.asarray(gauss, dtype=float).copy()
    if flip_azimuth and pars.ndim == 2:
        if pars.shape[1] > 1:
            pars[:, 1] = -pars[:, 1]
        if pars.shape[1] > 5:
            pars[:, 5] = np.pi - pars[:, 5]
    return pars


def style_figure_for_export(fig) -> None:
    for ax in fig.axes:
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(top=False, right=False)
        if hasattr(ax, "outline"):
            try:
                ax.outline.set_visible(False)
            except Exception:
                pass


def save_figure(fig, base_path: Path, *, dpi: int = 220) -> list[Path]:
    style_figure_for_export(fig)
    base_path.parent.mkdir(parents=True, exist_ok=True)
    outputs = save_figure_outputs(fig, base_path, dpi=dpi)
    plt.close(fig)
    return outputs


def save_figure_outputs(fig, base_path: Path, *, dpi: int = 220) -> list[Path]:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    png = base_path.with_suffix(".png")
    svg = base_path.with_suffix(".svg")
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    return [png, svg]


def add_zero_axis_lines(ax, *, x_zero: bool = True, y_zero: bool = True) -> None:
    """Draw grey x=0 and y=0 data-axis lines for EV/RF-size panels."""

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    if x_zero:
        x0, x1 = min(float(x0), 0.0), max(float(x1), 0.0)
    if y_zero:
        y0, y1 = min(float(y0), 0.0), max(float(y1), 0.0)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    if y_zero:
        ax.axhline(0, color="0.45", linewidth=0.9, alpha=0.85, zorder=0)
    if x_zero:
        ax.axvline(0, color="0.45", linewidth=0.9, alpha=0.85, zorder=0)


def make_two_tick_square_axes(ax, *, x_zero: bool = True, y_zero: bool = True) -> None:
    """Use only two ticks per axis and make EV/RF-size panels square."""

    add_zero_axis_lines(ax, x_zero=x_zero, y_zero=y_zero)
    ax.xaxis.set_major_locator(LinearLocator(2))
    ax.yaxis.set_major_locator(LinearLocator(2))
    try:
        ax.set_box_aspect(1)
    except Exception:
        pass
    ax.grid(False)


# =============================================================================
# Plot: RF size histogram
# =============================================================================


def plot_rf_size_histogram(
    result: FanciResult,
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    fov_name: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
) -> int:
    """Plot RF-size distribution for good RFs.

    RF size is defined from the fitted 2D Gaussian as:
        (abs(sigma_x) + abs(sigma_y)) / 2
    where gaussPars columns are [amp, x0, sigma_x, y0, sigma_y, theta, ...].
    """

    gauss = np.asarray(result.require("gaussPars"), dtype=float)
    n_cells = result.require("maps").shape[0]
    n = min(n_cells, gauss.shape[0])
    gauss = gauss[:n, :]

    has_valid_width = (
        (gauss.shape[1] >= 5)
        & np.isfinite(gauss[:, 2])
        & np.isfinite(gauss[:, 4])
        & (np.abs(gauss[:, 2]) > 0)
        & (np.abs(gauss[:, 4]) > 0)
    )
    good = good_rf_mask(
        result,
        min_ev=min_ev,
        max_p_value=max_p_value,
        min_peak_to_noise=min_peak_to_noise,
    )[:n]

    use = good & has_valid_width
    rf_size = (np.abs(gauss[use, 2]) + np.abs(gauss[use, 4])) / 2.0
    rf_indices = np.flatnonzero(use)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Save the values used for the histogram so the plot can be audited.
    csv_path = output_dir / f"{fov_name}_goodRF_size_values.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.DictWriter(
            file_obj,
            fieldnames=[
                "rf_result_index_0based",
                "rf_result_index_1based",
                "source_roi0",
                "sigma_x_deg",
                "sigma_y_deg",
                "rf_size_deg",
            ],
        )
        writer.writeheader()
        roi_index0 = np.ravel(np.asarray(result.get("roiIndex0", np.arange(n)), dtype=float))
        for rf_idx, size_deg in zip(rf_indices, rf_size):
            source_roi0 = roi_index0[rf_idx] if rf_idx < roi_index0.size else np.nan
            writer.writerow(
                {
                    "rf_result_index_0based": int(rf_idx),
                    "rf_result_index_1based": int(rf_idx + 1),
                    "source_roi0": int(source_roi0) if np.isfinite(source_roi0) else "",
                    "sigma_x_deg": float(gauss[rf_idx, 2]),
                    "sigma_y_deg": float(gauss[rf_idx, 4]),
                    "rf_size_deg": float(size_deg),
                }
            )

    fig, ax = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
    if rf_size.size > 0:
        max_size = float(np.nanmax(rf_size))
        if np.isfinite(max_size) and max_size > 0:
            bin_width = 2.0
            bins = np.arange(0.0, math.ceil(max_size / bin_width) * bin_width + bin_width, bin_width)
            if bins.size < 2:
                bins = 10
        else:
            bins = 10
        ax.hist(rf_size, bins=bins, color="0.75", edgecolor="0.45", linewidth=0.8)
        mean_size = float(np.nanmean(rf_size))
        ax.axvline(mean_size, linestyle="-", linewidth=2.2, color="black", label="mean")
        stats_text = (
            f"n = {rf_size.size}\n"
            f"mean = {mean_size:.2f} deg\n"
            f"median = {np.nanmedian(rf_size):.2f} deg"
        )
    else:
        ax.text(
            0.5,
            0.5,
            "No good RF neurons with valid Gaussian widths",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )
        stats_text = "n = 0"

    ax.text(
        0.98,
        0.95,
        stats_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )
    ax.set_xlabel("RF size (deg) = (|sigma_x| + |sigma_y|) / 2")
    ax.set_ylabel("Number of neurons")
    ax.set_title(
        f"{subject} {date} session {session} {fov_name}: good-RF size histogram\n"
        f"Good RF: EV > {min_ev:g}, p < {max_p_value:g}, P/N > {min_peak_to_noise:g}"
    )
    make_two_tick_square_axes(ax, x_zero=True, y_zero=True)

    save_figure(fig, output_dir / f"{fov_name}_goodRF_size_histogram", dpi=250)
    return 1


# =============================================================================
# Plot 1: Gaussian visual-space coverage
# =============================================================================


def plot_gaussian_coverage(
    result: FanciResult,
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    fov_name: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    flip_azimuth: bool = False,
    n_azimuth_grid: int = 500,
    n_elevation_grid: int = 250,
    gaussian_sigma_radius: float = 2.0,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    gauss_pars = gaussian_display_pars(result.require("gaussPars"), flip_azimuth=flip_azimuth)
    n_neurons = gauss_pars.shape[0]
    ev = result_column(result, "explVars", n_neurons)
    p_values = result_column(result, "pValues", n_neurons)
    peak_to_noise = result_column(result, "peakToNoise", n_neurons)

    has_gaussian = (
        np.all(np.isfinite(gauss_pars[:, :6]), axis=1)
        & (gauss_pars[:, 2] > 0)
        & (gauss_pars[:, 4] > 0)
    )
    significant = has_gaussian & np.isfinite(ev) & np.isfinite(p_values) & (ev > min_ev) & (p_values < max_p_value)
    peak_pass = significant & np.isfinite(peak_to_noise) & (peak_to_noise > min_peak_to_noise)
    non_peak_significant = significant & ~peak_pass
    all_significant = peak_pass | non_peak_significant

    edges = get_edges(result)
    azimuth_limits = tuple(np.sort(np.asarray(edges[:2], dtype=float)))
    elevation_limits = tuple(np.sort(np.asarray(edges[2:4], dtype=float)))

    azimuth_axis = np.linspace(azimuth_limits[0], azimuth_limits[1], n_azimuth_grid)
    elevation_axis = np.linspace(elevation_limits[0], elevation_limits[1], n_elevation_grid)
    azimuth_grid, elevation_grid = np.meshgrid(azimuth_axis, elevation_axis)

    peak_ids = np.flatnonzero(peak_pass)
    non_peak_ids = np.flatnonzero(non_peak_significant)
    all_ids = np.flatnonzero(all_significant)

    peak_coverage = gaussian_coverage_map(gauss_pars, peak_ids, azimuth_grid, elevation_grid, gaussian_sigma_radius)
    all_coverage = gaussian_coverage_map(gauss_pars, all_ids, azimuth_grid, elevation_grid, gaussian_sigma_radius)

    shared_max = np.nanmax([np.nanmax(peak_coverage), np.nanmax(all_coverage)])
    if not np.isfinite(shared_max) or shared_max <= 0:
        shared_max = 1.0

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), constrained_layout=True)
    plot_coverage_panel(
        axes[0],
        peak_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        shared_max,
        f"Peak/noise pass only (n = {peak_ids.size})",
        gauss_pars,
        peak_ids,
        np.empty(0, dtype=int),
        gaussian_sigma_radius,
    )
    plot_coverage_panel(
        axes[1],
        all_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        shared_max,
        f"All significant RFs (n = {all_ids.size}; P/N n = {peak_ids.size}, non-P/N n = {non_peak_ids.size})",
        gauss_pars,
        peak_ids,
        non_peak_ids,
        gaussian_sigma_radius,
    )
    fig.suptitle(
        f"{subject} {date} session {session} {fov_name}: Gaussian RF coverage | "
        f"EV > {min_ev:.3f}, p < {max_p_value:.3f}, P/N > {min_peak_to_noise:.2f}",
        fontweight="bold",
    )

    save_figure(fig, output_dir / f"{subject}_{date}_session{session}_{fov_name}_RF_gaussian_coverage_heatmap", dpi=300)
    return 1


def plot_coverage_panel(
    ax,
    coverage: np.ndarray,
    azimuth_axis: np.ndarray,
    elevation_axis: np.ndarray,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    shared_max: float,
    title: str,
    gauss_pars: np.ndarray,
    peak_ids: np.ndarray,
    non_peak_ids: np.ndarray,
    sigma_radius: float,
    cax=None,
) -> None:
    im = ax.imshow(
        coverage,
        origin="lower",
        extent=(azimuth_axis[0], azimuth_axis[-1], elevation_axis[0], elevation_axis[-1]),
        aspect="equal",
        cmap="viridis",
        vmin=0,
        vmax=shared_max,
    )
    ax.axvline(0, color="#333333", linewidth=0.75)
    ax.axhline(0, color="#333333", linewidth=0.75)
    ax.set_xlim(azimuth_limits)
    ax.set_ylim(elevation_limits)
    ax.set_xlabel("Azimuth (deg)")
    ax.set_ylabel("Elevation (deg)")
    ax.set_title(title)
    if cax is None:
        cb = plt.colorbar(im, ax=ax)
    else:
        cb = plt.colorbar(im, cax=cax)
    cb.set_label("Summed normalized Gaussian RF density")

    if non_peak_ids.size:
        plot_gaussian_ellipses(
            ax,
            gauss_pars,
            non_peak_ids,
            sigma_radius,
            color="#dddddd",
            linestyle="--",
            linewidth=0.8,
            azimuth_limits=azimuth_limits,
            elevation_limits=elevation_limits,
        )
    if peak_ids.size:
        plot_gaussian_ellipses(
            ax,
            gauss_pars,
            peak_ids,
            sigma_radius,
            color="#ffffff",
            linestyle="-",
            linewidth=1.1,
            azimuth_limits=azimuth_limits,
            elevation_limits=elevation_limits,
        )


def gaussian_coverage_map(
    gauss_pars: np.ndarray,
    cell_ids: np.ndarray,
    azimuth_grid: np.ndarray,
    elevation_grid: np.ndarray,
    sigma_radius: float,
) -> np.ndarray:
    coverage = np.zeros_like(azimuth_grid, dtype=float)
    for cell_id in cell_ids:
        pars = gauss_pars[int(cell_id), :]
        if pars.size < 6 or not np.all(np.isfinite(pars[:6])):
            continue
        x0, sx, y0, sy, theta = pars[1], abs(pars[2]), pars[3], abs(pars[4]), pars[5]
        if sx <= 0 or sy <= 0:
            continue
        dx = azimuth_grid - x0
        dy = elevation_grid - y0
        x_rot = dx * np.cos(theta) + dy * np.sin(theta)
        y_rot = -dx * np.sin(theta) + dy * np.cos(theta)
        radius_squared = (x_rot / sx) ** 2 + (y_rot / sy) ** 2
        gaussian = np.exp(-0.5 * radius_squared)
        if np.isfinite(sigma_radius) and sigma_radius > 0:
            gaussian[radius_squared > sigma_radius**2] = 0
        coverage += gaussian
    return coverage


def plot_gaussian_ellipses(
    ax,
    gauss_pars: np.ndarray,
    cell_ids: np.ndarray,
    sigma_radius: float,
    *,
    color: str,
    linestyle: str,
    linewidth: float,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
) -> None:
    theta_axis = np.linspace(-np.pi, np.pi, 100)
    for cell_id in cell_ids:
        pars = gauss_pars[int(cell_id), :]
        if pars.size < 6 or not np.all(np.isfinite(pars[:6])):
            continue
        x0, sx, y0, sy, theta = pars[1], abs(pars[2]), pars[3], abs(pars[4]), pars[5]
        if sx <= 0 or sy <= 0:
            continue
        x = sx * np.cos(theta_axis) * sigma_radius
        y = sy * np.sin(theta_axis) * sigma_radius
        x_rot = x0 + x * np.cos(theta) - y * np.sin(theta)
        y_rot = y0 + x * np.sin(theta) + y * np.cos(theta)
        outside = (
            (x_rot < azimuth_limits[0])
            | (x_rot > azimuth_limits[1])
            | (y_rot < elevation_limits[0])
            | (y_rot > elevation_limits[1])
        )
        x_rot[outside] = np.nan
        y_rot[outside] = np.nan
        ax.plot(x_rot, y_rot, color=color, linestyle=linestyle, linewidth=linewidth)


# =============================================================================
# Plot 2: Per-neuron RF map metadata panels
# =============================================================================


def plot_per_neuron_metadata(
    result: FanciResult,
    output_dir: Path,
    *,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    max_neuron_plots: int | None,
    selected_per_neuron_keys: set[tuple[str, str, str, str, str, int]] | None,
    group: str,
    subject: str,
    date: str,
    session: str,
    fov_name: str,
    flip_azimuth: bool,
    selected_good_rf_collection_dir: Path | None,
) -> int:
    """Plot per-neuron RF maps, split into goodRF PASS/FAIL folders.

    The per-neuron plot set is still restricted to neurons that pass the
    RF-inclusion test (RF EV > min_ev and p < max_p_value).  Within that
    significant set, plots are separated by the peak-to-noise criterion:

      * PASS_goodRFs_peakNoise: significant RFs with peak/noise > threshold
      * FAIL_significantRF_peakNoise: significant RFs failing peak/noise
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    pass_dir = output_dir / "PASS_goodRFs_peakNoise"
    fail_dir = output_dir / "FAIL_significantRF_peakNoise"
    pass_dir.mkdir(parents=True, exist_ok=True)
    fail_dir.mkdir(parents=True, exist_ok=True)

    maps = result.require("maps")
    n_neurons = maps.shape[0]
    ev = result_column(result, "explVars", n_neurons)
    p_values = result_column(result, "pValues", n_neurons)
    full_ev = result_column(result, "explVarsFullModel", n_neurons)
    peak_to_noise = result_column(result, "peakToNoise", n_neurons)

    significant = np.flatnonzero(
        np.isfinite(ev)
        & np.isfinite(p_values)
        & (ev > min_ev)
        & (p_values < max_p_value)
    )
    if selected_per_neuron_keys is not None:
        significant = np.asarray(
            [
                int(cell_id)
                for cell_id in significant
                if per_neuron_selection_key(group, subject, date, session, fov_name, int(cell_id))
                in selected_per_neuron_keys
            ],
            dtype=int,
        )
    elif max_neuron_plots is not None:
        significant = significant[:max_neuron_plots]

    edges = get_edges(result)
    timestamps = np.ravel(result.get("timestamps", np.arange(maps.shape[3], dtype=float)))
    rf_types = ("ON", "OFF", "ON+OFF")

    index_rows = []
    n_plotted = 0
    n_pass = 0
    n_fail = 0
    for cell_id in significant:
        plot_map, mx_time, best_subfield = plot_map_for_neuron(
            maps,
            result,
            int(cell_id),
            flip_azimuth=flip_azimuth,
        )
        if plot_map is None:
            continue

        passes_ev = (
            np.isfinite(ev[cell_id])
            and np.isfinite(p_values[cell_id])
            and ev[cell_id] > min_ev
            and p_values[cell_id] < max_p_value
        )
        passes_peak = np.isfinite(peak_to_noise[cell_id]) and peak_to_noise[cell_id] > min_peak_to_noise
        target_dir = pass_dir if passes_peak else fail_dir
        if passes_peak:
            n_pass += 1
        else:
            n_fail += 1

        mx = np.nanmax(np.abs(plot_map))
        if not np.isfinite(mx) or mx == 0:
            mx = 1.0
        display_map = plot_map

        fig, ax = plt.subplots(figsize=(9, 7.25), constrained_layout=True)
        extent = rf_image_extent(edges)
        im = ax.imshow(
            display_map,
            extent=extent,
            origin="upper",
            cmap=RF_PER_NEURON_CMAP,
            norm=TwoSlopeNorm(vmin=-mx, vcenter=0.0, vmax=mx),
        )
        ax.set_aspect("equal")
        ax.set_ylim(rf_elevation_limits(edges))
        ax.set_xlabel("Azimuth (deg)")
        ax.set_ylabel("Elevation (deg)")
        cb = plt.colorbar(im, ax=ax)
        cb.set_label("Signed RF amplitude")
        plot_gaussian_on_rf(ax, result, int(cell_id), flip_azimuth=flip_azimuth)

        rf_type = rf_types[best_subfield - 1] if best_subfield in (1, 2, 3) else "n/a"
        delay_text = "delay n/a"
        if 0 <= mx_time < timestamps.size:
            delay_text = f"delay {float(timestamps[mx_time]):.3f} s"

        ax.set_title(
            "\n".join(
                [
                    f"RF result ROI {cell_id + 1} | source ROI0 {source_roi_text(result, cell_id)} | {rf_type} | {delay_text}",
                    f"EV test: {pass_fail(passes_ev)} | RF EV {ev[cell_id]:.3f} > {min_ev:.3f}, p {p_values[cell_id]:.3f} < {max_p_value:.3f} | full EV {full_ev[cell_id]:.3f}",
                    f"Peak/noise test: {pass_fail(passes_peak)} | P/N {peak_to_noise[cell_id]:.2f} > {min_peak_to_noise:.2f}",
                    "RF-map colors: negative light gray, zero white, strong positive #00ff00, peak positive #126d36",
                ]
            )
        )

        base_name = (
            f"RFresult_{cell_id + 1:04d}_sourceROI0_{source_roi_text(result, cell_id)}_"
            f"EV_{pass_fail(passes_ev)}_peakToNoise_{pass_fail(passes_peak)}"
        )
        style_figure_for_export(fig)
        save_figure_outputs(fig, target_dir / base_name, dpi=220)
        if selected_good_rf_collection_dir is not None and selected_per_neuron_keys is not None:
            collection_name = (
                f"{subject}_{date}_session{session}_{fov_name}_"
                f"RFresult_{cell_id + 1:04d}_sourceROI0_{source_roi_text(result, cell_id)}"
            )
            save_figure_outputs(fig, selected_good_rf_collection_dir / collection_name, dpi=220)
        plt.close(fig)
        index_rows.append(
            {
                "rf_result_index_0based": int(cell_id),
                "rf_result_index_1based": int(cell_id + 1),
                "source_roi_index_0based": source_roi_text(result, int(cell_id)),
                "rf_ev": float(ev[cell_id]) if np.isfinite(ev[cell_id]) else np.nan,
                "p_value": float(p_values[cell_id]) if np.isfinite(p_values[cell_id]) else np.nan,
                "full_ev": float(full_ev[cell_id]) if np.isfinite(full_ev[cell_id]) else np.nan,
                "peak_to_noise": float(peak_to_noise[cell_id]) if np.isfinite(peak_to_noise[cell_id]) else np.nan,
                "passes_ev_and_p": bool(passes_ev),
                "passes_peak_to_noise": bool(passes_peak),
                "selected_by_group_diverse_sampler": selected_per_neuron_keys is not None,
                "rf_map_azimuth_flipped": bool(flip_azimuth),
                "folder": target_dir.name,
            }
        )
        n_plotted += 1

    write_metadata(index_rows, output_dir / "per_neuron_rf_map_index.csv")
    print(
        f"[PER-NEURON RF MAPS] {result.results_dir.name}: "
        f"{n_pass} peak/noise PASS -> {pass_dir.name}; "
        f"{n_fail} peak/noise FAIL -> {fail_dir.name}"
    )
    return n_plotted


def plot_map_for_neuron(
    maps: np.ndarray,
    result: FanciResult,
    cell_id: int,
    *,
    flip_azimuth: bool = False,
) -> tuple[np.ndarray | None, int, int]:
    selected_rf, best_subfield, _, _ = selected_rf_cube_for_neuron(result, cell_id)
    if selected_rf is None:
        return None, -1, -1

    optimal = result.get("optimalDelays")
    mx_time = -1
    if optimal is not None and cell_id < np.ravel(optimal).size and np.isfinite(np.ravel(optimal)[cell_id]):
        mx_time = int(round(float(np.ravel(optimal)[cell_id]))) - 1
    if mx_time < 0 or mx_time >= selected_rf.shape[2]:
        peak_per_delay = np.nanmax(selected_rf, axis=(0, 1))
        mx_time = int(np.nanargmax(peak_per_delay))

    plot_map = selected_rf[:, :, mx_time]
    if flip_azimuth:
        plot_map = np.fliplr(plot_map)
    return plot_map, mx_time, best_subfield


def selected_rf_cube_for_neuron(result: FanciResult, cell_id: int) -> tuple[np.ndarray | None, int, np.ndarray, np.ndarray]:
    maps = result.require("maps")
    rf = np.asarray(maps[cell_id, :, :, :, :], dtype=float).copy()
    if rf.ndim != 4 or np.all(np.isnan(rf)):
        return None, -1, np.full(2, np.nan), np.full(3, np.nan)

    # Match the MATLAB plotting convention: flip OFF sign for visual comparison.
    rf[:, :, :, 1] = -rf[:, :, :, 1]
    rf_delay_mean = np.nanmean(rf, axis=2)
    rf_mean = np.full((rf_delay_mean.shape[0], rf_delay_mean.shape[1], 3), np.nan, dtype=float)
    rf_mean[:, :, :2] = rf_delay_mean
    signs = np.full(2, np.nan)
    subs = np.full(3, np.nan)

    for sub in range(2):
        r = rf_mean[:, :, sub]
        if np.all(np.isnan(r)):
            continue
        idx = np.nanargmax(np.abs(r))
        signs[sub] = np.sign(np.ravel(r)[idx])
        if signs[sub] == 0:
            signs[sub] = 1
        subs[sub] = np.ravel(np.abs(r))[idx]

    rf_mean[:, :, 2] = (rf_mean[:, :, 0] * signs[0] + rf_mean[:, :, 1] * signs[1]) / 2.0
    subs[2] = np.nanmax(rf_mean[:, :, 2])
    best_subfield = int(np.nanargmax(subs)) + 1
    if subs[2] > 0.7 * np.nanmax(subs):
        best_subfield = 3

    if best_subfield < 3:
        selected_rf = rf[:, :, :, best_subfield - 1] * signs[best_subfield - 1]
    else:
        selected_rf = (rf[:, :, :, 0] * signs[0] + rf[:, :, :, 1] * signs[1]) / 2.0

    return selected_rf, best_subfield, signs, subs


def plot_gaussian_on_rf(ax, result: FanciResult, cell_id: int, *, flip_azimuth: bool = False) -> None:
    gauss = result.get("gaussPars")
    edges = get_edges(result)
    if gauss is None or cell_id >= gauss.shape[0]:
        return
    pars = gaussian_display_pars(gauss[[cell_id], :], flip_azimuth=flip_azimuth)[0, :]
    if pars.size < 6 or not np.all(np.isfinite(pars[:6])):
        return

    ellipse_x = np.linspace(-np.pi, np.pi, 100)
    x = pars[2] * np.cos(ellipse_x) * 2.0
    y = pars[4] * np.sin(ellipse_x) * 2.0
    x_rot = pars[1] + x * np.cos(pars[5]) - y * np.sin(pars[5])
    y_rot = pars[3] + x * np.sin(pars[5]) + y * np.cos(pars[5])
    az_limits = tuple(np.sort(edges[:2]))
    el_limits = tuple(np.sort(edges[2:4]))
    outside = (x_rot < az_limits[0]) | (x_rot > az_limits[1]) | (y_rot < el_limits[0]) | (y_rot > el_limits[1])
    x_rot[outside] = np.nan
    y_rot[outside] = np.nan
    ax.plot(x_rot, y_rot, color="black", linewidth=1.5)


# =============================================================================
# New plot 1: Circular-shift EV distribution, p<0.05 vs p>=0.05
# =============================================================================


def plot_shift_ev_histogram(
    result: FanciResult,
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    fov_name: str,
    max_p_value: float,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    n = result.require("maps").shape[0]
    p_values = result_column(result, "pValues", n)
    ev_values, ev_label = ev_for_shift_plots(result)

    finite = np.isfinite(ev_values) & np.isfinite(p_values)
    sig = finite & (p_values < max_p_value)
    nonsig = finite & (p_values >= max_p_value)
    if sig.sum() == 0 and nonsig.sum() == 0:
        return 0

    combined = ev_values[finite]
    lo, hi = np.nanpercentile(combined, [0, 99.5])
    lo = min(float(lo), float(np.nanmin(combined)))
    hi = max(float(hi), float(np.nanmax(combined)))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = -0.01, 0.1
    bins = np.linspace(lo, hi, 40)

    fig, ax = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
    if sig.sum() > 0:
        ax.hist(
            ev_values[sig],
            bins=bins,
            weights=np.full(sig.sum(), 100.0 / sig.sum()),
            histtype="stepfilled",
            color="0.72",
            edgecolor="0.42",
            alpha=0.45,
            label=f"p < {max_p_value:g} (n={sig.sum()})",
        )
        ax.hist(
            ev_values[sig],
            bins=bins,
            weights=np.full(sig.sum(), 100.0 / sig.sum()),
            histtype="step",
            color="0.35",
            linewidth=1.5,
        )
    if nonsig.sum() > 0:
        ax.hist(
            ev_values[nonsig],
            bins=bins,
            weights=np.full(nonsig.sum(), 100.0 / nonsig.sum()),
            histtype="stepfilled",
            color="0.82",
            edgecolor="0.55",
            alpha=0.35,
            label=f"p >= {max_p_value:g} (n={nonsig.sum()})",
        )
        ax.hist(
            ev_values[nonsig],
            bins=bins,
            weights=np.full(nonsig.sum(), 100.0 / nonsig.sum()),
            histtype="step",
            color="0.55",
            linewidth=1.5,
            linestyle="--",
        )
    mean_ev = float(np.nanmean(ev_values[finite]))
    ax.axvline(mean_ev, color="black", linewidth=2.2, label="mean")
    ax.set_xlabel(ev_label)
    ax.set_ylabel("Neurons (% within group)")
    ax.set_title(f"{subject} {date} session {session} {fov_name}: circular-shift RF EV distribution")
    make_two_tick_square_axes(ax, x_zero=True, y_zero=True)

    save_figure(fig, output_dir / f"{fov_name}_circularShift_EV_hist_pSig_vs_nonSig", dpi=250)
    return 1


# =============================================================================
# New plot 2: Predicted versus actual preprocessed calcium traces
# =============================================================================


def plot_predicted_vs_actual_traces(
    result: FanciResult,
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    fov_name: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    traces = result.get("tracesPreprocessed")
    times = result.get("tracesPreprocessed_times")
    preds = result.get("predictionsFullModel")

    if traces is None or times is None or preds is None:
        print(
            f"[WARN] {fov_name}: predicted-vs-actual skipped. Need tracesPreprocessed, "
            "tracesPreprocessed_times, and predictionsFullModel from _FANCi_rf.*.npy "
            "outputs or MAT fallback."
        )
        return 0

    traces = orient_time_by_cell(np.asarray(traces, dtype=float), n_cells=result.require("maps").shape[0])
    preds = align_predictions_to_traces(np.asarray(preds, dtype=float), traces)
    times = np.ravel(np.asarray(times, dtype=float))
    if times.size != traces.shape[0]:
        times = np.arange(traces.shape[0], dtype=float)

    top_ids = top_good_by_peak_noise(
        result,
        min_ev=min_ev,
        max_p_value=max_p_value,
        min_peak_to_noise=min_peak_to_noise,
        n_top=5,
    )
    median_id = neuron_closest_to_median_shift_ev(result)
    ids = list(map(int, top_ids))
    labels = ["top_peak_noise"] * len(ids)
    if median_id is not None and int(median_id) not in ids:
        ids.append(int(median_id))
        labels.append("median_shift_ev")
    elif median_id is not None:
        labels[ids.index(int(median_id))] += "_and_median_shift_ev"

    if not ids:
        return 0

    n_saved = 0
    for cell_id, label in zip(ids, labels):
        if cell_id >= traces.shape[1] or cell_id >= preds.shape[1]:
            continue
        finite = np.isfinite(traces[:, cell_id]) & np.isfinite(preds[:, cell_id]) & np.isfinite(times)
        if finite.sum() == 0:
            continue
        idx = np.flatnonzero(finite)[:: max(1, int(TRACE_PLOT_STRIDE))]
        ev_values, ev_label = ev_for_shift_plots(result)
        peak_to_noise = result_column(result, "peakToNoise", result.require("maps").shape[0])
        p_values = result_column(result, "pValues", result.require("maps").shape[0])

        fig, ax = plt.subplots(figsize=(13, 4.8), constrained_layout=True)
        ax.plot(times[idx], traces[idx, cell_id], linewidth=0.85, label="actual preprocessed calcium trace")
        ax.plot(times[idx], preds[idx, cell_id], linewidth=1.1, label="full model prediction")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("z-scored preprocessed trace / prediction")
        ax.set_title(
            f"{subject} {date} session {session} {fov_name}: predicted vs actual | "
            f"RF result ROI {cell_id + 1}, source ROI0 {source_roi_text(result, cell_id)}\n"
            f"{label}; P/N={peak_to_noise[cell_id]:.3g}, p={p_values[cell_id]:.3g}, {ev_label}={ev_values[cell_id]:.4g}"
        )
        ax.legend(frameon=False, loc="upper right")
        ax.grid(True, alpha=0.2)
        save_figure(
            fig,
            output_dir / f"{fov_name}_RFresult_{cell_id + 1:04d}_sourceROI0_{source_roi_text(result, cell_id)}_{label}_predicted_vs_actual",
            dpi=220,
        )
        n_saved += 1
    return n_saved


def align_predictions_to_traces(preds: np.ndarray, traces: np.ndarray) -> np.ndarray:
    preds = orient_time_by_cell(preds, n_cells=traces.shape[1])
    if preds.shape[0] == traces.shape[0]:
        return preds
    finite_rows = np.any(np.isfinite(preds), axis=1)
    if finite_rows.sum() == traces.shape[0]:
        return preds[finite_rows, :]
    n = min(preds.shape[0], traces.shape[0])
    return preds[:n, :]


def neuron_closest_to_median_shift_ev(result: FanciResult) -> int | None:
    ev_values, _ = ev_for_shift_plots(result)
    finite_ids = np.flatnonzero(np.isfinite(ev_values))
    if finite_ids.size == 0:
        return None
    median_ev = np.nanmedian(ev_values[finite_ids])
    return int(finite_ids[np.nanargmin(np.abs(ev_values[finite_ids] - median_ev))])


# =============================================================================
# New plot 3: Peak time courses at Gaussian RF center
# =============================================================================


def plot_peak_time_courses(
    result: FanciResult,
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    fov_name: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    top_ids = top_good_by_peak_noise(
        result,
        min_ev=min_ev,
        max_p_value=max_p_value,
        min_peak_to_noise=min_peak_to_noise,
        n_top=5,
    )
    if top_ids.size == 0:
        return 0

    maps = result.require("maps")
    gauss = result.require("gaussPars")
    edges = get_edges(result)
    lags = np.ravel(result.get("timestamps", np.arange(maps.shape[3], dtype=float)))
    peak_to_noise = result_column(result, "peakToNoise", maps.shape[0])
    p_values = result_column(result, "pValues", maps.shape[0])
    ev_values, ev_label = ev_for_shift_plots(result)

    n_rows, n_cols, n_lags = maps.shape[1], maps.shape[2], maps.shape[3]
    x_centers = bin_centers(edges[0], edges[1], n_cols)
    y_centers = bin_centers(edges[2], edges[3], n_rows)

    fig, axes = plt.subplots(top_ids.size, 1, figsize=(9.5, 2.65 * top_ids.size), sharex=True, constrained_layout=True)
    if top_ids.size == 1:
        axes = np.asarray([axes])

    rows = []
    for ax, cell_id in zip(axes, top_ids):
        pars = gauss[int(cell_id), :]
        x0, y0 = float(pars[1]), float(pars[3])
        row = int(np.nanargmin(np.abs(y_centers - y0)))
        col = int(np.nanargmin(np.abs(x_centers - x0)))
        selected_rf, best_subfield, signs, _ = selected_rf_cube_for_neuron(result, int(cell_id))
        if selected_rf is None:
            continue
        selected_profile = selected_rf[row, col, :]

        rf = np.asarray(maps[int(cell_id), :, :, :, :], dtype=float).copy()
        rf[:, :, :, 1] = -rf[:, :, :, 1]
        on_profile = rf[row, col, :, 0]
        off_profile = rf[row, col, :, 1]

        ax.plot(lags[:n_lags], on_profile, linewidth=0.9, alpha=0.45, label="ON at center")
        ax.plot(lags[:n_lags], off_profile, linewidth=0.9, alpha=0.45, label="OFF sign-flipped at center")
        ax.plot(lags[:n_lags], selected_profile, linewidth=2.0, label=f"selected profile ({subfield_name(best_subfield)})")
        ax.axhline(0, color="black", linewidth=0.7, alpha=0.45)
        ax.set_ylabel("RF amplitude")
        ax.set_title(
            f"RF result ROI {int(cell_id) + 1}, source ROI0 {source_roi_text(result, int(cell_id))}: "
            f"center=({x0:.1f} az, {y0:.1f} el), grid row/col=({row},{col}), "
            f"P/N={peak_to_noise[int(cell_id)]:.2f}, p={p_values[int(cell_id)]:.3g}, {ev_label}={ev_values[int(cell_id)]:.4g}",
            fontsize=10,
        )
        ax.grid(True, alpha=0.2)
        rows.append(
            {
                "rf_result_index_0based": int(cell_id),
                "rf_result_index_1based": int(cell_id) + 1,
                "source_roi_index_0based": int(source_roi_value(result, int(cell_id))),
                "source_roi_index_matlab": int(source_roi_matlab_value(result, int(cell_id))),
                "gauss_azimuth": x0,
                "gauss_elevation": y0,
                "grid_row_0based": row,
                "grid_col_0based": col,
                "grid_azimuth": float(x_centers[col]),
                "grid_elevation": float(y_centers[row]),
                "best_subfield": subfield_name(best_subfield),
                "peak_to_noise": float(peak_to_noise[int(cell_id)]),
                "p_value": float(p_values[int(cell_id)]),
                "shift_ev": float(ev_values[int(cell_id)]),
            }
        )

    axes[-1].set_xlabel("Temporal lag (s)")
    axes[0].legend(frameon=False, ncol=3, loc="upper right", fontsize=8)
    fig.suptitle(
        f"{subject} {date} session {session} {fov_name}: top 5 peak/noise RF temporal profiles at Gaussian center",
        fontweight="bold",
    )
    save_figure(fig, output_dir / f"{fov_name}_top5_peakNoise_RF_center_temporal_profiles", dpi=250)
    write_metadata(rows, output_dir / f"{fov_name}_top5_peakNoise_RF_center_temporal_profiles.csv")
    return 1


# =============================================================================
# FOV overlays
# =============================================================================


def plot_neutral_roi_overlays(
    raw_folder: Path,
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    fov_name: str,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    plane_id = plane_id_from_fov_name(fov_name, fallback=0)
    plane_folder = raw_folder / f"plane{plane_id}"
    if not plane_folder.is_dir():
        print(f"[WARN] {fov_name}: neutral ROI overlay skipped; missing {plane_folder}")
        return 0
    required = [plane_folder / "stat.npy", plane_folder / "iscell.npy", plane_folder / "reg_outputs.npy"]
    if not all(path.exists() for path in required):
        print(f"[WARN] {fov_name}: neutral ROI overlay skipped; missing one of {required}")
        return 0

    stat = np.load(plane_folder / "stat.npy", allow_pickle=True)
    iscell = np.load(plane_folder / "iscell.npy", allow_pickle=True)
    reg_outputs = load_plane_reg_outputs(plane_folder)
    fov = pick_fov_image(reg_outputs)

    all_indices = np.arange(len(stat), dtype=int)
    cell_indices = all_indices[iscell[:, 0].astype(bool)]

    plot_rois(
        fov=fov,
        stat=stat,
        roi_indices=all_indices,
        title=f"{subject} {date} session {session} {fov_name}: all Suite2p ROIs ({all_indices.size} total)",
        output_base=output_dir / f"{fov_name}_all_suite2p_ROIs_on_FOV",
    )
    plot_rois(
        fov=fov,
        stat=stat,
        roi_indices=cell_indices,
        title=f"{subject} {date} session {session} {fov_name}: Suite2p cell ROIs ({cell_indices.size} cells)",
        output_base=output_dir / f"{fov_name}_suite2p_cell_ROIs_on_FOV",
    )
    return 2


def plot_rf_roi_overlay(
    result: FanciResult,
    raw_folder: Path,
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    fov_name: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    flip_azimuth: bool = False,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    plane_id = plane_id_from_fov_name(fov_name, fallback=0)
    plane_folder = raw_folder / f"plane{plane_id}"
    if not plane_folder.is_dir():
        print(f"[WARN] {fov_name}: RF ROI overlay skipped; missing {plane_folder}")
        return 0
    required = [plane_folder / "stat.npy", plane_folder / "reg_outputs.npy"]
    if not all(path.exists() for path in required):
        print(f"[WARN] {fov_name}: RF ROI overlay skipped; missing one of {required}")
        return 0

    stat = np.load(plane_folder / "stat.npy", allow_pickle=True)
    reg_outputs = load_plane_reg_outputs(plane_folder)
    fov = pick_fov_image(reg_outputs)

    gauss_pars = gaussian_display_pars(result.require("gaussPars"), flip_azimuth=flip_azimuth)
    n_neurons = gauss_pars.shape[0]
    ev = result_column(result, "explVars", n_neurons)
    p_values = result_column(result, "pValues", n_neurons)
    peak_to_noise = result_column(result, "peakToNoise", n_neurons)
    roi_indices = np.ravel(np.asarray(result.get("roiIndex0", np.arange(n_neurons)), dtype=int))

    azimuth = gauss_pars[:, 1]
    elevation = gauss_pars[:, 3]
    has_gaussian = np.all(np.isfinite(gauss_pars[:, :6]), axis=1)
    significant = np.isfinite(ev) & np.isfinite(p_values) & (ev > min_ev) & (p_values < max_p_value) & has_gaussian
    peak_pass = significant & np.isfinite(peak_to_noise) & (peak_to_noise > min_peak_to_noise)

    edges = get_edges(result)
    azimuth_limits = tuple(np.sort(edges[:2]))
    elevation_limits = tuple(np.sort(edges[2:4]))
    bivariate_azimuth_limits, bivariate_elevation_limits = fov_overlay_bivariate_plot_limits(
        azimuth_limits,
        elevation_limits,
    )

    rf_colormap = make_rf_bivariate_colormap(bivariate_azimuth_limits, bivariate_elevation_limits)

    fig = plt.figure(figsize=(11.5, 8.5), constrained_layout=False)
    ax_fov = fig.add_axes([0.06, 0.07, 0.60, 0.82])
    ax_wheel = fig.add_axes([0.705, 0.425, 0.245, 0.185])

    fov_artist = ax_fov.imshow(fov, cmap="gray", origin="upper")
    fov_artist.set_rasterized(True)
    draw_rf_roi_masks(
        ax=ax_fov,
        stat=stat,
        image_shape=fov.shape,
        source_roi_indices=roi_indices,
        azimuth=azimuth,
        elevation=elevation,
        significant=significant,
        peak_pass=peak_pass,
        azimuth_limits=bivariate_azimuth_limits,
        elevation_limits=bivariate_elevation_limits,
        rf_colormap=rf_colormap,
    )
    ax_fov.set_title(
        f"{fov_name} FOV overlay: pycolorbar RF center color ({int(significant.sum())} significant, {int(peak_pass.sum())} P/N pass)",
        fontsize=15,
    )
    ax_fov.set_axis_off()
    plot_2d_colorwheel(
        ax_wheel,
        bivariate_azimuth_limits,
        bivariate_elevation_limits,
        rf_colormap=rf_colormap,
        azimuth=azimuth,
        elevation=elevation,
        significant=significant,
        peak_pass=peak_pass,
    )
    fig.suptitle(
        f"{subject} {date} session {session} {fov_name}: significant RF ROIs on Suite2p FOV | "
        f"pycolorbar={PYCOLORBAR_BIVARIATE_CMAP_NAME}; filled ROIs = P/N > {min_peak_to_noise:g}; contours = significant only",
        fontsize=17,
    )
    save_figure(fig, output_dir / f"{fov_name}_significant_ROIs_on_FOV_2D_colorwheel", dpi=220)

    fig = plt.figure(figsize=(11.5, 8.5), constrained_layout=False)
    ax_fov = fig.add_axes([0.06, 0.07, 0.60, 0.82])
    ax_wheel = fig.add_axes([0.705, 0.425, 0.245, 0.185])

    fov_artist = ax_fov.imshow(fov, cmap="gray", origin="upper")
    fov_artist.set_rasterized(True)
    draw_rf_roi_masks(
        ax=ax_fov,
        stat=stat,
        image_shape=fov.shape,
        source_roi_indices=roi_indices,
        azimuth=azimuth,
        elevation=elevation,
        significant=peak_pass,
        peak_pass=peak_pass,
        azimuth_limits=bivariate_azimuth_limits,
        elevation_limits=bivariate_elevation_limits,
        rf_colormap=rf_colormap,
    )
    ax_fov.set_title(
        f"{fov_name} FOV overlay: pycolorbar RF center color, good RFs only (n={int(peak_pass.sum())})",
        fontsize=15,
    )
    ax_fov.set_axis_off()
    plot_2d_colorwheel(
        ax_wheel,
        bivariate_azimuth_limits,
        bivariate_elevation_limits,
        rf_colormap=rf_colormap,
        azimuth=azimuth,
        elevation=elevation,
        significant=peak_pass,
        peak_pass=peak_pass,
    )
    fig.suptitle(
        f"{subject} {date} session {session} {fov_name}: good RF ROIs on Suite2p FOV | "
        f"pycolorbar={PYCOLORBAR_BIVARIATE_CMAP_NAME}; P/N > {min_peak_to_noise:g}",
        fontsize=17,
    )
    save_figure(fig, output_dir / f"{fov_name}_good_RFs_on_FOV_2D_colorwheel", dpi=220)

    plot_scalar_rf_roi_overlay(
        fov=fov,
        stat=stat,
        source_roi_indices=roi_indices,
        values=elevation,
        significant=significant,
        peak_pass=peak_pass,
        value_limits=good_value_limits(elevation, peak_pass, fallback=elevation_limits),
        value_label="Gaussian elevation (deg)",
        dimension_name="elevation",
        output_base=output_dir / f"{fov_name}_good_RFs_on_FOV_elevation_jet",
        title=(
            f"{subject} {date} session {session} {fov_name}: Good RF ROIs on Suite2p FOV | "
            f"elevation color; P/N > {min_peak_to_noise:g}"
        ),
        fov_title=f"{fov_name} FOV overlay: Good RF elevation, jet (n={int(peak_pass.sum())})",
    )
    plot_scalar_rf_roi_overlay(
        fov=fov,
        stat=stat,
        source_roi_indices=roi_indices,
        values=azimuth,
        significant=significant,
        peak_pass=peak_pass,
        value_limits=good_value_limits(azimuth, peak_pass, fallback=azimuth_limits),
        value_label="Gaussian azimuth (deg)",
        dimension_name="azimuth",
        output_base=output_dir / f"{fov_name}_good_RFs_on_FOV_azimuth_jet",
        title=(
            f"{subject} {date} session {session} {fov_name}: Good RF ROIs on Suite2p FOV | "
            f"azimuth color; P/N > {min_peak_to_noise:g}"
        ),
        fov_title=f"{fov_name} FOV overlay: Good RF azimuth, jet (n={int(peak_pass.sum())})",
    )

    rows = []
    for rf_idx in np.flatnonzero(significant):
        roi_idx = int(roi_indices[rf_idx])
        if roi_idx >= len(stat):
            continue
        roi = stat_entry(stat[roi_idx])
        y_med, x_med = roi_center(roi)
        rows.append(
            {
                "rf_result_index_1based": int(rf_idx + 1),
                "fov_name": fov_name,
                "plane_id": int(plane_id),
                "suite2p_roi_index_0based": roi_idx,
                "suite2p_roi_index_matlab": roi_idx + 1,
                "x_med": float(x_med),
                "y_med": float(y_med),
                    "azimuth_deg": float(azimuth[rf_idx]),
                    "elevation_deg": float(elevation[rf_idx]),
                    "azimuth_flipped": bool(flip_azimuth),
                    "explained_variance": float(ev[rf_idx]),
                "p_value": float(p_values[rf_idx]),
                "peak_to_noise": float(peak_to_noise[rf_idx]),
                "passes_peak_to_noise": bool(peak_pass[rf_idx]),
            }
        )
    write_metadata(rows, output_dir / f"{fov_name}_significant_roi_fov_overlay_metadata.csv")
    return 4
# =============================================================================
# Combined session-level plots across all FOVs. Existing per-FOV plots are left
# unchanged; these outputs are written under plots_session / "combined".
# =============================================================================


def plot_combined_session_outputs(
    records: list[dict],
    *,
    output_dir: Path,
    raw_suite2p_folder: Path | None,
    subject: str,
    date: str,
    session: str,
    group: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    max_neuron_plots: int | None,
    selected_per_neuron_keys: set[tuple[str, str, str, str, str, int]] | None,
    skip_fov_overlays: bool,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    counts: dict[str, Any] = {
        "combined_gaussian_coverage": 0,
        "combined_cell_roi_overlay": 0,
        "combined_rf_roi_overlay": 0,
        "combined_rf_size_lineplot": 0,
        "combined_shift_ev_lineplot": 0,
        "combined_per_neuron_rf_maps": 0,
        "errors": [],
    }

    if not records:
        counts["errors"].append("No FOV records were loaded for combined plotting.")
        return counts

    try:
        counts["combined_gaussian_coverage"] = plot_combined_gaussian_coverage(
            records,
            output_dir=output_dir / "RF_gaussian_coverage",
            subject=subject,
            date=date,
            session=session,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
        )
    except Exception as exc:
        message = f"combined Gaussian coverage failed: {exc}"
        counts["errors"].append(message)
        print(f"[WARN] {message}")

    if not skip_fov_overlays and raw_suite2p_folder is not None:
        if MAKE_NEUTRAL_ROI_OVERLAYS:
            try:
                counts["combined_cell_roi_overlay"] = plot_combined_cell_roi_overlay(
                    records,
                    raw_suite2p_folder=raw_suite2p_folder,
                    output_dir=output_dir / "FOV_overlays",
                    subject=subject,
                    date=date,
                    session=session,
                )
            except Exception as exc:
                message = f"combined Suite2p cell ROI overlay failed: {exc}"
                counts["errors"].append(message)
                print(f"[WARN] {message}")

        try:
            counts["combined_rf_roi_overlay"] = plot_combined_rf_roi_overlay(
                records,
                raw_suite2p_folder=raw_suite2p_folder,
                output_dir=output_dir / "FOV_overlays",
                subject=subject,
                date=date,
                session=session,
                min_ev=min_ev,
                max_p_value=max_p_value,
                min_peak_to_noise=min_peak_to_noise,
            )
        except Exception as exc:
            message = f"combined RF/FOV overlay failed: {exc}"
            counts["errors"].append(message)
            print(f"[WARN] {message}")
    else:
        print("[WARN] combined FOV overlay plots skipped; raw Suite2p folder unavailable or --skip-fov-overlays was set.")

    try:
        counts["combined_rf_size_lineplot"] = plot_combined_rf_size_lineplot(
            records,
            output_dir=output_dir / "RF_size_histogram",
            subject=subject,
            date=date,
            session=session,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
        )
    except Exception as exc:
        message = f"combined RF-size line plot failed: {exc}"
        counts["errors"].append(message)
        print(f"[WARN] {message}")

    try:
        counts["combined_shift_ev_lineplot"] = plot_combined_shift_ev_lineplot(
            records,
            output_dir=output_dir / "shift_ev_diagnostics",
            subject=subject,
            date=date,
            session=session,
            max_p_value=max_p_value,
        )
    except Exception as exc:
        message = f"combined shift-EV line plot failed: {exc}"
        counts["errors"].append(message)
        print(f"[WARN] {message}")

    try:
        counts["combined_per_neuron_rf_maps"] = plot_combined_per_neuron_rf_maps(
            records,
            output_dir=output_dir / "per_neuron_metadata",
            group=group,
            subject=subject,
            date=date,
            session=session,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
            max_neuron_plots=max_neuron_plots,
            selected_per_neuron_keys=selected_per_neuron_keys,
        )
    except Exception as exc:
        message = f"combined per-neuron RF maps failed: {exc}"
        counts["errors"].append(message)
        print(f"[WARN] {message}")

    (output_dir / "combined_plot_summary.json").write_text(json.dumps(counts, indent=2), encoding="utf-8")
    return counts


def plot_combined_gaussian_coverage(
    records: list[dict],
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    n_azimuth_grid: int = 500,
    n_elevation_grid: int = 250,
    gaussian_sigma_radius: float = 2.0,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    all_gauss: list[np.ndarray] = []
    peak_ids: list[int] = []
    non_peak_ids: list[int] = []

    for record in records:
        fov_name = record["fov_name"]
        result = record["result"]
        gauss_pars = gaussian_display_pars(
            result.require("gaussPars"),
            flip_azimuth=bool(record.get("flip_azimuth", False)),
        )
        n_neurons = min(gauss_pars.shape[0], result.require("maps").shape[0])
        gauss_pars = gauss_pars[:n_neurons, :]
        ev = result_column(result, "explVars", n_neurons)
        p_values = result_column(result, "pValues", n_neurons)
        peak_to_noise = result_column(result, "peakToNoise", n_neurons)
        roi_index0 = np.ravel(np.asarray(result.get("roiIndex0", np.arange(n_neurons)), dtype=float))

        if gauss_pars.shape[1] < 6:
            continue

        has_gaussian = (
            np.all(np.isfinite(gauss_pars[:, :6]), axis=1)
            & (gauss_pars[:, 2] > 0)
            & (gauss_pars[:, 4] > 0)
        )
        significant = (
            has_gaussian
            & np.isfinite(ev)
            & np.isfinite(p_values)
            & (ev > min_ev)
            & (p_values < max_p_value)
        )
        peak_pass = significant & np.isfinite(peak_to_noise) & (peak_to_noise > min_peak_to_noise)

        for rf_idx in np.flatnonzero(significant):
            combined_idx = len(all_gauss)
            all_gauss.append(gauss_pars[rf_idx, :])
            if peak_pass[rf_idx]:
                peak_ids.append(combined_idx)
            else:
                non_peak_ids.append(combined_idx)
            source_roi0 = roi_index0[rf_idx] if rf_idx < roi_index0.size else np.nan
            rows.append(
                {
                    "fov": fov_name,
                    "rf_result_index_0based": int(rf_idx),
                    "rf_result_index_1based": int(rf_idx + 1),
                    "source_roi0": int(source_roi0) if np.isfinite(source_roi0) else "",
                    "azimuth_deg": float(gauss_pars[rf_idx, 1]),
                    "elevation_deg": float(gauss_pars[rf_idx, 3]),
                    "sigma_x_deg": float(gauss_pars[rf_idx, 2]),
                    "sigma_y_deg": float(gauss_pars[rf_idx, 4]),
                    "theta_rad": float(gauss_pars[rf_idx, 5]),
                    "explained_variance": float(ev[rf_idx]) if np.isfinite(ev[rf_idx]) else np.nan,
                    "p_value": float(p_values[rf_idx]) if np.isfinite(p_values[rf_idx]) else np.nan,
                    "peak_to_noise": float(peak_to_noise[rf_idx]) if np.isfinite(peak_to_noise[rf_idx]) else np.nan,
                    "passes_peak_to_noise": bool(peak_pass[rf_idx]),
                }
            )

    write_metadata(rows, output_dir / "combined_RF_gaussian_coverage_values.csv")

    if all_gauss:
        gauss_pars = np.vstack(all_gauss)
    else:
        gauss_pars = np.empty((0, 6), dtype=float)
    peak_ids_array = np.asarray(peak_ids, dtype=int)
    non_peak_ids_array = np.asarray(non_peak_ids, dtype=int)
    all_ids = np.arange(gauss_pars.shape[0], dtype=int)

    azimuth_limits, elevation_limits = combined_visual_limits(records)
    azimuth_axis = np.linspace(azimuth_limits[0], azimuth_limits[1], n_azimuth_grid)
    elevation_axis = np.linspace(elevation_limits[0], elevation_limits[1], n_elevation_grid)
    azimuth_grid, elevation_grid = np.meshgrid(azimuth_axis, elevation_axis)

    peak_coverage = gaussian_coverage_map(
        gauss_pars,
        peak_ids_array,
        azimuth_grid,
        elevation_grid,
        gaussian_sigma_radius,
    )
    all_coverage = gaussian_coverage_map(
        gauss_pars,
        all_ids,
        azimuth_grid,
        elevation_grid,
        gaussian_sigma_radius,
    )
    shared_max = np.nanmax([np.nanmax(peak_coverage), np.nanmax(all_coverage)])
    if not np.isfinite(shared_max) or shared_max <= 0:
        shared_max = 1.0

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), constrained_layout=True)
    plot_coverage_panel(
        axes[0],
        peak_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        shared_max,
        f"Peak/noise pass only (n = {peak_ids_array.size})",
        gauss_pars,
        peak_ids_array,
        np.empty(0, dtype=int),
        gaussian_sigma_radius,
    )
    plot_coverage_panel(
        axes[1],
        all_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        shared_max,
        (
            f"All significant RFs (n = {all_ids.size}; "
            f"P/N n = {peak_ids_array.size}, non-P/N n = {non_peak_ids_array.size})"
        ),
        gauss_pars,
        peak_ids_array,
        non_peak_ids_array,
        gaussian_sigma_radius,
    )
    fig.suptitle(
        f"{subject} {date} session {session}: combined Gaussian RF coverage | "
        f"EV > {min_ev:.3f}, p < {max_p_value:.3f}, P/N > {min_peak_to_noise:.2f}",
        fontweight="bold",
    )
    save_figure(fig, output_dir / "combined_RF_gaussian_coverage_heatmap", dpi=300)
    return 1


def combined_visual_limits(records: list[dict]) -> tuple[tuple[float, float], tuple[float, float]]:
    azimuth_limits_all = []
    elevation_limits_all = []
    for record in records:
        edges = get_edges(record["result"])
        azimuth_limits_all.extend(np.sort(edges[:2]).astype(float).tolist())
        elevation_limits_all.extend(np.sort(edges[2:4]).astype(float).tolist())
    if not azimuth_limits_all or not elevation_limits_all:
        return (-135.0, 135.0), (-40.0, 40.0)
    return (
        (float(np.nanmin(azimuth_limits_all)), float(np.nanmax(azimuth_limits_all))),
        (float(np.nanmin(elevation_limits_all)), float(np.nanmax(elevation_limits_all))),
    )


def rf_summary_vectors(
    result: FanciResult,
    *,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    flip_azimuth: bool = False,
) -> dict[str, np.ndarray]:
    maps = result.require("maps")
    gauss = gaussian_display_pars(result.require("gaussPars"), flip_azimuth=flip_azimuth)
    n = min(int(maps.shape[0]), int(gauss.shape[0]))
    ev = result_column(result, "explVars", maps.shape[0])[:n]
    p_values = result_column(result, "pValues", maps.shape[0])[:n]
    peak_to_noise = result_column(result, "peakToNoise", maps.shape[0])[:n]
    roi_indices = np.ravel(np.asarray(result.get("roiIndex0", np.arange(maps.shape[0])), dtype=int))[:n]
    gauss = gauss[:n, :]
    has_gaussian = np.zeros(n, dtype=bool)
    if gauss.shape[1] >= 6:
        has_gaussian = np.all(np.isfinite(gauss[:, :6]), axis=1)
    significant = np.isfinite(ev) & np.isfinite(p_values) & (ev > min_ev) & (p_values < max_p_value) & has_gaussian
    peak_pass = significant & np.isfinite(peak_to_noise) & (peak_to_noise > min_peak_to_noise)
    return {
        "ev": ev,
        "p_values": p_values,
        "peak_to_noise": peak_to_noise,
        "roi_indices": roi_indices,
        "gauss": gauss,
        "azimuth": gauss[:, 1] if gauss.shape[1] > 1 else np.full(n, np.nan),
        "elevation": gauss[:, 3] if gauss.shape[1] > 3 else np.full(n, np.nan),
        "significant": significant,
        "peak_pass": peak_pass,
    }


def plot_combined_rf_size_lineplot(
    records: list[dict],
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    rf_sizes: list[float] = []

    for record in records:
        result = record["result"]
        fov_name = record["fov_name"]
        gauss = np.asarray(result.require("gaussPars"), dtype=float)
        n_cells = result.require("maps").shape[0]
        n = min(n_cells, gauss.shape[0])
        gauss = gauss[:n, :]
        good = good_rf_mask(
            result,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
        )[:n]
        has_valid_width = (
            (gauss.shape[1] >= 5)
            & np.isfinite(gauss[:, 2])
            & np.isfinite(gauss[:, 4])
            & (np.abs(gauss[:, 2]) > 0)
            & (np.abs(gauss[:, 4]) > 0)
        )
        use = good & has_valid_width
        roi_index0 = np.ravel(np.asarray(result.get("roiIndex0", np.arange(n)), dtype=float))
        for rf_idx in np.flatnonzero(use):
            size_deg = float((abs(gauss[rf_idx, 2]) + abs(gauss[rf_idx, 4])) / 2.0)
            rf_sizes.append(size_deg)
            source_roi0 = roi_index0[rf_idx] if rf_idx < roi_index0.size else np.nan
            rows.append(
                {
                    "fov": fov_name,
                    "rf_result_index_0based": int(rf_idx),
                    "rf_result_index_1based": int(rf_idx + 1),
                    "source_roi0": int(source_roi0) if np.isfinite(source_roi0) else "",
                    "sigma_x_deg": float(gauss[rf_idx, 2]),
                    "sigma_y_deg": float(gauss[rf_idx, 4]),
                    "rf_size_deg": size_deg,
                }
            )

    write_metadata(rows, output_dir / "combined_goodRF_size_values.csv")
    rf_size = np.asarray(rf_sizes, dtype=float)

    fig, ax = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
    if rf_size.size > 0:
        max_size = float(np.nanmax(rf_size))
        bin_width = 2.0
        if np.isfinite(max_size) and max_size > 0:
            bins = np.arange(0.0, math.ceil(max_size / bin_width) * bin_width + bin_width, bin_width)
            if bins.size < 2:
                bins = np.linspace(0.0, max_size + bin_width, 10)
        else:
            bins = np.linspace(0.0, 1.0, 10)
        counts, edges = np.histogram(rf_size[np.isfinite(rf_size)], bins=bins)
        centers = (edges[:-1] + edges[1:]) / 2.0
        ax.plot(centers, counts, color="0.55", linewidth=1.8)
        mean_size = float(np.nanmean(rf_size))
        ax.axvline(mean_size, linestyle="-", linewidth=2.2, color="black", label="mean")
        stats_text = (
            f"n = {rf_size.size}\n"
            f"mean = {mean_size:.2f} deg\n"
            f"median = {np.nanmedian(rf_size):.2f} deg"
        )
    else:
        ax.text(
            0.5,
            0.5,
            "No good RF neurons with valid Gaussian widths across FOVs",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )
        stats_text = "n = 0"
    ax.text(
        0.98,
        0.95,
        stats_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )
    ax.set_xlabel("RF size (deg) = (|sigma_x| + |sigma_y|) / 2")
    ax.set_ylabel("Number of neurons")
    ax.set_title(
        f"{subject} {date} session {session}: combined good-RF size line plot\n"
        f"Good RF: EV > {min_ev:g}, p < {max_p_value:g}, P/N > {min_peak_to_noise:g}"
    )
    make_two_tick_square_axes(ax, x_zero=True, y_zero=True)
    save_figure(fig, output_dir / "combined_goodRF_size_lineplot", dpi=250)
    return 1


def plot_combined_shift_ev_lineplot(
    records: list[dict],
    output_dir: Path,
    *,
    subject: str,
    date: str,
    session: str,
    max_p_value: float,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    values: list[float] = []
    pvals: list[float] = []
    rows: list[dict] = []
    ev_label = "RF-only EV from final circular-shift test"

    for record in records:
        fov_name = record["fov_name"]
        result = record["result"]
        n = result.require("maps").shape[0]
        p_values = result_column(result, "pValues", n)
        ev_values, label = ev_for_shift_plots(result)
        ev_label = label
        for rf_idx, (ev_value, p_value) in enumerate(zip(ev_values, p_values)):
            values.append(float(ev_value) if np.isfinite(ev_value) else np.nan)
            pvals.append(float(p_value) if np.isfinite(p_value) else np.nan)
            rows.append(
                {
                    "fov": fov_name,
                    "rf_result_index_0based": int(rf_idx),
                    "rf_result_index_1based": int(rf_idx + 1),
                    "ev_value": float(ev_value) if np.isfinite(ev_value) else np.nan,
                    "p_value": float(p_value) if np.isfinite(p_value) else np.nan,
                    "p_lt_threshold": bool(np.isfinite(p_value) and p_value < max_p_value),
                }
            )

    write_metadata(rows, output_dir / "combined_circularShift_EV_values.csv")
    ev_values = np.asarray(values, dtype=float)
    p_values = np.asarray(pvals, dtype=float)
    finite = np.isfinite(ev_values) & np.isfinite(p_values)
    sig = finite & (p_values < max_p_value)
    nonsig = finite & (p_values >= max_p_value)
    if sig.sum() == 0 and nonsig.sum() == 0:
        return 0

    combined = ev_values[finite]
    lo, hi = np.nanpercentile(combined, [0, 99.5])
    lo = min(float(lo), float(np.nanmin(combined)))
    hi = max(float(hi), float(np.nanmax(combined)))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = -0.01, 0.1
    bins = np.linspace(lo, hi, 40)
    centers = (bins[:-1] + bins[1:]) / 2.0

    fig, ax = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
    if sig.sum() > 0:
        counts, _ = np.histogram(ev_values[sig], bins=bins)
        ax.plot(
            centers,
            counts.astype(float) * 100.0 / float(sig.sum()),
            color="0.45",
            linestyle="-",
            linewidth=1.8,
            label=f"p < {max_p_value:g} (n={sig.sum()})",
        )
    if nonsig.sum() > 0:
        counts, _ = np.histogram(ev_values[nonsig], bins=bins)
        ax.plot(
            centers,
            counts.astype(float) * 100.0 / float(nonsig.sum()),
            color="0.60",
            linestyle="--",
            linewidth=1.8,
            label=f"p >= {max_p_value:g} (n={nonsig.sum()})",
        )
    mean_ev = float(np.nanmean(ev_values[finite]))
    ax.axvline(mean_ev, color="black", linewidth=2.2, label="mean")
    ax.set_xlabel(ev_label)
    ax.set_ylabel("Neurons (% within group)")
    ax.set_title(f"{subject} {date} session {session}: combined circular-shift RF EV distribution")
    make_two_tick_square_axes(ax, x_zero=True, y_zero=True)
    save_figure(fig, output_dir / "combined_circularShift_EV_line_pSig_vs_nonSig", dpi=250)
    return 1


def plot_combined_per_neuron_rf_maps(
    records: list[dict],
    output_dir: Path,
    *,
    group: str,
    subject: str,
    date: str,
    session: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    max_neuron_plots: int | None,
    selected_per_neuron_keys: set[tuple[str, str, str, str, str, int]] | None,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    pass_dir = output_dir / "PASS_goodRFs_peakNoise"
    fail_dir = output_dir / "FAIL_significantRF_peakNoise"
    pass_dir.mkdir(parents=True, exist_ok=True)
    fail_dir.mkdir(parents=True, exist_ok=True)

    entries: list[dict] = []
    for record in records:
        result = record["result"]
        fov_name = record["fov_name"]
        maps = result.require("maps")
        n_neurons = maps.shape[0]
        ev = result_column(result, "explVars", n_neurons)
        p_values = result_column(result, "pValues", n_neurons)
        peak_to_noise = result_column(result, "peakToNoise", n_neurons)
        significant = np.flatnonzero(
            np.isfinite(ev)
            & np.isfinite(p_values)
            & (ev > min_ev)
            & (p_values < max_p_value)
        )
        for cell_id in significant:
            if (
                selected_per_neuron_keys is not None
                and per_neuron_selection_key(group, subject, date, session, fov_name, int(cell_id))
                not in selected_per_neuron_keys
            ):
                continue
            passes_peak = bool(np.isfinite(peak_to_noise[cell_id]) and peak_to_noise[cell_id] > min_peak_to_noise)
            entries.append(
                {
                    "record": record,
                    "cell_id": int(cell_id),
                    "passes_peak": passes_peak,
                    "peak_to_noise": float(peak_to_noise[cell_id]) if np.isfinite(peak_to_noise[cell_id]) else np.nan,
                }
            )

    if selected_per_neuron_keys is None and max_neuron_plots is not None:
        entries = entries[: int(max_neuron_plots)]

    rows: list[dict] = []
    n_plotted = 0
    for entry in entries:
        record = entry["record"]
        result = record["result"]
        fov_name = record["fov_name"]
        cell_id = int(entry["cell_id"])
        maps = result.require("maps")
        n_neurons = maps.shape[0]
        ev = result_column(result, "explVars", n_neurons)
        p_values = result_column(result, "pValues", n_neurons)
        full_ev = result_column(result, "explVarsFullModel", n_neurons)
        peak_to_noise = result_column(result, "peakToNoise", n_neurons)

        flip_azimuth = bool(record.get("flip_azimuth", False))
        plot_map, mx_time, best_subfield = plot_map_for_neuron(
            maps,
            result,
            cell_id,
            flip_azimuth=flip_azimuth,
        )
        if plot_map is None:
            continue

        passes_ev = bool(
            np.isfinite(ev[cell_id])
            and np.isfinite(p_values[cell_id])
            and ev[cell_id] > min_ev
            and p_values[cell_id] < max_p_value
        )
        passes_peak = bool(np.isfinite(peak_to_noise[cell_id]) and peak_to_noise[cell_id] > min_peak_to_noise)
        target_dir = pass_dir if passes_peak else fail_dir

        mx = np.nanmax(np.abs(plot_map))
        if not np.isfinite(mx) or mx == 0:
            mx = 1.0
        display_map = plot_map

        edges = get_edges(result)
        timestamps = np.ravel(result.get("timestamps", np.arange(maps.shape[3], dtype=float)))
        fig, ax = plt.subplots(figsize=(9, 7.25), constrained_layout=True)
        im = ax.imshow(
            display_map,
            extent=rf_image_extent(edges),
            origin="upper",
            cmap=RF_PER_NEURON_CMAP,
            norm=TwoSlopeNorm(vmin=-mx, vcenter=0.0, vmax=mx),
        )
        ax.set_aspect("equal")
        ax.set_ylim(rf_elevation_limits(edges))
        ax.set_xlabel("Azimuth (deg)")
        ax.set_ylabel("Elevation (deg)")
        cb = plt.colorbar(im, ax=ax)
        cb.set_label("Signed RF amplitude")
        plot_gaussian_on_rf(ax, result, cell_id, flip_azimuth=flip_azimuth)

        delay_text = "delay n/a"
        if 0 <= mx_time < timestamps.size:
            delay_text = f"delay {float(timestamps[mx_time]):.3f} s"
        ax.set_title(
            "\n".join(
                [
                    f"{fov_name} | RF result ROI {cell_id + 1} | source ROI0 {source_roi_text(result, cell_id)} | {subfield_name(best_subfield)} | {delay_text}",
                    f"EV test: {pass_fail(passes_ev)} | RF EV {ev[cell_id]:.3f} > {min_ev:.3f}, p {p_values[cell_id]:.3f} < {max_p_value:.3f} | full EV {full_ev[cell_id]:.3f}",
                    f"Peak/noise test: {pass_fail(passes_peak)} | P/N {peak_to_noise[cell_id]:.2f} > {min_peak_to_noise:.2f}",
                    "RF-map colors: negative light gray, zero white, strong positive #00ff00, peak positive #126d36",
                ]
            )
        )

        base_name = (
            f"{fov_name}_RFresult_{cell_id + 1:04d}_sourceROI0_{source_roi_text(result, cell_id)}_"
            f"EV_{pass_fail(passes_ev)}_peakToNoise_{pass_fail(passes_peak)}_signed_gray_white_darkgreen"
        )
        save_figure(fig, target_dir / base_name, dpi=220)
        rows.append(
            {
                "fov": fov_name,
                "rf_result_index_0based": int(cell_id),
                "rf_result_index_1based": int(cell_id + 1),
                "source_roi_index_0based": source_roi_text(result, cell_id),
                "rf_ev": float(ev[cell_id]) if np.isfinite(ev[cell_id]) else np.nan,
                "p_value": float(p_values[cell_id]) if np.isfinite(p_values[cell_id]) else np.nan,
                "full_ev": float(full_ev[cell_id]) if np.isfinite(full_ev[cell_id]) else np.nan,
                "peak_to_noise": float(peak_to_noise[cell_id]) if np.isfinite(peak_to_noise[cell_id]) else np.nan,
                "passes_ev_and_p": bool(passes_ev),
                "passes_peak_to_noise": bool(passes_peak),
                "selected_by_group_diverse_sampler": selected_per_neuron_keys is not None,
                "rf_map_azimuth_flipped": bool(flip_azimuth),
                "folder": target_dir.name,
            }
        )
        n_plotted += 1

    write_metadata(rows, output_dir / "combined_per_neuron_rf_map_index.csv")
    print(
        f"[COMBINED PER-NEURON RF MAPS] {sum(row['passes_peak_to_noise'] for row in rows)} peak/noise PASS; "
        f"{sum(not row['passes_peak_to_noise'] for row in rows)} peak/noise FAIL"
    )
    return n_plotted


def plot_combined_cell_roi_overlay(
    records: list[dict],
    *,
    raw_suite2p_folder: Path,
    output_dir: Path,
    subject: str,
    date: str,
    session: str,
) -> int:
    """Plot concatenated FOV images with all Suite2p cell ROIs overlaid.

    The filled version is drawn as one transparent RGBA image layer rather than
    per-pixel vector markers, keeping SVG output compact.  A separate contour
    version is also saved for clean outline inspection.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    payloads: list[dict] = []
    metadata_rows: list[dict] = []

    for record in records:
        fov_name = record["fov_name"]
        plane_id = plane_id_from_fov_name(fov_name, fallback=0)
        plane_folder = raw_suite2p_folder / f"plane{plane_id}"
        required = [plane_folder / "stat.npy", plane_folder / "iscell.npy", plane_folder / "reg_outputs.npy"]
        if not plane_folder.is_dir() or not all(path.exists() for path in required):
            print(f"[WARN] {fov_name}: combined cell ROI overlay skipped for this FOV; missing {plane_folder} or required files")
            continue

        stat = np.load(plane_folder / "stat.npy", allow_pickle=True)
        iscell = np.load(plane_folder / "iscell.npy", allow_pickle=True)
        reg_outputs = load_plane_reg_outputs(plane_folder)
        fov = pick_fov_image_combined(reg_outputs)

        if iscell.ndim < 2 or iscell.shape[0] != len(stat):
            print(f"[WARN] {fov_name}: combined cell ROI overlay skipped for this FOV; iscell/stat length mismatch")
            continue

        cell_indices = np.flatnonzero(iscell[:, 0].astype(bool))
        payloads.append(
            {
                "fov_name": fov_name,
                "plane_id": int(plane_id),
                "stat": stat,
                "fov": fov,
                "cell_indices": cell_indices.astype(int, copy=False),
            }
        )

    if not payloads:
        print("[WARN] Combined Suite2p cell ROI overlay skipped; no FOV images/cell masks could be loaded.")
        return 0

    mosaic, offsets, panel_shapes = concatenate_fov_images(
        [payload["fov"] for payload in payloads],
        gap_pixels=COMBINED_FOV_GAP_PIXELS,
    )
    total_cells = int(sum(payload["cell_indices"].size for payload in payloads))

    for payload, offset_x, shape in zip(payloads, offsets, panel_shapes):
        metadata_rows.append(
            {
                "fov_name": payload["fov_name"],
                "plane_id": int(payload["plane_id"]),
                "n_cell_rois": int(payload["cell_indices"].size),
                "x_offset_pixels": int(offset_x),
                "fov_height_pixels": int(shape[0]),
                "fov_width_pixels": int(shape[1]),
            }
        )

    def render_cell_roi_overlay(
        *,
        mode: str,
        output_name: str,
        mosaic_image: np.ndarray,
        roi_offsets: list,
        shapes: list[tuple[int, int]],
        layout_label: str,
        draw_separators: bool,
    ) -> None:
        fig_width = max(12.0, min(22.0, mosaic_image.shape[1] / max(1, mosaic_image.shape[0]) * 9.0))
        fig, ax = plt.subplots(figsize=(fig_width, 9.0), constrained_layout=True)
        fov_artist = ax.imshow(mosaic_image, cmap="gray", origin="upper", vmin=0.0, vmax=1.0)
        fov_artist.set_rasterized(True)
        ax.set_axis_off()

        if mode == "filled":
            roi_overlay = combined_cell_roi_rgba_overlay(
                payloads,
                roi_offsets,
                mosaic_shape=mosaic_image.shape,
                color="#00ff00",
                alpha=COMBINED_CELL_ROI_FILLED_ALPHA,
            )
            overlay_artist = ax.imshow(roi_overlay, origin="upper", interpolation="nearest", zorder=8)
            overlay_artist.set_rasterized(True)

        for payload, offset, shape in zip(payloads, roi_offsets, shapes):
            offset_x, offset_y = offset_xy(offset)
            fov_name = payload["fov_name"]
            cell_indices = payload["cell_indices"]
            stat = payload["stat"]

            if mode == "contour":
                for roi_idx in cell_indices:
                    draw_single_roi_contour_offset(
                        ax,
                        stat,
                        int(roi_idx),
                        int(offset_x),
                        int(offset_y),
                        "#00ff00",
                        linewidth=COMBINED_CELL_ROI_CONTOUR_LINEWIDTH,
                        zorder=8,
                    )
            elif mode != "filled":
                raise ValueError(f"Unknown combined cell ROI overlay mode: {mode}")

            ax.text(
                float(offset_x + shape[1] / 2.0) if draw_separators else float(offset_x + 8),
                float(offset_y - 5) if draw_separators else float(offset_y + 8),
                f"{fov_name}\n{cell_indices.size} cells",
                ha="center",
                va="bottom" if draw_separators else "top",
                color="black" if draw_separators else "white",
                fontsize=11,
                clip_on=False,
                bbox=None if draw_separators else {"facecolor": "black", "alpha": 0.55, "edgecolor": "none", "pad": 2},
            )
            if draw_separators and offset_y == 0 and offset_x > 0:
                ax.axvline(offset_x - 0.5, color="black", linewidth=1.0, alpha=0.5)

        mode_label = (
            f"filled cell ROIs, raster mask, alpha={COMBINED_CELL_ROI_FILLED_ALPHA:.2f}"
            if mode == "filled"
            else "cell ROI contours only"
        )
        ax.set_title(
            f"{subject} {date} session {session}: combined Suite2p cell ROIs "
            f"({total_cells} cells across {len(payloads)} FOVs) | {layout_label}; {mode_label}",
            fontsize=16,
            pad=12,
        )
        save_figure(fig, output_dir / output_name, dpi=ROI_CONTOUR_OVERLAY_DPI)

    render_cell_roi_overlay(
        mode="filled",
        output_name=f"combined_suite2p_cell_ROIs_on_concatenated_FOVs_filled_raster_alpha{int(round(COMBINED_CELL_ROI_FILLED_ALPHA * 100)):03d}",
        mosaic_image=mosaic,
        roi_offsets=offsets,
        shapes=panel_shapes,
        layout_label="horizontal concatenation",
        draw_separators=True,
    )
    render_cell_roi_overlay(
        mode="contour",
        output_name="combined_suite2p_cell_ROIs_on_concatenated_FOVs_contours",
        mosaic_image=mosaic,
        roi_offsets=offsets,
        shapes=panel_shapes,
        layout_label="horizontal concatenation",
        draw_separators=True,
    )
    write_metadata(metadata_rows, output_dir / "combined_suite2p_cell_roi_overlay_metadata.csv")

    n_written = 2
    try:
        spatial_mosaic, spatial_offsets, spatial_shapes, spatial_metadata_rows = spatial_fov_image_mosaic(
            payloads,
            raw_suite2p_folder=raw_suite2p_folder,
        )
    except Exception as exc:
        print(f"[WARN] Combined spatial Suite2p cell ROI overlays skipped: {exc}")
    else:
        render_cell_roi_overlay(
            mode="filled",
            output_name=f"combined_suite2p_cell_ROIs_on_spatial_FOVs_filled_raster_alpha{int(round(COMBINED_CELL_ROI_FILLED_ALPHA * 100)):03d}",
            mosaic_image=spatial_mosaic,
            roi_offsets=spatial_offsets,
            shapes=spatial_shapes,
            layout_label="spatial FOV layout",
            draw_separators=False,
        )
        render_cell_roi_overlay(
            mode="contour",
            output_name="combined_suite2p_cell_ROIs_on_spatial_FOVs_contours",
            mosaic_image=spatial_mosaic,
            roi_offsets=spatial_offsets,
            shapes=spatial_shapes,
            layout_label="spatial FOV layout",
            draw_separators=False,
        )
        write_metadata(spatial_metadata_rows, output_dir / "combined_suite2p_cell_roi_spatial_overlay_metadata.csv")
        n_written += 2

    return n_written


def plot_combined_rf_roi_overlay(
    records: list[dict],
    *,
    raw_suite2p_folder: Path,
    output_dir: Path,
    subject: str,
    date: str,
    session: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    azimuth_limits, elevation_limits = combined_visual_limits(records)
    bivariate_azimuth_limits, bivariate_elevation_limits = fov_overlay_bivariate_plot_limits(
        azimuth_limits,
        elevation_limits,
    )
    rf_colormap = make_rf_bivariate_colormap(bivariate_azimuth_limits, bivariate_elevation_limits)

    payloads: list[dict] = []
    metadata_rows: list[dict] = []
    all_azimuth: list[np.ndarray] = []
    all_elevation: list[np.ndarray] = []
    all_significant: list[np.ndarray] = []
    all_peak_pass: list[np.ndarray] = []

    for record in records:
        fov_name = record["fov_name"]
        result = record["result"]
        plane_id = plane_id_from_fov_name(fov_name, fallback=0)
        plane_folder = raw_suite2p_folder / f"plane{plane_id}"
        required = [plane_folder / "stat.npy", plane_folder / "reg_outputs.npy"]
        if not plane_folder.is_dir() or not all(path.exists() for path in required):
            print(f"[WARN] {fov_name}: combined RF overlay skipped for this FOV; missing {plane_folder} or required files")
            continue

        stat = np.load(plane_folder / "stat.npy", allow_pickle=True)
        reg_outputs = load_plane_reg_outputs(plane_folder)
        fov = pick_fov_image_combined(reg_outputs)
        summary = rf_summary_vectors(
            result,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
            flip_azimuth=bool(record.get("flip_azimuth", False)),
        )

        payloads.append(
            {
                "fov_name": fov_name,
                "plane_id": plane_id,
                "stat": stat,
                "fov": fov,
                "flip_azimuth": bool(record.get("flip_azimuth", False)),
                **summary,
            }
        )
        all_azimuth.append(summary["azimuth"])
        all_elevation.append(summary["elevation"])
        all_significant.append(summary["significant"])
        all_peak_pass.append(summary["peak_pass"])

    if not payloads:
        print("[WARN] Combined RF/FOV overlay skipped; no FOV images could be loaded.")
        return 0

    mosaic, offsets, panel_shapes = concatenate_fov_images([payload["fov"] for payload in payloads], gap_pixels=COMBINED_FOV_GAP_PIXELS)

    fig_width = max(12.0, min(22.0, mosaic.shape[1] / max(1, mosaic.shape[0]) * 9.0))
    fig = plt.figure(figsize=(fig_width, 8.8), constrained_layout=False)
    ax_fov = fig.add_axes([0.04, 0.06, 0.67, 0.84])
    ax_wheel = fig.add_axes([0.76, 0.40, 0.20, 0.22])

    fov_artist = ax_fov.imshow(mosaic, cmap="gray", origin="upper", vmin=0.0, vmax=1.0)
    fov_artist.set_rasterized(True)
    ax_fov.set_axis_off()

    for payload, offset_x, shape in zip(payloads, offsets, panel_shapes):
        fov_name = payload["fov_name"]
        draw_combined_rf_roi_masks(
            ax=ax_fov,
            stat=payload["stat"],
            x_offset=int(offset_x),
            mosaic_shape=mosaic.shape,
            source_roi_indices=payload["roi_indices"],
            azimuth=payload["azimuth"],
            elevation=payload["elevation"],
            significant=payload["significant"],
            peak_pass=payload["peak_pass"],
            azimuth_limits=bivariate_azimuth_limits,
            elevation_limits=bivariate_elevation_limits,
            rf_colormap=rf_colormap,
        )
        ax_fov.text(
            float(offset_x + shape[1] / 2.0),
            -5,
            fov_name,
            ha="center",
            va="bottom",
            color="black",
            fontsize=11,
            clip_on=False,
        )
        if offset_x > 0:
            ax_fov.axvline(offset_x - 0.5, color="black", linewidth=1.0, alpha=0.5)

        for rf_idx in np.flatnonzero(payload["significant"]):
            roi_idx = int(payload["roi_indices"][rf_idx])
            if roi_idx < 0 or roi_idx >= len(payload["stat"]):
                continue
            roi = stat_entry(payload["stat"][roi_idx])
            y_med, x_med = roi_center(roi)
            metadata_rows.append(
                {
                    "fov_name": fov_name,
                    "plane_id": int(payload["plane_id"]),
                    "rf_result_index_1based": int(rf_idx + 1),
                    "suite2p_roi_index_0based": roi_idx,
                    "suite2p_roi_index_matlab": roi_idx + 1,
                    "x_med_concatenated": float(x_med + offset_x),
                    "x_med_local": float(x_med),
                    "y_med": float(y_med),
                    "azimuth_deg": float(payload["azimuth"][rf_idx]),
                    "elevation_deg": float(payload["elevation"][rf_idx]),
                    "azimuth_flipped": bool(payload.get("flip_azimuth", False)),
                    "explained_variance": float(payload["ev"][rf_idx]),
                    "p_value": float(payload["p_values"][rf_idx]),
                    "peak_to_noise": float(payload["peak_to_noise"][rf_idx]),
                    "passes_peak_to_noise": bool(payload["peak_pass"][rf_idx]),
                }
            )

    if all_azimuth:
        azimuth = np.concatenate(all_azimuth)
        elevation = np.concatenate(all_elevation)
        significant = np.concatenate(all_significant)
        peak_pass = np.concatenate(all_peak_pass)
    else:
        azimuth = elevation = np.asarray([], dtype=float)
        significant = peak_pass = np.asarray([], dtype=bool)

    plot_2d_colorwheel(
        ax_wheel,
        bivariate_azimuth_limits,
        bivariate_elevation_limits,
        rf_colormap=rf_colormap,
        azimuth=azimuth,
        elevation=elevation,
        significant=significant,
        peak_pass=peak_pass,
    )
    ax_fov.set_title(
        f"Combined FOV overlay: concatenated Suite2p images | {int(significant.sum())} significant, {int(peak_pass.sum())} P/N pass",
        fontsize=15,
    )
    fig.suptitle(
        f"{subject} {date} session {session}: all-FOV RF ROI overlay | filled = good RF, open/contour = significant but P/N fail",
        fontsize=17,
    )
    save_figure(fig, output_dir / "combined_significant_ROIs_on_concatenated_FOVs_2D_colorwheel", dpi=220)

    fig_width = max(12.0, min(22.0, mosaic.shape[1] / max(1, mosaic.shape[0]) * 9.0))
    fig = plt.figure(figsize=(fig_width, 8.8), constrained_layout=False)
    ax_fov = fig.add_axes([0.04, 0.06, 0.67, 0.84])
    ax_wheel = fig.add_axes([0.76, 0.40, 0.20, 0.22])

    fov_artist = ax_fov.imshow(mosaic, cmap="gray", origin="upper", vmin=0.0, vmax=1.0)
    fov_artist.set_rasterized(True)
    ax_fov.set_axis_off()

    for payload, offset_x, shape in zip(payloads, offsets, panel_shapes):
        fov_name = payload["fov_name"]
        draw_combined_rf_roi_masks(
            ax=ax_fov,
            stat=payload["stat"],
            x_offset=int(offset_x),
            mosaic_shape=mosaic.shape,
            source_roi_indices=payload["roi_indices"],
            azimuth=payload["azimuth"],
            elevation=payload["elevation"],
            significant=payload["peak_pass"],
            peak_pass=payload["peak_pass"],
            azimuth_limits=bivariate_azimuth_limits,
            elevation_limits=bivariate_elevation_limits,
            rf_colormap=rf_colormap,
        )
        ax_fov.text(
            float(offset_x + shape[1] / 2.0),
            -5,
            fov_name,
            ha="center",
            va="bottom",
            color="black",
            fontsize=11,
            clip_on=False,
        )
        if offset_x > 0:
            ax_fov.axvline(offset_x - 0.5, color="black", linewidth=1.0, alpha=0.5)

    plot_2d_colorwheel(
        ax_wheel,
        bivariate_azimuth_limits,
        bivariate_elevation_limits,
        rf_colormap=rf_colormap,
        azimuth=azimuth,
        elevation=elevation,
        significant=peak_pass,
        peak_pass=peak_pass,
    )
    ax_fov.set_title(
        f"Combined FOV overlay: concatenated Suite2p images | good RFs only (n={int(peak_pass.sum())})",
        fontsize=15,
    )
    fig.suptitle(
        f"{subject} {date} session {session}: all-FOV good RF ROI overlay | pycolorbar={PYCOLORBAR_BIVARIATE_CMAP_NAME}",
        fontsize=17,
    )
    save_figure(fig, output_dir / "combined_good_RFs_on_concatenated_FOVs_2D_colorwheel", dpi=220)

    plot_combined_scalar_rf_roi_overlay(
        payloads=payloads,
        offsets=offsets,
        panel_shapes=panel_shapes,
        mosaic=mosaic,
        value_key="elevation",
        value_limits=good_value_limits(elevation, peak_pass, fallback=elevation_limits),
        value_label="Gaussian elevation (deg)",
        dimension_name="elevation",
        output_base=output_dir / "combined_good_RFs_on_concatenated_FOVs_elevation_jet",
        title=(
            f"{subject} {date} session {session}: combined Good RF ROIs on concatenated Suite2p FOVs | "
            f"elevation color; P/N > {min_peak_to_noise:g}"
        ),
        fov_title=f"Combined FOV overlay: Good RF elevation, jet (n={int(peak_pass.sum())})",
    )
    plot_combined_scalar_rf_roi_overlay(
        payloads=payloads,
        offsets=offsets,
        panel_shapes=panel_shapes,
        mosaic=mosaic,
        value_key="azimuth",
        value_limits=good_value_limits(azimuth, peak_pass, fallback=azimuth_limits),
        value_label="Gaussian azimuth (deg)",
        dimension_name="azimuth",
        output_base=output_dir / "combined_good_RFs_on_concatenated_FOVs_azimuth_jet",
        title=(
            f"{subject} {date} session {session}: combined Good RF ROIs on concatenated Suite2p FOVs | "
            f"azimuth color; P/N > {min_peak_to_noise:g}"
        ),
        fov_title=f"Combined FOV overlay: Good RF azimuth, jet (n={int(peak_pass.sum())})",
    )

    n_written = 4
    try:
        spatial_mosaic, spatial_offsets, spatial_panel_shapes, spatial_layout_rows = spatial_fov_image_mosaic(
            payloads,
            raw_suite2p_folder=raw_suite2p_folder,
        )
    except Exception as exc:
        print(f"[WARN] Combined spatial RF/FOV overlay plots skipped: {exc}")
    else:
        spatial_roi_metadata_rows: list[dict] = []

        fig_width = max(12.0, min(22.0, spatial_mosaic.shape[1] / max(1, spatial_mosaic.shape[0]) * 9.0))
        fig = plt.figure(figsize=(fig_width, 8.8), constrained_layout=False)
        ax_fov = fig.add_axes([0.04, 0.06, 0.67, 0.84])
        ax_wheel = fig.add_axes([0.76, 0.40, 0.20, 0.22])

        fov_artist = ax_fov.imshow(spatial_mosaic, cmap="gray", origin="upper", vmin=0.0, vmax=1.0)
        fov_artist.set_rasterized(True)
        ax_fov.set_axis_off()

        for payload, offset, shape in zip(payloads, spatial_offsets, spatial_panel_shapes):
            offset_x, offset_y = offset_xy(offset)
            fov_name = payload["fov_name"]
            draw_combined_rf_roi_masks(
                ax=ax_fov,
                stat=payload["stat"],
                x_offset=int(offset_x),
                y_offset=int(offset_y),
                mosaic_shape=spatial_mosaic.shape,
                source_roi_indices=payload["roi_indices"],
                azimuth=payload["azimuth"],
                elevation=payload["elevation"],
                significant=payload["significant"],
                peak_pass=payload["peak_pass"],
                azimuth_limits=bivariate_azimuth_limits,
                elevation_limits=bivariate_elevation_limits,
                rf_colormap=rf_colormap,
            )
            ax_fov.text(
                float(offset_x + 8),
                float(offset_y + 8),
                fov_name,
                ha="left",
                va="top",
                color="white",
                fontsize=11,
                clip_on=False,
                bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none", "pad": 2},
            )
            for rf_idx in np.flatnonzero(payload["significant"]):
                roi_idx = int(payload["roi_indices"][rf_idx])
                if roi_idx < 0 or roi_idx >= len(payload["stat"]):
                    continue
                roi = stat_entry(payload["stat"][roi_idx])
                y_med, x_med = roi_center(roi)
                spatial_roi_metadata_rows.append(
                    {
                        "fov_name": fov_name,
                        "plane_id": int(payload["plane_id"]),
                        "rf_result_index_1based": int(rf_idx + 1),
                        "suite2p_roi_index_0based": roi_idx,
                        "suite2p_roi_index_matlab": roi_idx + 1,
                        "x_med_spatial": float(x_med + offset_x),
                        "y_med_spatial": float(y_med + offset_y),
                        "x_med_local": float(x_med),
                        "y_med_local": float(y_med),
                        "azimuth_deg": float(payload["azimuth"][rf_idx]),
                        "elevation_deg": float(payload["elevation"][rf_idx]),
                        "azimuth_flipped": bool(payload.get("flip_azimuth", False)),
                        "explained_variance": float(payload["ev"][rf_idx]),
                        "p_value": float(payload["p_values"][rf_idx]),
                        "peak_to_noise": float(payload["peak_to_noise"][rf_idx]),
                        "passes_peak_to_noise": bool(payload["peak_pass"][rf_idx]),
                    }
                )

        plot_2d_colorwheel(
            ax_wheel,
            bivariate_azimuth_limits,
            bivariate_elevation_limits,
            rf_colormap=rf_colormap,
            azimuth=azimuth,
            elevation=elevation,
            significant=significant,
            peak_pass=peak_pass,
        )
        ax_fov.set_title(
            f"Combined FOV overlay: spatial Suite2p images | {int(significant.sum())} significant, {int(peak_pass.sum())} P/N pass",
            fontsize=15,
        )
        fig.suptitle(
            f"{subject} {date} session {session}: all-FOV RF ROI overlay in spatial layout | filled = good RF, open/contour = significant but P/N fail",
            fontsize=17,
        )
        save_figure(fig, output_dir / "combined_significant_ROIs_on_spatial_FOVs_2D_colorwheel", dpi=220)

        fig_width = max(12.0, min(22.0, spatial_mosaic.shape[1] / max(1, spatial_mosaic.shape[0]) * 9.0))
        fig = plt.figure(figsize=(fig_width, 8.8), constrained_layout=False)
        ax_fov = fig.add_axes([0.04, 0.06, 0.67, 0.84])
        ax_wheel = fig.add_axes([0.76, 0.40, 0.20, 0.22])

        fov_artist = ax_fov.imshow(spatial_mosaic, cmap="gray", origin="upper", vmin=0.0, vmax=1.0)
        fov_artist.set_rasterized(True)
        ax_fov.set_axis_off()

        for payload, offset, shape in zip(payloads, spatial_offsets, spatial_panel_shapes):
            offset_x, offset_y = offset_xy(offset)
            fov_name = payload["fov_name"]
            draw_combined_rf_roi_masks(
                ax=ax_fov,
                stat=payload["stat"],
                x_offset=int(offset_x),
                y_offset=int(offset_y),
                mosaic_shape=spatial_mosaic.shape,
                source_roi_indices=payload["roi_indices"],
                azimuth=payload["azimuth"],
                elevation=payload["elevation"],
                significant=payload["peak_pass"],
                peak_pass=payload["peak_pass"],
                azimuth_limits=bivariate_azimuth_limits,
                elevation_limits=bivariate_elevation_limits,
                rf_colormap=rf_colormap,
            )
            ax_fov.text(
                float(offset_x + 8),
                float(offset_y + 8),
                fov_name,
                ha="left",
                va="top",
                color="white",
                fontsize=11,
                clip_on=False,
                bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none", "pad": 2},
            )

        plot_2d_colorwheel(
            ax_wheel,
            bivariate_azimuth_limits,
            bivariate_elevation_limits,
            rf_colormap=rf_colormap,
            azimuth=azimuth,
            elevation=elevation,
            significant=peak_pass,
            peak_pass=peak_pass,
        )
        ax_fov.set_title(
            f"Combined FOV overlay: spatial Suite2p images | good RFs only (n={int(peak_pass.sum())})",
            fontsize=15,
        )
        fig.suptitle(
            f"{subject} {date} session {session}: all-FOV good RF ROI overlay in spatial layout | pycolorbar={PYCOLORBAR_BIVARIATE_CMAP_NAME}",
            fontsize=17,
        )
        save_figure(fig, output_dir / "combined_good_RFs_on_spatial_FOVs_2D_colorwheel", dpi=220)

        plot_combined_scalar_rf_roi_overlay(
            payloads=payloads,
            offsets=spatial_offsets,
            panel_shapes=spatial_panel_shapes,
            mosaic=spatial_mosaic,
            value_key="elevation",
            value_limits=good_value_limits(elevation, peak_pass, fallback=elevation_limits),
            value_label="Gaussian elevation (deg)",
            dimension_name="elevation",
            output_base=output_dir / "combined_good_RFs_on_spatial_FOVs_elevation_jet",
            title=(
                f"{subject} {date} session {session}: combined Good RF ROIs on spatial Suite2p FOVs | "
                f"elevation color; P/N > {min_peak_to_noise:g}"
            ),
            fov_title=f"Combined FOV overlay: spatial Good RF elevation, jet (n={int(peak_pass.sum())})",
            draw_separators=False,
        )
        plot_combined_scalar_rf_roi_overlay(
            payloads=payloads,
            offsets=spatial_offsets,
            panel_shapes=spatial_panel_shapes,
            mosaic=spatial_mosaic,
            value_key="azimuth",
            value_limits=good_value_limits(azimuth, peak_pass, fallback=azimuth_limits),
            value_label="Gaussian azimuth (deg)",
            dimension_name="azimuth",
            output_base=output_dir / "combined_good_RFs_on_spatial_FOVs_azimuth_jet",
            title=(
                f"{subject} {date} session {session}: combined Good RF ROIs on spatial Suite2p FOVs | "
                f"azimuth color; P/N > {min_peak_to_noise:g}"
            ),
            fov_title=f"Combined FOV overlay: spatial Good RF azimuth, jet (n={int(peak_pass.sum())})",
            draw_separators=False,
        )
        write_metadata(spatial_layout_rows, output_dir / "combined_fov_spatial_layout_metadata.csv")
        write_metadata(spatial_roi_metadata_rows, output_dir / "combined_significant_roi_fov_spatial_overlay_metadata.csv")
        n_written += 4

    write_metadata(metadata_rows, output_dir / "combined_significant_roi_fov_overlay_metadata.csv")
    return n_written


def combined_cell_roi_rgba_overlay(
    payloads: list[dict],
    offsets: list,
    *,
    mosaic_shape: tuple[int, int],
    color: str,
    alpha: float,
) -> np.ndarray:
    """Return one RGBA image containing all filled cell ROIs.

    This keeps filled cell-ROI SVGs compact by saving the ROI layer as a single
    transparent raster image rather than hundreds of thousands of vector
    marker elements.
    """

    height, width = int(mosaic_shape[0]), int(mosaic_shape[1])
    overlay = np.zeros((height, width, 4), dtype=np.float32)
    rgb = np.asarray(matplotlib.colors.to_rgb(color), dtype=np.float32)
    alpha = float(np.clip(alpha, 0.0, 1.0))

    for payload, offset in zip(payloads, offsets):
        offset_x, offset_y = offset_xy(offset)
        stat = payload["stat"]
        for roi_idx in np.asarray(payload["cell_indices"], dtype=int):
            if roi_idx < 0 or roi_idx >= len(stat):
                continue
            roi = stat_entry(stat[int(roi_idx)])
            ypix = np.asarray(roi["ypix"], dtype=int) + int(offset_y)
            xpix = np.asarray(roi["xpix"], dtype=int) + int(offset_x)
            valid = (ypix >= 0) & (ypix < height) & (xpix >= 0) & (xpix < width)
            if not np.any(valid):
                continue
            overlay[ypix[valid], xpix[valid], :3] = rgb
            overlay[ypix[valid], xpix[valid], 3] = np.maximum(overlay[ypix[valid], xpix[valid], 3], alpha)

    return overlay


def pick_fov_image_combined(reg_outputs: dict) -> np.ndarray:
    for key in ("meanImgE", "meanImg", "refImg"):
        if key in reg_outputs:
            image = np.asarray(reg_outputs[key], dtype=float)
            if image.ndim == 2:
                return normalize_fov_black_floor(image)
    raise KeyError("Could not find meanImgE, meanImg, or refImg in reg_outputs.npy")


def normalize_fov_black_floor(image: np.ndarray) -> np.ndarray:
    finite = np.isfinite(image)
    if not finite.any():
        return np.zeros_like(image, dtype=float)
    lo = float(np.nanpercentile(image[finite], FOV_BLACK_FLOOR_PERCENTILE))
    hi = float(np.nanpercentile(image[finite], 99.8))
    if not np.isfinite(hi) or hi <= lo:
        hi = float(np.nanmax(image[finite]))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(image, dtype=float)
    return np.clip((image - lo) / (hi - lo), 0.0, 1.0)


def concatenate_fov_images(images: list[np.ndarray], *, gap_pixels: int) -> tuple[np.ndarray, list[int], list[tuple[int, int]]]:
    shapes = [(int(image.shape[0]), int(image.shape[1])) for image in images]
    max_height = max(height for height, _ in shapes)
    total_width = sum(width for _, width in shapes) + gap_pixels * max(0, len(images) - 1)
    mosaic = np.zeros((max_height, total_width), dtype=float)
    offsets: list[int] = []
    x0 = 0
    for image, (height, width) in zip(images, shapes):
        offsets.append(x0)
        mosaic[:height, x0 : x0 + width] = image
        x0 += width + gap_pixels
    return mosaic, offsets, shapes


def spatial_fov_image_mosaic(
    payloads: list[dict],
    *,
    raw_suite2p_folder: Path,
) -> tuple[np.ndarray, list[tuple[int, int]], list[tuple[int, int]], list[dict]]:
    combined_stat_path = raw_suite2p_folder / "combined" / "stat.npy"
    if not combined_stat_path.exists():
        raise FileNotFoundError(f"Spatial FOV layout requires combined Suite2p stat.npy: {combined_stat_path}")

    combined_stat = np.load(combined_stat_path, allow_pickle=True)
    combined_x, combined_y, combined_iplane = stat_centers_and_iplanes(combined_stat)
    if not combined_x.size or not np.isfinite(combined_iplane).any():
        raise ValueError(f"Spatial FOV layout requires finite iplane coordinates in {combined_stat_path}")

    pieces: list[dict] = []
    metadata_rows: list[dict] = []
    for payload in payloads:
        image = np.asarray(payload["fov"], dtype=float)
        height, width = int(image.shape[0]), int(image.shape[1])
        plane_id = int(payload["plane_id"])
        plane_x, plane_y, _ = stat_centers_and_iplanes(payload["stat"])
        mask = np.isfinite(combined_iplane) & (combined_iplane == float(plane_id))
        this_combined_x = combined_x[mask]
        this_combined_y = combined_y[mask]
        n = min(this_combined_x.size, plane_x.size)
        if n < 10:
            raise ValueError(
                f"Spatial FOV layout for {payload['fov_name']} needs at least 10 matched combined/stat ROIs; found {n}"
            )

        ddx = this_combined_x[:n] - plane_x[:n]
        ddy = this_combined_y[:n] - plane_y[:n]
        ddx = ddx[np.isfinite(ddx)]
        ddy = ddy[np.isfinite(ddy)]
        if not ddx.size or not ddy.size:
            raise ValueError(f"Spatial FOV layout for {payload['fov_name']} has no finite combined/stat offsets")
        dx = float(np.nanmedian(ddx))
        dy = float(np.nanmedian(ddy))

        pieces.append(
            {
                "payload": payload,
                "image": image,
                "x": float(dx),
                "y": float(dy),
                "width": width,
                "height": height,
                "source": "combined_stat_iplane",
            }
        )

    min_x = min(0.0, min(piece["x"] for piece in pieces))
    min_y = min(0.0, min(piece["y"] for piece in pieces))
    shift_x = -min_x
    shift_y = -min_y
    mosaic_width = int(np.ceil(max(piece["x"] + shift_x + piece["width"] for piece in pieces)))
    mosaic_height = int(np.ceil(max(piece["y"] + shift_y + piece["height"] for piece in pieces)))
    mosaic = np.zeros((mosaic_height, mosaic_width), dtype=float)
    offsets: list[tuple[int, int]] = []
    shapes: list[tuple[int, int]] = []

    for piece in pieces:
        x0 = int(round(piece["x"] + shift_x))
        y0 = int(round(piece["y"] + shift_y))
        image = piece["image"]
        height, width = int(piece["height"]), int(piece["width"])
        mosaic[y0 : y0 + height, x0 : x0 + width] = np.maximum(
            mosaic[y0 : y0 + height, x0 : x0 + width],
            image,
        )
        offsets.append((x0, y0))
        shapes.append((height, width))
        metadata_rows.append(
            {
                "fov_name": piece["payload"]["fov_name"],
                "plane_id": int(piece["payload"]["plane_id"]),
                "x_offset_pixels": x0,
                "y_offset_pixels": y0,
                "raw_x_offset_pixels": float(piece["x"]),
                "raw_y_offset_pixels": float(piece["y"]),
                "fov_height_pixels": height,
                "fov_width_pixels": width,
                "offset_source": piece["source"],
            }
        )

    return mosaic, offsets, shapes, metadata_rows


def stat_centers_and_iplanes(stat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    iplanes: list[float] = []
    for entry in np.asarray(stat, dtype=object).ravel():
        try:
            roi = stat_entry(entry)
        except Exception:
            xs.append(np.nan)
            ys.append(np.nan)
            iplanes.append(np.nan)
            continue

        try:
            y_med, x_med = roi_center(roi)
        except Exception:
            y_med, x_med = np.nan, np.nan
        try:
            iplane = float(np.asarray(roi.get("iplane", np.nan), dtype=float).squeeze())
        except Exception:
            iplane = np.nan
        xs.append(float(x_med))
        ys.append(float(y_med))
        iplanes.append(float(iplane))
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float), np.asarray(iplanes, dtype=float)


def offset_xy(offset) -> tuple[int, int]:
    if isinstance(offset, (tuple, list, np.ndarray)):
        if len(offset) >= 2:
            return int(offset[0]), int(offset[1])
        if len(offset) == 1:
            return int(offset[0]), 0
    return int(offset), 0


def draw_combined_rf_roi_masks(
    *,
    ax,
    stat: np.ndarray,
    x_offset: int,
    y_offset: int = 0,
    mosaic_shape: tuple[int, int],
    source_roi_indices: np.ndarray,
    azimuth: np.ndarray,
    elevation: np.ndarray,
    significant: np.ndarray,
    peak_pass: np.ndarray,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    rf_colormap,
) -> None:
    non_peak = significant & ~peak_pass
    for rf_idx in np.flatnonzero(non_peak):
        draw_single_roi_contour_offset(
            ax,
            stat,
            int(source_roi_indices[rf_idx]),
            int(x_offset),
            int(y_offset),
            rf_position_to_rgb(azimuth[rf_idx], elevation[rf_idx], azimuth_limits, elevation_limits, rf_colormap=rf_colormap),
            linewidth=NON_PEAK_CONTOUR_LINEWIDTH,
            zorder=7,
        )
    roi_overlay = rf_roi_rgba_overlay(
        stat=stat,
        source_roi_indices=source_roi_indices,
        azimuth=azimuth,
        elevation=elevation,
        use_mask=peak_pass,
        azimuth_limits=azimuth_limits,
        elevation_limits=elevation_limits,
        rf_colormap=rf_colormap,
        image_shape=mosaic_shape,
        x_offset=int(x_offset),
        y_offset=int(y_offset),
        alpha=PEAK_ROI_ALPHA,
    )
    if np.any(roi_overlay[:, :, 3] > 0):
        overlay_artist = ax.imshow(roi_overlay, origin="upper", interpolation="nearest", zorder=9)
        overlay_artist.set_rasterized(True)


def plot_combined_scalar_rf_roi_overlay(
    *,
    payloads: list[dict],
    offsets: list,
    panel_shapes: list[tuple[int, int]],
    mosaic: np.ndarray,
    value_key: str,
    value_limits: tuple[float, float],
    value_label: str,
    dimension_name: str,
    output_base: Path,
    title: str,
    fov_title: str,
    draw_separators: bool = True,
) -> None:
    fig_width = max(12.0, min(22.0, mosaic.shape[1] / max(1, mosaic.shape[0]) * 9.0))
    fig = plt.figure(figsize=(fig_width, 8.8), constrained_layout=False)
    ax_fov = fig.add_axes([0.04, 0.06, 0.70, 0.84])
    ax_cbar = fig.add_axes([0.82, 0.24, 0.035, 0.48])

    fov_artist = ax_fov.imshow(mosaic, cmap="gray", origin="upper", vmin=0.0, vmax=1.0)
    fov_artist.set_rasterized(True)
    ax_fov.set_axis_off()

    roi_overlay = combined_scalar_roi_rgba_overlay(
        payloads=payloads,
        offsets=offsets,
        value_key=value_key,
        value_limits=value_limits,
        mosaic_shape=mosaic.shape,
        alpha=PEAK_ROI_ALPHA,
    )
    if np.any(roi_overlay[:, :, 3] > 0):
        overlay_artist = ax_fov.imshow(roi_overlay, origin="upper", interpolation="nearest", zorder=9)
        overlay_artist.set_rasterized(True)

    for payload, offset, shape in zip(payloads, offsets, panel_shapes):
        offset_x, offset_y = offset_xy(offset)
        fov_name = payload["fov_name"]
        ax_fov.text(
            float(offset_x + shape[1] / 2.0) if draw_separators else float(offset_x + 8),
            float(offset_y - 5) if draw_separators else float(offset_y + 8),
            fov_name,
            ha="center" if draw_separators else "left",
            va="bottom" if draw_separators else "top",
            color="black" if draw_separators else "white",
            fontsize=11,
            clip_on=False,
            bbox=None if draw_separators else {"facecolor": "black", "alpha": 0.55, "edgecolor": "none", "pad": 2},
        )
        if draw_separators and offset_y == 0 and offset_x > 0:
            ax_fov.axvline(offset_x - 0.5, color="black", linewidth=1.0, alpha=0.5)

    norm = scalar_value_norm(value_limits)
    colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.get_cmap("jet")), cax=ax_cbar)
    colorbar.set_label(value_label, fontsize=11)
    colorbar.ax.tick_params(labelsize=9)
    colorbar.ax.set_title(dimension_name.title(), fontsize=11, pad=7)

    ax_fov.set_title(fov_title, fontsize=15)
    fig.suptitle(title, fontsize=17)
    save_figure(fig, output_base, dpi=220)


def combined_scalar_roi_rgba_overlay(
    *,
    payloads: list[dict],
    offsets: list,
    value_key: str,
    value_limits: tuple[float, float],
    mosaic_shape: tuple[int, int],
    alpha: float,
) -> np.ndarray:
    height, width = int(mosaic_shape[0]), int(mosaic_shape[1])
    overlay = np.zeros((height, width, 4), dtype=np.float32)

    for payload, offset in zip(payloads, offsets):
        offset_x, offset_y = offset_xy(offset)
        fov_overlay = scalar_roi_rgba_overlay(
            stat=payload["stat"],
            source_roi_indices=payload["roi_indices"],
            values=payload[value_key],
            use_mask=payload["peak_pass"],
            value_limits=value_limits,
            image_shape=mosaic_shape,
            x_offset=int(offset_x),
            y_offset=int(offset_y),
            alpha=alpha,
        )
        use_pixels = fov_overlay[:, :, 3] > 0
        if np.any(use_pixels):
            overlay[use_pixels, :] = fov_overlay[use_pixels, :]
    return overlay


def draw_single_roi_contour_offset(
    ax,
    stat: np.ndarray,
    roi_idx: int,
    x_offset: int,
    y_offset_or_color,
    color=None,
    *,
    linewidth: float,
    zorder: int,
) -> None:
    if color is None:
        y_offset = 0
        color = y_offset_or_color
    else:
        y_offset = int(y_offset_or_color)
    if roi_idx >= len(stat) or roi_idx < 0:
        return
    roi = stat_entry(stat[roi_idx])
    ypix = np.asarray(roi["ypix"], dtype=int) + int(y_offset)
    xpix = np.asarray(roi["xpix"], dtype=int) + int(x_offset)
    if ypix.size == 0 or xpix.size == 0:
        return
    y0, y1 = int(ypix.min()) - 1, int(ypix.max()) + 1
    x0, x1 = int(xpix.min()) - 1, int(xpix.max()) + 1
    local_mask = np.zeros((y1 - y0 + 1, x1 - x0 + 1), dtype=float)
    local_mask[ypix - y0, xpix - x0] = 1.0
    contour_set = ax.contour(
        np.arange(x0, x1 + 1),
        np.arange(y0, y1 + 1),
        local_mask,
        levels=[0.5],
        colors=[color],
        linewidths=linewidth,
        alpha=1.0,
        zorder=zorder,
    )
    rasterize_contour_set(contour_set)


def draw_single_roi_mask_offset(
    ax,
    stat: np.ndarray,
    roi_idx: int,
    x_offset: int,
    color,
    *,
    alpha: float,
    zorder: int,
    y_offset: int = 0,
) -> None:
    if roi_idx >= len(stat) or roi_idx < 0:
        return
    roi = stat_entry(stat[roi_idx])
    ypix = np.asarray(roi["ypix"], dtype=int) + int(y_offset)
    xpix = np.asarray(roi["xpix"], dtype=int) + int(x_offset)
    artist = ax.scatter(xpix, ypix, s=ROI_MARKER_SIZE, c=[color], marker="s", linewidths=0, alpha=alpha, zorder=zorder)
    artist.set_rasterized(True)


def plot_scalar_rf_roi_overlay(
    *,
    fov: np.ndarray,
    stat: np.ndarray,
    source_roi_indices: np.ndarray,
    values: np.ndarray,
    significant: np.ndarray,
    peak_pass: np.ndarray,
    value_limits: tuple[float, float],
    value_label: str,
    dimension_name: str,
    output_base: Path,
    title: str,
    fov_title: str,
) -> None:
    fig = plt.figure(figsize=(11.5, 8.5), constrained_layout=False)
    ax_fov = fig.add_axes([0.06, 0.07, 0.66, 0.82])
    ax_cbar = fig.add_axes([0.79, 0.24, 0.035, 0.48])

    fov_artist = ax_fov.imshow(fov, cmap="gray", origin="upper")
    fov_artist.set_rasterized(True)
    draw_rf_roi_masks_by_scalar(
        ax=ax_fov,
        stat=stat,
        image_shape=fov.shape,
        source_roi_indices=source_roi_indices,
        values=values,
        significant=significant,
        peak_pass=peak_pass,
        value_limits=value_limits,
    )
    ax_fov.set_title(fov_title, fontsize=15)
    ax_fov.set_axis_off()

    norm = scalar_value_norm(value_limits)
    colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.get_cmap("jet")), cax=ax_cbar)
    colorbar.set_label(value_label, fontsize=11)
    colorbar.ax.tick_params(labelsize=9)
    colorbar.ax.set_title(dimension_name.title(), fontsize=11, pad=7)
    add_scalar_colorbar_reference_dots(
        colorbar.ax,
        values=values,
        use_mask=peak_pass,
        value_limits=value_limits,
    )

    fig.suptitle(title, fontsize=17)
    save_figure(fig, output_base, dpi=220)


def draw_rf_roi_masks_by_scalar(
    *,
    ax,
    stat: np.ndarray,
    image_shape: tuple[int, int],
    source_roi_indices: np.ndarray,
    values: np.ndarray,
    significant: np.ndarray,
    peak_pass: np.ndarray,
    value_limits: tuple[float, float],
) -> None:
    roi_overlay = scalar_roi_rgba_overlay(
        stat=stat,
        source_roi_indices=source_roi_indices,
        values=values,
        use_mask=peak_pass,
        value_limits=value_limits,
        image_shape=image_shape,
        x_offset=0,
        alpha=PEAK_ROI_ALPHA,
    )
    if np.any(roi_overlay[:, :, 3] > 0):
        overlay_artist = ax.imshow(roi_overlay, origin="upper", interpolation="nearest", zorder=9)
        overlay_artist.set_rasterized(True)



def draw_rf_roi_masks(
    *,
    ax,
    stat: np.ndarray,
    image_shape: tuple[int, int],
    source_roi_indices: np.ndarray,
    azimuth: np.ndarray,
    elevation: np.ndarray,
    significant: np.ndarray,
    peak_pass: np.ndarray,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    rf_colormap,
) -> None:
    non_peak = significant & ~peak_pass
    for rf_idx in np.flatnonzero(non_peak):
        draw_single_roi_contour(
            ax,
            stat,
            int(source_roi_indices[rf_idx]),
            rf_position_to_rgb(azimuth[rf_idx], elevation[rf_idx], azimuth_limits, elevation_limits, rf_colormap=rf_colormap),
            linewidth=NON_PEAK_CONTOUR_LINEWIDTH,
            zorder=7,
        )
    roi_overlay = rf_roi_rgba_overlay(
        stat=stat,
        source_roi_indices=source_roi_indices,
        azimuth=azimuth,
        elevation=elevation,
        use_mask=peak_pass,
        azimuth_limits=azimuth_limits,
        elevation_limits=elevation_limits,
        rf_colormap=rf_colormap,
        image_shape=image_shape,
        x_offset=0,
        alpha=PEAK_ROI_ALPHA,
    )
    if np.any(roi_overlay[:, :, 3] > 0):
        overlay_artist = ax.imshow(roi_overlay, origin="upper", interpolation="nearest", zorder=9)
        overlay_artist.set_rasterized(True)


def draw_single_roi_contour(ax, stat: np.ndarray, roi_idx: int, color, *, linewidth: float, zorder: int) -> None:
    if roi_idx >= len(stat) or roi_idx < 0:
        return
    roi = stat_entry(stat[roi_idx])
    ypix = np.asarray(roi["ypix"], dtype=int)
    xpix = np.asarray(roi["xpix"], dtype=int)
    if ypix.size == 0 or xpix.size == 0:
        return
    y0, y1 = int(ypix.min()) - 1, int(ypix.max()) + 1
    x0, x1 = int(xpix.min()) - 1, int(xpix.max()) + 1
    local_mask = np.zeros((y1 - y0 + 1, x1 - x0 + 1), dtype=float)
    local_mask[ypix - y0, xpix - x0] = 1.0
    contour_set = ax.contour(
        np.arange(x0, x1 + 1),
        np.arange(y0, y1 + 1),
        local_mask,
        levels=[0.5],
        colors=[color],
        linewidths=linewidth,
        alpha=1.0,
        zorder=zorder,
    )
    rasterize_contour_set(contour_set)


def draw_single_roi_mask(ax, stat: np.ndarray, roi_idx: int, color, *, alpha: float, zorder: int) -> None:
    if roi_idx >= len(stat) or roi_idx < 0:
        return
    roi = stat_entry(stat[roi_idx])
    ypix = np.asarray(roi["ypix"], dtype=int)
    xpix = np.asarray(roi["xpix"], dtype=int)
    artist = ax.scatter(xpix, ypix, s=ROI_MARKER_SIZE, c=[color], marker="s", linewidths=0, alpha=alpha, zorder=zorder)
    artist.set_rasterized(True)


def plot_2d_colorwheel(
    ax,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    *,
    rf_colormap,
    azimuth: np.ndarray | None = None,
    elevation: np.ndarray | None = None,
    significant: np.ndarray | None = None,
    peak_pass: np.ndarray | None = None,
    force_equal_aspect: bool = False,
    draw_zero_lines: bool = True,
) -> None:
    azimuth_axis = np.linspace(azimuth_limits[0], azimuth_limits[1], BIVARIATE_CMAP_WHEEL_RESOLUTION)
    elevation_axis = np.linspace(elevation_limits[0], elevation_limits[1], BIVARIATE_CMAP_WHEEL_RESOLUTION)
    azimuth_grid, elevation_grid = np.meshgrid(azimuth_axis, elevation_axis)
    colorwheel = rf_position_to_rgb(
        azimuth_grid,
        elevation_grid,
        azimuth_limits,
        elevation_limits,
        rf_colormap=rf_colormap,
    )
    ax.imshow(
        colorwheel,
        origin="lower",
        extent=(azimuth_limits[0], azimuth_limits[1], elevation_limits[0], elevation_limits[1]),
        aspect="equal" if force_equal_aspect else "auto",
    )
    if draw_zero_lines:
        ax.axvline(0, color="black", linewidth=0.8, alpha=0.45)
        ax.axhline(0, color="black", linewidth=0.8, alpha=0.45)

    if azimuth is not None and elevation is not None and significant is not None and peak_pass is not None:
        azimuth = np.asarray(azimuth, dtype=float)
        elevation = np.asarray(elevation, dtype=float)
        significant = np.asarray(significant, dtype=bool)
        peak_pass = np.asarray(peak_pass, dtype=bool)

        finite_centers = np.isfinite(azimuth) & np.isfinite(elevation)
        in_bounds = (
            (azimuth >= azimuth_limits[0])
            & (azimuth <= azimuth_limits[1])
            & (elevation >= elevation_limits[0])
            & (elevation <= elevation_limits[1])
        )
        good = significant & peak_pass & finite_centers & in_bounds
        significant_non_good = significant & ~peak_pass & finite_centers & in_bounds

        if np.any(significant_non_good):
            ax.scatter(
                azimuth[significant_non_good],
                elevation[significant_non_good],
                s=22,
                marker="o",
                facecolors="none",
                edgecolors="black",
                linewidths=0.9,
                alpha=0.95,
                zorder=6,
            )
        if np.any(good):
            ax.scatter(
                azimuth[good],
                elevation[good],
                s=18,
                marker="o",
                c="black",
                edgecolors="black",
                linewidths=0.4,
                alpha=0.95,
                zorder=7,
            )
    ax.set_xlim(azimuth_limits)
    ax.set_ylim(elevation_limits)
    ax.set_xlabel("Gaussian azimuth (deg)", fontsize=9, labelpad=2)
    ax.set_ylabel("Gaussian elevation (deg)", fontsize=9, labelpad=2)
    ax.set_title(f"pycolorbar {PYCOLORBAR_BIVARIATE_CMAP_NAME}", fontsize=11, pad=5)
    ax.tick_params(labelsize=8, pad=2)
    if force_equal_aspect:
        ax.set_aspect("equal", adjustable="box")


def fov_overlay_bivariate_plot_limits(
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
) -> tuple[tuple[float, float], tuple[float, float]]:
    plot_azimuth_limits = tuple(sorted((float(azimuth_limits[0]), float(azimuth_limits[1]))))
    if FOV_OVERLAY_BIVARIATE_AZIMUTH_LIMITS is not None:
        configured = tuple(sorted(float(value) for value in FOV_OVERLAY_BIVARIATE_AZIMUTH_LIMITS))
        if len(configured) != 2 or not all(np.isfinite(configured)) or configured[0] == configured[1]:
            raise ValueError(
                "FOV_OVERLAY_BIVARIATE_AZIMUTH_LIMITS must be None or two finite, distinct azimuth values"
            )
        plot_azimuth_limits = configured

    plot_elevation_limits = tuple(sorted((float(elevation_limits[0]), float(elevation_limits[1]))))
    return plot_azimuth_limits, plot_elevation_limits


def make_rf_bivariate_colormap(
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
):
    if BivariateColormap is None:
        raise ImportError(
            "pycolorbar is required for RF-center FOV overlay colors. "
            "Install it with: pip install pycolorbar"
        )
    return BivariateColormap.from_name(
        PYCOLORBAR_BIVARIATE_CMAP_NAME,
        n=BIVARIATE_CMAP_WHEEL_RESOLUTION,
        diagonal_tilt=PYCOLORBAR_TEULING_DIAGONAL_TILT,
        offdiag_tilt=PYCOLORBAR_TEULING_OFFDIAG_TILT,
    )


def sample_rf_bivariate_colormap(
    rf_colormap,
    x: np.ndarray,
    y: np.ndarray,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
) -> np.ndarray:
    norm_x = matplotlib.colors.Normalize(vmin=float(azimuth_limits[0]), vmax=float(azimuth_limits[1]), clip=True)
    norm_y = matplotlib.colors.Normalize(vmin=float(elevation_limits[0]), vmax=float(elevation_limits[1]), clip=True)
    color = rf_colormap(np.asarray(x, dtype=float), np.asarray(y, dtype=float), norm_x=norm_x, norm_y=norm_y)
    color = np.asarray(color, dtype=float).reshape(-1)
    if np.nanmax(color) > 1.0:
        color = color / 255.0
    if color.size and color.size % 4 == 0:
        color = color.reshape(np.asarray(x).shape + (4,))[..., :3]
    elif color.size and color.size % 3 == 0:
        color = color.reshape(np.asarray(x).shape + (3,))
    else:
        raise ValueError("pycolorbar returned a color array with an unexpected shape")
    return np.clip(color, 0.0, 1.0)


def rf_position_to_rgb(
    azimuth: np.ndarray | float,
    elevation: np.ndarray | float,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    *,
    rf_colormap,
) -> np.ndarray:
    """Map RF centers to RGB using pycolorbar.

    The ranges are the visual-field azimuth/elevation limits, so the same
    coordinates used in the Gaussian RF fits are passed directly to the 2D
    colormap.
    """
    azimuth = np.asarray(azimuth, dtype=float)
    elevation = np.asarray(elevation, dtype=float)
    az_b, el_b = np.broadcast_arrays(azimuth, elevation)
    rgb = sample_rf_bivariate_colormap(rf_colormap, az_b, el_b, azimuth_limits, elevation_limits)
    finite = np.isfinite(az_b) & np.isfinite(el_b)
    if not np.all(finite):
        rgb = np.asarray(rgb, dtype=float).copy()
        rgb[~finite, :] = np.array([0.5, 0.5, 0.5], dtype=float)
    return rgb


def scalar_value_norm(value_limits: tuple[float, float]) -> matplotlib.colors.Normalize:
    lo, hi = sorted((float(value_limits[0]), float(value_limits[1])))
    if not np.isfinite(lo) or not np.isfinite(hi):
        lo, hi = 0.0, 1.0
    if lo == hi:
        pad = 1.0 if lo == 0.0 else abs(lo) * 0.05
        lo -= pad
        hi += pad
    return matplotlib.colors.Normalize(vmin=lo, vmax=hi, clip=True)


def good_value_limits(values: np.ndarray, good_mask: np.ndarray, *, fallback: tuple[float, float]) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    good_mask = np.asarray(good_mask, dtype=bool)
    good_values = values[good_mask & np.isfinite(values)]
    if good_values.size == 0:
        return tuple(sorted((float(fallback[0]), float(fallback[1]))))
    return (float(np.nanmin(good_values)), float(np.nanmax(good_values)))


def add_scalar_colorbar_reference_dots(
    ax,
    *,
    values: np.ndarray,
    use_mask: np.ndarray,
    value_limits: tuple[float, float],
) -> None:
    values = np.asarray(values, dtype=float)
    use_mask = np.asarray(use_mask, dtype=bool)
    lo, hi = sorted((float(value_limits[0]), float(value_limits[1])))
    show_values = values[use_mask & np.isfinite(values) & (values >= lo) & (values <= hi)]
    if show_values.size == 0:
        return
    transform = blended_transform_factory(ax.transAxes, ax.transData)
    ax.scatter(
        np.full(show_values.shape, 0.5, dtype=float),
        show_values,
        s=18,
        marker="o",
        facecolors="black",
        edgecolors="white",
        linewidths=0.35,
        alpha=0.8,
        transform=transform,
        clip_on=True,
        zorder=10,
    )


def scalar_value_to_rgb(
    value: np.ndarray | float,
    value_limits: tuple[float, float],
    *,
    cmap_name: str = "jet",
) -> np.ndarray:
    values = np.asarray(value, dtype=float)
    rgb = np.empty(values.shape + (3,), dtype=float)
    flat_values = values.ravel()
    flat_rgb = rgb.reshape(-1, 3)
    norm = scalar_value_norm(value_limits)
    cmap = plt.get_cmap(cmap_name)
    for i, this_value in enumerate(flat_values):
        if np.isfinite(this_value):
            flat_rgb[i, :] = np.asarray(cmap(norm(float(this_value)))[:3], dtype=float)
        else:
            flat_rgb[i, :] = np.array([0.5, 0.5, 0.5], dtype=float)
    return rgb


def rasterize_contour_set(contour_set) -> None:
    if hasattr(contour_set, "set_rasterized"):
        contour_set.set_rasterized(True)
    collections = getattr(contour_set, "collections", [])
    for collection in collections:
        collection.set_rasterized(True)


def roi_rgba_overlay(
    *,
    stat: np.ndarray,
    roi_indices: np.ndarray,
    image_shape: tuple[int, int],
    color: str,
    alpha: float,
    x_offset: int = 0,
    y_offset: int = 0,
) -> np.ndarray:
    height, width = int(image_shape[0]), int(image_shape[1])
    overlay = np.zeros((height, width, 4), dtype=np.float32)
    rgb = np.asarray(matplotlib.colors.to_rgb(color), dtype=np.float32)
    alpha = float(np.clip(alpha, 0.0, 1.0))

    for roi_idx in np.asarray(roi_indices, dtype=int):
        if roi_idx < 0 or roi_idx >= len(stat):
            continue
        roi = stat_entry(stat[int(roi_idx)])
        ypix = np.asarray(roi["ypix"], dtype=int) + int(y_offset)
        xpix = np.asarray(roi["xpix"], dtype=int) + int(x_offset)
        valid = (ypix >= 0) & (ypix < height) & (xpix >= 0) & (xpix < width)
        if not np.any(valid):
            continue
        overlay[ypix[valid], xpix[valid], :3] = rgb
        overlay[ypix[valid], xpix[valid], 3] = np.maximum(overlay[ypix[valid], xpix[valid], 3], alpha)
    return overlay


def rf_roi_rgba_overlay(
    *,
    stat: np.ndarray,
    source_roi_indices: np.ndarray,
    azimuth: np.ndarray,
    elevation: np.ndarray,
    use_mask: np.ndarray,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    rf_colormap,
    image_shape: tuple[int, int],
    x_offset: int,
    alpha: float,
    y_offset: int = 0,
) -> np.ndarray:
    height, width = int(image_shape[0]), int(image_shape[1])
    overlay = np.zeros((height, width, 4), dtype=np.float32)
    alpha = float(np.clip(alpha, 0.0, 1.0))

    for rf_idx in np.flatnonzero(np.asarray(use_mask, dtype=bool)):
        if rf_idx >= len(source_roi_indices):
            continue
        roi_idx = int(source_roi_indices[rf_idx])
        if roi_idx < 0 or roi_idx >= len(stat):
            continue
        roi = stat_entry(stat[roi_idx])
        ypix = np.asarray(roi["ypix"], dtype=int) + int(y_offset)
        xpix = np.asarray(roi["xpix"], dtype=int) + int(x_offset)
        valid = (ypix >= 0) & (ypix < height) & (xpix >= 0) & (xpix < width)
        if not np.any(valid):
            continue
        rgb = np.asarray(
            rf_position_to_rgb(
                azimuth[rf_idx],
                elevation[rf_idx],
                azimuth_limits,
                elevation_limits,
                rf_colormap=rf_colormap,
            ),
            dtype=np.float32,
        ).reshape(3)
        overlay[ypix[valid], xpix[valid], :3] = rgb
        overlay[ypix[valid], xpix[valid], 3] = np.maximum(overlay[ypix[valid], xpix[valid], 3], alpha)
    return overlay


def scalar_roi_rgba_overlay(
    *,
    stat: np.ndarray,
    source_roi_indices: np.ndarray,
    values: np.ndarray,
    use_mask: np.ndarray,
    value_limits: tuple[float, float],
    image_shape: tuple[int, int],
    x_offset: int,
    alpha: float,
    y_offset: int = 0,
) -> np.ndarray:
    height, width = int(image_shape[0]), int(image_shape[1])
    overlay = np.zeros((height, width, 4), dtype=np.float32)
    alpha = float(np.clip(alpha, 0.0, 1.0))

    for rf_idx in np.flatnonzero(np.asarray(use_mask, dtype=bool)):
        if rf_idx >= len(source_roi_indices):
            continue
        roi_idx = int(source_roi_indices[rf_idx])
        if roi_idx < 0 or roi_idx >= len(stat):
            continue
        roi = stat_entry(stat[roi_idx])
        ypix = np.asarray(roi["ypix"], dtype=int) + int(y_offset)
        xpix = np.asarray(roi["xpix"], dtype=int) + int(x_offset)
        valid = (ypix >= 0) & (ypix < height) & (xpix >= 0) & (xpix < width)
        if not np.any(valid):
            continue
        rgb = np.asarray(scalar_value_to_rgb(values[rf_idx], value_limits), dtype=np.float32).reshape(3)
        overlay[ypix[valid], xpix[valid], :3] = rgb
        overlay[ypix[valid], xpix[valid], 3] = np.maximum(overlay[ypix[valid], xpix[valid], 3], alpha)
    return overlay


def plot_rois(*, fov: np.ndarray, stat: np.ndarray, roi_indices: np.ndarray, title: str, output_base: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 11), constrained_layout=True)
    fov_artist = ax.imshow(fov, cmap="gray", origin="upper")
    fov_artist.set_rasterized(True)
    roi_overlay = roi_rgba_overlay(
        stat=stat,
        roi_indices=roi_indices,
        image_shape=fov.shape,
        color="#00ff00",
        alpha=NEUTRAL_ROI_ALPHA,
    )
    overlay_artist = ax.imshow(roi_overlay, origin="upper", interpolation="nearest", zorder=8)
    overlay_artist.set_rasterized(True)
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_axis_off()
    save_figure(fig, output_base, dpi=ROI_CONTOUR_OVERLAY_DPI)


# =============================================================================
# File/path utilities
# =============================================================================


def discover_fovs(alf_session: Path, results_session: Path) -> list[str]:
    names = set()
    if alf_session.is_dir():
        names.update(path.name for path in alf_session.iterdir() if path.is_dir() and path.name.upper().startswith("FOV_"))
    if results_session.is_dir():
        names.update(
            path.name
            for path in results_session.iterdir()
            if path.is_dir() and path.name.upper().startswith("FOV_") and (path / "_FANCi_rf.maps.npy").exists()
        )
    return sorted(names)


def find_raw_suite2p_folder(*, subject: str, date: str, session: str) -> Path:
    """Find the selected session's Suite2p folder without scanning the full NAS."""

    subject_root = REMOTE_2P_ROOT / subject
    date_label = normalize_date_label(date)
    compact_date = date_label.replace("-", "")
    date_roots = [
        subject_root / compact_date,
        subject_root / date_label,
        subject_root / "Processed" / compact_date,
        subject_root / "Processed" / date_label,
    ]

    search_roots: list[Path] = []
    for date_root in date_roots:
        if not date_root.is_dir():
            continue
        search_roots.extend([date_root / str(session), date_root / "Processed", date_root / "processed"])
        for child in date_root.iterdir():
            if child.is_dir() and session_name_matches(child.name, str(session)):
                search_roots.extend([child, child / "suite2p", child / "Processed", child / "processed"])
        search_roots.append(date_root)

    run_folders: list[Path] = []
    seen = set()
    for root in unique_paths(search_roots):
        for run_folder in suite2p_run_folders_under(root):
            key = str(run_folder).lower()
            if key not in seen:
                run_folders.append(run_folder)
                seen.add(key)

    if not run_folders:
        raise FileNotFoundError(
            "Could not find reg_outputs.npy for "
            f"{subject} {date_label} session {session} under {subject_root}"
        )

    run_folders.sort(key=lambda path: selected_run_preference_key(path, session=str(session)))
    if len(run_folders) > 1:
        top_rank = selected_run_preference_key(run_folders[0], session=str(session))[:3]
        tied = [path for path in run_folders if selected_run_preference_key(path, session=str(session))[:3] == top_rank]
        if len(tied) > 1:
            details = "\n  ".join(str(path) for path in tied)
            raise FileNotFoundError(
                "Found multiple equally plausible Suite2p run folders for "
                f"{subject} {date_label} session {session}; refusing to guess:\n  {details}"
            )
    return run_folders[0]


def unique_paths(paths: list[Path]) -> list[Path]:
    seen = set()
    out = []
    for path in paths:
        key = str(path).lower()
        if key not in seen:
            out.append(path)
            seen.add(key)
    return out


def suite2p_run_folders_under(root: Path) -> list[Path]:
    if not root.is_dir():
        return []
    reg_paths = [
        path
        for path in root.rglob("reg_outputs.npy")
        if re.fullmatch(r"plane\d+", path.parent.name, flags=re.IGNORECASE)
    ]
    return sorted({path.parent.parent for path in reg_paths}, key=run_folder_preference_key)


def session_name_matches(name: str, session: str) -> bool:
    if name == session:
        return True
    if not re.fullmatch(r"\d+(?:_\d+)*", name):
        return False
    tokens = [token for token in re.split(r"\D+", name) if token]
    return session in tokens


def selected_run_preference_key(path: Path, *, session: str) -> tuple[int, int, int, str]:
    parts = [part.lower() for part in path.parts]
    processed_rank = 0 if any(part == "processed" for part in parts) else 1
    exact_session_rank = 2
    for part in parts:
        if part == session:
            exact_session_rank = min(exact_session_rank, 0)
        elif session_name_matches(part, session):
            exact_session_rank = min(exact_session_rank, 1)
    return (exact_session_rank, processed_rank, len(path.parts), str(path).lower())


def normalize_date_label(value: str) -> str:
    text = str(value).strip()
    match = re.fullmatch(r"((?:19|20)\d{2})-?(\d{2})-?(\d{2})", text)
    if not match:
        return text
    return f"{match.group(1)}-{match.group(2)}-{match.group(3)}"


def normalize_fov_name(value: str) -> str:
    text = str(value).strip()
    match = re.fullmatch(r"(?i)FOV[_-]?(\d+)", text)
    if match:
        return f"FOV_{int(match.group(1)):02d}"
    return text


def run_folder_preference_key(path: Path) -> tuple[int, int, str]:
    parts = [part.lower() for part in path.parts]
    has_processed = any(part == "processed" for part in parts)
    return (0 if has_processed else 1, len(path.parts), str(path).lower())


def load_plane_reg_outputs(plane_folder: Path) -> dict:
    reg_path = plane_folder / "reg_outputs.npy"
    if not reg_path.exists():
        raise FileNotFoundError(f"Missing per-plane registration output: {reg_path}")
    return np.load(reg_path, allow_pickle=True).item()


def pick_fov_image(reg_outputs: dict) -> np.ndarray:
    for key in ("meanImgE", "meanImg", "refImg"):
        if key in reg_outputs:
            image = np.asarray(reg_outputs[key], dtype=float)
            if image.ndim == 2:
                return normalize_fov_black_floor(image)
    raise KeyError("Could not find meanImgE, meanImg, or refImg in reg_outputs.npy")


def normalize_fov(image: np.ndarray) -> np.ndarray:
    finite = np.isfinite(image)
    if not finite.any():
        return np.zeros_like(image, dtype=float)
    lo, hi = np.percentile(image[finite], [1, 99.8])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = np.nanmin(image), np.nanmax(image)
    if hi <= lo:
        return np.zeros_like(image, dtype=float)
    return np.clip((image - lo) / (hi - lo), 0, 1)


def plane_id_from_fov_name(name: str, *, fallback: int) -> int:
    text = str(name)
    upper = text.upper()
    digits = "".join(ch for ch in text if ch.isdigit())
    if digits:
        if upper.startswith("FOV_"):
            return int(digits)
        return max(0, int(digits) - 1)
    return int(fallback)


def get_edges(result: FanciResult) -> np.ndarray:
    edges = np.ravel(np.asarray(result.get("edges", np.asarray([-135.0, 135.0, 40.0, -40.0])), dtype=float))
    if edges.size < 4:
        return np.asarray([-135.0, 135.0, 40.0, -40.0], dtype=float)
    return edges[:4]


def rf_image_extent(edges: np.ndarray) -> tuple[float, float, float, float]:
    """Return Matplotlib extent for saved RF maps.

    The RF MATLAB code saves edges as [left, right, top, bottom].
    Matplotlib imshow wants (left, right, bottom, top).  With origin="upper",
    row 0 of the 7 x 27 RF map is drawn at the top/elevated side of the
    screen, matching the Gaussian fit and the MATLAB fitting grid.
    """
    edges = np.ravel(np.asarray(edges, dtype=float))
    return (float(edges[0]), float(edges[1]), float(edges[3]), float(edges[2]))


def rf_elevation_limits(edges: np.ndarray) -> tuple[float, float]:
    edges = np.ravel(np.asarray(edges, dtype=float))
    return (float(min(edges[2], edges[3])), float(max(edges[2], edges[3])))


def bin_centers(start: float, stop: float, n: int) -> np.ndarray:
    e = np.linspace(float(start), float(stop), n + 1)
    return (e[:-1] + e[1:]) / 2.0


def source_roi_value(result: FanciResult, cell_id: int) -> int:
    arr = np.ravel(np.asarray(result.get("roiIndex0", np.arange(result.require("maps").shape[0])), dtype=int))
    if cell_id < arr.size:
        return int(arr[cell_id])
    return int(cell_id)


def source_roi_matlab_value(result: FanciResult, cell_id: int) -> int:
    arr = np.ravel(np.asarray(result.get("roiIndexMatlab", np.arange(result.require("maps").shape[0]) + 1), dtype=int))
    if cell_id < arr.size:
        return int(arr[cell_id])
    return int(cell_id + 1)


def source_roi_text(result: FanciResult, cell_id: int) -> str:
    try:
        return str(source_roi_value(result, cell_id))
    except Exception:
        return "NA"


def stat_entry(entry) -> dict:
    if isinstance(entry, dict):
        return entry
    if hasattr(entry, "item"):
        return entry.item()
    raise TypeError(f"Unsupported stat.npy entry type: {type(entry)!r}")


def roi_center(roi: dict) -> tuple[float, float]:
    if "med" in roi and len(roi["med"]) >= 2:
        return float(roi["med"][0]), float(roi["med"][1])
    return float(np.mean(roi["ypix"])), float(np.mean(roi["xpix"]))


def write_metadata(rows: list[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        output_path.write_text("", encoding="utf-8")
        return
    with output_path.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.DictWriter(file_obj, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_fov_summary(result: FanciResult, output_path: Path) -> None:
    maps = result.require("maps")
    n = maps.shape[0]
    ev = result_column(result, "explVars", n)
    p = result_column(result, "pValues", n)
    pn = result_column(result, "peakToNoise", n)
    ev_shift, ev_label = ev_for_shift_plots(result)
    data = {
        "results_dir": str(result.results_dir),
        "mat_path": str(result.mat_path) if result.mat_path is not None else None,
        "n_rois": int(n),
        "maps_shape": list(map(int, maps.shape)),
        "n_p_lt_0_05": int(np.sum(np.isfinite(p) & (p < 0.05))),
        "n_ev_gt_min": int(np.sum(np.isfinite(ev) & (ev > MIN_EV))),
        "n_peak_noise_gt_min": int(np.sum(np.isfinite(pn) & (pn > MIN_PEAK_TO_NOISE))),
        "ev_for_shift_label": ev_label,
        "ev_for_shift_min": float(np.nanmin(ev_shift)) if np.isfinite(ev_shift).any() else None,
        "ev_for_shift_max": float(np.nanmax(ev_shift)) if np.isfinite(ev_shift).any() else None,
        "loaded_arrays": sorted(result.arrays.keys()),
    }
    output_path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def session_aggregate_rows(
    records: list[dict],
    *,
    group: str,
    subject: str,
    date: str,
    session: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
) -> list[dict]:
    rows: list[dict] = []
    session_label = f"{subject} {date} s{session}"
    session_total_neurons = int(sum(record["result"].require("maps").shape[0] for record in records))
    for record in records:
        result = record["result"]
        fov_name = record["fov_name"]
        flip_azimuth = bool(record.get("flip_azimuth", False))
        summary = rf_summary_vectors(
            result,
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
            flip_azimuth=flip_azimuth,
        )
        gauss = summary["gauss"]
        n = gauss.shape[0]
        fov_total_neurons = int(result.require("maps").shape[0])
        if gauss.shape[1] < 6:
            continue
        shift_ev, shift_ev_label = ev_for_shift_plots(result)
        shift_ev = result_column_from_array(shift_ev, n)
        roi_index0 = np.ravel(np.asarray(result.get("roiIndex0", np.arange(n)), dtype=float))
        edges = get_edges(result)
        azimuth_limits = tuple(np.sort(edges[:2]).astype(float))
        elevation_limits = tuple(np.sort(edges[2:4]).astype(float))
        valid_width = np.isfinite(gauss[:, 2]) & np.isfinite(gauss[:, 4]) & (np.abs(gauss[:, 2]) > 0) & (np.abs(gauss[:, 4]) > 0)
        use = summary["significant"] & valid_width

        for rf_idx in np.flatnonzero(use):
            source_roi0 = roi_index0[rf_idx] if rf_idx < roi_index0.size else np.nan
            is_good = bool(summary["peak_pass"][rf_idx])
            rows.append(
                {
                    "group": group,
                    "session_label": session_label,
                    "subject": subject,
                    "date": date,
                    "session": session,
                    "fov": fov_name,
                    "fov_total_neurons": fov_total_neurons,
                    "session_total_neurons": session_total_neurons,
                    "rf_result_index_0based": int(rf_idx),
                    "rf_result_index_1based": int(rf_idx + 1),
                    "source_roi0": int(source_roi0) if np.isfinite(source_roi0) else "",
                    "azimuth_deg": float(gauss[rf_idx, 1]),
                    "elevation_deg": float(gauss[rf_idx, 3]),
                    "sigma_x_deg": float(gauss[rf_idx, 2]),
                    "sigma_y_deg": float(gauss[rf_idx, 4]),
                    "theta_rad": float(gauss[rf_idx, 5]),
                    "rf_width_deg": float((abs(gauss[rf_idx, 2]) + abs(gauss[rf_idx, 4])) / 2.0),
                    "explained_variance": float(summary["ev"][rf_idx]) if np.isfinite(summary["ev"][rf_idx]) else np.nan,
                    "shift_ev": float(shift_ev[rf_idx]) if np.isfinite(shift_ev[rf_idx]) else np.nan,
                    "shift_ev_label": shift_ev_label,
                    "p_value": float(summary["p_values"][rf_idx]) if np.isfinite(summary["p_values"][rf_idx]) else np.nan,
                    "peak_to_noise": float(summary["peak_to_noise"][rf_idx]) if np.isfinite(summary["peak_to_noise"][rf_idx]) else np.nan,
                    "is_valid": True,
                    "is_good": is_good,
                    "is_valid_not_good": bool(not is_good),
                    "azimuth_flipped": flip_azimuth,
                    "azimuth_min": float(azimuth_limits[0]),
                    "azimuth_max": float(azimuth_limits[1]),
                    "elevation_min": float(elevation_limits[0]),
                    "elevation_max": float(elevation_limits[1]),
                }
            )
    return rows


def session_ev_aggregate_rows(
    records: list[dict],
    *,
    group: str,
    subject: str,
    date: str,
    session: str,
    max_p_value: float,
) -> list[dict]:
    rows: list[dict] = []
    session_label = f"{subject} {date} s{session}"
    for record in records:
        result = record["result"]
        fov_name = record["fov_name"]
        n = result.require("maps").shape[0]
        p_values = result_column(result, "pValues", n)
        ev_values, ev_label = ev_for_shift_plots(result)
        ev_values = result_column_from_array(ev_values, n)
        for rf_idx, (ev_value, p_value) in enumerate(zip(ev_values, p_values)):
            finite_ev = np.isfinite(ev_value)
            finite_p = np.isfinite(p_value)
            rows.append(
                {
                    "group": group,
                    "session_label": session_label,
                    "subject": subject,
                    "date": date,
                    "session": session,
                    "fov": fov_name,
                    "rf_result_index_0based": int(rf_idx),
                    "rf_result_index_1based": int(rf_idx + 1),
                    "ev_value": float(ev_value) if finite_ev else np.nan,
                    "ev_label": ev_label,
                    "p_value": float(p_value) if finite_p else np.nan,
                    "p_lt_threshold": bool(finite_p and p_value < max_p_value),
                }
            )
    return rows


def result_column_from_array(values: np.ndarray, n_expected: int) -> np.ndarray:
    out = np.full(n_expected, np.nan, dtype=float)
    incoming = np.ravel(np.asarray(values, dtype=float))
    n = min(n_expected, incoming.size)
    out[:n] = incoming[:n]
    return out


def plot_group_aggregates(
    group: str,
    rows: list[dict],
    ev_rows: list[dict],
    output_dir: Path,
    *,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_metadata(rows, output_dir / f"{group}_aggregate_RF_values.csv")
    write_metadata(ev_rows, output_dir / f"{group}_aggregate_session_EV_values.csv")
    plot_aggregate_center_map(group, rows, output_dir, good_only=False, min_peak_to_noise=min_peak_to_noise)
    plot_aggregate_center_map(group, rows, output_dir, good_only=True, min_peak_to_noise=min_peak_to_noise)
    plot_aggregate_gaussian_density(group, rows, output_dir, min_ev=min_ev, max_p_value=max_p_value, min_peak_to_noise=min_peak_to_noise)
    plot_aggregate_vertical_center_density_summary(
        group,
        rows,
        output_dir,
        min_ev=min_ev,
        max_p_value=max_p_value,
        min_peak_to_noise=min_peak_to_noise,
    )
    plot_aggregate_vertical_center_density_contour_summary(
        group,
        rows,
        output_dir,
        min_ev=min_ev,
        max_p_value=max_p_value,
        min_peak_to_noise=min_peak_to_noise,
    )
    plot_aggregate_ev_diagnostics(group, ev_rows, output_dir, max_p_value=max_p_value)
    plot_aggregate_rf_widths(group, rows, output_dir, min_peak_to_noise=min_peak_to_noise)
    print(f"[DONE] {group} aggregate plots: {output_dir}")



def plot_subject_debug_aggregates(
    sessions: list[dict[str, Any]],
    group_rows: dict[str, list[dict]],
    plots_root: Path,
    *,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
) -> None:
    """Make the group vertical summary plot again after aggregating within each subject.

    Outputs are written to:
        <plots_root>/<group>/<subject>/debug/
    """

    sessions_by_group_subject: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for info in sessions:
        key = (str(info.get("group", "")), str(info.get("subject", "")))
        sessions_by_group_subject.setdefault(key, []).append(info)

    for group, rows in sorted(group_rows.items()):
        rows_by_subject: dict[str, list[dict]] = {}
        for row in rows:
            rows_by_subject.setdefault(str(row.get("subject", "")), []).append(row)

        for subject, subject_rows in sorted(rows_by_subject.items()):
            if not subject:
                continue
            output_dir = Path(plots_root) / group / subject / "debug"
            output_dir.mkdir(parents=True, exist_ok=True)
            write_metadata(subject_rows, output_dir / f"{subject}_aggregate_RF_values.csv")

            flip_summary = format_subject_flip_summary(
                sessions_by_group_subject.get((group, subject), [])
            )
            made = plot_aggregate_vertical_center_density_summary(
                subject,
                subject_rows,
                output_dir,
                min_ev=min_ev,
                max_p_value=max_p_value,
                min_peak_to_noise=min_peak_to_noise,
                display_label=f"{group} / {subject}",
                title_extra=flip_summary,
                output_stem=f"{subject}_aggregate_good_RF_center_density_vertical_summary",
            )
            if made:
                print(f"[DONE] {group}/{subject} subject debug aggregate plots: {output_dir}")
            else:
                print(f"[INFO] {group}/{subject} subject debug aggregate skipped: no good RF rows")


def format_subject_flip_summary(sessions: list[dict[str, Any]]) -> str:
    """Return a compact title line listing flipped FOVs for every session."""

    if not sessions:
        return "Flipped FOVs by session: none found for this subject"

    entries: list[str] = []
    for info in sorted(
        sessions,
        key=lambda item: (
            normalize_date_label(str(item.get("date", ""))),
            str(item.get("session", "")),
        ),
    ):
        date = normalize_date_label(str(info.get("date", "")))
        session = str(info.get("session", ""))
        flips = sorted(normalize_fov_name(name) for name in set(info.get("flip_fovs", set())))
        flipped_text = ", ".join(flips) if flips else "none"
        entries.append(f"{date} s{session}: {flipped_text}")

    return "Flipped FOVs by session: " + "; ".join(entries)


def wrap_title_text(text: str, *, width: int = 125) -> str:
    """Wrap long figure-title annotations without losing explicit newlines."""

    wrapped_lines: list[str] = []
    for line in str(text).splitlines():
        wrapped = textwrap.wrap(line, width=width, break_long_words=False, break_on_hyphens=False)
        wrapped_lines.extend(wrapped if wrapped else [""])
    return "\n".join(wrapped_lines)

def aggregate_visual_limits(rows: list[dict]) -> tuple[tuple[float, float], tuple[float, float]]:
    if not rows:
        return (-135.0, 135.0), (-40.0, 40.0)
    return (
        (float(np.nanmin([row["azimuth_min"] for row in rows])), float(np.nanmax([row["azimuth_max"] for row in rows]))),
        (float(np.nanmin([row["elevation_min"] for row in rows])), float(np.nanmax([row["elevation_max"] for row in rows]))),
    )


def plot_aggregate_center_map(
    group: str,
    rows: list[dict],
    output_dir: Path,
    *,
    good_only: bool,
    min_peak_to_noise: float,
) -> int:
    use_rows = [row for row in rows if row["is_good"] or not good_only]
    azimuth_limits, elevation_limits = aggregate_visual_limits(rows)
    rf_colormap = make_rf_bivariate_colormap(azimuth_limits, elevation_limits)
    azimuth = np.asarray([row["azimuth_deg"] for row in use_rows], dtype=float)
    elevation = np.asarray([row["elevation_deg"] for row in use_rows], dtype=float)
    significant = np.ones(len(use_rows), dtype=bool)
    peak_pass = np.asarray([row["is_good"] for row in use_rows], dtype=bool)
    if good_only:
        peak_pass[:] = True

    fig, ax = plt.subplots(figsize=(7.5, 6.5), constrained_layout=True)
    plot_2d_colorwheel(
        ax,
        azimuth_limits,
        elevation_limits,
        rf_colormap=rf_colormap,
        azimuth=azimuth,
        elevation=elevation,
        significant=significant,
        peak_pass=peak_pass,
        force_equal_aspect=True,
        draw_zero_lines=False,
    )
    title_suffix = f"good RF centers only (n={len(use_rows)})" if good_only else f"valid and good RF centers (n={len(use_rows)}, good={int(peak_pass.sum())})"
    ax.set_title(f"{group}: aggregate {title_suffix}\nGood RF: P/N > {min_peak_to_noise:g}")
    name = f"{group}_aggregate_RF_centers_good_only" if good_only else f"{group}_aggregate_RF_centers_valid_and_good"
    save_figure(fig, output_dir / name, dpi=300)
    return 1


def plot_aggregate_gaussian_density(
    group: str,
    rows: list[dict],
    output_dir: Path,
    *,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    n_azimuth_grid: int = 500,
    n_elevation_grid: int = 250,
    gaussian_sigma_radius: float = 2.0,
) -> int:
    valid_rows = [row for row in rows if row["is_valid"]]
    good_ids = np.asarray([idx for idx, row in enumerate(valid_rows) if row["is_good"]], dtype=int)
    all_ids = np.arange(len(valid_rows), dtype=int)
    gauss = np.asarray(
        [[1.0, row["azimuth_deg"], row["sigma_x_deg"], row["elevation_deg"], row["sigma_y_deg"], row["theta_rad"]] for row in valid_rows],
        dtype=float,
    )
    azimuth_limits, elevation_limits = aggregate_visual_limits(rows)
    azimuth_axis = np.linspace(azimuth_limits[0], azimuth_limits[1], n_azimuth_grid)
    elevation_axis = np.linspace(elevation_limits[0], elevation_limits[1], n_elevation_grid)
    azimuth_grid, elevation_grid = np.meshgrid(azimuth_axis, elevation_axis)
    good_coverage = gaussian_coverage_map(gauss, good_ids, azimuth_grid, elevation_grid, gaussian_sigma_radius)
    valid_coverage = gaussian_coverage_map(gauss, all_ids, azimuth_grid, elevation_grid, gaussian_sigma_radius)
    good_vmax = coverage_vmax(good_coverage)
    valid_vmax = coverage_vmax(valid_coverage)

    # Clean density plot: keep the original selected_RFplot_Final.py output unchanged.
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), constrained_layout=True)
    plot_density_panel(axes[0], good_coverage, azimuth_axis, elevation_axis, azimuth_limits, elevation_limits, good_vmax, f"Good RFs (n={good_ids.size})")
    plot_density_panel(axes[1], valid_coverage, azimuth_axis, elevation_axis, azimuth_limits, elevation_limits, valid_vmax, f"All valid RFs (n={all_ids.size})")
    fig.suptitle(
        f"{group}: aggregate Gaussian RF density | EV > {min_ev:g}, p < {max_p_value:g}, good P/N > {min_peak_to_noise:g}",
        fontweight="bold",
    )
    save_figure(fig, output_dir / f"{group}_aggregate_gaussian_density", dpi=300)

    # Second density plot: same heatmaps, with only good-fit Gaussian ellipses overlaid.
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), constrained_layout=True)
    plot_density_panel(
        axes[0],
        good_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        good_vmax,
        f"Good RFs (n={good_ids.size}) with good-fit Gaussian ellipses",
    )
    plot_density_panel(
        axes[1],
        valid_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        valid_vmax,
        f"All valid RFs (n={all_ids.size}) with good-fit Gaussian ellipses",
    )
    if good_ids.size:
        for ax in axes:
            plot_gaussian_ellipses(
                ax,
                gauss,
                good_ids,
                gaussian_sigma_radius,
                color="#ffffff",
                linestyle="-",
                linewidth=1.1,
                azimuth_limits=azimuth_limits,
                elevation_limits=elevation_limits,
            )
    fig.suptitle(
        f"{group}: aggregate Gaussian RF density + good Gaussian ellipses | "
        f"EV > {min_ev:g}, p < {max_p_value:g}, good P/N > {min_peak_to_noise:g}",
        fontweight="bold",
    )
    save_figure(fig, output_dir / f"{group}_aggregate_gaussian_density_with_good_ellipses", dpi=300)
    return 2



def plot_aggregate_vertical_center_density_summary(
    group: str,
    rows: list[dict],
    output_dir: Path,
    *,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    n_azimuth_grid: int = 500,
    n_elevation_grid: int = 250,
    gaussian_sigma_radius: float = 2.0,
    display_label: str | None = None,
    title_extra: str | None = None,
    output_stem: str | None = None,
) -> int:
    """Stack center color maps and good-only density maps without legends/colorbars."""

    valid_rows = [row for row in rows if row["is_valid"]]
    good_rows = [row for row in valid_rows if row["is_good"]]
    if not good_rows:
        return 0

    azimuth_limits, elevation_limits = aggregate_visual_limits(rows)
    rf_colormap = make_rf_bivariate_colormap(azimuth_limits, elevation_limits)

    gauss = np.asarray(
        [
            [
                1.0,
                row["azimuth_deg"],
                row["sigma_x_deg"],
                row["elevation_deg"],
                row["sigma_y_deg"],
                row["theta_rad"],
            ]
            for row in good_rows
        ],
        dtype=float,
    )
    good_ids = np.arange(len(good_rows), dtype=int)

    azimuth_axis = np.linspace(azimuth_limits[0], azimuth_limits[1], n_azimuth_grid)
    elevation_axis = np.linspace(elevation_limits[0], elevation_limits[1], n_elevation_grid)
    azimuth_grid, elevation_grid = np.meshgrid(azimuth_axis, elevation_axis)
    good_coverage = gaussian_coverage_map(gauss, good_ids, azimuth_grid, elevation_grid, gaussian_sigma_radius)
    good_vmax = coverage_vmax(good_coverage)

    azimuth = np.asarray([row["azimuth_deg"] for row in good_rows], dtype=float)
    elevation = np.asarray([row["elevation_deg"] for row in good_rows], dtype=float)
    significant = np.ones(len(good_rows), dtype=bool)
    peak_pass = np.ones(len(good_rows), dtype=bool)

    fig, axes = plt.subplots(
        4,
        1,
        figsize=(7.6, 10.2),
        constrained_layout=False,
        sharex=True,
    )
    ax_center_clean, ax_density_clean, ax_center_ellipses, ax_density_ellipses = axes

    for ax in (ax_center_clean, ax_center_ellipses):
        plot_2d_colorwheel(
            ax,
            azimuth_limits,
            elevation_limits,
            rf_colormap=rf_colormap,
            azimuth=azimuth,
            elevation=elevation,
            significant=significant,
            peak_pass=peak_pass,
            force_equal_aspect=True,
            draw_zero_lines=False,
        )

    plot_density_panel_no_colorbar(
        ax_density_clean,
        good_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        good_vmax,
        f"Good RF Gaussian density (n={len(good_rows)})",
    )
    plot_density_panel_no_colorbar(
        ax_density_ellipses,
        good_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        good_vmax,
        f"Good RF Gaussian density + good ellipses (n={len(good_rows)})",
    )
    plot_gaussian_ellipses(
        ax_density_ellipses,
        gauss,
        good_ids,
        gaussian_sigma_radius,
        color="#ffffff",
        linestyle="-",
        linewidth=1.1,
        azimuth_limits=azimuth_limits,
        elevation_limits=elevation_limits,
    )

    ax_center_clean.set_title(f"Good RF center color map (n={len(good_rows)})", pad=4)
    ax_center_ellipses.set_title(f"Good RF center color map (n={len(good_rows)})", pad=4)

    for ax in axes:
        ax.set_xlim(azimuth_limits)
        ax.set_ylim(elevation_limits)
        ax.set_aspect("equal", adjustable="box")
        ax.set_ylabel("Elevation (deg)")

    # Keep the stack compact: only the bottom panel carries the x-axis label.
    for ax in axes[:-1]:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)
    axes[-1].set_xlabel("Azimuth (deg)")

    label = display_label if display_label is not None else group
    title_lines = [
        (
            f"{label}: aggregate good-RF center maps and Gaussian densities | "
            f"EV > {min_ev:g}, p < {max_p_value:g}, good P/N > {min_peak_to_noise:g}"
        )
    ]
    if title_extra:
        title_lines.append(wrap_title_text(title_extra))
    fig.suptitle(
        "\n".join(title_lines),
        fontweight="bold",
        y=0.992,
        fontsize=10 if title_extra else None,
    )
    top_margin = 0.895 if title_extra else 0.94
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.055, top=top_margin, hspace=0.24)
    stem = output_stem if output_stem is not None else f"{group}_aggregate_good_RF_center_density_vertical_summary"
    save_figure(fig, output_dir / stem, dpi=300)
    return 1



def plot_aggregate_vertical_center_density_contour_summary(
    group: str,
    rows: list[dict],
    output_dir: Path,
    *,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    n_azimuth_grid: int = 500,
    n_elevation_grid: int = 250,
    gaussian_sigma_radius: float = 2.0,
) -> int:
    """Duplicate the vertical group summary with the FOV-overlay azimuth crop and contour panels.

    This intentionally leaves the original
    <group>_aggregate_good_RF_center_density_vertical_summary output unchanged.
    The duplicate uses FOV_OVERLAY_BIVARIATE_AZIMUTH_LIMITS exactly through
    fov_overlay_bivariate_plot_limits(), matching the 2D overlay colorwheel logic.
    """

    valid_rows = [row for row in rows if row["is_valid"]]
    good_rows = [row for row in valid_rows if row["is_good"]]
    if not good_rows:
        return 0

    full_azimuth_limits, full_elevation_limits = aggregate_visual_limits(rows)
    azimuth_limits, elevation_limits = fov_overlay_bivariate_plot_limits(
        full_azimuth_limits,
        full_elevation_limits,
    )
    rf_colormap = make_rf_bivariate_colormap(azimuth_limits, elevation_limits)

    gauss = np.asarray(
        [
            [
                1.0,
                row["azimuth_deg"],
                row["sigma_x_deg"],
                row["elevation_deg"],
                row["sigma_y_deg"],
                row["theta_rad"],
            ]
            for row in good_rows
        ],
        dtype=float,
    )
    good_ids = np.arange(len(good_rows), dtype=int)

    # Build the density grid directly on the cropped azimuth domain so the
    # heatmap, contour, and contourf panels all use the same display extent.
    azimuth_axis = np.linspace(azimuth_limits[0], azimuth_limits[1], n_azimuth_grid)
    elevation_axis = np.linspace(elevation_limits[0], elevation_limits[1], n_elevation_grid)
    azimuth_grid, elevation_grid = np.meshgrid(azimuth_axis, elevation_axis)
    good_coverage = gaussian_coverage_map(gauss, good_ids, azimuth_grid, elevation_grid, gaussian_sigma_radius)
    good_vmax = coverage_vmax(good_coverage)

    azimuth = np.asarray([row["azimuth_deg"] for row in good_rows], dtype=float)
    elevation = np.asarray([row["elevation_deg"] for row in good_rows], dtype=float)
    significant = np.ones(len(good_rows), dtype=bool)
    peak_pass = np.ones(len(good_rows), dtype=bool)

    fig, axes = plt.subplots(
        7,
        1,
        figsize=(7.6, 17.6),
        constrained_layout=False,
        sharex=True,
    )
    (
        ax_center_clean,
        ax_density_clean,
        ax_center_contours,
        ax_density_contours,
        ax_density_contourf,
        ax_center_90pct,
        ax_density_90pct,
    ) = axes

    for ax in (ax_center_clean, ax_center_contours, ax_center_90pct):
        plot_2d_colorwheel(
            ax,
            azimuth_limits,
            elevation_limits,
            rf_colormap=rf_colormap,
            azimuth=azimuth,
            elevation=elevation,
            significant=significant,
            peak_pass=peak_pass,
            force_equal_aspect=True,
            draw_zero_lines=False,
        )

    plot_density_panel_no_colorbar(
        ax_density_clean,
        good_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        good_vmax,
        f"Good RF Gaussian density (n={len(good_rows)})",
    )
    plot_density_contour_panel_no_colorbar(
        ax_density_contours,
        good_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        good_vmax,
        f"Good RF Gaussian density contours (n={len(good_rows)})",
        filled=False,
    )
    plot_density_contour_panel_no_colorbar(
        ax_density_contourf,
        good_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        good_vmax,
        f"Good RF Gaussian density filled contours (n={len(good_rows)})",
        filled=True,
    )
    plot_density_panel_no_colorbar(
        ax_density_90pct,
        good_coverage,
        azimuth_axis,
        elevation_axis,
        azimuth_limits,
        elevation_limits,
        good_vmax,
        f"Good RF Gaussian density + 66% density contour (n={len(good_rows)})",
    )
    plot_density_mass_contour(
        ax_density_90pct,
        good_coverage,
        azimuth_axis,
        elevation_axis,
        mass_fraction=0.666,
        color="red",
        linewidth=1.8,
    )

    ax_center_clean.set_title(f"Good RF center color map (n={len(good_rows)})", pad=4)
    ax_center_contours.set_title(f"Good RF center color map (n={len(good_rows)})", pad=4)
    ax_center_90pct.set_title(f"Good RF center color map (n={len(good_rows)})", pad=4)

    for ax in axes:
        ax.set_xlim(azimuth_limits)
        ax.set_ylim(elevation_limits)
        ax.set_aspect("equal", adjustable="box")
        ax.set_ylabel("Elevation (deg)")

    for ax in axes[:-1]:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)
    axes[-1].set_xlabel("Azimuth (deg)")

    if FOV_OVERLAY_BIVARIATE_AZIMUTH_LIMITS is None:
        crop_text = "Azimuth crop: none; full aggregate azimuth limits shown"
    else:
        crop_text = f"Azimuth crop: {azimuth_limits[0]:g} to {azimuth_limits[1]:g} deg"

    fig.suptitle(
        (
            f"{group}: aggregate good-RF center maps and Gaussian density contours | "
            f"EV > {min_ev:g}, p < {max_p_value:g}, good P/N > {min_peak_to_noise:g}\n"
            f"{crop_text}; contours replace the Gaussian-ellipse density panel; bottom density panel adds the 66% mass contour"
        ),
        fontweight="bold",
        y=0.995,
        fontsize=10,
    )
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.035, top=0.94, hspace=0.24)
    save_figure(
        fig,
        output_dir / f"{group}_aggregate_good_RF_center_density_vertical_summary_azimuth_crop_contours",
        dpi=300,
    )
    return 1



def density_mass_contour_level(coverage: np.ndarray, *, mass_fraction: float = 0.666) -> float | None:
    """Return an HPD-style density threshold enclosing the requested mass fraction.

    The returned level is chosen so that the sum of all grid-cell density values
    at or above that level contains approximately `mass_fraction` of the total
    positive density. With a uniform grid this is the discrete analogue of a
    contour enclosing 66% of the density mass.
    """

    values = np.asarray(coverage, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size == 0:
        return None

    total = float(np.sum(values))
    if not np.isfinite(total) or total <= 0:
        return None

    target = float(np.clip(mass_fraction, 0.0, 1.0)) * total
    ordered = np.sort(values)[::-1]
    cumulative = np.cumsum(ordered)
    index = int(np.searchsorted(cumulative, target, side="left"))
    index = min(max(index, 0), ordered.size - 1)
    level = float(ordered[index])
    if not np.isfinite(level) or level <= 0:
        return None
    return level


def plot_density_mass_contour(
    ax,
    coverage: np.ndarray,
    azimuth_axis: np.ndarray,
    elevation_axis: np.ndarray,
    *,
    mass_fraction: float = 0.666,
    color: str = "red",
    linewidth: float = 1.8,
) -> None:
    """Overlay one contour line enclosing the requested density mass fraction."""

    level = density_mass_contour_level(coverage, mass_fraction=mass_fraction)
    if level is None:
        ax.text(
            0.5,
            0.08,
            "66.6% density contour unavailable",
            ha="center",
            va="bottom",
            transform=ax.transAxes,
            color=color,
        )
        return

    azimuth_grid, elevation_grid = np.meshgrid(azimuth_axis, elevation_axis)
    coverage = np.asarray(coverage, dtype=float)
    finite_values = coverage[np.isfinite(coverage)]
    if finite_values.size == 0:
        return
    min_value = float(np.nanmin(finite_values))
    max_value = float(np.nanmax(finite_values))
    if not (min_value < level < max_value):
        if max_value > min_value:
            eps = np.finfo(float).eps * max(abs(max_value), 1.0) * 16.0
            level = min(max(level, min_value + eps), max_value - eps)
        else:
            ax.text(
                0.5,
                0.08,
                "90% density contour unavailable",
                ha="center",
                va="bottom",
                transform=ax.transAxes,
                color=color,
            )
            return

    ax.contour(
        azimuth_grid,
        elevation_grid,
        coverage,
        levels=[level],
        colors=[color],
        linewidths=linewidth,
    )

def plot_density_contour_panel_no_colorbar(
    ax,
    coverage: np.ndarray,
    azimuth_axis: np.ndarray,
    elevation_axis: np.ndarray,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    panel_vmax: float,
    title: str,
    *,
    filled: bool,
    n_levels: int = 10,
) -> None:
    """Draw MATLAB-style contour or contourf panels for an RF-density map."""

    azimuth_grid, elevation_grid = np.meshgrid(azimuth_axis, elevation_axis)
    coverage = np.asarray(coverage, dtype=float)
    finite_values = coverage[np.isfinite(coverage)]
    if finite_values.size == 0 or float(np.nanmax(finite_values)) <= 0:
        ax.text(
            0.5,
            0.5,
            "No positive density",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    else:
        vmax = float(panel_vmax)
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = float(np.nanmax(finite_values))
        levels = np.linspace(0.0, vmax, int(n_levels) + 1)
        if filled:
            ax.contourf(
                azimuth_grid,
                elevation_grid,
                coverage,
                levels=levels,
                cmap="viridis",
                vmin=0,
                vmax=vmax,
            )
        else:
            line_levels = levels[1:] if levels.size > 1 else levels
            ax.contour(
                azimuth_grid,
                elevation_grid,
                coverage,
                levels=line_levels,
                linewidths=0.9,
            )

    ax.set_xlim(azimuth_limits)
    ax.set_ylim(elevation_limits)
    ax.set_xlabel("Azimuth (deg)")
    ax.set_ylabel("Elevation (deg)")
    ax.set_title(title, pad=4)

def plot_density_panel_no_colorbar(
    ax,
    coverage: np.ndarray,
    azimuth_axis: np.ndarray,
    elevation_axis: np.ndarray,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    panel_vmax: float,
    title: str,
) -> None:
    ax.imshow(
        coverage,
        origin="lower",
        extent=(azimuth_axis[0], azimuth_axis[-1], elevation_axis[0], elevation_axis[-1]),
        aspect="equal",
        cmap="viridis",
        vmin=0,
        vmax=panel_vmax,
    )
    ax.set_xlim(azimuth_limits)
    ax.set_ylim(elevation_limits)
    ax.set_xlabel("Azimuth (deg)")
    ax.set_ylabel("Elevation (deg)")
    ax.set_title(title, pad=4)


def coverage_vmax(coverage: np.ndarray) -> float:
    vmax = float(np.nanmax(coverage)) if np.size(coverage) else 1.0
    if not np.isfinite(vmax) or vmax <= 0:
        return 1.0
    return vmax


def plot_density_panel(
    ax,
    coverage: np.ndarray,
    azimuth_axis: np.ndarray,
    elevation_axis: np.ndarray,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    panel_vmax: float,
    title: str,
    *,
    cax=None,
) -> None:
    im = ax.imshow(
        coverage,
        origin="lower",
        extent=(azimuth_axis[0], azimuth_axis[-1], elevation_axis[0], elevation_axis[-1]),
        aspect="equal",
        cmap="viridis",
        vmin=0,
        vmax=panel_vmax,
    )
    ax.set_xlim(azimuth_limits)
    ax.set_ylim(elevation_limits)
    ax.set_xlabel("Azimuth (deg)")
    ax.set_ylabel("Elevation (deg)")
    ax.set_title(title)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("Summed normalized Gaussian RF density")


def plot_aggregate_ev_diagnostics(
    group: str,
    ev_rows: list[dict],
    output_dir: Path,
    *,
    max_p_value: float,
) -> int:
    """Overlay session EV curves in grey and binwise session-average curves in black."""

    values = np.asarray([row["ev_value"] for row in ev_rows], dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0

    bins = np.linspace(AGGREGATE_EV_X_LIMITS[0], AGGREGATE_EV_X_LIMITS[1], 41)
    centers = (bins[:-1] + bins[1:]) / 2.0
    session_labels = sorted({row["session_label"] for row in ev_rows})
    sig_curves = []
    nonsig_curves = []
    ev_label = next((str(row.get("ev_label", "")) for row in ev_rows if row.get("ev_label")), "RF-only EV from final circular-shift test")

    fig, ax = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
    for label in session_labels:
        session_rows = [row for row in ev_rows if row["session_label"] == label]
        sig_values = [row["ev_value"] for row in session_rows if bool(row["p_lt_threshold"])]
        nonsig_values = [
            row["ev_value"]
            for row in session_rows
            if np.isfinite(row.get("p_value", np.nan)) and not bool(row["p_lt_threshold"])
        ]

        sig_curve = histogram_percent(sig_values, bins)
        if sig_curve is not None:
            sig_curves.append(sig_curve)
            ax.plot(centers, sig_curve, color="0.68", linestyle="-", linewidth=1.15, alpha=0.75)

        nonsig_curve = histogram_percent(nonsig_values, bins)
        if nonsig_curve is not None:
            nonsig_curves.append(nonsig_curve)
            ax.plot(centers, nonsig_curve, color="0.80", linestyle="--", linewidth=1.05, alpha=0.75)

    if sig_curves:
        ax.plot(
            centers,
            np.nanmean(np.vstack(sig_curves), axis=0),
            color="black",
            linestyle="-",
            linewidth=2.4,
        )
    if nonsig_curves:
        ax.plot(
            centers,
            np.nanmean(np.vstack(nonsig_curves), axis=0),
            color="black",
            linestyle="--",
            linewidth=2.4,
        )

    ax.set_xlabel(ev_label)
    ax.set_ylabel("Neurons (% within session p-value group)")
    ax.set_title(f"{group}: aggregate circular-shift RF EV session overlays")
    make_two_tick_square_axes(ax, x_zero=True, y_zero=True)
    ax.set_xlim(AGGREGATE_EV_X_LIMITS)
    ax.set_xticks(list(AGGREGATE_EV_X_LIMITS))
    save_figure(fig, output_dir / f"{group}_aggregate_EV_session_overlays", dpi=250)
    return 1


def plot_aggregate_rf_widths(
    group: str,
    rows: list[dict],
    output_dir: Path,
    *,
    min_peak_to_noise: float,
) -> int:
    """Overlay per-session good-RF width curves in grey and the average in black."""

    good_rows = [row for row in rows if row["is_good"] and np.isfinite(row["rf_width_deg"])]
    widths = np.asarray([row["rf_width_deg"] for row in good_rows], dtype=float)
    if widths.size == 0:
        return 0

    bins = np.arange(0.0, math.ceil(float(np.nanmax(widths)) / 2.0) * 2.0 + 2.0, 2.0)
    if bins.size < 2:
        bins = np.linspace(0.0, float(np.nanmax(widths)) + 2.0, 10)
    centers = (bins[:-1] + bins[1:]) / 2.0
    session_labels = sorted({row["session_label"] for row in rows})
    curves = []

    fig, ax = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
    for label in session_labels:
        session_widths = [row["rf_width_deg"] for row in good_rows if row["session_label"] == label]
        arr = np.asarray(session_widths, dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            continue
        curve, _ = np.histogram(arr, bins=bins)
        curve = curve.astype(float)
        curves.append(curve)
        ax.plot(centers, curve, color="0.72", linewidth=1.25, alpha=0.75)

    if curves:
        ax.plot(centers, np.nanmean(np.vstack(curves), axis=0), color="black", linewidth=2.4)

    ax.set_xlabel("RF width (deg) = (|sigma_x| + |sigma_y|) / 2")
    ax.set_ylabel("Number of good RFs")
    ax.set_title(f"{group}: aggregate good-RF width session overlays | P/N > {min_peak_to_noise:g}")
    make_two_tick_square_axes(ax, x_zero=True, y_zero=True)
    save_figure(fig, output_dir / f"{group}_aggregate_good_RF_widths", dpi=250)
    return 1


def group_total_neurons_from_rows(rows: list[dict]) -> int:
    """Return the pooled neuron denominator across unique sessions in aggregate rows."""

    totals_by_session: dict[tuple[str, str, str], int] = {}
    for row in rows:
        key = (str(row.get("subject", "")), str(row.get("date", "")), str(row.get("session", "")))
        try:
            total = int(row.get("session_total_neurons", 0))
        except Exception:
            total = 0
        if total > 0:
            totals_by_session[key] = max(totals_by_session.get(key, 0), total)
    total_neurons = int(sum(totals_by_session.values()))
    return max(total_neurons, 1)


def aggregate_bins(values: np.ndarray, *, n_bins: int) -> np.ndarray:
    """Bins with explicit padding so EV curves are not clipped at plot edges."""

    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        lo, hi = -0.01, 0.1
    else:
        lo = float(np.nanmin(arr))
        hi = float(np.nanmax(arr))
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            lo, hi = -0.01, 0.1
    pad = max((hi - lo) * 0.04, 0.005)
    return np.linspace(lo - pad, hi + pad, int(n_bins) + 1)


def histogram_percent(values: list[float] | np.ndarray, bins: np.ndarray) -> np.ndarray | None:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    counts, _ = np.histogram(arr, bins=bins)
    return counts.astype(float) * 100.0 / float(arr.size)


def session_colors(labels: list[str]) -> dict[str, Any]:
    cmap = plt.get_cmap("tab20", max(1, len(labels)))
    return {label: cmap(idx) for idx, label in enumerate(labels)}


def subfield_name(best_subfield: int) -> str:
    return {1: "ON", 2: "OFF", 3: "ON+OFF"}.get(int(best_subfield), "unknown")


def pass_fail(value: bool) -> str:
    return "PASS" if bool(value) else "FAIL"


if __name__ == "__main__":
    raise SystemExit(main())
