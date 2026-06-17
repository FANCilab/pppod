"""Plot FANCi ALF RF-mapping outputs from main_mapRFALF_debug_curated_outputs_verbose.m.

Hardcoded defaults at the top mirror the current RF script:
    ALF_ROOT / SUBJECT / DATE / SESSION / alf / FOV_XX
    RESULTS_ROOT / SUBJECT / DATE / SESSION / alf / FOV_XX

All figures are written to:
    PLOTS_ROOT / SUBJECT / DATE / SESSION / alf / FOV_XX

Every plot is saved as both PNG and SVG. Per-neuron RF maps are split into peak/noise PASS and FAIL folders. RF-size histograms are made from valid Gaussian sigma parameters.

RF map image display note:
    MATLAB saves RF edges as [left, right, top, bottom], e.g. [-135 135 40 -40].
    Matplotlib imshow expects extent=(left, right, bottom, top).  The helper
    rf_image_extent() below performs this conversion so the RF pixels and
    Gaussian fit overlay use the same elevation convention.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

try:
    import h5py
except Exception:  # pragma: no cover - only used when h5py is absent on user machine
    h5py = None

try:
    from pycolormap_2d import (
        ColorMap2DBremm,
        ColorMap2DCubeDiagonal,
        ColorMap2DSchumann,
        ColorMap2DSteiger,
        ColorMap2DTeuling2,
        ColorMap2DZiegler,
    )
except Exception:  # pragma: no cover - optional plotting dependency on user machine
    ColorMap2DBremm = None
    ColorMap2DCubeDiagonal = None
    ColorMap2DSchumann = None
    ColorMap2DSteiger = None
    ColorMap2DTeuling2 = None
    ColorMap2DZiegler = None


# =============================================================================
# User-editable session configuration
# =============================================================================

SUBJECT = "NM023"
DATE = "2026-02-09"
SESSION = "2"

ALF_ROOT = Path(r"D:\Data\ALF2P")
RESULTS_ROOT = Path(r"D:\Results\ALF2P")
PLOTS_ROOT = Path(r"D:\Plots")
#RAW_2P_ROOT = Path(r"\\10.233.25.135\FANCiNAS1\Data\2P")
RAW_2P_ROOT = Path(r"D:\Data\2P")
# Set to None to process every FOV_* found under the ALF session.
FOV_NAMES: list[str] | None = None

MIN_EV = 0.01
MAX_P_VALUE = 0.05
MIN_PEAK_TO_NOISE = 7.7
MAX_NEURON_RF_PLOTS: int | None = None  # None = plot all significant RF neurons

MAKE_GAUSSIAN_COVERAGE = True
MAKE_PER_NEURON_RF_MAPS = True
MAKE_FOV_OVERLAYS = True
MAKE_NEUTRAL_ROI_OVERLAYS = True
MAKE_SHIFT_EV_HISTOGRAM = True
MAKE_PREDICTED_VS_ACTUAL = True
MAKE_PEAK_TIME_COURSES = True
MAKE_RF_SIZE_HISTOGRAM = True

# Trace panel downsampling for very dense plots. 1 keeps every point.
TRACE_PLOT_STRIDE = 1

# pycolormap-2d is used for RF-center colors in the FOV overlay plots and
# their 2D color legend. Install with: pip install pycolormap-2d
PYCOLORMAP_2D_NAME = "Bremm"
PYCOLORMAP_2D_WHEEL_RESOLUTION = 301


# =============================================================================
# Styling
# =============================================================================

PEAK_ROI_ALPHA = 1.0
ROI_MARKER_SIZE = 3
NON_PEAK_CONTOUR_LINEWIDTH = 1.1
NEUTRAL_ROI_ALPHA = 0.1

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


# =============================================================================
# Main orchestration
# =============================================================================


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)

    subject = args.subject
    date = args.date
    session = args.session
    alf_session = args.alf_root / subject / date / session / "alf"
    results_session = args.results_root / subject / date / session / "alf"
    plots_session = args.plots_root / subject / date / session / "alf"

    if args.fovs:
        fov_names = args.fovs
    elif FOV_NAMES is not None:
        fov_names = FOV_NAMES
    else:
        fov_names = discover_fovs(alf_session, results_session)

    if not fov_names:
        raise FileNotFoundError(
            f"No FOV folders found under ALF={alf_session} or results={results_session}"
        )

    print("\n========== FANCi RF plotting ==========")
    print(f"subject/date/session: {subject} / {date} / {session}")
    print(f"ALF session:          {alf_session}")
    print(f"Results session:      {results_session}")
    print(f"Plots session:        {plots_session}")
    print(f"FOVs:                 {', '.join(fov_names)}")

    raw_suite2p = None
    if not args.skip_fov_overlays and (MAKE_FOV_OVERLAYS or MAKE_NEUTRAL_ROI_OVERLAYS):
        try:
            raw_suite2p = find_raw_suite2p_folder(
                args.raw_2p_root,
                subject=subject,
                date=date,
                session=session,
            )
            print(f"Raw Suite2p folder:   {raw_suite2p}")
        except FileNotFoundError as exc:
            raw_suite2p = None
            print(f"[WARN] FOV overlay plots will be skipped: {exc}")

    summary_rows = []
    for fov_name in fov_names:
        print(f"\n=== Plotting {fov_name} ===")
        fov_results_dir = results_session / fov_name
        fov_alf_dir = alf_session / fov_name
        fov_output_dir = plots_session / fov_name
        fov_output_dir.mkdir(parents=True, exist_ok=True)

        result = load_fanci_result(fov_results_dir)
        counts = plot_one_fov(
            result,
            fov_name=fov_name,
            fov_alf_dir=fov_alf_dir,
            output_dir=fov_output_dir,
            raw_suite2p_folder=raw_suite2p,
            subject=subject,
            date=date,
            session=session,
            min_ev=args.min_ev,
            max_p_value=args.max_p_value,
            min_peak_to_noise=args.min_peak_to_noise,
            max_neuron_plots=args.max_neuron_plots,
            skip_fov_overlays=args.skip_fov_overlays,
        )
        summary_rows.append({"fov": fov_name, **counts, "output_dir": str(fov_output_dir)})
        print(f"[DONE] {fov_name}: {counts}")

    summary_path = plots_session / "plot_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(
            {
                "subject": subject,
                "date": date,
                "session": session,
                "alf_session": str(alf_session),
                "results_session": str(results_session),
                "plots_session": str(plots_session),
                "raw_suite2p_folder": str(raw_suite2p) if raw_suite2p is not None else None,
                "fovs": summary_rows,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"\nSummary written to: {summary_path}")
    return 0


def plot_one_fov(
    result: FanciResult,
    *,
    fov_name: str,
    fov_alf_dir: Path,
    output_dir: Path,
    raw_suite2p_folder: Path | None,
    subject: str,
    date: str,
    session: str,
    min_ev: float,
    max_p_value: float,
    min_peak_to_noise: float,
    max_neuron_plots: int | None,
    skip_fov_overlays: bool,
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
        )

    if MAKE_PER_NEURON_RF_MAPS:
        counts["per_neuron_rf_maps"] = plot_per_neuron_metadata(
            result,
            output_dir / "per_neuron_metadata",
            min_ev=min_ev,
            max_p_value=max_p_value,
            min_peak_to_noise=min_peak_to_noise,
            max_neuron_plots=max_neuron_plots,
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
    return FanciResult(results_dir, arrays, mat_path if mat_path.exists() else None)


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


def save_figure(fig, base_path: Path, *, dpi: int = 220) -> list[Path]:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    png = base_path.with_suffix(".png")
    svg = base_path.with_suffix(".svg")
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return [png, svg]


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

    fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
    if rf_size.size > 0:
        max_size = float(np.nanmax(rf_size))
        if np.isfinite(max_size) and max_size > 0:
            bin_width = 2.0
            bins = np.arange(0.0, math.ceil(max_size / bin_width) * bin_width + bin_width, bin_width)
            if bins.size < 2:
                bins = 10
        else:
            bins = 10
        ax.hist(rf_size, bins=bins, edgecolor="black", linewidth=0.8)
        ax.axvline(np.nanmedian(rf_size), linestyle="--", linewidth=1.3, color="black")
        stats_text = (
            f"n = {rf_size.size}\n"
            f"median = {np.nanmedian(rf_size):.2f} deg\n"
            f"mean = {np.nanmean(rf_size):.2f} deg"
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
    ax.grid(True, axis="y", alpha=0.25)

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
    n_azimuth_grid: int = 500,
    n_elevation_grid: int = 250,
    gaussian_sigma_radius: float = 2.0,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    gauss_pars = result.require("gaussPars")
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
    cb = plt.colorbar(im, ax=ax)
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
    if max_neuron_plots is not None:
        significant = significant[:max_neuron_plots]

    edges = get_edges(result)
    timestamps = np.ravel(result.get("timestamps", np.arange(maps.shape[3], dtype=float)))
    rf_types = ("ON", "OFF", "ON+OFF")

    index_rows = []
    n_plotted = 0
    n_pass = 0
    n_fail = 0
    for cell_id in significant:
        plot_map, mx_time, best_subfield = plot_map_for_neuron(maps, result, int(cell_id))
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

        cmap = RF_ON_CMAP if best_subfield == 1 else RF_OFF_CMAP if best_subfield == 2 else RF_BWR_CMAP
        mx = np.nanmax(np.abs(plot_map))
        if not np.isfinite(mx) or mx == 0:
            mx = 1.0

        fig, ax = plt.subplots(figsize=(9, 7.25), constrained_layout=True)
        extent = rf_image_extent(edges)
        im = ax.imshow(plot_map, extent=extent, origin="upper", cmap=cmap, vmin=-mx, vmax=mx)
        ax.set_aspect("equal")
        ax.set_ylim(rf_elevation_limits(edges))
        ax.set_xlabel("Azimuth (deg)")
        ax.set_ylabel("Elevation (deg)")
        plt.colorbar(im, ax=ax)
        plot_gaussian_on_rf(ax, result, int(cell_id))

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
                ]
            )
        )

        base_name = (
            f"RFresult_{cell_id + 1:04d}_sourceROI0_{source_roi_text(result, cell_id)}_"
            f"EV_{pass_fail(passes_ev)}_peakToNoise_{pass_fail(passes_peak)}"
        )
        save_figure(fig, target_dir / base_name, dpi=220)
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

    return selected_rf[:, :, mx_time], mx_time, best_subfield


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


def plot_gaussian_on_rf(ax, result: FanciResult, cell_id: int) -> None:
    gauss = result.get("gaussPars")
    edges = get_edges(result)
    if gauss is None or cell_id >= gauss.shape[0]:
        return
    pars = gauss[cell_id, :]
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

    fig, ax = plt.subplots(figsize=(9.5, 6.2), constrained_layout=True)
    if sig.sum() > 0:
        ax.hist(
            ev_values[sig],
            bins=bins,
            weights=np.full(sig.sum(), 100.0 / sig.sum()),
            histtype="stepfilled",
            alpha=0.45,
            label=f"p < {max_p_value:g} (n={sig.sum()})",
        )
        ax.hist(
            ev_values[sig],
            bins=bins,
            weights=np.full(sig.sum(), 100.0 / sig.sum()),
            histtype="step",
            linewidth=1.8,
        )
    if nonsig.sum() > 0:
        ax.hist(
            ev_values[nonsig],
            bins=bins,
            weights=np.full(nonsig.sum(), 100.0 / nonsig.sum()),
            histtype="stepfilled",
            alpha=0.35,
            label=f"p >= {max_p_value:g} (n={nonsig.sum()})",
        )
        ax.hist(
            ev_values[nonsig],
            bins=bins,
            weights=np.full(nonsig.sum(), 100.0 / nonsig.sum()),
            histtype="step",
            linewidth=1.8,
        )
    ax.axvline(0, color="black", linewidth=0.8, alpha=0.45)
    ax.set_xlabel(ev_label)
    ax.set_ylabel("Neurons (% within group)")
    ax.set_title(f"{subject} {date} session {session} {fov_name}: circular-shift RF EV distribution")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.2)

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

    gauss_pars = result.require("gaussPars")
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

    rf_colormap = make_pycolormap_2d(azimuth_limits, elevation_limits)

    fig = plt.figure(figsize=(11.5, 8.5), constrained_layout=False)
    ax_fov = fig.add_axes([0.06, 0.07, 0.60, 0.82])
    ax_wheel = fig.add_axes([0.705, 0.425, 0.245, 0.185])

    ax_fov.imshow(fov, cmap="gray", origin="upper")
    draw_rf_roi_masks(
        ax=ax_fov,
        stat=stat,
        source_roi_indices=roi_indices,
        azimuth=azimuth,
        elevation=elevation,
        significant=significant,
        peak_pass=peak_pass,
        azimuth_limits=azimuth_limits,
        elevation_limits=elevation_limits,
        rf_colormap=rf_colormap,
    )
    ax_fov.set_title(
        f"{fov_name} FOV overlay: pycolormap-2d RF center color ({int(significant.sum())} significant, {int(peak_pass.sum())} P/N pass)",
        fontsize=15,
    )
    ax_fov.set_axis_off()
    plot_2d_colorwheel(
        ax_wheel,
        azimuth_limits,
        elevation_limits,
        rf_colormap=rf_colormap,
        azimuth=azimuth,
        elevation=elevation,
        significant=significant,
        peak_pass=peak_pass,
    )
    fig.suptitle(
        f"{subject} {date} session {session} {fov_name}: significant RF ROIs on Suite2p FOV | "
        f"pycolormap-2d={PYCOLORMAP_2D_NAME}; filled ROIs = P/N > {min_peak_to_noise:g}; contours = significant only",
        fontsize=17,
    )
    save_figure(fig, output_dir / f"{fov_name}_significant_ROIs_on_FOV_2D_colorwheel", dpi=220)

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
                "explained_variance": float(ev[rf_idx]),
                "p_value": float(p_values[rf_idx]),
                "peak_to_noise": float(peak_to_noise[rf_idx]),
                "passes_peak_to_noise": bool(peak_pass[rf_idx]),
            }
        )
    write_metadata(rows, output_dir / f"{fov_name}_significant_roi_fov_overlay_metadata.csv")
    return 1


def draw_rf_roi_masks(
    *,
    ax,
    stat: np.ndarray,
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
    for rf_idx in np.flatnonzero(peak_pass):
        draw_single_roi_mask(
            ax,
            stat,
            int(source_roi_indices[rf_idx]),
            rf_position_to_rgb(azimuth[rf_idx], elevation[rf_idx], azimuth_limits, elevation_limits, rf_colormap=rf_colormap),
            alpha=PEAK_ROI_ALPHA,
            zorder=9,
        )


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
    ax.contour(
        np.arange(x0, x1 + 1),
        np.arange(y0, y1 + 1),
        local_mask,
        levels=[0.5],
        colors=[color],
        linewidths=linewidth,
        alpha=1.0,
        zorder=zorder,
    )


def draw_single_roi_mask(ax, stat: np.ndarray, roi_idx: int, color, *, alpha: float, zorder: int) -> None:
    if roi_idx >= len(stat) or roi_idx < 0:
        return
    roi = stat_entry(stat[roi_idx])
    ypix = np.asarray(roi["ypix"], dtype=int)
    xpix = np.asarray(roi["xpix"], dtype=int)
    ax.scatter(xpix, ypix, s=ROI_MARKER_SIZE, c=[color], marker="s", linewidths=0, alpha=alpha, zorder=zorder)


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
) -> None:
    azimuth_axis = np.linspace(azimuth_limits[0], azimuth_limits[1], PYCOLORMAP_2D_WHEEL_RESOLUTION)
    elevation_axis = np.linspace(elevation_limits[0], elevation_limits[1], PYCOLORMAP_2D_WHEEL_RESOLUTION)
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
        aspect="auto",
    )
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
    ax.set_title(f"pycolormap-2d {PYCOLORMAP_2D_NAME}", fontsize=11, pad=5)
    ax.tick_params(labelsize=8, pad=2)


def make_pycolormap_2d(
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
):
    if ColorMap2DBremm is None:
        raise ImportError(
            "pycolormap-2d is required for RF-center FOV overlay colors. "
            "Install it with: pip install pycolormap-2d"
        )
    cmap_classes = {
        "Bremm": ColorMap2DBremm,
        "CubeDiagonal": ColorMap2DCubeDiagonal,
        "Schumann": ColorMap2DSchumann,
        "Steiger": ColorMap2DSteiger,
        "Teuling2": ColorMap2DTeuling2,
        "Ziegler": ColorMap2DZiegler,
    }
    cmap_class = cmap_classes.get(PYCOLORMAP_2D_NAME)
    if cmap_class is None:
        valid = ", ".join(sorted(cmap_classes))
        raise ValueError(f"Unknown PYCOLORMAP_2D_NAME={PYCOLORMAP_2D_NAME!r}. Valid options: {valid}")
    return cmap_class(
        range_x=(float(azimuth_limits[0]), float(azimuth_limits[1])),
        range_y=(float(elevation_limits[0]), float(elevation_limits[1])),
    )


def sample_pycolormap_2d(rf_colormap, x: float, y: float) -> np.ndarray:
    if hasattr(rf_colormap, "sample"):
        color = rf_colormap.sample(float(x), float(y))
    else:
        color = rf_colormap(float(x), float(y))
    color = np.asarray(color, dtype=float).reshape(-1)
    if color.size < 3:
        raise ValueError("pycolormap-2d returned a color with fewer than 3 channels")
    color = color[:3]
    if np.nanmax(color) > 1.0:
        color = color / 255.0
    return np.clip(color, 0.0, 1.0)


def rf_position_to_rgb(
    azimuth: np.ndarray | float,
    elevation: np.ndarray | float,
    azimuth_limits: tuple[float, float],
    elevation_limits: tuple[float, float],
    *,
    rf_colormap,
) -> np.ndarray:
    """Map RF centers to RGB using pycolormap-2d.

    The ranges are the visual-field azimuth/elevation limits, so the same
    coordinates used in the Gaussian RF fits are passed directly to the 2D
    colormap. This replaces the older hand-rolled HSV mapping.
    """
    azimuth = np.asarray(azimuth, dtype=float)
    elevation = np.asarray(elevation, dtype=float)
    az_b, el_b = np.broadcast_arrays(azimuth, elevation)
    rgb = np.empty(az_b.shape + (3,), dtype=float)

    flat_rgb = rgb.reshape(-1, 3)
    flat_az = az_b.ravel()
    flat_el = el_b.ravel()
    for i, (az, el) in enumerate(zip(flat_az, flat_el)):
        if np.isfinite(az) and np.isfinite(el):
            flat_rgb[i, :] = sample_pycolormap_2d(rf_colormap, float(az), float(el))
        else:
            flat_rgb[i, :] = np.array([0.5, 0.5, 0.5], dtype=float)
    return rgb

def plot_rois(*, fov: np.ndarray, stat: np.ndarray, roi_indices: np.ndarray, title: str, output_base: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 11), constrained_layout=True)
    ax.imshow(fov, cmap="gray", origin="upper")
    for roi_idx in roi_indices:
        roi = stat_entry(stat[int(roi_idx)])
        ypix = np.asarray(roi["ypix"], dtype=int)
        xpix = np.asarray(roi["xpix"], dtype=int)
        ax.scatter(xpix, ypix, s=0.5, c="#00ff00", edgecolors="#00ff00", marker="s", linewidths=0, alpha=NEUTRAL_ROI_ALPHA)
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_axis_off()
    save_figure(fig, output_base, dpi=220)


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


def find_raw_suite2p_folder(raw_2p_root: Path, *, subject: str, date: str, session: str) -> Path:
    compact_date = date.replace("-", "")
    candidates = [
        raw_2p_root / subject / compact_date / session / "processed",
        raw_2p_root / subject / compact_date / session / "Processed",
        raw_2p_root / subject / compact_date / session,
        raw_2p_root / subject / "Processed" / compact_date / session,
        raw_2p_root / subject / "processed" / compact_date / session,
        raw_2p_root / subject / date / session / "processed",
        raw_2p_root / subject / date / session / "Processed",
        raw_2p_root / subject / date / session,
    ]
    for candidate in candidates:
        if candidate.is_dir() and any(candidate.glob("plane*/reg_outputs.npy")):
            return candidate
    tried = "\n  ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"Could not find local raw Suite2p folder with plane*/reg_outputs.npy. Tried:\n  {tried}")


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
                return normalize_fov(image)
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


def subfield_name(best_subfield: int) -> str:
    return {1: "ON", 2: "OFF", 3: "ON+OFF"}.get(int(best_subfield), "unknown")


def pass_fail(value: bool) -> str:
    return "PASS" if bool(value) else "FAIL"


# =============================================================================
# CLI
# =============================================================================


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot latest FANCi RF mapping outputs.")
    parser.add_argument("--subject", default=SUBJECT)
    parser.add_argument("--date", default=DATE)
    parser.add_argument("--session", default=SESSION)
    parser.add_argument("--alf-root", type=Path, default=ALF_ROOT)
    parser.add_argument("--results-root", type=Path, default=RESULTS_ROOT)
    parser.add_argument("--plots-root", type=Path, default=PLOTS_ROOT)
    parser.add_argument("--raw-2p-root", type=Path, default=RAW_2P_ROOT)
    parser.add_argument("--fovs", nargs="*", default=None, help="Explicit FOV names, e.g. FOV_00 FOV_01")
    parser.add_argument("--min-ev", type=float, default=MIN_EV)
    parser.add_argument("--max-p-value", type=float, default=MAX_P_VALUE)
    parser.add_argument("--min-peak-to-noise", type=float, default=MIN_PEAK_TO_NOISE)
    parser.add_argument("--max-neuron-plots", type=int, default=MAX_NEURON_RF_PLOTS)
    parser.add_argument("--skip-fov-overlays", action="store_true")
    return parser


if __name__ == "__main__":
    raise SystemExit(main())
