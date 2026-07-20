# ALF exports from `ONEpy`

This note documents what the current Python exporters in `ONEpy/` write to disk.  It is based on the exporter code plus the local IBL documentation in `IBL documentation/`: `[PUBLIC] IBL DatasetTypes.xlsx`, `[PUBLIC] IBL dataset types.docx`, and `Open_Neurophysiology_Environment_Filename_Convention (5).pdf`.

## Minimal ALF structure for these RF analyses

ALF is a filename convention, not a file format.  The IBL filename convention uses names of the form:

```text
object.attribute.extension
```

The convention document states that files sharing an object name must share the same leading dimension, and that timing arrays whose attribute is `times` or ends in `times` are in seconds relative to experiment start.  Collections are represented as subdirectories when the same dataset type occurs multiple times in one experiment.

For this repository's sparse-noise RF workflow, the minimal exported session is therefore:

```text
<dataoutroot>/<subject>/<date>/<session>/
  alf/
    _ibl_passiveRFM.times.npy
    _ibl_sparseNoise.times.npy
    _ibl_sparseNoise.xy.npy
    _ibl_wheel.position.npy
    _ibl_wheel.timestamps.npy
    _ibl_wheel.velocity.npy
    FOV_00/
      mpci.*.npy
      mpciROIs.*.npy
      mpciROITypes.names.tsv
      mpciROIs.masks.sparse_npz
      mpciMeanImage.images.npy
    FOV_01/
      ...
  raw_passive_data/
    _iblrig_RFMapStim.raw.bin
```

The exporters currently use documented IBL dataset names and namespaces (`_ibl_`, `_iblrig_`, and standard `mpci*` names).  They do not currently emit datasets under the lab `_FANCi_` namespace.

## `single_session_calcium_export`

Defined in `calcium_export.py`.  This exporter writes one ALF collection per Suite2p plane/FOV:

```text
<dataoutroot>/<subject>/<date>/<session>/alf/FOV_NN/
```

`FOV_NN` comes from plane-level `db.npy["iroi"]` when present, otherwise from the sorted `plane*` folder order.  The exporter currently accepts `stim_type="sparse_noise"`; `stim_type="gratings"` raises `NotImplementedError`.

Dimension symbols used below:

- `nFrames`: exported frames from the selected sparse-noise raw session.
- `nROIs`: number of Suite2p ROIs in that plane.
- `H`, `W`: mean-image height and width from `ops["meanImg"]`.
- `nROItypes`: currently 2, because the exporter writes `0 = nonCell` and `1 = neuron`.

### `mpci.times.npy`

- IBL dimension: `[nframes]`
- Exporter shape: `(nFrames,)`
- dtype: `float64`
- Contents: time of each imaging frame.  The IBL spreadsheet notes possible finer offsets in `mpciStack.timeshift.npy`; this exporter does not write that file.
- Exporter details: times are converted to session-relative seconds.  With `fix_interleave=True`, the exporter corrects mROI FOV timing by offsetting each FOV within each DI3/None scan cycle according to cumulative FOV height.

### `mpci.badFrames.npy`

- IBL dimension: `[nframes]`
- Exporter shape: `(nFrames,)`
- dtype: `bool`
- Contents: bad frames identified by processing software; the IBL spreadsheet describes these as frames where registration was considered poor.
- Exporter details: starts from `ops["badframes"]` when present, otherwise all false.  Timing-derived bad frames are OR'ed in.  In the supported terminal missing-`None` timing edge case, the final frame is preserved and marked bad.

### `mpci.ROIActivityF.npy`

- IBL dimension: `[nFrames, nROIs]`
- Exporter shape: `(nFrames, nROIs)`
- dtype: `float32`
- Contents: mean activity of all pixels in each ROI.
- Exporter details: Suite2p `F.npy` is loaded as `[nROIs, nFrames]`, sliced to the selected session, then transposed before saving.

### `mpci.ROINeuropilActivityF.npy`

- IBL dimension: `[nFrames, nROIs]`
- Exporter shape: `(nFrames, nROIs)`
- dtype: `float32`
- Contents: mean activity of neuropil pixels neighboring each ROI.
- Exporter details: Suite2p `Fneu.npy` is sliced and transposed in the same way as `F.npy`.

### `mpci.ROIActivityDeconvolved.npy`

- IBL dimension: `[nFrames, nROIs]`
- Exporter shape: `(nFrames, nROIs)`
- dtype: `float32`
- Contents: neuropil-subtracted deconvolved activity of each ROI using standard parameters.
- Exporter details: Suite2p `spks.npy` is sliced and transposed in the same way as `F.npy`.

### `mpciROIs.cellClassifier.npy`

- IBL dimension: `[nROIs]`
- Exporter shape: `(nROIs,)`
- dtype: `float64`
- Contents: floating-point cell classifier score for each ROI, documented as ranging between 0 and 1.
- Exporter details: saved from `iscell[:, 1]` when available.  If `iscell.npy` has only one column, the exporter writes the 0/1 labels as the classifier values and warns.

### `mpciROIs.mpciROITypes.npy`

- IBL dimension: `[nROIs]`
- Exporter shape: `(nROIs,)`
- dtype: `int16`
- Contents: numerical enumeration of ROI type for each ROI.
- Exporter details: saved from `iscell[:, 0]`, so current values are Suite2p cell labels: `0` for non-cell and `1` for neuron/cell.

### `mpciROITypes.names.tsv`

- IBL dimension: `[nROItypes, 2]`
- Exporter shape: 2 data rows plus a header row.
- file type: TSV
- Contents: human-readable definitions of ROI types.  The IBL spreadsheet describes the first column as `roi_values` and the second as a string label column.
- Exporter details: writes columns `roi_values` and `roi_names`, with rows `(0, nonCell)` and `(1, neuron)`.  Note the local header `roi_names` differs from the spreadsheet wording, which says `roi_labels`.

### `mpciROIs.stackPos.npy`

- IBL dimension: `[nROIs, 3]`
- Exporter shape: `(nROIs, 3)`
- dtype: `int64`
- Contents: Y, X, and Z/plane pixel coordinates of each ROI centroid.
- Exporter details: uses Suite2p `stat[*]["med"]` for rounded `(Y, X)` and the Suite2p plane index for `Z`.

### `mpciROIs.masks.sparse_npz`

- IBL dimension: `[nROIs, H, W]`
- Exporter shape: sparse COO array with shape `(nROIs, H, W)`
- dtype: `float32` data values
- Contents: floating-point mask of each ROI, saved as a sparse array.
- Exporter details: builds the sparse mask from Suite2p `stat[*]["ypix"]`, `stat[*]["xpix"]`, and `stat[*]["lam"]`.  The first sparse axis is ROI index, followed by image Y and X.

### `mpciMeanImage.images.npy`

- IBL dimension: `[H, W, nPlanes, nChannels]`
- Exporter shape: `(H, W, 1, 1)`
- dtype: `float64`
- Contents: mean image for each plane and channel.  The IBL spreadsheet notes that channel 0 should be the calcium-sensitive channel.
- Exporter details: saved from Suite2p `ops["meanImg"]` after adding singleton plane and channel axes.

## `single_session_sparse_noise_export`

Defined in `sparse_noise_export.py`.  This exporter writes session-level passive RF mapping files to:

```text
<dataoutroot>/<subject>/<date>/<session>/alf/
<dataoutroot>/<subject>/<date>/<session>/raw_passive_data/
```

Dimension symbols used below:

- `nFrames`: number of sparse-noise movie frames, which must match `times_grid.csv`.
- `nRows`, `nCols`: sparse-noise grid size inferred from the raw binary; current accepted shapes are `(10, 30)` and `(7, 27)`.
- `nSparseNoise`: number of non-gray sparse-noise squares across all frames.

### `_ibl_passiveRFM.times.npy`

- IBL dimension: `[nFrames]`
- Exporter shape: `(nFrames,)`
- dtype: `float64`
- Contents: passive RF mapping frame times.
- Exporter details: loaded from `times_grid.csv`, flattened, validated against the decoded movie frame count, and converted to session-relative seconds.

### `_iblrig_RFMapStim.raw.bin`

- IBL dimension: `[nFrames, nx, ny]`
- Exporter payload: flat raw binary bytes, with expected size `nFrames * nRows * nCols`
- dtype: `uint8`
- collection: `raw_passive_data`
- Contents: raw RF mapping matrix.
- Exporter details: preserves the original `SparseNoise_Log.bin` bytes as unsigned 8-bit data.  The decoded stimulus movie is used only to extract sparse-noise event times and positions.  Because this is a `.bin` file, the shape is not carried in the file itself; the exporter validates the byte count against the decoded movie shape before writing.

### `_ibl_sparseNoise.times.npy`

- IBL dimension: `[nSparseNoise]`
- Exporter shape: `(nSparseNoise,)`
- dtype: saved as `float64`; the IBL spreadsheet does not list a dtype.
- Contents: times when sparse-noise stimulus squares appeared.
- Exporter details: one time is emitted for every non-zero sparse-noise square.  If multiple squares are non-zero in one frame, that frame time appears multiple times.

### `_ibl_sparseNoise.xy.npy`

- IBL dimension: `[nSparseNoise, 2]`
- Exporter shape: `(nSparseNoise, 2)`
- dtype: saved as `float64`; the IBL spreadsheet does not list a dtype.
- Contents: x/y screen coordinates of sparse-noise stimulus squares.  The IBL spreadsheet explicitly leaves the unit unresolved with `(WHAT UNIT?)`.
- Exporter details: this exporter uses degrees from `stimulus_edges_deg`, defaulting to `(left=-135, right=135, bottom=-40, top=40)`.  Coordinates are bin centers computed from the inferred stimulus grid.

## `single_session_wheel_export`

Defined in `wheel_export.py`.  This exporter writes session-level wheel files to:

```text
<dataoutroot>/<subject>/<date>/<session>/alf/
```

Dimension symbols used below:

- `nWheelSamples`: number of rows in `encoder_log.csv`.

### `_ibl_wheel.position.npy`

- IBL dimension: `[nWheelSamples]`
- Exporter shape: `(nWheelSamples,)`
- dtype: `float64`
- Contents: absolute wheel rotation in radians, using the mathematical convention.
- Exporter details: treats column 2 of `encoder_log.csv` as a signed 16-bit wrapped encoder position, unwraps it, and converts counts to radians using `2*pi / encoder_counts_per_revolution` with a default of `65536` counts/revolution.

### `_ibl_wheel.timestamps.npy`

- IBL dimension in the spreadsheet: `[nWheelSamples, 2]`
- Exporter shape: `(nWheelSamples,)`
- dtype: `float64`
- Contents: wheel sample times in seconds from session start, non-evenly spaced.
- Exporter details: uses column 1 of `encoder_log.csv`, converted to session-relative seconds.  The exporter validates that timestamps are strictly increasing.  Note that the current exporter writes a 1D vector even though the local IBL spreadsheet lists `[nWheelSamples, 2]`.

### `_ibl_wheel.velocity.npy`

- IBL dimension: `[nWheelSamples]`
- Exporter shape: `(nWheelSamples,)`
- dtype: `float64`
- Contents: tangential wheel velocity in rad/s, with positive values documented as counter-clockwise.
- Exporter details: computed as the first difference of exported position divided by the first difference of timestamps.  The first value is filled from the first finite velocity when available; for fewer than two samples, the exporter writes zeros.

## Files not written by these exporters

The IBL multi-photon calcium section also defines related datasets such as `mpci.mpciFrameQC.npy`, `mpciFrameQC.names.tsv`, `mpciROIs.neuropilMasks.sparse_npz`, `mpciROIs.mlapdv*`, `mpciROIs.brainLocationIds*`, `mpciROIs.uuids`, and `_suite2p_ROIData.raw`.  These are documented IBL dataset types, but the current `ONEpy` exporters do not write them.

Likewise, the wheel section defines wheel-movement datasets such as `_ibl_wheelMoves.intervals` and `_ibl_wheelMoves.peakAmplitude`; the current wheel exporter writes only position, timestamps, and velocity.
