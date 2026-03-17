Grating MATLAB workflow

This package analyzes neuronal responses to grating stimuli from flattened or
canonical-format inputs and saves per-neuron summary plots.

Input format (preferred)
------------------------
data.tc   : [nStim, nRepeats, nTimes, nNeurons]
data.amp  : [nStim, nRepeats, nNeurons]
where nStim = nDir * nSize * nTf * nSf and stimuli are linearized in
canonical order [dir, size, tf, sf].

Required metadata
-----------------
data.t
.data.directions
.data.sizes
.data.tfs
.data.sfs
Optional:
data.activeParams

data.activeParams can contain any subset of:
{'direction', 'size', 'tf', 'sf'}

Orientation analysis
--------------------
Whenever direction is active, the package also performs orientation analysis.
Orientation is defined as mod(direction + 90, 180).
Opposite directions map to the same orientation, so orientation grouping pools
those trials together.

Main outputs
------------
results.pVals
results.isResponsive
results.rVals
results.direction
results.orientation
results.size
results.tf
results.sf
results.meta

Plots saved by the workflow
---------------------------
- odd/even responsiveness scatter plots
- time-course plots for each analyzed parameter, including orientation
- combined 1D tuning figures for all active/derived parameters
- pairwise 2D response matrices for all active/derived parameter pairs
- population summary-distribution figure

New helper function
-------------------
plot_results_distributions(results)
plot_results_distributions(results, savePath)

This plots density-style line distributions across neurons for available
summary statistics such as preferred orientation, OSI, preferred direction,
DSI, preferred size, preferred TF, and preferred SF.
