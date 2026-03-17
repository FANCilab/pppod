Grating MATLAB Workflow
=======================

Files
-----
run_grating_workflow.m              Example entry-point script
demo_synthetic_dataset.m            Synthetic demo dataset and example run
analyze_grating_experiment.m        Main workflow function
validate_grating_inputs.m           Input checks
infer_grating_input_layout.m        Detects supported input layout and sizes
canonicalize_grating_dimensions.m   Restores canonical stimulus dimensions
make_output_folders.m               Creates output folder tree
get_active_parameter_info.m         Metadata for active stimulus parameters
compute_visual_responsiveness.m     Odd-vs-even repeat correlation analysis
plot_even_odd_scatter.m             Scatter plot for responsiveness summary
analyze_parameter_tuning.m          Mean-based analysis for one parameter
reshape_tc_by_parameter.m           Groups time courses by one stimulus parameter
reshape_amp_by_parameter.m          Groups amplitudes by one stimulus parameter
plot_parameter_timecourses.m        One-row time-course plotting with shared y-axis
plot_parameter_tuning.m             Combined active-parameter tuning figure (2 rows)
compute_pairwise_response_matrix.m  Mean response matrix for one parameter pair
plot_pairwise_response_matrices.m   Pairwise 2D response-matrix plotting
compute_preferred_value.m           Preferred stimulus value from mean response
compute_direction_selectivity.m     Computes DSI for direction tuning
find_opposite_direction_idx.m       Finds the sampled direction nearest to 180 deg opposite
save_figure_compat.m                Figure saving helper

Preferred data layout
---------------------
Provide inputs in the following preferred format:

Time-course matrix:
    data.tc  : [nStim, nRep, nTime, nNeuron]

Amplitude matrix:
    data.amp : [nStim, nRep, nNeuron]

where:
    nStim = nDir * nSize * nTf * nSf

The stimulus index must follow the canonical unsqueezed stimulus ordering:
    [direction, size, tf, sf]

That means these reshapes must be valid:
    tcCanon  = reshape(data.tc,  [nDir, nSize, nTf, nSf, nRep, nTime, nNeuron])
    ampCanon = reshape(data.amp, [nDir, nSize, nTf, nSf, nRep, nNeuron])

The workflow uses the parameter vectors:
    data.directions, data.sizes, data.tfs, data.sfs
and optionally:
    data.activeParams

to infer how stimuli should be restored into canonical dimensions.

Backward-compatible layouts
---------------------------
The package also accepts:

1) Full canonical arrays:
       tc  = [dir, size, tf, sf, rep, time, neuron]
       amp = [dir, size, tf, sf, rep, neuron]

2) Legacy squeezed canonical arrays where inactive stimulus dimensions were
   removed before repeat/time/neuron.

Usage
-----
1) Put all .m files on your MATLAB path.
2) Edit run_grating_workflow.m and assign your arrays and vectors.
3) Run:
       run_grating_workflow

Quick test
----------
If you want to verify that the package runs end-to-end before using your own
data, run:
       demo_synthetic_dataset

Outputs
-------
A folder tree is created under targetFolder:
    responsiveness/
    direction_timecourses/
    size_timecourses/
    tf_timecourses/
    sf_timecourses/
    combined_tuning/
    pairwise_response_matrices/

and a MAT-file:
    grating_analysis_results.mat

Plotting updates in this version
--------------------------------
- Time-course plots use one row of subplots for each active parameter.
- Time-course subplots share the same y-axis limits.
- Means are used instead of medians across trials.
- All active 1D tuning plots are combined in one figure per responsive neuron:
    row 1 = single trials + mean
    row 2 = mean +/- standard deviation
- Pairwise 2D response matrices are plotted for all combinations of active
  stimulus parameters.

Notes
-----
- Responsiveness is defined by p < opts.alpha from the odd-vs-even Pearson
  correlation across stimulus conditions.
- Only responsive neurons are plotted for tuning and pairwise matrix figures.
- Preferred values and DSI are computed from mean responses across pooled trials.
- If the exact opposite direction is not present, DSI uses the nearest sampled
  direction to preferred+180 degrees.
- Blank responses are accepted in the input struct but are not used by the
  requested workflow.
