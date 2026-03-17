function results = analyze_grating_experiment(data, targetFolder, opts)
% ANALYZE_GRATING_EXPERIMENT Analyze neuronal responses to grating stimuli.
%
% RESULTS = ANALYZE_GRATING_EXPERIMENT(DATA, TARGETFOLDER, OPTS)
%
% Preferred inputs
%   data.tc          : time-course array of size [nStim, nRep, nTime, nNeuron]
%                      where nStim = nDir * nSize * nTf * nSf.
%                      Stimuli must be linearized according to canonical
%                      stimulus order [dir, size, tf, sf].
%   data.amp         : amplitude array of size [nStim, nRep, nNeuron]
%                      using the same stimulus order.
%   data.blankTc     : blank time-course array (optional, not used here)
%   data.blankAmp    : blank amplitude array   (optional, not used here)
%   data.t           : 1 x nTime vector of time stamps
%   data.directions  : 1 x nDir vector
%   data.sizes       : 1 x nSize vector
%   data.tfs         : 1 x nTf vector
%   data.sfs         : 1 x nSf vector
%   data.activeParams: optional cell array, e.g. {'direction','size'}
%
% Backward compatibility
%   The workflow also accepts:
%   - full canonical arrays [dir size tf sf rep time neuron] and
%     [dir size tf sf rep neuron]
%   - legacy squeezed canonical arrays where inactive stimulus dimensions
%     were removed before repeat/time/neuron.
%
% Outputs
%   results.pVals         : odd-vs-even correlation p-values
%   results.isResponsive  : logical vector of responsive neurons
%   results.rVals         : odd-vs-even correlation coefficients
%   results.direction     : mean tuning outputs for direction
%   results.size          : mean tuning outputs for size
%   results.tf            : mean tuning outputs for temporal frequency
%   results.sf            : mean tuning outputs for spatial frequency
%
% Notes
%   - Internally, data are reshaped to canonical form with dimensions:
%       tc  = [dir, size, tf, sf, rep, time, neuron]
%       amp = [dir, size, tf, sf, rep, neuron]
%   - Blank responses are accepted and saved into results.meta, but are not
%     used by the requested analyses.

    if nargin < 2 || isempty(targetFolder)
        error('You must provide a targetFolder.');
    end

    if nargin < 3
        opts = struct();
    end

    if ~isfield(opts, 'alpha') || isempty(opts.alpha)
        opts.alpha = 0.05;
    end
    if ~isfield(opts, 'saveExt') || isempty(opts.saveExt)
        opts.saveExt = 'png';
    end
    if ~isfield(opts, 'visible') || isempty(opts.visible)
        opts.visible = 'off';
    end

    validate_grating_inputs(data);
    out = make_output_folders(targetFolder);
    canon = canonicalize_grating_dimensions(data);

    [pVals, isResponsive, rVals] = compute_visual_responsiveness( ...
        canon.amp6, out.responsiveness, opts);

    results = struct();
    results.pVals = pVals;
    results.isResponsive = isResponsive;
    results.rVals = rVals;
    results.activeParamNames = canon.activeParams;

    results.direction = local_init_param_result(numel(canon.directions), canon.nNeurons, 'direction');
    results.size = local_init_param_result(numel(canon.sizes), canon.nNeurons, 'size');
    results.tf = local_init_param_result(numel(canon.tfs), canon.nNeurons, 'tf');
    results.sf = local_init_param_result(numel(canon.sfs), canon.nNeurons, 'sf');

    paramInfos = get_active_parameter_info(canon);

    for iParam = 1:numel(paramInfos)
        paramInfo = paramInfos(iParam);
        paramResult = analyze_parameter_tuning( ...
            canon, paramInfo.name, paramInfo.values, paramInfo.dim, ...
            isResponsive, out.(paramInfo.timecourseFolderField), opts);

        results.(paramInfo.resultField) = paramResult;
    end

    for iNeuron = 1:canon.nNeurons
        if ~isResponsive(iNeuron)
            continue
        end

        if ~isempty(paramInfos)
            plot_parameter_tuning(canon, paramInfos, iNeuron, out.combinedTuning, opts);
        end

        if numel(paramInfos) >= 2
            plot_pairwise_response_matrices(canon, paramInfos, iNeuron, out.pairwiseMatrices, opts);
        end
    end

    results.meta = struct();
    results.meta.targetFolder = targetFolder;
    results.meta.alpha = opts.alpha;
    results.meta.saveExt = opts.saveExt;
    results.meta.visible = opts.visible;
    results.meta.activeParams = canon.activeParams;
    results.meta.inputFormat = canon.inputFormat;
    results.meta.nDirections = numel(canon.directions);
    results.meta.nSizes = numel(canon.sizes);
    results.meta.nTFs = numel(canon.tfs);
    results.meta.nSFs = numel(canon.sfs);
    results.meta.nRepeats = canon.nRep;
    results.meta.nTimes = canon.nTime;
    results.meta.nNeurons = canon.nNeurons;
    results.meta.blankDataProvided = struct( ...
        'blankTc', ~isempty(canon.blankTc), ...
        'blankAmp', ~isempty(canon.blankAmp));

    save(fullfile(targetFolder, 'grating_analysis_results.mat'), 'results', '-v7.3');
end

function out = local_init_param_result(nLevels, nNeurons, paramName)
% LOCAL_INIT_PARAM_RESULT Initialize results struct for one parameter.

    out = struct();
    out.meanResponses = nan(nLevels, nNeurons);
    out.stdResponses = nan(nLevels, nNeurons);
    out.preferredValue = nan(1, nNeurons);

    switch lower(paramName)
        case 'direction'
            out.preferredDirection = nan(1, nNeurons);
            out.DSI = nan(1, nNeurons);
        case 'size'
            out.preferredSize = nan(1, nNeurons);
        case 'tf'
            out.preferredTF = nan(1, nNeurons);
        case 'sf'
            out.preferredSF = nan(1, nNeurons);
    end
end
