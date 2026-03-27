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
% Outputs include derived orientation analysis whenever direction is active.

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
    [pVals, isResponsive, rVals] = compute_visual_responsiveness(canon.amp6, out.responsiveness, opts);
    canon.isResponsive = isResponsive;

    results = struct();
    results.pVals = pVals;
    results.isResponsive = isResponsive;
    results.rVals = rVals;

    paramInfos = get_active_parameter_info(canon);
    results.activeParamNames = {paramInfos.name};

    % initialize supported outputs
    results.direction = local_init_param_result(numel(canon.directions), canon.nNeurons, 'direction');
    results.orientation = local_init_param_result(numel(canon.orientations), canon.nNeurons, 'orientation');
    results.size = local_init_param_result(numel(canon.sizes), canon.nNeurons, 'size');
    results.tf = local_init_param_result(numel(canon.tfs), canon.nNeurons, 'tf');
    results.sf = local_init_param_result(numel(canon.sfs), canon.nNeurons, 'sf');
    
    % best 1D tuning for each active par
    results.oneDTuning = struct();
    for iPar = 1:numel(canon.activeParams)
        pName = canon.activeParams{iPar};
        results.oneDTuning.(pName) = compute_best_condition_1d_tuning(canon, pName, opts, targetFolder);
    end

    % Pairwise response analyses for all active-parameter pairs
    results.pairwise = struct();
    if isfield(canon, 'activeParams') && numel(canon.activeParams) >= 2
        for iPar = 1:numel(canon.activeParams)-1
            for jPar = iPar+1:numel(canon.activeParams)
                p1 = canon.activeParams{iPar};
                p2 = canon.activeParams{jPar};

                % pairFolder = fullfile(targetFolder, [p1 '_vs_' p2]);
                pairMat = compute_plot_pairwise_responses(canon, p1, p2, targetFolder, opts);

                results.pairwise.([p1 '_vs_' p2]) = pairMat;
                results.([p1 '_vs_' p2]) = pairMat;
            end
        end
    end

    for iParam = 1:numel(paramInfos)
        paramInfo = paramInfos(iParam);
        paramResult = analyze_parameter_tuning(canon, paramInfo, isResponsive, out.(paramInfo.timecourseFolderField), opts);
        results.(paramInfo.resultField) = paramResult;
    end

    % for iNeuron = 1:canon.nNeurons
    %     if ~isResponsive(iNeuron)
    %         continue
    %     end
    % 
    %     if ~isempty(paramInfos)
    %         plot_parameter_tuning(canon, paramInfos, iNeuron, out.combinedTuning, opts);
    %     end
    % 
    %     if numel(paramInfos) >= 2
    %         plot_pairwise_response_matrices(canon, paramInfos, iNeuron, out.pairwiseMatrices, opts);
    %     end
    % end

    
    results.meta = struct();
    results.meta.targetFolder = targetFolder;
    results.meta.alpha = opts.alpha;
    results.meta.saveExt = opts.saveExt;
    results.meta.visible = opts.visible;
    results.meta.activeParams = canon.activeParams;
    results.meta.analysisParams = {paramInfos.name};
    results.meta.inputFormat = canon.inputFormat;
    results.meta.directions = canon.directions;
    results.meta.orientations = canon.orientations;
    results.meta.sizes = canon.sizes;
    results.meta.tfs = canon.tfs;
    results.meta.sfs = canon.sfs;
    results.meta.nDirections = numel(canon.directions);
    results.meta.nOrientations = numel(canon.orientations);
    results.meta.nSizes = numel(canon.sizes);
    results.meta.nTFs = numel(canon.tfs);
    results.meta.nSFs = numel(canon.sfs);
    results.meta.nRepeats = canon.nRep;
    results.meta.nTimes = canon.nTime;
    results.meta.nNeurons = canon.nNeurons;
    results.meta.blankDataProvided = struct( ...
        'blankTc', ~isempty(canon.blankTc), ...
        'blankAmp', ~isempty(canon.blankAmp));

    distPath = out.summaryDistributions;
    figDist = plot_results_distributions(results, distPath);
    if ishghandle(figDist)
        close(figDist);
    end

    save(fullfile(targetFolder, 'grating_analysis_results.mat'), 'results', '-v7.3');
end

function out = local_init_param_result(nLevels, nNeurons, paramName)
    out = struct();
    out.meanResponses = nan(nLevels, nNeurons);
    out.stdResponses = nan(nLevels, nNeurons);
    out.preferredValue = nan(1, nNeurons);

    switch lower(paramName)
        case 'direction'
            out.preferredDirection = nan(1, nNeurons);
            out.DSI = nan(1, nNeurons);
        case 'orientation'
            out.preferredOrientation = nan(1, nNeurons);
            out.OSI = nan(1, nNeurons);
        case 'size'
            out.preferredSize = nan(1, nNeurons);
        case 'tf'
            out.preferredTF = nan(1, nNeurons);
        case 'sf'
            out.preferredSF = nan(1, nNeurons);
    end
end
