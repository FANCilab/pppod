function oneDTuning = compute_best_condition_1d_tuning(canon, paramName, opts, targetFolder)
%COMPUTE_BEST_CONDITION_1D_TUNING Compute 1D tuning after fixing other active
%parameters at their best-response values.
%
% oneDTuning = compute_best_condition_1d_tuning(canon, paramName, opts, targetFolder)
%
% Inputs
%   canon        struct with canonical fields:
%                tc7  [dir size tf sf rep time neuron]
%                amp6 [dir size tf sf rep neuron]
%                t, directions, sizes, tfs, sfs, activeParams
%   paramName    active parameter of interest: 'direction','size','tf','sf'
%   opts         plotting options struct (optional). Supported fields:
%                  .visible  = 'on'/'off'   (default 'off')
%                  .saveExt  = file extension (default 'png')
%                  .lineWidthThin = scalar   (default 0.5)
%                  .lineWidthThick = scalar  (default 2)
%   targetFolder folder where figures are saved (optional)
%
% Output
%   oneDTuning   struct with fields:
%                  .parameterName
%                  .parameterValues
%                  .meanResponses   [nValues x nNeurons]
%                  .bestValues      struct with one 1xnNeurons vector per active param
%                  .bestIndices     struct with one 1xnNeurons vector per active param
%                  .bestCombinationResponses [1 x nNeurons]
%
% Notes
%   For each neuron, the function:
%   1) averages amp6 across repeats
%   2) finds the active-parameter combination with the largest average response
%   3) fixes all other active parameters at those best values
%   4) computes a 1D tuning curve for the selected parameter using only the
%      trials matching the fixed best values of the other active parameters
%   5) saves a tuning-curve plot and time-course plot for that neuron
%
% The time-course plot uses one row of subplots, one subplot per value of
% paramName, with thin single-trial traces and a thick mean trace.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if nargin < 4 || isempty(targetFolder)
        targetFolder = pwd;
    end

    opts = local_apply_defaults(opts);
    local_validate_inputs(canon, paramName);

    paramInfos = local_get_param_infos(canon);
    activeParamNames = local_get_active_param_names(canon);
    activeMask = ismember({paramInfos.name}, activeParamNames);

    targetInfo = paramInfos(strcmp({paramInfos.name}, paramName));
    if isempty(targetInfo)
        error('Unknown parameter name: %s', paramName);
    end
    targetInfo = targetInfo(1);

    nNeurons = size(canon.amp6, 6);
    nValues = numel(targetInfo.values);

    oneDTuning = struct();
    oneDTuning.parameterName = paramName;
    oneDTuning.parameterValues = targetInfo.values(:).';
    oneDTuning.meanResponses = nan(nValues, nNeurons);
    oneDTuning.bestCombinationResponses = nan(1, nNeurons);
    oneDTuning.bestValues = struct();
    oneDTuning.bestIndices = struct();

    for iInfo = 1:numel(paramInfos)
        if activeMask(iInfo)
            oneDTuning.bestValues.(paramInfos(iInfo).name) = nan(1, nNeurons);
            oneDTuning.bestIndices.(paramInfos(iInfo).name) = nan(1, nNeurons);
        end
    end

    tuningFolder = fullfile(targetFolder, [paramName '_best_condition_tuning']);
    tcFolder = fullfile(targetFolder, [paramName '_best_condition_timecourses']);
    local_mkdir_if_needed(tuningFolder);
    local_mkdir_if_needed(tcFolder);

    ampMean = mean(canon.amp6, 5, 'omitnan'); % [dir size tf sf neuron]

    for iNeuron = 1:nNeurons
        ampMeanNeuron = ampMean(:,:,:,:,iNeuron);

        % Find best active-parameter combination for this neuron.
        ampForBest = ampMeanNeuron;
        inactiveDims = find(~activeMask);
        for iDim = 1:numel(inactiveDims)
            ampForBest = mean(ampForBest, inactiveDims(iDim), 'omitnan');
        end
        ampForBest = squeeze(ampForBest);

        if all(isnan(ampForBest(:)))
            continue
        end

        [maxResp, linIdx] = max(ampForBest(:));
        oneDTuning.bestCombinationResponses(iNeuron) = maxResp;

        activeSizes = [paramInfos(activeMask).nValues];
        activeSubs = cell(1, numel(activeSizes));
        [activeSubs{:}] = ind2sub(activeSizes, linIdx);
        activeSubs = cellfun(@double, activeSubs);

        % Build full index set into canonical [dir size tf sf] dims.
        fixedIndices = ones(1, numel(paramInfos));
        activeCounter = 0;
        for iInfo = 1:numel(paramInfos)
            if activeMask(iInfo)
                activeCounter = activeCounter + 1;
                fixedIndices(iInfo) = activeSubs(activeCounter);
                oneDTuning.bestIndices.(paramInfos(iInfo).name)(iNeuron) = fixedIndices(iInfo);
                oneDTuning.bestValues.(paramInfos(iInfo).name)(iNeuron) = paramInfos(iInfo).values(fixedIndices(iInfo));
            end
        end

        % For the selected parameter, vary its full axis while fixing the other
        % active parameters at best values.
        ampByValue = local_extract_amp_by_value(canon.amp6, fixedIndices, targetInfo.dim, iNeuron, targetInfo.nValues);
        oneDTuning.meanResponses(:, iNeuron) = mean(ampByValue, 2, 'omitnan');

        tcByValue = local_extract_tc_by_value(canon.tc7, fixedIndices, targetInfo.dim, iNeuron, targetInfo.nValues);

        if ~canon.isResponsive(iNeuron)
            continue
        end
        local_plot_tuning(paramName, oneDTuning.parameterValues, ampByValue, iNeuron, tuningFolder, opts);
        local_plot_timecourses(paramName, oneDTuning.parameterValues, tcByValue, canon.t, iNeuron, tcFolder, opts);
    end
end

function opts = local_apply_defaults(opts)
    if ~isfield(opts, 'visible') || isempty(opts.visible)
        opts.visible = 'off';
    end
    if ~isfield(opts, 'saveExt') || isempty(opts.saveExt)
        opts.saveExt = 'png';
    end
    if ~isfield(opts, 'lineWidthThin') || isempty(opts.lineWidthThin)
        opts.lineWidthThin = 0.5;
    end
    if ~isfield(opts, 'lineWidthThick') || isempty(opts.lineWidthThick)
        opts.lineWidthThick = 2;
    end
end

function local_validate_inputs(canon, paramName)
    requiredFields = {'tc7','amp6','t','directions','sizes','tfs','sfs'};
    for iField = 1:numel(requiredFields)
        if ~isfield(canon, requiredFields{iField})
            error('canon is missing required field "%s".', requiredFields{iField});
        end
    end

    validParams = {'direction','size','tf','sf'};
    if ~ischar(paramName) && ~isstring(paramName)
        error('paramName must be a character vector or string.');
    end
    paramName = char(paramName);
    if ~ismember(paramName, validParams)
        error('paramName must be one of: direction, size, tf, sf.');
    end

    paramInfos = local_get_param_infos(canon);
    targetInfo = paramInfos(strcmp({paramInfos.name}, paramName));
    if isempty(targetInfo) || targetInfo.nValues < 2
        error('The requested parameter "%s" is not variable in this canon structure.', paramName);
    end

    activeNames = local_get_active_param_names(canon);
    if ~ismember(paramName, activeNames)
        error('The requested parameter "%s" is not listed in canon.activeParams.', paramName);
    end
end

function paramInfos = local_get_param_infos(canon)
    paramInfos = struct( ...
        'name', {'direction','size','tf','sf'}, ...
        'values', {local_row(canon.directions), local_row(canon.sizes), local_row(canon.tfs), local_row(canon.sfs)}, ...
        'dim', {1, 2, 3, 4}, ...
        'nValues', {numel(canon.directions), numel(canon.sizes), numel(canon.tfs), numel(canon.sfs)} );
end

function names = local_get_active_param_names(canon)
    if isfield(canon, 'activeParams') && ~isempty(canon.activeParams)
        raw = canon.activeParams;
        if iscell(raw)
            names = cellfun(@char, raw, 'UniformOutput', false);
        else
            names = cell(1, numel(raw));
            for i = 1:numel(raw)
                names{i} = char(raw{i});
            end
        end
    else
        infos = local_get_param_infos(canon);
        names = {infos([infos.nValues] > 1).name};
    end
    names = cellfun(@char, names, 'UniformOutput', false);
end

function ampByValue = local_extract_amp_by_value(amp6, fixedIndices, varyingDim, iNeuron, nParam)
    idx = {fixedIndices(1), fixedIndices(2), fixedIndices(3), fixedIndices(4), ':', iNeuron};
    idx{varyingDim} = ':';
    A = amp6(idx{:}); % [param x repeats] with possible singleton dims
    A = squeeze(A);
    if isvector(A)
        A = A(:);
    end
    if size(A, 1) ~= nParam && size(A, 2) == nParam
        A = A.';
    end
    if ndims(A) ~= 2
        A = reshape(A, size(A, 1), []);
    end
    ampByValue = A;
end

function tcByValue = local_extract_tc_by_value(tc7, fixedIndices, varyingDim, iNeuron, nParam)
    idx = {fixedIndices(1), fixedIndices(2), fixedIndices(3), fixedIndices(4), ':', ':', iNeuron};
    idx{varyingDim} = ':';
    A = tc7(idx{:}); % [param x repeats x time] with possible singleton dims
    A = squeeze(A);

    if ismatrix(A)
        % Could be [param x time] when nRep==1, or [rep x time] when nParam==1.
        nTime = size(tc7, 6);
        nParam = nParam;
        if size(A, 2) == nTime && size(A, 1) == nParam
            A = reshape(A, [nParam, 1, nTime]);
        elseif size(A, 1) == nTime && size(A, 2) == nParam
            A = permute(reshape(A, [nTime, nParam]), [2 3 1]);
        else
            A = reshape(A, [nParam, 1, nTime]);
        end
    end

    if ndims(A) ~= 3
        sz = size(A);
        A = reshape(A, [sz(1), sz(2), prod(sz(3:end))]);
    end

    % Enforce [param x repeat x time]
    nTime = size(tc7, 6);
    if size(A, 3) ~= nTime && size(A, 2) == nTime
        A = permute(A, [1 3 2]);
    end
    tcByValue = A;
end

function local_plot_tuning(paramName, values, ampByValue, iNeuron, folder, opts)
    fig = local_make_figure(opts.visible);
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on');

    nValues = numel(values);
    for iValue = 1:nValues
        y = ampByValue(iValue, :);
        x = repmat(values(iValue), size(y));
        plot(ax, x, y, 'o', 'MarkerFaceColor', 'none');
    end
    plot(ax, values, mean(ampByValue, 2, 'omitnan'), '-', 'LineWidth', opts.lineWidthThick);

    xlabel(ax, local_pretty_label(paramName));
    ylabel(ax, 'Response amplitude');
    title(ax, sprintf('Neuron %d - best-condition %s tuning', iNeuron, local_pretty_label(paramName)));
    local_set_x_axis(ax, values);
    grid(ax, 'on');

    local_save_figure(fig, fullfile(folder, sprintf('neuron_%04d_best_condition_tuning.%s', iNeuron, opts.saveExt)), opts);
end

function local_plot_timecourses(paramName, values, tcByValue, t, iNeuron, folder, opts)
    nValues = numel(values);
    fig = local_make_figure(opts.visible);
    tl = tiledlayout(fig, 1, nValues, 'Padding', 'compact', 'TileSpacing', 'compact');

    yMin = inf;
    yMax = -inf;
    for iValue = 1:nValues
        trials = squeeze(tcByValue(iValue, :, :));
        if isvector(trials)
            trials = trials(:).';
        end
        yMin = min([yMin; trials(:)], [], 'omitnan');
        yMax = max([yMax; trials(:)], [], 'omitnan');
    end
    if ~isfinite(yMin) || ~isfinite(yMax)
        yMin = 0; yMax = 1;
    elseif yMin == yMax
        yMin = yMin - 1; yMax = yMax + 1;
    end

    for iValue = 1:nValues
        ax = nexttile(tl);
        hold(ax, 'on');
        trials = squeeze(tcByValue(iValue, :, :)); % [repeats x time]
        if isvector(trials)
            trials = trials(:).';
        end
        plot(ax, local_row(t), trials.', 'LineWidth', opts.lineWidthThin);
        plot(ax, local_row(t), mean(trials, 1, 'omitnan'), 'LineWidth', opts.lineWidthThick);
        xline(ax, 0, '--');
        title(ax, sprintf('%.4g', values(iValue)));
        xlabel(ax, 'Time from stimulus onset');
        ylabel(ax, 'Response');
        ylim(ax, [yMin yMax]);
        grid(ax, 'on');
    end

    sgtitle(tl, sprintf('Neuron %d - best-condition %s timecourses', iNeuron, local_pretty_label(paramName)));
    local_save_figure(fig, fullfile(folder, sprintf('neuron_%04d_best_condition_timecourses.%s', iNeuron, opts.saveExt)), opts);
end

function local_set_x_axis(ax, values)
    values = local_row(values);
    values = sort(values);
    if numel(values) > 1 && all(diff(values) > 0)
        set(ax, 'XTick', values, 'XLim', [values(1) values(end)]);
    else
        set(ax, 'XTick', values);
    end
end

function fig = local_make_figure(visibleSetting)
    if nargin < 1 || isempty(visibleSetting)
        visibleSetting = 'off';
    end
    fig = figure('Visible', visibleSetting, 'Color', 'w');
end

function local_save_figure(fig, filePath, opts)
    [folder, ~, ~] = fileparts(filePath);
    local_mkdir_if_needed(folder);
    try
        exportgraphics(fig, filePath, 'Resolution', 150);
    catch
        saveas(fig, filePath);
    end
    if strcmpi(opts.visible, 'off') && isgraphics(fig)
        close(fig);
    end
end

function local_mkdir_if_needed(folder)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end

function x = local_row(x)
    x = double(x(:)).';
end

function label = local_pretty_label(name)
    switch lower(char(name))
        case 'direction'
            label = 'Direction';
        case 'size'
            label = 'Size';
        case 'tf'
            label = 'Temporal frequency';
        case 'sf'
            label = 'Spatial frequency';
        case 'parameter'
            label = 'Parameter value';
        otherwise
            label = char(name);
    end
end

