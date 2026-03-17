function [figDistributions, figTunings] = plot_results_distributions(results, saveFolder, opts)
%PLOT_RESULTS_DISTRIBUTIONS Plot population summary distributions and tuning curves.
%
% [figDistributions, figTunings] = plot_results_distributions(results, saveFolder, opts)
%
% This function operates on the output structure "results" from the grating
% analysis pipeline.
%
% Figure 1
%   One subplot per summary metric:
%   - preferredValues for each active parameter listed in results.activeParNames
%   - results.direction.DSI (if present)
%   - results.orientation.OSI (if present)
%   NaNs are excluded. Histograms are normalised by the number of valid
%   neurons so the y-axis shows density per bin.
%
% Figure 2
%   One subplot per active parameter:
%   - thin lines: one mean tuning curve per neuron
%   - thick line : across-neuron average tuning curve
%   X-values are taken from results.meta.<parameter values>.
%
% Inputs
%   results    : results structure
%   saveFolder : optional folder where figures will be saved
%   opts       : optional struct with fields
%       .visible : 'on' or 'off' (default 'off')
%       .saveExt : image extension, e.g. 'png' (default 'png')
%
% Outputs
%   figDistributions : handle to distribution figure (or [])
%   figTunings       : handle to tuning-curve figure (or [])

if nargin < 2 || isempty(saveFolder)
    saveFolder = '';
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'visible') || isempty(opts.visible)
    opts.visible = 'off';
end
if ~isfield(opts, 'saveExt') || isempty(opts.saveExt)
    opts.saveExt = 'png';
end

figDistributions = [];
figTunings = [];

activeNames = local_get_active_param_names(results);

%% ---------------- Figure 1: distributions ----------------
distInfos = local_collect_distribution_infos(results, activeNames);

if ~isempty(distInfos)
    nPlots = numel(distInfos);
    figDistributions = local_make_figure(opts.visible);
    tl1 = tiledlayout(figDistributions, local_n_rows(nPlots), local_n_cols(nPlots), ...
        'Padding', 'compact', 'TileSpacing', 'compact');

    for iPlot = 1:nPlots
        ax = nexttile(tl1); %#ok<LAXES>
        hold(ax, 'on');

        values = distInfos(iPlot).values(:);
        values = values(isfinite(values));

        if isempty(values)
            title(ax, distInfos(iPlot).title, 'Interpreter', 'none');
            xlabel(ax, distInfos(iPlot).xLabel, 'Interpreter', 'none');
            ylabel(ax, 'Density');
            box(ax, 'off');
            continue;
        end

        [xLine, yLine] = local_density_histogram(values, distInfos(iPlot).xValues);
        stairs(ax, xLine, yLine, 'LineWidth', 1, 'Color','k');

        xlabel(ax, distInfos(iPlot).xLabel, 'Interpreter', 'none');
        ylabel(ax, 'Density');
        ylim([0 0.6])
        title(ax, distInfos(iPlot).title, 'Interpreter', 'none');
        box(ax, 'off');

        if ~isempty(distInfos(iPlot).xValues)
            xlim(ax, [min(distInfos(iPlot).xValues) max(distInfos(iPlot).xValues)]);
        end
    end

    sgtitle(figDistributions, 'Population summary distributions');

    if ~isempty(saveFolder)
        local_safe_export(figDistributions, fullfile(saveFolder, ['population_summary_distributions.' opts.saveExt]));
    end
end

%% ---------------- Figure 2: population tuning curves ----------------
tuningInfos = local_collect_tuning_infos(results, activeNames);

if ~isempty(tuningInfos)
    nPlots = numel(tuningInfos);
    figTunings = local_make_figure(opts.visible);
    tl2=tiledlayout(figTunings, 1, nPlots, 'Padding', 'compact', 'TileSpacing', 'compact');

    for iPlot = 1:nPlots
        ax = nexttile(tl2); %#ok<LAXES>
        hold(ax, 'on');

        xValues = tuningInfos(iPlot).xValues(:).';
        Y = tuningInfos(iPlot).Y;

        if isempty(xValues) || isempty(Y)
            title(ax, tuningInfos(iPlot).title, 'Interpreter', 'none');
            xlabel(ax, tuningInfos(iPlot).xLabel, 'Interpreter', 'none');
            ylabel(ax, 'Mean response');
            box(ax, 'off');
            continue;
        end

        % Ensure Y is nValues x nNeurons
        if size(Y,1) ~= numel(xValues) && size(Y,2) == numel(xValues)
            Y = Y.';
        end

        % Sort x-axis for cleaner plotting
        [xValues, sortIdx] = sort(xValues, 'ascend');
        Y = Y(sortIdx, :);

        % Plot one thin line per neuron
        for iNeuron = 1:size(Y, 2)
            y = Y(:, iNeuron);
            if any(isfinite(y))
                plot(ax, xValues, y, 'LineWidth', 0.5, 'Color',[0.7 0.7 0.7]);
            end
        end

        % Plot population average as thick line
        meanY = mean(Y, 2, 'omitnan');
        plot(ax, xValues, meanY, 'LineWidth', 1, 'Color','k');

        xlabel(ax, tuningInfos(iPlot).xLabel, 'Interpreter', 'none');
        ylabel(ax, 'Mean response');
        title(ax, tuningInfos(iPlot).title, 'Interpreter', 'none');
        xlim(ax, [min(xValues) max(xValues)]);
        box(ax, 'off');
    end

    sgtitle(figTunings, 'Population tuning curves');

    if ~isempty(saveFolder)
        local_safe_export(figTunings, fullfile(saveFolder, ['population_tuning_curves.' opts.saveExt]));
    end
end

end

%% =====================================================================
function activeNames = local_get_active_param_names(results)

activeNames = {};

if isfield(results, 'activeParNames') && ~isempty(results.activeParNames)
    activeNames = local_cellstr(results.activeParNames);
elseif isfield(results, 'activeParamNames') && ~isempty(results.activeParamNames)
    activeNames = local_cellstr(results.activeParamNames);
elseif isfield(results, 'meta') && isfield(results.meta, 'activeParams') && ~isempty(results.meta.activeParams)
    activeNames = local_cellstr(results.meta.activeParams);
end

% Keep only parameter names that actually exist as fields in results
keep = false(size(activeNames));
for i = 1:numel(activeNames)
    keep(i) = isfield(results, activeNames{i});
end
activeNames = activeNames(keep);

end

function infos = local_collect_distribution_infos(results, activeNames)

infos = struct('title', {}, 'xLabel', {}, 'values', {}, 'xValues', {});

% Preferred values for each active parameter
for i = 1:numel(activeNames)
    parName = activeNames{i};
    if ~isfield(results, parName) || ~isstruct(results.(parName))
        continue;
    end

    parStruct = results.(parName);

    if isfield(parStruct, 'preferredValues')
        vals = parStruct.preferredValues;
    elseif isfield(parStruct, 'preferredValue')
        vals = parStruct.preferredValue;
    else
        vals = [];
    end

    if isempty(vals)
        continue;
    end

    info.title = ['Preferred ' local_pretty_name(parName)];
    info.xLabel = local_pretty_name(parName);
    info.values = vals(:);
    info.xValues = local_get_parameter_values(results, parName);
    infos(end+1) = info; %#ok<AGROW>
end

% Add DSI and OSI if present
if isfield(results, 'direction') && isstruct(results.direction) && isfield(results.direction, 'DSI')
    info.title = 'DSI';
    info.xLabel = 'DSI';
    info.values = results.direction.DSI(:);
    info.xValues = 0:0.1:1; %[]; 
    infos(end+1) = info; %#ok<AGROW>
end

if isfield(results, 'orientation') && isstruct(results.orientation) && isfield(results.orientation, 'OSI')
    info.title = 'OSI';
    info.xLabel = 'OSI';
    info.values = results.orientation.OSI(:);
    info.xValues = 0:0.1:1; %[]; 
    infos(end+1) = info; %#ok<AGROW>
end

end

function infos = local_collect_tuning_infos(results, activeNames)

infos = struct('title', {}, 'xLabel', {}, 'xValues', {}, 'Y', {});

for i = 1:numel(activeNames)
    parName = activeNames{i};
    if ~isfield(results, parName) || ~isstruct(results.(parName))
        continue;
    end

    parStruct = results.(parName);

    if ~isfield(parStruct, 'meanResponses') || isempty(parStruct.meanResponses)
        continue;
    end

    xValues = local_get_parameter_values(results, parName);
    Y = parStruct.meanResponses;

    if isempty(xValues) || isempty(Y)
        continue;
    end


    switch parName
        case 'direction'
            prefPos = 1+ parStruct.preferredValue/median(diff(xValues));
            prefPos(isnan(prefPos)) =0;
            central = 180/median(diff(xValues)) + 1;
            for iN =1:size(Y, 2)
            Y(:, iN) = circshift(Y(:, iN),central-prefPos(iN), 1);
            % Y = circshift(Y,central-prefPos(iN), 2);

            end
            xValues = xValues - 180;
        case 'orientation'
                 prefPos = 1+ parStruct.preferredValue/median(diff(xValues));
            prefPos(isnan(prefPos)) =0;
            central = 90/median(diff(xValues)) + 1;
            for iN =1:size(Y, 2)
            Y(:, iN) = circshift(Y(:, iN),central-prefPos(iN), 1);
            % Y = circshift(Y,central-prefPos(iN), 2);

            end
            xValues = xValues - 90;
            
    end

    info.title = local_pretty_name(parName);
    info.xLabel = local_pretty_name(parName);
    info.xValues = xValues(:).';
    info.Y = Y;
    infos(end+1) = info; %#ok<AGROW>
end

end

function values = local_get_parameter_values(results, parName)

values = [];

if ~isfield(results, 'meta') || ~isstruct(results.meta)
    return
end

meta = results.meta;

switch lower(parName)
    case 'direction'
        candidateFields = {'direction', 'directions'};
    case 'orientation'
        candidateFields = {'orientation', 'orientations'};
    case 'size'
        candidateFields = {'size', 'sizes'};
    case 'tf'
        candidateFields = {'tf', 'tfs'};
    case 'sf'
        candidateFields = {'sf', 'sfs'};
    otherwise
        candidateFields = {parName, [parName 's']};
end

for i = 1:numel(candidateFields)
    if isfield(meta, candidateFields{i})
        values = meta.(candidateFields{i});
        return
    end
end

end

function [xLine, yLine] = local_density_histogram(values, xSupport)

values = values(:);
values = values(isfinite(values));

if isempty(values)
    xLine = [];
    yLine = [];
    return
end

if nargin < 2
    xSupport = [];
end

% For discrete preferred values, use bins centered on the known parameter values
if ~isempty(xSupport)
    xSupport = xSupport(:).';
    xSupport = sort(unique(xSupport));

    if numel(xSupport) == 1
        xLine = xSupport;
        yLine = 1;
        return
    end

    dx = diff(xSupport);
    edges = [xSupport(1) - dx(1)/2, (xSupport(1:end-1) + xSupport(2:end))/2, xSupport(end) + dx(end)/2];
    counts = histcounts(values, edges);
    yLine = counts / sum(counts);
    xLine = xSupport;
    return
end

% Otherwise use automatic histogram bins
nBins = min(20, max(5, ceil(sqrt(numel(values)))));
[counts, edges] = histcounts(values, nBins);
centers = edges(1:end-1) + diff(edges)/2;

if sum(counts) > 0
    counts = counts / sum(counts);
end

xLine = centers;
yLine = counts;

end

function out = local_cellstr(x)

if ischar(x)
    out = {x};
elseif isstring(x)
    out = cellstr(x(:));
elseif iscell(x)
    out = cell(size(x));
    for i = 1:numel(x)
        if isstring(x{i}) || ischar(x{i})
            out{i} = char(x{i});
        else
            out{i} = char(string(x{i}));
        end
    end
    out = out(:).';
else
    out = cellstr(string(x(:)));
end

out = out(:).';
for i = 1:numel(out)
    out{i} = strtrim(out{i});
end

end

function txt = local_pretty_name(name)

switch lower(char(name))
    case 'sf'
        txt = 'SF';
    case 'tf'
        txt = 'TF';
    case 'osi'
        txt = 'OSI';
    case 'dsi'
        txt = 'DSI';
    otherwise
        txt = char(name);
        if ~isempty(txt)
            txt(1) = upper(txt(1));
        end
    end

end

function n = local_n_rows(nPlots)
n = ceil(sqrt(nPlots));
end

function n = local_n_cols(nPlots)
n = ceil(nPlots / local_n_rows(nPlots));
end

function fig = local_make_figure(visibleState)

if nargin < 1 || isempty(visibleState)
    visibleState = 'off';
end

oldDefault = get(groot, 'DefaultFigureVisible');
cleanupObj = onCleanup(@() set(groot, 'DefaultFigureVisible', oldDefault)); %#ok<NASGU>
set(groot, 'DefaultFigureVisible', visibleState);
fig = figure('Visible', visibleState);

end

function local_safe_export(fig, filePath)

[folderPath, ~, ~] = fileparts(filePath);
if ~isempty(folderPath) && ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

try
    exportgraphics(fig, filePath, 'Resolution', 150);
catch
    saveas(fig, filePath);
end

end
