function pairwiseMeanResponses = compute_plot_pairwise_responses(canon, paramName1, paramName2, targetFolder, opts)
%COMPUTE_PLOT_PAIRWISE_RESPONSES Compute and plot pairwise response summaries.
%
% pairwiseMeanResponses = compute_plot_pairwise_responses(canon, paramName1, paramName2, targetFolder)
% pairwiseMeanResponses = compute_plot_pairwise_responses(canon, paramName1, paramName2, targetFolder, opts)
%
% Inputs
%   canon       : canonical struct created inside analyze_grating_experiment
%   paramName1  : name of first parameter ('direction','size','tf','sf')
%   paramName2  : name of second parameter ('direction','size','tf','sf')
%   targetFolder: folder where plots are saved
%   opts        : optional struct with fields:
%                   .visible  -> 'on' or 'off' (default 'off')
%                   .saveExt  -> file extension, e.g. 'png' (default 'png')
%
% Output
%   pairwiseMeanResponses : [nParam1 x nParam2 x nNeurons] array of mean
%                           amplitude responses, averaged across repeats
%                           and any remaining stimulus dimensions.
%
% Notes
% - The output array follows the order [paramName1 x paramName2 x neuron].
% - The time-course figure for each neuron is displayed as a matrix of
%   subplots with rows corresponding to paramName2 and columns to paramName1.
% - Within each subplot, thin lines are single trials and the thick line is
%   the across-trial mean time course.

if nargin < 5 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'visible') || isempty(opts.visible)
    opts.visible = 'off';
end
if ~isfield(opts, 'saveExt') || isempty(opts.saveExt)
    opts.saveExt = 'png';
end

paramName1 = char(string(paramName1));
paramName2 = char(string(paramName2));

if strcmpi(paramName1, paramName2)
    error('paramName1 and paramName2 must be different.');
end

[dim1, values1, label1] = local_get_param_info(canon, paramName1);
[dim2, values2, label2] = local_get_param_info(canon, paramName2);

n1 = numel(values1);
n2 = numel(values2);
nNeurons = double(canon.nNeurons);

pairName = sprintf('%s_vs_%s', local_safe_name(paramName1), local_safe_name(paramName2));
outFolder = fullfile(targetFolder, pairName);
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

pairwiseMeanResponses = nan(n1, n2, nNeurons);

[values1Plot, order1] = sort(values1(:).');
[values2Plot, order2] = sort(values2(:).');

for iNeuron = 1:nNeurons
    ampMat = local_compute_pairwise_amp_matrix(canon.amp6, iNeuron, dim1, dim2);
    ampMatPlot = ampMat(order1, order2);
    pairwiseMeanResponses(:,:,iNeuron) = ampMat;

    if ~canon.isResponsive(iNeuron)
        continue
    end
    figAmp = local_make_figure(opts.visible);
    axAmp = axes('Parent', figAmp);
    imagesc(axAmp, values1Plot, values2Plot, ampMatPlot.');
    axis(axAmp, 'tight');
    set(axAmp, 'YDir', 'normal');
    xlabel(axAmp, label1, 'Interpreter', 'none');
    ylabel(axAmp, label2, 'Interpreter', 'none');
    title(axAmp, sprintf('Neuron %d: %s vs %s mean response', iNeuron, label1, label2), ...
        'Interpreter', 'none');
    colorbar(axAmp);

    local_save_figure(figAmp, fullfile(outFolder, ...
        sprintf('neuron_%04d_%s_matrix.%s', iNeuron, pairName, opts.saveExt)), opts);
    close(figAmp);

    tcCell = local_collect_pairwise_timecourses(canon.tc7, iNeuron, dim1, dim2, order1, order2);

    figTc = local_make_figure(opts.visible);

    yMin = inf;
    yMax = -inf;
    for i2 = 1:n2
        for i1 = 1:n1
            trials = tcCell{i2, i1};
            if ~isempty(trials)
                yMin = min(yMin, min(trials(:)));
                yMax = max(yMax, max(trials(:)));
            end
        end
    end
    if ~isfinite(yMin) || ~isfinite(yMax) || yMin == yMax
        yMin = -1;
        yMax = 1;
    end

    for i2 = 1:n2
        for i1 = 1:n1
            ax = subplot(n2, n1, (i2-1)*n1 + i1, 'Parent', figTc);
            trials = tcCell{i2, i1};

            if ~isempty(trials)
                plot(ax, canon.t(:), trials.', 'LineWidth', 0.5, 'Color',[0.7 0.7 0.7]);
                hold(ax, 'on');
                meanTrace = mean(trials, 1, 'omitnan');
                plot(ax, canon.t(:), meanTrace(:), 'LineWidth', 2.0, 'Color','k');
            end

            xline(ax, 0, '--');
            ylim(ax, [yMin, yMax]);

            if i2 == n2
                xlabel(ax, sprintf('%s = %g', label1, values1Plot(i1)), 'Interpreter', 'none');
            else
                set(ax, 'XTickLabel', []);
            end

            if i1 == 1
                ylabel(ax, sprintf('%s = %g', label2, values2Plot(i2)), 'Interpreter', 'none');
            else
                set(ax, 'YTickLabel', []);
            end

            title(ax, sprintf('%g / %g', values1Plot(i1), values2Plot(i2)));
            box(ax, 'off');

        end
    end

    sgtitle(figTc, sprintf('Neuron %d: %s vs %s time courses', iNeuron, label1, label2), ...
        'Interpreter', 'none');

    local_save_figure(figTc, fullfile(outFolder, ...
        sprintf('neuron_%04d_%s_timecourses.%s', iNeuron, pairName, opts.saveExt)), opts);
    close(figTc);
end
end

function ampMat = local_compute_pairwise_amp_matrix(amp6, iNeuron, dim1, dim2)
A = amp6(:,:,:,:,:,iNeuron);
otherDims = setdiff(1:5, [dim1 dim2], 'stable');
A = permute(A, [dim1, dim2, otherDims]);

sz = size(A);
n1 = sz(1);
n2 = sz(2);
nRest = prod(sz(3:end));

A = reshape(A, [n1, n2, nRest]);
ampMat = mean(A, 3, 'omitnan');
end

function tcCell = local_collect_pairwise_timecourses(tc7, iNeuron, dim1, dim2, order1, order2)
A = tc7(:,:,:,:,:,:,iNeuron);
otherDims = setdiff(1:5, [dim1 dim2], 'stable');
A = permute(A, [dim2, dim1, otherDims, 6]);

sz = size(A);
n2 = sz(1);
n1 = sz(2);
nTime = sz(end);
nTrials = prod(sz(3:end-1));

A = reshape(A, [n2, n1, nTrials, nTime]);
A = A(order2, order1, :, :);

tcCell = cell(n2, n1);
for i2 = 1:n2
    for i1 = 1:n1
        trials = squeeze(A(i2, i1, :, :));
        if isvector(trials)
            trials = reshape(trials, 1, []);
        end
        tcCell{i2, i1} = trials;
    end
end
end

function [dimIdx, values, label] = local_get_param_info(canon, paramName)
switch lower(strtrim(paramName))
    case {'direction', 'dir'}
        dimIdx = 1;
        values = canon.directions(:).';
        label = 'direction';
    case {'size'}
        dimIdx = 2;
        values = canon.sizes(:).';
        label = 'size';
    case {'tf', 'temporalfrequency', 'temporal_frequency'}
        dimIdx = 3;
        values = canon.tfs(:).';
        label = 'tf';
    case {'sf', 'spatialfrequency', 'spatial_frequency'}
        dimIdx = 4;
        values = canon.sfs(:).';
        label = 'sf';
    otherwise
        error('Unsupported parameter name: %s', paramName);
end

if isempty(values)
    error('No values found in canon for parameter "%s".', paramName);
end
end

function fig = local_make_figure(visibleSetting)
if nargin < 1 || isempty(visibleSetting)
    visibleSetting = 'off';
end
fig = figure('Visible', visibleSetting, 'Color', 'w');
end

function local_save_figure(fig, filePath, opts)
[folderPath, ~, ~] = fileparts(filePath);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

try
    if strcmp(opts.saveExt, 'svg') || strcmp(opts.saveExt, 'pdf')
        exportgraphics(fig, filePath, 'ContentType','vector');
    else
        exportgraphics(fig, filePath, 'Resolution', 150);

    end
catch
    saveas(fig, filePath);
end
end

function out = local_safe_name(txt)
out = regexprep(lower(char(txt)), '\s+', '_');
out = regexprep(out, '[^a-z0-9_]', '');
end
