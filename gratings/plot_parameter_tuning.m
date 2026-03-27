function plot_parameter_tuning(canon, paramInfos, iNeuron, saveFolder, opts)
% PLOT_PARAMETER_TUNING Plot all active 1D tuning curves for one neuron.
%
% The output figure has 2 rows:
%   Row 1: single-trial responses (open circles) plus mean responses
%   Row 2: mean +/- standard deviation across trials
%
% Each active parameter occupies one column.
% The top row and bottom row use separate shared y-limits.

    nParams = numel(paramInfos);
    if nParams == 0
        return
    end

    ampByParam = cell(1, nParams);
    meanResp = cell(1, nParams);
    stdResp = cell(1, nParams);
    sortedValues = cell(1, nParams);
    allTopY = [];
    allBottomY = [];

    for iParam = 1:nParams
        currentAmp = group_trials_by_parameter_amp(canon, paramInfos(iParam), iNeuron);
        currentValues = paramInfos(iParam).values(:)';
        [currentValues, sortIdx] = sort(currentValues);
        currentAmp = currentAmp(sortIdx, :);

        ampByParam{iParam} = currentAmp;
        sortedValues{iParam} = currentValues;
        meanResp{iParam} = mean(currentAmp, 2, 'omitnan');
        stdResp{iParam} = std(currentAmp, 0, 2, 'omitnan');

        finiteTrials = currentAmp(isfinite(currentAmp));
        finiteMeans = meanResp{iParam}(isfinite(meanResp{iParam}));
        finiteUpper = meanResp{iParam} + stdResp{iParam};
        finiteLower = meanResp{iParam} - stdResp{iParam};

        allTopY = [allTopY; finiteTrials(:); finiteMeans(:)]; %#ok<AGROW>
        allBottomY = [allBottomY; finiteUpper(isfinite(finiteUpper)); finiteLower(isfinite(finiteLower))]; %#ok<AGROW>
    end

    yLimitsTop = local_compute_y_limits(allTopY);
    yLimitsBottom = local_compute_y_limits(allBottomY);
    figWidth = max(320 * nParams, 500);
  
    fig = local_make_figure(opts.visible);
    set(fig, 'Position', [100 100 figWidth 650]);

    for iParam = 1:nParams
        values = sortedValues{iParam};

        axTop = subplot(2, nParams, iParam, 'Parent', fig); %#ok<LAXES>
        hold(axTop, 'on');

        currentAmp = ampByParam{iParam};
        for iLevel = 1:numel(values)
            y = currentAmp(iLevel, :);
            valid = isfinite(y);
            if any(valid)
                x = repmat(values(iLevel), [1, sum(valid)]);
                plot(axTop, x, y(valid), 'o', ...
                    'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', [0.55 0.55 0.55], ...
                    'LineWidth', 1.0, ...
                    'MarkerSize', 4);
            end
        end

        plot(axTop, values, meanResp{iParam}, 'k-', 'LineWidth', 1.5);
        plot(axTop, values, meanResp{iParam}, 'ko', ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 6, ...
            'LineWidth', 1.0);

        ylim(axTop, yLimitsTop);
        local_set_x_axis(axTop, values);
        title(axTop, sprintf('%s | trials + mean', paramInfos(iParam).prettyName), ...
            'Interpreter', 'none');
        if iParam == 1
            ylabel(axTop, 'Response amplitude');
        end
        box(axTop, 'off');

        axBottom = subplot(2, nParams, nParams + iParam, 'Parent', fig); %#ok<LAXES>
        hold(axBottom, 'on');

        errorbar(axBottom, values, meanResp{iParam}, stdResp{iParam}, 'k-o', ...
            'LineWidth', 1.5, ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 6);

        ylim(axBottom, yLimitsBottom);
        local_set_x_axis(axBottom, values);
        xlabel(axBottom, paramInfos(iParam).prettyName, 'Interpreter', 'none');
        if iParam == 1
            ylabel(axBottom, 'Response amplitude');
        end
        title(axBottom, sprintf('%s | mean +/- SD', paramInfos(iParam).prettyName), ...
            'Interpreter', 'none');
        box(axBottom, 'off');
    end

    add_sgtitle_compat(fig, sprintf('Neuron %d - active parameter tuning', iNeuron));

    fileName = fullfile(saveFolder, sprintf('neuron_%04d_combined_tuning.%s', ...
        iNeuron, opts.saveExt));
    save_figure_compat(fig, fileName);
    close(fig);
end

function local_set_x_axis(ax, values)
    values = sort(values(:)');
    if isempty(values)
        return
    end

    if numel(values) == 1
        pad = max(1e-6, 0.05 * max(abs(values(1)), 1));
        xlim(ax, [values(1) - pad, values(1) + pad]);
    else
        xMin = min(values);
        xMax = max(values);
        if xMin == xMax
            pad = max(1e-6, 0.05 * max(abs(xMin), 1));
        else
            pad = 0.05 * (xMax - xMin);
        end
        xlim(ax, [xMin - pad, xMax + pad]);
    end

    set(ax, 'XTick', values);
end

function yLimits = local_compute_y_limits(allY)
    finiteVals = allY(isfinite(allY));
    if isempty(finiteVals)
        yLimits = [0 1];
        return
    end

    yMin = min(finiteVals);
    yMax = max(finiteVals);

    if yMin == yMax
        pad = max(1e-6, 0.05 * max(abs(yMin), 1));
        yLimits = [yMin - pad, yMax + pad];
    else
        pad = 0.05 * (yMax - yMin);
        yLimits = [yMin - pad, yMax + pad];
    end
end

function fig = local_make_figure(visibleSetting)
    if nargin < 1 || isempty(visibleSetting)
        visibleSetting = 'off';
    end
    fig = figure('Visible', visibleSetting, 'Color', 'w');
end