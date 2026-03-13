function plot_parameter_tuning(canon, paramInfos, iNeuron, saveFolder, opts)
% PLOT_PARAMETER_TUNING Plot all active 1D tuning curves for one neuron.
%
% The output figure has 2 rows:
%   Row 1: single-trial responses (open circles) plus mean responses
%   Row 2: mean +/- standard deviation across trials
%
% Each active parameter occupies one column.

    nParams = numel(paramInfos);
    if nParams == 0
        return
    end

    ampByParam = cell(1, nParams);
    meanResp = cell(1, nParams);
    stdResp = cell(1, nParams);
    allY = [];

    for iParam = 1:nParams
        ampByParam{iParam} = reshape_amp_by_parameter(canon.amp6, paramInfos(iParam).dim, iNeuron);
        meanResp{iParam} = mean(ampByParam{iParam}, 2, 'omitnan');
        stdResp{iParam} = std(ampByParam{iParam}, 0, 2, 'omitnan');

        finiteTrials = ampByParam{iParam}(isfinite(ampByParam{iParam}));
        finiteUpper = meanResp{iParam} + stdResp{iParam};
        finiteLower = meanResp{iParam} - stdResp{iParam};
        allY = [allY; finiteTrials(:); finiteUpper(isfinite(finiteUpper)); finiteLower(isfinite(finiteLower))]; %#ok<AGROW>
    end

    yLimits = local_compute_y_limits(allY);
    figWidth = max(320 * nParams, 500);
    fig = figure('Visible', opts.visible, 'Color', 'w', ...
        'Position', [100 100 figWidth 650]);

    for iParam = 1:nParams
        values = paramInfos(iParam).values(:)';

        axTop = subplot(2, nParams, iParam); %#ok<LAXES>
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

        ylim(axTop, yLimits);
        local_set_x_axis(axTop, values);
        title(axTop, sprintf('%s | trials + mean', paramInfos(iParam).prettyName), ...
            'Interpreter', 'none');
        if iParam == 1
            ylabel(axTop, 'Response amplitude');
        end
        box(axTop, 'off');

        axBottom = subplot(2, nParams, nParams + iParam); %#ok<LAXES>
        hold(axBottom, 'on');

        errorbar(axBottom, values, meanResp{iParam}, stdResp{iParam}, 'k-o', ...
            'LineWidth', 1.5, ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 6);

        ylim(axBottom, yLimits);
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
% LOCAL_SET_X_AXIS Configure x-limits and tick labels for numeric parameter values.

    values = values(:)';
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
% LOCAL_COMPUTE_Y_LIMITS Compute shared y-limits for tuning plots.

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
