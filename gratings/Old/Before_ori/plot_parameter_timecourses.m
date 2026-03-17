function plot_parameter_timecourses(tcByParam, t, paramValues, paramName, iNeuron, saveFolder, opts)
% PLOT_PARAMETER_TIMECOURSES Plot single-trial and mean traces by parameter value.
%
% tcByParam has size [nLevels, nTrials, nTime].
% All subplots are arranged in one row and share the same y-axis limits.
% The thick line is the mean across trials at each time sample.

    nLevels = size(tcByParam, 1);
    if nLevels == 0
        return
    end

    figWidth = max(320 * nLevels, 500);
    fig = figure('Visible', opts.visible, 'Color', 'w', ...
        'Position', [100 100 figWidth 300]);

    yLimits = local_compute_y_limits(tcByParam);
    axesHandles = [];

    for iLevel = 1:nLevels
        ax = subplot(1, nLevels, iLevel); %#ok<LAXES>
        axesHandles = [axesHandles, ax]; %#ok<AGROW>
        hold(ax, 'on');

        trials = squeeze(tcByParam(iLevel, :, :));
        if isvector(trials)
            trials = reshape(trials, 1, []);
        end

        plot(ax, t, trials', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);

        meanTrace = mean(trials, 1, 'omitnan');
        plot(ax, t, meanTrace, 'k-', 'LineWidth', 2.0);

        if exist('xline', 'file') == 2
            xline(ax, 0, '--', 'Color', [0.2 0.2 0.2]);
        else
            line(ax, [0 0], yLimits, 'LineStyle', '--', 'Color', [0.2 0.2 0.2]);
        end

        ylim(ax, yLimits);
        title(ax, sprintf('%s = %.6g', paramName, paramValues(iLevel)), ...
            'Interpreter', 'none');
        xlabel(ax, 'Time from stimulus onset');
        if iLevel == 1
            ylabel(ax, 'Response');
        end
        box(ax, 'off');
    end

    if exist('linkaxes', 'file') == 2 && nLevels > 1
        linkaxes(axesHandles, 'y');
    end

    add_sgtitle_compat(fig, sprintf('Neuron %d - %s time courses', iNeuron, paramName));

    fileName = fullfile(saveFolder, sprintf('neuron_%04d_%s_timecourses.%s', ...
        iNeuron, paramName, opts.saveExt));
    save_figure_compat(fig, fileName);
    close(fig);
end

function yLimits = local_compute_y_limits(tcByParam)
% LOCAL_COMPUTE_Y_LIMITS Shared y-limits across all time-course subplots.

    finiteVals = tcByParam(isfinite(tcByParam));
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
