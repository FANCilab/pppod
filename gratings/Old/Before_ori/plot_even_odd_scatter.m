function plot_even_odd_scatter(x, y, R, P, isResponsive, iNeuron, saveFolder, opts)
% PLOT_EVEN_ODD_SCATTER Save odd-vs-even mean response scatter for one neuron.

    valid = isfinite(x) & isfinite(y);
    xv = x(valid);
    yv = y(valid);

    fig = figure('Visible', opts.visible, 'Color', 'w');
    ax = axes('Parent', fig); %#ok<LAXES>
    hold(ax, 'on');

    if ~isempty(xv)
        plot(ax, xv, yv, 'o', ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.25 0.25 0.25], ...
            'LineWidth', 1.0);

        minVal = min([xv; yv]);
        maxVal = max([xv; yv]);
        if minVal == maxVal
            minVal = minVal - 1;
            maxVal = maxVal + 1;
        end
        plot(ax, [minVal, maxVal], [minVal, maxVal], '--', ...
            'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
        xlim(ax, [minVal, maxVal]);
        ylim(ax, [minVal, maxVal]);
    end

    xlabel(ax, 'Odd-repeat mean response');
    ylabel(ax, 'Even-repeat mean response');
    axis(ax, 'square');
    box(ax, 'off');

    if isResponsive
        statusText = 'responsive';
    else
        statusText = 'not responsive';
    end

    title(ax, sprintf('Neuron %d | r = %.3f | p = %.3g | %s', ...
        iNeuron, R, P, statusText), 'Interpreter', 'none');

    fileName = fullfile(saveFolder, sprintf('neuron_%04d_even_odd_scatter.%s', ...
        iNeuron, opts.saveExt));
    save_figure_compat(fig, fileName);
    close(fig);
end
