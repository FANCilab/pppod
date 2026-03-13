function plot_pairwise_response_matrices(canon, paramInfos, iNeuron, saveFolder, opts)
% PLOT_PAIRWISE_RESPONSE_MATRICES Plot 2D mean response matrices for all active pairs.
%
% For each pair of active stimulus parameters, the response matrix is the
% mean across repeats and across all remaining stimulus dimensions.

    nParams = numel(paramInfos);
    if nParams < 2
        return
    end

    pairIdx = nchoosek(1:nParams, 2);
    nPairs = size(pairIdx, 1);

    matrixList = cell(1, nPairs);
    allVals = [];

    for iPair = 1:nPairs
        p1 = pairIdx(iPair, 1);
        p2 = pairIdx(iPair, 2);

        matrixList{iPair} = compute_pairwise_response_matrix( ...
            canon.amp6, paramInfos(p1).dim, paramInfos(p2).dim, iNeuron);

        finiteVals = matrixList{iPair}(isfinite(matrixList{iPair}));
        allVals = [allVals; finiteVals(:)]; %#ok<AGROW>
    end

    cLimits = local_compute_clim(allVals);
    nCols = ceil(sqrt(nPairs));
    nRows = ceil(nPairs / nCols);

    figWidth = max(320 * nCols, 500);
    figHeight = max(280 * nRows, 320);
    fig = figure('Visible', opts.visible, 'Color', 'w', ...
        'Position', [100 100 figWidth figHeight]);
    colormap(parula);

    for iPair = 1:nPairs
        p1 = pairIdx(iPair, 1);
        p2 = pairIdx(iPair, 2);

        ax = subplot(nRows, nCols, iPair); %#ok<LAXES>
        imagesc(ax, matrixList{iPair});
        set(ax, 'YDir', 'normal');
        if ~isempty(cLimits)
            caxis(ax, cLimits);
        end

        valuesX = paramInfos(p2).values(:)';
        valuesY = paramInfos(p1).values(:)';
        set(ax, 'XTick', 1:numel(valuesX), 'XTickLabel', local_tick_labels(valuesX));
        set(ax, 'YTick', 1:numel(valuesY), 'YTickLabel', local_tick_labels(valuesY));

        if exist('xtickangle', 'file') == 2
            xtickangle(ax, 45);
        end

        xlabel(ax, paramInfos(p2).prettyName, 'Interpreter', 'none');
        ylabel(ax, paramInfos(p1).prettyName, 'Interpreter', 'none');
        title(ax, sprintf('%s vs %s', paramInfos(p1).prettyName, paramInfos(p2).prettyName), ...
            'Interpreter', 'none');
        box(ax, 'off');
        colorbar(ax);
    end

    add_sgtitle_compat(fig, sprintf('Neuron %d - pairwise mean response matrices', iNeuron));

    fileName = fullfile(saveFolder, sprintf('neuron_%04d_pairwise_response_matrices.%s', ...
        iNeuron, opts.saveExt));
    save_figure_compat(fig, fileName);
    close(fig);
end

function cLimits = local_compute_clim(allVals)
% LOCAL_COMPUTE_CLIM Shared color limits across all pairwise matrices.

    finiteVals = allVals(isfinite(allVals));
    if isempty(finiteVals)
        cLimits = [];
        return
    end

    cMin = min(finiteVals);
    cMax = max(finiteVals);

    if cMin == cMax
        pad = max(1e-6, 0.05 * max(abs(cMin), 1));
        cLimits = [cMin - pad, cMax + pad];
    else
        cLimits = [cMin, cMax];
    end
end

function labels = local_tick_labels(values)
% LOCAL_TICK_LABELS Format numeric axis tick labels.

    labels = arrayfun(@(x) sprintf('%.3g', x), values, 'UniformOutput', false);
end
