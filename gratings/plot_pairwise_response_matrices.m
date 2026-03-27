function plot_pairwise_response_matrices(canon, paramInfos, iNeuron, saveFolder, opts)
% PLOT_PAIRWISE_RESPONSE_MATRICES Plot 2D mean response matrices for all active pairs.
%
% For each pair of active or derived stimulus parameters, the response matrix
% is the mean across repeats and across all remaining stimulus dimensions.

    nParams = numel(paramInfos);
    if nParams < 2
        return
    end

    pairIdx = nchoosek(1:nParams, 2);
    nPairs = size(pairIdx, 1);

    matrixList = cell(1, nPairs);
    xValuesList = cell(1, nPairs);
    yValuesList = cell(1, nPairs);
    allVals = [];

    for iPair = 1:nPairs
        p1 = pairIdx(iPair, 1);
        p2 = pairIdx(iPair, 2);

        matrixList{iPair} = compute_pairwise_response_matrix(canon, paramInfos(p1), paramInfos(p2), iNeuron);
        [yValuesList{iPair}, yOrder] = sort(paramInfos(p1).values(:)');
        [xValuesList{iPair}, xOrder] = sort(paramInfos(p2).values(:)');
        matrixList{iPair} = matrixList{iPair}(yOrder, xOrder);

        finiteVals = matrixList{iPair}(isfinite(matrixList{iPair}));
        allVals = [allVals; finiteVals(:)]; %#ok<AGROW>
    end

    cLimits = local_compute_clim(allVals);
    nCols = ceil(sqrt(nPairs));
    nRows = ceil(nPairs / nCols);

    figWidth = max(320 * nCols, 500);
    figHeight = max(280 * nRows, 320);

     fig = local_make_figure(opts.visible);

    set(fig, 'Position', [100 100 figWidth figHeight]);
    colormap(parula);

    for iPair = 1:nPairs
        p1 = pairIdx(iPair, 1);
        p2 = pairIdx(iPair, 2);

        ax = subplot(nRows, nCols, iPair, 'Parent', fig); %#ok<LAXES>
        imagesc(ax, matrixList{iPair});
        set(ax, 'YDir', 'normal');
        if ~isempty(cLimits)
            caxis(ax, cLimits);
        end

        valuesX = xValuesList{iPair};
        valuesY = yValuesList{iPair};
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
    labels = arrayfun(@(x) sprintf('%.3g', x), values, 'UniformOutput', false);
end

function fig = local_make_figure(visibleSetting)
    if nargin < 1 || isempty(visibleSetting)
        visibleSetting = 'off';
    end
    fig = figure('Visible', visibleSetting, 'Color', 'w');
end
