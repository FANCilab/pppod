function [pVals, isResponsive, rVals] = compute_visual_responsiveness(amp6, saveFolder, opts)
% COMPUTE_VISUAL_RESPONSIVENESS Odd-vs-even repeat correlation per neuron.
%
% For each neuron:
%   1) average responses over odd repeats
%   2) average responses over even repeats
%   3) correlate stimulus-wise odd and even averages
%   4) classify as responsive when p < opts.alpha

    [~, ~, ~, ~, nRep, nNeurons] = size(amp6);

    oddIdx = 1:2:nRep;
    evenIdx = 2:2:nRep;

    pVals = nan(1, nNeurons);
    rVals = nan(1, nNeurons);
    isResponsive = false(1, nNeurons);

    for iNeuron = 1:nNeurons
        oddMean = mean(amp6(:,:,:,:,oddIdx,iNeuron), 5, 'omitnan');
        evenMean = mean(amp6(:,:,:,:,evenIdx,iNeuron), 5, 'omitnan');

        x = oddMean(:);
        y = evenMean(:);
        valid = isfinite(x) & isfinite(y);

        R = NaN;
        P = NaN;

        if numel(evenIdx) >= 1 && numel(oddIdx) >= 1 && sum(valid) >= 3
            xValid = x(valid);
            yValid = y(valid);
            if std(xValid) > 0 && std(yValid) > 0
                [Rmat, Pmat] = corrcoef(xValid, yValid, 'Rows', 'complete');
                R = Rmat(1, 2);
                P = Pmat(1, 2);
            end
        end

        rVals(iNeuron) = R;
        pVals(iNeuron) = P;
        isResponsive(iNeuron) = ~isnan(P) && (P < opts.alpha);

        plot_even_odd_scatter(x, y, R, P, isResponsive(iNeuron), iNeuron, saveFolder, opts);
    end
end
