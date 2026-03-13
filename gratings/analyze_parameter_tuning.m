function out = analyze_parameter_tuning(canon, paramName, paramValues, paramDim, ...
    isResponsive, saveFolderTc, opts)
% ANALYZE_PARAMETER_TUNING Generic mean-based analysis for one stimulus parameter.
%
% For each responsive neuron this function:
%   - groups time courses by the requested parameter and saves the summary plot
%   - groups trial amplitudes by the requested parameter
%   - computes mean and standard deviation across pooled trials
%   - finds the preferred parameter value from the maximum mean response
%   - computes DSI for direction tuning

    nLevels = numel(paramValues);
    nNeurons = canon.nNeurons;

    out = struct();
    out.meanResponses = nan(nLevels, nNeurons);
    out.stdResponses = nan(nLevels, nNeurons);
    out.preferredValue = nan(1, nNeurons);

    switch lower(paramName)
        case 'direction'
            out.preferredDirection = nan(1, nNeurons);
            out.DSI = nan(1, nNeurons);
        case 'size'
            out.preferredSize = nan(1, nNeurons);
        case 'tf'
            out.preferredTF = nan(1, nNeurons);
        case 'sf'
            out.preferredSF = nan(1, nNeurons);
    end

    for iNeuron = 1:nNeurons
        if ~isResponsive(iNeuron)
            continue
        end

        tcByParam = reshape_tc_by_parameter(canon.tc7, paramDim, iNeuron);
        plot_parameter_timecourses( ...
            tcByParam, canon.t, paramValues, paramName, iNeuron, saveFolderTc, opts);

        ampByParam = reshape_amp_by_parameter(canon.amp6, paramDim, iNeuron);
        meanResp = mean(ampByParam, 2, 'omitnan');
        stdResp = std(ampByParam, 0, 2, 'omitnan');

        out.meanResponses(:, iNeuron) = meanResp;
        out.stdResponses(:, iNeuron) = stdResp;

        [prefVal, prefIdx] = compute_preferred_value(meanResp, paramValues);
        out.preferredValue(iNeuron) = prefVal;

        switch lower(paramName)
            case 'direction'
                out.preferredDirection(iNeuron) = prefVal;
                out.DSI(iNeuron) = compute_direction_selectivity(meanResp, paramValues, prefIdx);
            case 'size'
                out.preferredSize(iNeuron) = prefVal;
            case 'tf'
                out.preferredTF(iNeuron) = prefVal;
            case 'sf'
                out.preferredSF(iNeuron) = prefVal;
        end
    end
end
