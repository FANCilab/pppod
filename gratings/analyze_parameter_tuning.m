function out = analyze_parameter_tuning(canon, paramInfo, isResponsive, saveFolderTc, opts)
% ANALYZE_PARAMETER_TUNING Generic mean-based analysis for one stimulus parameter.
%
% For each responsive neuron this function:
%   - groups time courses by the requested parameter and saves the summary plot
%   - groups trial amplitudes by the requested parameter
%   - computes mean and standard deviation across pooled trials
%   - finds the preferred parameter value from the maximum mean response
%   - computes DSI for direction tuning and OSI for orientation tuning

    nLevels = numel(paramInfo.values);
    nNeurons = canon.nNeurons;

    out = struct();
    out.meanResponses = nan(nLevels, nNeurons);
    out.stdResponses = nan(nLevels, nNeurons);
    out.preferredValue = nan(1, nNeurons);

    switch lower(paramInfo.name)
        case 'direction'
            out.preferredDirection = nan(1, nNeurons);
            out.DSI = nan(1, nNeurons);
        case 'orientation'
            out.preferredOrientation = nan(1, nNeurons);
            out.OSI = nan(1, nNeurons);
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

        tcByParam = group_trials_by_parameter_tc(canon, paramInfo, iNeuron);
        % plot_parameter_timecourses( ...
        %     tcByParam, canon.t, paramInfo.values, paramInfo.name, iNeuron, saveFolderTc, opts);

        ampByParam = group_trials_by_parameter_amp(canon, paramInfo, iNeuron);
        meanResp = mean(ampByParam, 2, 'omitnan');
        stdResp = std(ampByParam, 0, 2, 'omitnan');

        out.meanResponses(:, iNeuron) = meanResp;
        out.stdResponses(:, iNeuron) = stdResp;

        [prefVal, prefIdx] = compute_preferred_value(meanResp, paramInfo.values);
        out.preferredValue(iNeuron) = prefVal;

        switch lower(paramInfo.name)
            case 'direction'
                out.preferredDirection(iNeuron) = prefVal;
                out.DSI(iNeuron) = compute_direction_selectivity(meanResp, paramInfo.values, prefIdx);
            case 'orientation'
                out.preferredOrientation(iNeuron) = prefVal;
                out.OSI(iNeuron) = compute_orientation_selectivity(meanResp, paramInfo.values);
            case 'size'
                out.preferredSize(iNeuron) = prefVal;
            case 'tf'
                out.preferredTF(iNeuron) = prefVal;
            case 'sf'
                out.preferredSF(iNeuron) = prefVal;
        end
    end
end
