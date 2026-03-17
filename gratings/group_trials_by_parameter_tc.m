function tcByParam = group_trials_by_parameter_tc(canon, paramInfo, iNeuron)
% GROUP_TRIALS_BY_PARAMETER_TC Return [nLevels x nTrialsMax x nTime] groups.

    [tcTrials, labels] = flatten_neuron_trials_tc(canon, iNeuron);
    paramVals = labels.(paramInfo.name);
    levels = paramInfo.values(:)';

    counts = zeros(1, numel(levels));
    for iLevel = 1:numel(levels)
        counts(iLevel) = sum(paramVals == levels(iLevel));
    end
    maxCount = max([counts, 1]);

    tcByParam = nan(numel(levels), maxCount, canon.nTime);
    for iLevel = 1:numel(levels)
        idx = find(paramVals == levels(iLevel));
        if ~isempty(idx)
            tcByParam(iLevel, 1:numel(idx), :) = tcTrials(idx, :);
        end
    end
end
