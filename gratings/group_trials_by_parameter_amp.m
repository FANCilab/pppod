function ampByParam = group_trials_by_parameter_amp(canon, paramInfo, iNeuron)
% GROUP_TRIALS_BY_PARAMETER_AMP Return [nLevels x nTrialsMax] amplitude groups.

    [ampTrials, labels] = flatten_neuron_trials_amp(canon, iNeuron);
    paramVals = labels.(paramInfo.name);
    levels = paramInfo.values(:)';

    counts = zeros(1, numel(levels));
    for iLevel = 1:numel(levels)
        counts(iLevel) = sum(paramVals == levels(iLevel));
    end
    maxCount = max([counts, 1]);

    ampByParam = nan(numel(levels), maxCount);
    for iLevel = 1:numel(levels)
        idx = find(paramVals == levels(iLevel));
        if ~isempty(idx)
            ampByParam(iLevel, 1:numel(idx)) = ampTrials(idx);
        end
    end
end
