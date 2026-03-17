function meanMatrix = compute_pairwise_response_matrix(canon, paramInfo1, paramInfo2, iNeuron)
% COMPUTE_PAIRWISE_RESPONSE_MATRIX Mean response matrix for one parameter pair.
%
% Derived parameters such as orientation are supported by flattening trial
% responses and grouping by the requested pair of parameter values.

    [ampTrials, labels] = flatten_neuron_trials_amp(canon, iNeuron);
    values1 = paramInfo1.values(:)';
    values2 = paramInfo2.values(:)';
    vals1 = labels.(paramInfo1.name);
    vals2 = labels.(paramInfo2.name);

    meanMatrix = nan(numel(values1), numel(values2));
    for i1 = 1:numel(values1)
        for i2 = 1:numel(values2)
            idx = (vals1 == values1(i1)) & (vals2 == values2(i2));
            if any(idx)
                meanMatrix(i1, i2) = mean(ampTrials(idx), 'omitnan');
            end
        end
    end
end
