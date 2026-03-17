function tcByParam = reshape_tc_by_parameter(tc7, paramDim, iNeuron)
% RESHAPE_TC_BY_PARAMETER Return [nLevels, nTrials, nTime] for one neuron.
%
% tc7 canonical dimension order:
%   [dir, size, tf, sf, repeat, time, neuron]
%
% paramDim should be one of the stimulus dimensions:
%   1 = direction, 2 = size, 3 = tf, 4 = sf

    if ~ismember(paramDim, 1:4)
        error('paramDim must be 1..4 for one of the stimulus dimensions.');
    end

    A = tc7(:,:,:,:,:,:,iNeuron);      % [dir, size, tf, sf, rep, time]
    otherDims = setdiff(1:5, paramDim, 'stable');
    A = permute(A, [paramDim, otherDims, 6]);

    sz = size(A);
    nLevels = sz(1);
    nTime = sz(end);
    if numel(sz) <= 2
        nTrials = 1;
    else
        nTrials = prod(sz(2:end-1));
    end

    tcByParam = reshape(A, [nLevels, nTrials, nTime]);
end
