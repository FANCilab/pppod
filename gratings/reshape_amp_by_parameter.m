function ampByParam = reshape_amp_by_parameter(amp6, paramDim, iNeuron)
% RESHAPE_AMP_BY_PARAMETER Return [nLevels, nTrials] for one neuron.
%
% amp6 canonical dimension order:
%   [dir, size, tf, sf, repeat, neuron]
%
% paramDim should be one of the stimulus dimensions:
%   1 = direction, 2 = size, 3 = tf, 4 = sf

    if ~ismember(paramDim, 1:4)
        error('paramDim must be 1..4 for one of the stimulus dimensions.');
    end

    A = amp6(:,:,:,:,:,iNeuron);       % [dir, size, tf, sf, rep]
    otherDims = setdiff(1:5, paramDim, 'stable');
    A = permute(A, [paramDim, otherDims]);

    sz = size(A);
    nLevels = sz(1);
    if numel(sz) == 1
        nTrials = 1;
    else
        nTrials = prod(sz(2:end));
    end

    ampByParam = reshape(A, [nLevels, nTrials]);
end
