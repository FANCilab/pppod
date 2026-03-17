function [ampTrials, labels] = flatten_neuron_trials_amp(canon, iNeuron)
% FLATTEN_NEURON_TRIALS_AMP Flatten amplitude responses and stimulus labels.
%
% Returns:
%   ampTrials : [nTrials x 1]
%   labels    : struct with per-trial parameter values for direction,
%               orientation, size, tf, sf, and repeat index.

    A = canon.amp6(:,:,:,:,:,iNeuron);
    ampTrials = A(:);

    [dirIdx, sizeIdx, tfIdx, sfIdx, repIdx] = ndgrid( ...
        1:numel(canon.directions), ...
        1:numel(canon.sizes), ...
        1:numel(canon.tfs), ...
        1:numel(canon.sfs), ...
        1:canon.nRep);

    labels = struct();
    labels.direction = canon.directions(dirIdx(:));
    labels.orientation = canon.orientationPerDirection(dirIdx(:));
    labels.size = canon.sizes(sizeIdx(:));
    labels.tf = canon.tfs(tfIdx(:));
    labels.sf = canon.sfs(sfIdx(:));
    labels.repeat = repIdx(:);
end
