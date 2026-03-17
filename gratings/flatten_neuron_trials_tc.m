function [tcTrials, labels] = flatten_neuron_trials_tc(canon, iNeuron)
% FLATTEN_NEURON_TRIALS_TC Flatten time-course responses and stimulus labels.
%
% Returns:
%   tcTrials : [nTrials x nTime]
%   labels   : struct with per-trial parameter values for direction,
%              orientation, size, tf, sf, and repeat index.

    A = canon.tc7(:,:,:,:,:,:,iNeuron);
    tcTrials = reshape(A, [], canon.nTime);

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
