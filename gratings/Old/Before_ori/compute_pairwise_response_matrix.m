function meanMatrix = compute_pairwise_response_matrix(amp6, dim1, dim2, iNeuron)
% COMPUTE_PAIRWISE_RESPONSE_MATRIX Mean response matrix for one parameter pair.
%
% amp6 canonical dimension order:
%   [dir, size, tf, sf, repeat, neuron]
%
% The output is the mean across repeats and all stimulus dimensions other
% than dim1 and dim2.

    if dim1 == dim2
        error('dim1 and dim2 must be different stimulus dimensions.');
    end

    A = amp6(:,:,:,:,:,iNeuron);   % [dir, size, tf, sf, rep]
    keepDims = [dim1, dim2];
    otherDims = setdiff(1:5, keepDims, 'stable');
    A = permute(A, [keepDims, otherDims]);

    sz = size(A);
    sz(end+1:3) = 1;
    n1 = sz(1);
    n2 = sz(2);
    nPool = prod(sz(3:end));

    A = reshape(A, [n1, n2, nPool]);
    meanMatrix = mean(A, 3, 'omitnan');
end
