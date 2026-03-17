function canon = canonicalize_grating_dimensions(data)
% CANONICALIZE_GRATING_DIMENSIONS Restore canonical stimulus dimensions.
%
% Preferred input layout:
%   data.tc  : [nStim, nRep, nTime, nNeuron]
%   data.amp : [nStim, nRep, nNeuron]
% where nStim = nDir * nSize * nTf * nSf and the stimulus linearization
% follows canonical order [dir, size, tf, sf].
%
% Supported output shapes:
%   canon.tc7  = [dir, size, tf, sf, repeat, time, neuron]
%   canon.amp6 = [dir, size, tf, sf, repeat, neuron]

    info = infer_grating_input_layout(data);

    nDir = numel(data.directions);
    nSize = numel(data.sizes);
    nTf = numel(data.tfs);
    nSf = numel(data.sfs);

    canon = struct();
    canon.tc7 = reshape(data.tc, [nDir, nSize, nTf, nSf, info.nRep, info.nTime, info.nNeurons]);
    canon.amp6 = reshape(data.amp, [nDir, nSize, nTf, nSf, info.nRep, info.nNeurons]);

    canon.blankTc = [];
    canon.blankAmp = [];
    if isfield(data, 'blankTc')
        canon.blankTc = data.blankTc;
    end
    if isfield(data, 'blankAmp')
        canon.blankAmp = data.blankAmp;
    end

    canon.t = data.t(:)';
    canon.directions = data.directions(:)';
    canon.sizes = data.sizes(:)';
    canon.tfs = data.tfs(:)';
    canon.sfs = data.sfs(:)';
    canon.activeParams = info.activeParams;
    canon.inputFormat = info.format;

    canon.nRep = info.nRep;
    canon.nTime = info.nTime;
    canon.nNeurons = info.nNeurons;
end
