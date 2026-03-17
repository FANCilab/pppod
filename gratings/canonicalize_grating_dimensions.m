function canon = canonicalize_grating_dimensions(data)
% CANONICALIZE_GRATING_DIMENSIONS Convert supported input layouts to canonical form.
%
% Canonical outputs:
%   tc7  : [nDir, nSize, nTf, nSf, nRep, nTime, nNeuron]
%   amp6 : [nDir, nSize, nTf, nSf, nRep, nNeuron]
%
% Supported input layouts:
%   1) Flattened stimulus layout (preferred):
%        data.tc  = [nStim, nRep, nTime, nNeuron]
%        data.amp = [nStim, nRep, nNeuron]
%      where nStim = nDir * nSize * nTf * nSf and the flattened stimulus
%      index follows canonical order [dir, size, tf, sf].
%
%   2) Full canonical layout:
%        data.tc  = [nDir, nSize, nTf, nSf, nRep, nTime, nNeuron]
%        data.amp = [nDir, nSize, nTf, nSf, nRep, nNeuron]
%
%   3) Legacy squeezed canonical layout. Singleton inactive stimulus
%      dimensions are inferred from the parameter vectors and reinserted.

    nDir = numel(data.directions);
    nSize = numel(data.sizes);
    nTf = numel(data.tfs);
    nSf = numel(data.sfs);
    nStimExpected = nDir * nSize * nTf * nSf;

    layout = infer_grating_input_layout(data);

    switch layout.format
        case 'flattened_stimulus'
            canon.tc7 = reshape(data.tc, [nDir, nSize, nTf, nSf, size(data.tc, 2), size(data.tc, 3), size(data.tc, 4)]);
            canon.amp6 = reshape(data.amp, [nDir, nSize, nTf, nSf, size(data.amp, 2), size(data.amp, 3)]);

        case 'canonical_full'
            canon.tc7 = data.tc;
            canon.amp6 = data.amp;

        case 'canonical_squeezed'
            tcSize = size(data.tc);
            ampSize = size(data.amp);

            tcSize = [tcSize, ones(1, max(0, 7 - numel(tcSize)))];
            ampSize = [ampSize, ones(1, max(0, 6 - numel(ampSize)))];

            activeCounts = [nDir > 1, nSize > 1, nTf > 1, nSf > 1];
            activeDims = find(activeCounts);

            nRep = ampSize(numel(activeDims) + 1);
            nNeuron = ampSize(numel(activeDims) + 2);
            nTime = tcSize(numel(activeDims) + 2);

            canonicalStimShape = [nDir, nSize, nTf, nSf];
            activeShape = canonicalStimShape(activeDims);
            if isempty(activeShape)
                activeShape = 1;
            end

            canon.tc7 = reshape(data.tc, [activeShape, nRep, nTime, nNeuron]);
            canon.amp6 = reshape(data.amp, [activeShape, nRep, nNeuron]);

            canon.tc7 = reshape(canon.tc7, [nDir, nSize, nTf, nSf, nRep, nTime, nNeuron]);
            canon.amp6 = reshape(canon.amp6, [nDir, nSize, nTf, nSf, nRep, nNeuron]);

        otherwise
            error('Unsupported input layout: %s', layout);
    end

    canon.t = data.t(:)';
    canon.directions = data.directions(:)';
    canon.sizes = data.sizes(:)';
    canon.tfs = data.tfs(:)';
    canon.sfs = data.sfs(:)';
    canon.activeParams = data.activeParams(:)';

    canon.nRep = size(canon.amp6, 5);
    canon.nTime = size(canon.tc7, 6);
    canon.nNeurons = size(canon.amp6, 6);
    canon.inputFormat = layout.format;
    canon.hasDirection = ismember('direction', canon.activeParams) && numel(canon.directions) > 1;
    canon.hasOrientation = canon.hasDirection;
    canon.orientationPerDirection = mod(canon.directions + 90, 180);
    canon.orientations = unique(canon.orientationPerDirection);
    [~, canon.orientationIndexPerDirection] = ismember(canon.orientationPerDirection, canon.orientations);

    if isfield(data, 'blankTc')
        canon.blankTc = data.blankTc;
    else
        canon.blankTc = [];
    end

    if isfield(data, 'blankAmp')
        canon.blankAmp = data.blankAmp;
    else
        canon.blankAmp = [];
    end
end
