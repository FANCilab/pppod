function info = infer_grating_input_layout(data)
% INFER_GRATING_INPUT_LAYOUT Infer input layout and dimensions.
%
% Preferred supported layout:
%   data.tc  : [nStim, nRep, nTime, nNeuron]
%   data.amp : [nStim, nRep, nNeuron]
% where nStim = nDir * nSize * nTf * nSf and stimulus order follows the
% canonical unsqueezed ordering [dir, size, tf, sf].
%
% For backward compatibility, this helper also accepts:
%   - full canonical arrays with explicit singleton dimensions
%   - legacy squeezed canonical arrays with inactive stimulus dimensions
%     removed before repeat/time/neuron.
%
% Output fields:
%   info.format      : 'flattened-stim', 'full-canonical', or 'legacy-squeezed'
%   info.activeParams: normalized active-parameter cell array
%   info.activeMask  : logical mask for {'direction','size','tf','sf'}
%   info.nStim       : total number of stimulus combinations
%   info.nRep        : number of repeats
%   info.nTime       : number of time samples
%   info.nNeurons    : number of neurons

    paramNames = {'direction', 'size', 'tf', 'sf'};
    paramCounts = [numel(data.directions), numel(data.sizes), ...
                   numel(data.tfs), numel(data.sfs)];
    nStim = prod(paramCounts);
    nTimeExpected = numel(data.t);

    if isfield(data, 'activeParams') && ~isempty(data.activeParams)
        if isstring(data.activeParams)
            activeParams = cellstr(data.activeParams);
        else
            activeParams = data.activeParams;
        end
        activeParams = cellfun(@char, activeParams, 'UniformOutput', false);
    else
        activeParams = paramNames(paramCounts > 1);
    end

    activeMask = ismember(paramNames, activeParams);
    activeCounts = paramCounts(activeMask);
    nActive = sum(activeMask);

    [flatAmpOK, flatNRepAmp, flatNNeuronAmp] = local_match_flat_amp(data.amp, nStim);
    [flatTcOK, flatNRepTc, flatNTimeTc, flatNNeuronTc] = ...
        local_match_flat_tc(data.tc, nStim, nTimeExpected);

    [fullAmpOK, fullNRepAmp, fullNNeuronAmp] = local_match_full_amp(data.amp, paramCounts);
    [fullTcOK, fullNRepTc, fullNTimeTc, fullNNeuronTc] = ...
        local_match_full_tc(data.tc, paramCounts, nTimeExpected);

    [legacyAmpOK, legacyNRepAmp, legacyNNeuronAmp] = ...
        local_match_legacy_amp(data.amp, activeCounts, nStim, nActive);
    [legacyTcOK, legacyNRepTc, legacyNTimeTc, legacyNNeuronTc] = ...
        local_match_legacy_tc(data.tc, activeCounts, nStim, nActive, nTimeExpected);

    flatOK = flatAmpOK && flatTcOK && ...
        flatNRepAmp == flatNRepTc && flatNNeuronAmp == flatNNeuronTc;

    fullOK = fullAmpOK && fullTcOK && ...
        fullNRepAmp == fullNRepTc && fullNNeuronAmp == fullNNeuronTc;

    legacyOK = legacyAmpOK && legacyTcOK && ...
        legacyNRepAmp == legacyNRepTc && legacyNNeuronAmp == legacyNNeuronTc;

    info = struct();
    info.activeParams = activeParams;
    info.activeMask = activeMask;
    info.paramCounts = paramCounts;
    info.nStim = nStim;

    if flatOK
        info.format = 'flattened-stim';
        info.nRep = flatNRepAmp;
        info.nTime = flatNTimeTc;
        info.nNeurons = flatNNeuronAmp;
    elseif fullOK
        info.format = 'full-canonical';
        info.nRep = fullNRepAmp;
        info.nTime = fullNTimeTc;
        info.nNeurons = fullNNeuronAmp;
    elseif legacyOK
        info.format = 'legacy-squeezed';
        info.nRep = legacyNRepAmp;
        info.nTime = legacyNTimeTc;
        info.nNeurons = legacyNNeuronAmp;
    else
        error(['Could not interpret input array layout. Supported layouts are: ', ...
               '(1) flattened stimulus inputs data.tc=[nStim nRep nTime nNeuron] ', ...
               'and data.amp=[nStim nRep nNeuron]; ', ...
               '(2) full canonical inputs [dir size tf sf rep time neuron] / ', ...
               '[dir size tf sf rep neuron]; ', ...
               '(3) legacy squeezed canonical inputs with inactive stimulus ', ...
               'dimensions removed before repeat/time/neuron.']);
    end
end

function [ok, nRep, nNeurons] = local_match_flat_amp(A, nStim)
    sz = size(A);
    sz(end+1:3) = 1;

    nRep = sz(2);
    nNeurons = sz(3);

    ok = sz(1) == nStim && numel(A) == nStim * nRep * nNeurons;
end

function [ok, nRep, nTime, nNeurons] = local_match_flat_tc(A, nStim, nTimeExpected)
    sz = size(A);
    sz(end+1:4) = 1;

    nRep = sz(2);
    nTime = sz(3);
    nNeurons = sz(4);

    ok = sz(1) == nStim && nTime == nTimeExpected && ...
        numel(A) == nStim * nRep * nTime * nNeurons;
end

function [ok, nRep, nNeurons] = local_match_full_amp(A, paramCounts)
    sz = size(A);
    sz(end+1:6) = 1;

    nRep = sz(5);
    nNeurons = sz(6);

    ok = isequal(sz(1:4), paramCounts) && ...
        numel(A) == prod(paramCounts) * nRep * nNeurons;
end

function [ok, nRep, nTime, nNeurons] = local_match_full_tc(A, paramCounts, nTimeExpected)
    sz = size(A);
    sz(end+1:7) = 1;

    nRep = sz(5);
    nTime = sz(6);
    nNeurons = sz(7);

    ok = isequal(sz(1:4), paramCounts) && nTime == nTimeExpected && ...
        numel(A) == prod(paramCounts) * nRep * nTime * nNeurons;
end

function [ok, nRep, nNeurons] = local_match_legacy_amp(A, activeCounts, nStim, nActive)
    repIdx = nActive + 1;
    neuronIdx = nActive + 2;

    sz = size(A);
    sz(end+1:neuronIdx) = 1;

    prefixOK = true;
    if nActive > 0
        prefixOK = isequal(sz(1:nActive), activeCounts);
    end

    nRep = sz(repIdx);
    nNeurons = sz(neuronIdx);

    ok = prefixOK && numel(A) == nStim * nRep * nNeurons;
end

function [ok, nRep, nTime, nNeurons] = local_match_legacy_tc(A, activeCounts, nStim, nActive, nTimeExpected)
    repIdx = nActive + 1;
    timeIdx = nActive + 2;
    neuronIdx = nActive + 3;

    sz = size(A);
    sz(end+1:neuronIdx) = 1;

    prefixOK = true;
    if nActive > 0
        prefixOK = isequal(sz(1:nActive), activeCounts);
    end

    nRep = sz(repIdx);
    nTime = sz(timeIdx);
    nNeurons = sz(neuronIdx);

    ok = prefixOK && nTime == nTimeExpected && ...
        numel(A) == nStim * nRep * nTime * nNeurons;
end
