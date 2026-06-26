function [responses, aveResponse, seResponse, kernelTimes, stimDur] = ...
getStimulusSweepsLFR(F, stimTimes, stimMatrix, FrameRate, stimSet)

% define parameters
[nT, nroi] =  size(F);

iSi = floor(min(stimTimes.onset(2:end)-stimTimes.offset(1:(end-1)))); %inter stimulus interval
stimDur= mean(stimTimes.offset-stimTimes.onset);

% iSi = 2;
% stimDur = 2;

[stimSeq, onsets] = find([stimMatrix(:,1), diff(stimMatrix, 1, 2)] > 0);

nStim = numel(unique(stimSeq));
% Fail if nRep is not the same for each angle
[uniqueStim, ~, stimGroup] = unique(stimSeq);
stimCounts = accumarray(stimGroup, 1);
if any(stimCounts ~= stimCounts(1))
    error('getStimulusSweepsLFR:UnequalStimRepetitions', ...
        'Stimuli do not all have the same number of repetitions. Counts per stimulus: %s', ...
        mat2str([uniqueStim, stimCounts]));
end
nRep = stimCounts(1);

if nargin < 5
    stimSet = 1:nStim;
else
    nStim = numel(stimSet);
end

stimON = zeros(nStim, nRep);

for iStim = 1:nStim
    stimON(iStim,:) = onsets(stimSeq == stimSet(iStim));
end


%% chunk up the F in sweeps

% F = [NaN(startpad, nroi);F ;  NaN(3*iSi, nroi)];
iSiFr = floor(iSi*FrameRate); %inter stimulus frames
F = [zeros(iSiFr, nroi);F ;  zeros(3*iSiFr, nroi)]; 
[nT, nroi] =  size(F);
t = 1:nT; t = t/FrameRate;
stimON =  (stimON + iSiFr)/FrameRate;

%response size is nroi x nStim x nRep x framesXsweep
% sweepON =  - iSi;
% sweepOFF = stimDur +iSi;

sweepON =  -2; 
sweepOFF =  4;
kernelTimes = linspace(sweepON, sweepOFF, range([sweepON, sweepOFF])*10);
frameXsweep = numel(kernelTimes);

responses = NaN(nroi, nStim, nRep, frameXsweep);
aveResponse = NaN(nroi, nStim, frameXsweep);
seResponse = NaN(nroi, nStim, frameXsweep);

for iroi = 1:nroi
    
    [resp, aveResp, seResp] = magicETA(t, F(:, iroi), stimON, kernelTimes);
%     [resp, aveResp, seResp] = magicETA2(t, F, stimON, kernelTimes);

    responses(iroi, :,:,:) = permute(resp, [3,1,2]);
    aveResponse(iroi, :,:) = aveResp';
    seResponse(iroi, :,:) = seResp';
    
end


end

