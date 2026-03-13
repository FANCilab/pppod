function stimMatrix=buildStimMatrix(stimSequence, stimTimes, timeAxis)

nSamples=length(timeAxis);
nStim=numel(unique(stimSequence));

stimMatrix=false(nStim, nSamples);

for iStim=1:length(stimSequence)
    ind=(timeAxis>=stimTimes.onset(iStim) & timeAxis<=stimTimes.offset(iStim));
    stimMatrix(stimSequence(iStim), ind)=true;
end
end