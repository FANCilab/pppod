function [signalCorrMat, noiseCorrMat, totCorrMat, globalSignal] = globalCorr(data)
% computes correlations catenating experiments with different stimulus types;
% data is a cell array, each cell contains responses from experiment(i),
% responses are  a matrix of dim nTrial(i)*nT(i)*nStim(i)*nN

% [~,nStimType] =  size(data); %nTrial*nT*nStim*nN
nStimType =  sum(~cellfun(@isempty,data)); %nTrial*nT*nStim*nN
startType = find(~cellfun(@isempty,data),1,'first');

globalSignal = [];
globalNoise = [];
globalResp= [];

for iTy = startType:(nStimType+startType-1)
    
    if ~isempty(data{1, iTy})
    response = data{1, iTy}; %nTrial*nT*nStim*nN
    
    [nTrial, nT, nStim, nN] = size(response);
    
    signal = shiftdim(mean(response, 1),1); %nT*nStim*nN

    noise = bsxfun(@minus, response, reshape(signal, 1, nT, nStim, nN));%nTrial*nT*nStim*nN
    
    response = reshape(permute(response, [2,1,3,4]), nT*nTrial*nStim, nN); %(nT*nTrial*nStim)*nN
    
    signal = reshape(signal, nT*nStim, nN); %(nT*nStim)*nN
    
    noise = reshape(permute(noise, [2,1,3,4]), nT*nTrial*nStim, nN); %(nT*nTrial*nStim)*nN

    globalResp = cat(1, globalResp, response); %nT*nN
    
    globalNoise = cat(1, globalNoise, noise); %nT*nN
    
    globalSignal = cat(1, globalSignal, signal); %nT*nN
    end
end

noiseCorrMat = corrcoef(globalNoise); % compute pairwise noise correlations
noiseCorrMat(1:nN+1:nN*nN) = NaN; % set diagonal elements to NaN

totCorrMat = corrcoef(globalResp); 
totCorrMat(1:nN+1:nN*nN) = NaN; 

signalCorrMat = corrcoef(globalSignal); % compute pairwise noise correlations
signalCorrMat(1:nN+1:nN*nN) = NaN; % set diagonal elements to NaN

end