function [resPeak, aveResPeak, seResPeak] = gratingOnResp(resp, kernelTimes, win)

% responses is (nroi, nStim, nRep, frameXsweep) and are baseline subtracted

respFrames = kernelTimes > win(1) & kernelTimes < win(2);

resPeak = mean(resp(:,:,:,respFrames), 4); % (nroi, nStim, nRep)

aveResPeak = mean(resPeak, 3);  % (nroi, nStim)

seResPeak = std(resPeak, 1, 3)/sqrt(size(resPeak, 3));  % (nroi, nStim)


end