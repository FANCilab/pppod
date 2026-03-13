function [ R, tKernel,x, rfContours, periT] = ...
    medianSTRF(F, fT, stimFrames, stimFrameTimes, stimTimes, respType, x0, doPlot)
if nargin <8
    doPlot = 0;
end

if nargin < 7
    x0 = [];
end

if nargin < 6
    respType = 'any';
end

% % periT = -0.6:0.1:1; %GCaMP6f
% periT = -0.6:0.1:2; %GCaMP6s
% periT = -0.6:0.1:0.8; %GCaMP6f
periT = -0.6:0.1:1; %GCaMP8s

nRep = numel(stimTimes.onset);
[w, h, nFrames] = size(stimFrames);

switch respType
    case 'on'
        igood = stimFrames(:,:,2:nFrames) == 1 & stimFrames(:,:,1:nFrames-1) == 0;
    case 'off'
        igood = stimFrames(:,:,2:nFrames) == -1 & stimFrames(:,:,1:nFrames-1) == 0;
    case 'onoff'
        igoodOn = stimFrames(:,:,2:nFrames) == 1 & stimFrames(:,:,1:nFrames-1) == 0;
        igoodOff = stimFrames(:,:,2:nFrames) == -1 & stimFrames(:,:,1:nFrames-1) == 0;
        igood = igoodOn | igoodOff;
    case 'any'
        igood = stimFrames(:,:,2:nFrames) ~= stimFrames(:,:,1:nFrames-1);
end

igood = reshape(igood, w*h, []);

onidx = cell(w,h);
for iRep = 1:nRep
    for ix = 1: w*h
    onidx{ix} = [onidx{ix}, stimFrameTimes(find(igood(ix, :))+1) + stimTimes.onset(iRep)];
    end
end

nt = numel(periT);
psth = zeros(w*h, nt);
for ix = 1:numel(onidx)
    
 ETAmat = magicETA(fT, F, onidx{ix}, periT);

 psth(ix, :) = nanmean(ETAmat, 1);
 
end
psth = my_conv(psth , 1);

psth = reshape(psth,  w,h, nt);

R = nanmean(psth(:, :, periT > 0),3) - nanmean(psth(:, :, periT <= 0), 3);

R = reshape(R, w, h);

R = imgaussfilt(R, 1.5, 'Padding', mean(R(:)));

% figure; imagesc(R); colormap jet

% [gfit, fitResp] = fitGauss2D(1:size(R, 2), ...
%     1:size(R, 1), R, x0, bounds, 1, doPlot);

[xPref, yPref, rfContours] = getRFContour(R, x0, doPlot);

x = [yPref, xPref];

mask = poly2mask(rfContours.xx, rfContours.yy, w,h);

psth = reshape(psth, [], nt);

tKernel = (mask(:)'*psth)/sum(mask(:));
end
