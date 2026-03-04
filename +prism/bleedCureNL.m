function  [curedR, bleedThrough, redList, redAll, globalThresh] = ...
    bleedCureNL(R, G, intercept, nBins, smoothing, doPlot)

if nargin <6 || isempty(doPlot)
    doPlot = false;
end

if nargin <5 || isempty(smoothing)
    smoothing = 0.000001; %0.000000
end
if nargin < 3 || isempty(intercept)
    intercept = 0;
end

if nargin <4 || isempty(nBins)
    nBins = 50;
end
[nPxX, nPxY, nImg] = size(G);

curedR = zeros(nPxX*nPxY, nImg);
redPx = zeros(nPxX*nPxY, nImg);
redList = cell(nImg, 1);
bleedThrough = nan(nImg, 2);
redAll= [];

for iM = 1:nImg

thisG = makeVec(G(:,:, iM)); thisR = makeVec(R(:,:, iM));


[~, fitY, lowY] = lowEnvelopeReg(thisG, thisR,nBins);

if intercept
    fitY = fitY';
    lowY = lowY';
else
    bleedThrough(iM, :) = robustfit(thisG, thisR);
    
    fitY = ([0, fitY])';
    lowY = ([bleedThrough(iM, 1), lowY])';
    
end


% lowY = gaussFilt(lowY, 2); 

g2r = fit(fitY, lowY, 'smoothingspline', 'SmoothingParam', smoothing);

if doPlot
    %%
figure; hold on
[xy, xbins, ybins] = cameraLucida.plot_Density2D (thisG, thisR, 30, 1, [0 8000 0 8000], 0, 0);
xy = xy/sum(xy(:));
% plot(thisG, thisR,'.'); plot(fitY, lowY); plot(g2r); axis image
imagesc(xbins, ybins, imgaussfilt(xy,2)); axis image; hold on;
caxis([0 0.000005])
% colormap(1-gray);
colormap(cat(1, [1 1 1], RedWhite(100),flip(Reds(100))));

plot(g2r);
xlim([0 5000])
ylim([0 5000])
xlabel('G')
ylabel('R')
formatAxes
%%
end

g2r = g2r(thisG);

curedR(:, iM) = thisR - g2r;

globalThresh(iM) = std(thisR);

redPx(:,iM) = curedR(:,iM) > globalThresh(iM) & curedR(:,iM)<3000;

redList{iM} = find(redPx(:,iM));

redAll = cat(1, redAll, redList{iM} + (iM-1)*nPxX*nPxY);

end



curedR = reshape(curedR, [nPxY, nPxX, nImg]);
redPx = reshape(redPx, [nPxY, nPxX, nImg]);


if doPlot
   rm =  cat(1, [1 1 1], cbrewer('seq', 'Reds',100));
   bm =  cat(1, [1 1 1], cbrewer('seq', 'Blues',100));
   gm =  cat(1, [1 1 1], cbrewer('seq', 'Greens',100));
%    gm =  Green(100);
%%
    figure;
    for iM = 1:nImg
        
    r1= subplot(nImg,4,1 + (iM-1)*4);
    imagesc(imgaussfilt(R(:,:,iM))); axis image; colormap(r1, rm)
    caxis(prctile(makeVec(R(:,:,iM)), [30, 98]));
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    xlim([170 470]);ylim([80 380]);
    
    g1= subplot(nImg,4,2+ (iM-1)*4);
    imagesc(imgaussfilt(G(:,:,iM))); axis image; colormap(g1, gm)
    caxis(prctile(makeVec(G(:,:,iM)), [15, 95]));
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    xlim([170 470]);ylim([80 380]);

    r2 = subplot(nImg,4,3+ (iM-1)*4);
    imagesc(imgaussfilt(curedR(:,:,iM))); axis image; colormap(r2, rm)
    caxis( [25 prctile(makeVec(R(:,:,iM)), [75])]);
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    xlim([170 470]);ylim([80 380]);

    r3 = subplot(nImg,4,4+ (iM-1)*4);
    imagesc(redPx(:,:,iM)); axis image; colormap(r3, flip(gray,1))
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    xlim([170 470]);ylim([80 380]);

    end
    
end
%%
end