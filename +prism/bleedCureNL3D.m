function  [curedR, bleedThrough, g2r] = ...
    bleedCureNL3D(R, G, intercept, nBins, model_type, doPlot)

if nargin <6 || isempty(doPlot)
    doPlot = false;
end

if nargin <5 || isempty(model_type)
   model_type = 'linear';
end
if nargin < 3 || isempty(intercept)
    intercept = 0;
end

if nargin <4 || isempty(nBins)
    nBins = 100;
end

[nPxX, nPxY, nImg] = size(G);

curedR = zeros(nPxX*nPxY, nImg);
redPx = zeros(nPxX*nPxY, nImg);
redList = cell(nImg, 1);
bleedThrough = nan(nImg, 2);
redAll= [];

G = makeVec(G); R = makeVec(R);

non_sat = G <40000 & R <40000 & G >2000;

thisG = G(non_sat); thisR= R(non_sat);


[~, fitY, lowY] = lowEnvelopeReg(thisG, thisR,nBins, 1);

if intercept
    fitY = fitY';
    lowY = lowY';
else
    bleedThrough = robustfit(thisG, thisR);
    fitY = ([0, fitY])';
    lowY = ([bleedThrough(1), lowY])';
    
end


% lowY = gaussFilt(lowY, 2); 
switch model_type
    case 'linear'
        g2r = fit(fitY, lowY, 'poly1');
    case 'spline'
        smoothing = 0.00001;
        g2r = fit(fitY, lowY, 'smoothingspline', 'SmoothingParam', smoothing);
end

if doPlot
    %%
figure; hold on
[xy, xbins, ybins] = plot_Density2D (thisG, thisR, 30, 1, [0 40000 0 40000], 0, 0);
xy = xy/sum(xy(:));
% plot(thisG, thisR,'.'); plot(fitY, lowY); plot(g2r); axis image
imagesc(xbins, ybins, imgaussfilt(xy,2)); axis image; hold on;
caxis([0 0.000005])
% colormap(1-gray);
colormap(cat(1, [1 1 1], RedWhite(100),flip(Reds(100))));

plot(g2r);
xlim([0 40000])
ylim([0 40000])
xlabel('G')
ylabel('R')
formatAxes
%%
end

g2r = g2r(G);

curedR = R - g2r;

% globalThresh = std(R);
% 
% redPx = curedR > globalThresh & curedR<3000;
% 
% redList = find(redPx);
% 
% redAll = redList;

curedR = reshape(curedR, [nPxY, nPxX, nImg]);

G = reshape(G, [nPxY, nPxX, nImg]);
R = reshape(R, [nPxY, nPxX, nImg]);

redPx = reshape(redPx, [nPxY, nPxX, nImg]);

if doPlot

   maxRc = max(curedR, [], 3); maxG = max(G, [], 3); maxR = max(R, [], 3);
   rm =  cat(1, [1 1 1], cbrewer('seq', 'Reds',100));
   bm =  cat(1, [1 1 1], cbrewer('seq', 'Blues',100));
   gm =  cat(1, [1 1 1], cbrewer('seq', 'Greens',100));
   gm =  Green(100);
%%
    figure;
    r1= subplot(1,3,1);
    imagesc(imgaussfilt(maxR)); axis image; colormap(r1, rm)
    caxis(prctile(makeVec(maxR), [30, 98]));
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    xlim([170 470]);ylim([80 380]);

    g1= subplot(1,3,2);
    imagesc(imgaussfilt(maxG)); axis image; colormap(g1, gm)
    caxis(prctile(makeVec(maxG), [15, 95]));
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    xlim([170 470]);ylim([80 380]);

    r2 = subplot(1,3,3);
    imagesc(imgaussfilt(maxRc)); axis image; colormap(r2, rm)
    caxis( [prctile(makeVec(maxRc), [30, 98])]);
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    xlim([170 470]);ylim([80 380]);


end
%
end