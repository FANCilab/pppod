function  [curedR, bleedThrough, redList, redAll, globalThresh] = bleedCure(R, G, doPlot)

if nargin <3
    doPlot = false;
end

[nPxY, nPxX, nImg] = size(G);

curedR = zeros(nPxX*nPxY, nImg);
redPx = zeros(nPxX*nPxY, nImg);
redList = cell(nImg, 1);
bleedThrough = nan(nImg, 2);
redAll= [];

for iM = 1:nImg

thisG = makeVec(G(:,:, iM)); thisR = makeVec(R(:,:, iM));

bleedThrough(iM, :) = robustfit(thisG, thisR);

if doPlot
    figure;
    plot(thisG, thisR, '.');
    hold on
    plot([0 max(thisG(:))], [0 max(thisG(:))]*bleedThrough(iM, 2) + bleedThrough(iM, 1), 'r-');
end
curedR(:, iM) = thisR - bleedThrough(iM, 2)*thisG - bleedThrough(iM, 1);

globalThresh(iM) = 2*std(thisR);

redPx(:,iM) = curedR(:,iM) > globalThresh(iM);

redList{iM} = find(redPx(:,iM));

redAll = cat(1, redAll, redList{iM} + (iM-1)*nPxX*nPxY);

end



curedR = reshape(curedR, [nPxY, nPxX, nImg]);
redPx = reshape(redPx, [nPxY, nPxX, nImg]);


if doPlot
   
    figure;
    for iM = 1:nImg
    subplot(nImg,4,1 + (iM-1)*4)
    imagesc(R(:,:,iM)); axis image; %colormap(flip(gray,1))
    caxis(prctile(makeVec(G(:,:,iM)), [5, 95]));
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    subplot(nImg,4,2+ (iM-1)*4)
    imagesc(G(:,:,iM)); axis image;%colormap(flip(gray,1))
    caxis(prctile(makeVec(G(:,:,iM)), [5, 95]));
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    subplot(nImg,4,3+ (iM-1)*4)
    imagesc(curedR(:,:,iM)); axis image;%colormap(flip(gray,1))
    caxis(prctile(makeVec(G(:,:,iM)), [5, 95]));
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    subplot(nImg,4,4+ (iM-1)*4)
    imagesc(redPx(:,:,iM)); axis image;%colormap(flip(gray,1))
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])

    end
    
end

end