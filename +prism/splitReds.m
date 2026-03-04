function [dsRedImg, mCherryImg, redList, cc] = splitReds(redCh, greenCh, waveL,laserPw, z, doPlot)


if nargin <5
    
   doPlot = false;
end

if nargin <4 || isempty(laserPw)
    
else
    nPhotons = bsxfun(@times, laserPw/100, waveL.*prism.lookUpPower(waveL)); % nphotons ~ pw*waveL
    nPhotons = nPhotons.^2;
%     
%     if z>510
%     laserSNR =  prism.lookUpWaveSNR(waveL, z);
%     else
%         laserSNR = ones(size(waveL));
%     end
    if z>510
     laserSNR = prism.laserPenetration(waveL, z);
     laserSNR = laserSNR/mean(laserSNR);
     else
     laserSNR = ones(size(waveL));
             
    end
% %     nPhotons = nPhotons.*laserSNR;
%     nPhotons = nPhotons./mean(nPhotons(:));
% nPhotons= ones(size(waveL));
%     redCh = bsxfun(@rdivide, redCh, reshape(nPhotons, [1 1 size(nPhotons)]));
%     greenCh = bsxfun(@rdivide, greenCh, reshape(nPhotons, [1 1 size(nPhotons)]));
end
nPhotons = ones(size(nPhotons));

[mCherry_exc] = prism.excitationSpectrum('mCherry', waveL);
[dsRed_exc] = prism.excitationSpectrum('dsRed2', waveL);

%bscope red channel 607/70
[mCherry_ems] = prism.emissionSpectrum('mCherry', 572:642);
[dsRed_ems] = prism.emissionSpectrum('dsRed2', 572:642);

[G6Ca_ems, pmt] = prism.emissionSpectrum('dsRed2', 572:642);
[G6free_ems] = prism.emissionSpectrum('dsRed2', 572:642);

mCherry = mCherry_exc.*nansum(mCherry_ems.*pmt);
dsRed = dsRed_exc.*nansum(dsRed_ems.*pmt);

% ch1 and ch2 are nPxX X nPxY X nPlanes X nWaves
nCh = 2; 

[nX, nY, nPlanes, nWave] = size(redCh);

for iPlane = 1:nPlanes
    for iW = 1:nWave
         rr= imgaussfilt(redCh(:,:,iPlane,iW), 0.1);
         gg= imgaussfilt(greenCh(:,:,iPlane,iW), 0.1);
        [curedR{iPlane}(:,:, iW), bleedThrough, ~, redAll{iW}, globalThresh(iPlane, iW)] = prism.bleedCure(rr, gg,0);
%         thStack = zstack.removeBackground(curedR{iPlane}(:,:, iW), 0.1, [8 8 1]);
        [bwBloob, redList{iPlane, iW}, cc{iW}] = zstack.forgrBlobs(curedR{iPlane}(:,:, iW), globalThresh(iPlane, iW), [10 10 1], 0);
%          [bwBloob, redList{iPlane, iW}, cc{iW}] = zstack.forgrBlobs_dev(curedR{iPlane}(:,:, iW), globalThresh(iPlane, iW), [10 10 1], 0);

    end
    
    
    
%     ccc = cat(1, cc{:});
    
%     for ic = 1:numel(ccc)
%         for iw = 1:nWave
%         pxList = ccc(ic).PixelList;
%         pxList = sub2ind([nX, nY], pxList(:,2), pxList(:,1));
%         imgR = squeeze(redCh(:,:,iPlane,iw));
%         imgG = squeeze(greenCh(:,:,iPlane,iw));
%         dataR(ic, iw) = mean(imgR(pxList));
%         dataG(ic, iw) = mean(imgG(pxList));
%         
%         end    
%     end
%     dataAll = cat(2, dataR, dataG);
%     waveAll = cat(2, waveL, waveL);
    
    
%     redPx{iPlane} = unique(cat(1, redAll{:}));
      redPx{iPlane} = unique(cat(1, redList{:}));

%     redPx{iPlane} = redList{waveL == 890};

    redPxAcross{iPlane} = [];
    
    for iC = 1:nWave
        
        redPxAcross{iPlane} = cat(1,redPxAcross{iPlane}, redPx{iPlane} + (iC-1)*nX*nY);
        
    end
    
    colorPx{iPlane} = reshape(curedR{iPlane}(redPxAcross{iPlane}), numel(redPx{iPlane}), nWave);
iPlane
end

allPx = cat(1, colorPx{:});

allPx = bsxfun(@rdivide, allPx, nPhotons.*laserSNR);



%%
mixing = [dsRed(:), mCherry(:)];
mixing = bsxfun(@rdivide, mixing, mean(mixing(:)));

% [sources, mixing] = prism.learnSources(allPx', mixing, waveL);


sources = prism.nnlsSeparation(mixing, allPx, waveL, 890, 0);


dsRedImg = nan(512,512, nPlanes);
mCherryImg = nan(512,512, nPlanes);


redPxAcrossPlanes = [];
for iPlane = 1:nPlanes
    redPxAcrossPlanes =cat(1, redPxAcrossPlanes, redPx{iPlane}+ (iPlane-1)*nX*nY);
end

dsRedImg(redPxAcrossPlanes) = sources(1,:);
mCherryImg(redPxAcrossPlanes) = sources(2,:);

% dsRedImg(redPx) = sources(:,1);
% mCherryImg(redPx) = sources(:,2);
for iPlane = 1:nPlanes

dsRedImg(:,:,iPlane) = nanGaussFilt(dsRedImg(:,:,iPlane), 3,1);
mCherryImg(:,:,iPlane) = nanGaussFilt(mCherryImg(:,:,iPlane), 3,1);

end

dsRedImg(isnan(dsRedImg(:))) = 0;
mCherryImg(isnan(mCherryImg(:))) = 0;

imgPxList = redPx{iPlane}+(iPlane-1)*nX*nY;

dsRedImgPx = dsRedImg(imgPxList);
mCherryImgPx = mCherryImg(imgPxList);

% [justCherry, edges] = histcounts(dsRedImgPx, 1000);
% [~, justCherryThrsh] = max(diff(gaussFilt(justCherry',10)));
% justCherryThrsh = edges(justCherryThrsh-1);
% justCherryPx = mCherryImgPx(dsRedImgPx < justCherryThrsh);
% justRedThrsh = prctile(justCherryPx, 20);
% 
% notCherryPx = mCherryImgPx < justRedThrsh & dsRedImgPx >justCherryThrsh;
% 
% mCherryImg(imgPxList(notCherryPx)) = 0;


if doPlot

for iPlane = 1:nPlanes
lims1 = prctile(dsRedImg(imgPxList ), [1 99]);
lims2 = prctile(mCherryImg(imgPxList ), [1 99]);

figure;
r3 = subplot(2,3,1);
img1 = curedR{iPlane}(:,:, find(waveL ==890 | waveL == 910, 1, 'first'));
imagesc(img1); axis image;caxis(r3,[0 3*globalThresh(iPlane, find(waveL ==890| waveL == 910, 1, 'first'))])

r4 = subplot(2,3,2);
try
img2   = curedR{iPlane}(:,:, find(waveL ==780, 1, 'first'));
imagesc(img2); axis image; caxis(r4,[0 3*globalThresh(iPlane, find(waveL ==780, 1, 'first'))])
catch
img2 = curedR{iPlane}(:,:, find(waveL ==1020, 1, 'first'));
imagesc(img2); axis image; caxis(r4,[0 3*globalThresh(iPlane, find(waveL ==1020, 1, 'first'))])
end
r1 = subplot(2,3,4);
imagesc(dsRedImg(:,:,iPlane)); axis image;caxis(r1,lims1)
r2 = subplot(2,3,5);
imagesc(mCherryImg(:,:,iPlane)); axis image; caxis(r2,lims2)
linkaxes([r1 r2 r3 r4], 'xy')

subplot(2,3,6)
plot(dsRedImgPx, mCherryImgPx, '.'); axis square; hold on
% plot([justCherryThrsh justCherryThrsh], [0 max(mCherryImg(:))], '--m')
% plot([0 max(dsRedImg(:))],[justRedThrsh justRedThrsh],  '--r')

subplot(2,3,3)

plot(img1(:), img2(:), '.'); axis square; hold on


end
end
% figure
% h = scatter(dsRedImgPx,mCherryImgPx);
% hold on
% gmPDF = @(x1,x2)reshape(pdf(gmm,[x1(:) x2(:)]),size(x1));
% g = gca;
% fcontour(gmPDF,[g.XLim g.YLim])



% figure;
% h=scatterhist(dsRedImg(redPx{iPlane}+(iPlane-1)*nX*nY),mCherryImg(redPx{iPlane}+(iPlane-1)*nX*nY),100); axis square


% figure; 
% plot(dsRedImg(redPxAcrossPlanes),mCherryImg(redPxAcrossPlanes), '.'); axis square
% axis image
% 
% figure;
% r1 = subplot(1,2,1);
% imagesc(squeeze(redCh(:,:, 1, 1))); axis image
% r2 = subplot(1,2,2);
% imagesc(squeeze(redCh(:,:, 1, 3))); axis image
% linkaxes([r1 r2], 'xy')
% 
% figure;
% r1 = subplot(1,2,1);
% imagesc(curedR(:,:, 1)); axis image
% r2 = subplot(1,2,2);
% imagesc(curedR(:,:, 3)); axis image
% linkaxes([r1 r2], 'xy')

end

% saveastiff(int16(dsRedImg), 'C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\Data\FR141\R_prism.tif');
% saveastiff(int16(mCherryImg), 'C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\Data\FR141\M_prism.tif');
