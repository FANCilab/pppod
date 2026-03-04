function [dsRedImg, mCherryImg, curedR] =splitColorsVolume(redCh, greenCh, waveL, z, forceZero, rectifySource)

if nargin < 6
    rectifySource = 0;
end

if nargin < 6
    	forceZero = 0;
end


addpath(genpath('C:\Users\Federico\Documents\MATLAB\FastICA_2.5\FastICA_25'));


[fpR, fpG] = prism.mixFP({'dsRed2', 'mCherry'}, waveL);


% ch1 and ch2 are nPxX X nPxY X nPlanes X nWaves

nCh = 2;

[nX, nY, nPlanes, nWave] = size(redCh);

allImgs = [];
for iPlane = 1:nPlanes
    for iW = 1:nWave
         rr= imgaussfilt(redCh(:,:,iPlane,iW), 0.1);
        gg= imgaussfilt(greenCh(:,:,iPlane,iW), 0.1);
        %     [curedR(:,:, iW), bleedThrough, redList(:, iW), redAll{iW}] = prism.bleedCure(redCh(:,:,:,iW), greenCh(:,:,:,iW),0);
        [curedR{iPlane}(:,:, iW), ~,redList(iPlane, iW),~,globalThresh(iPlane, iW)] = prism.bleedCure(rr, gg,0);
%         curedR{iPlane}(:,:, iW) = zstack.removeBackground(curedR{iPlane}(:,:, iW), 0.5, [10 10 1]);
%         [bwBlob{iPlane, iW}, redList{iPlane, iW}] = zstack.forgrBlobs_dev(curedR{iPlane}(:,:, iW), [10 10 1]);
   
    end
    
    redPx{iPlane} = unique(cat(1, redList{:}));


    redPxAcross{iPlane} = [];
    
    for iC = 1:nWave
        
        redPxAcross{iPlane} = cat(1,redPxAcross{iPlane}, redPx{iPlane} + (iC-1)*nX*nY);
        
    end
    
    colorPx{iPlane} = reshape(curedR{iPlane}(redPxAcross{iPlane}), numel(redPx{iPlane}), nWave);
    allImgs = cat(1, allImgs, reshape(curedR{iPlane}, [], nWave));
iPlane
end

 allPx = cat(1, colorPx{:});

%%


[~, mixing] = prism.learnSources_iterative2(allPx', fpR,forceZero, 1);
%         [~, mixing] = prism.learnSources_iterative(nonZeroPx', fpR,forceZero, 1);

sources = pinv(mixing)*allImgs';

%         [~, mixing] = prism.learnSources_iterative2(colorPx{iPlane}', fpR, 0, 1);
%     sources = pinv(mixing)*colorPx{iPlane}';
%     nonZeroPx = find(sum(colorPx{iPlane}>0, 2));

dsRedImg = nan( nX,nY, nPlanes);
mCherryImg = nan( nX,nY, nPlanes);

dsRedImg(:) = sources(1,:);
mCherryImg(:) = sources(2,:);




%%

for iPlane = 1:nPlanes

figure;

r3 = subplot(2,3,1);

img1 = curedR{iPlane}(:,:, find(waveL >=890 & waveL <=980, 1, 'first'));
imagesc(img1); axis image; colorbar
caxis(r3,[0 3*globalThresh(iPlane, find(waveL >=890 & waveL <=980, 1, 'first'))])

r4 = subplot(2,3,2);
try
    img2   = curedR{iPlane}(:,:, find(waveL ==780, 1, 'first'));
    imagesc(img2); axis image; colorbar;
    caxis(r4,[0 3*globalThresh(iPlane, find(waveL ==1020, 1, 'first'))])
catch
    img2 = curedR{iPlane}(:,:, find(waveL ==1020, 1, 'first'));
    imagesc(img2); axis image; colorbar;
    caxis(r4,[0 3*globalThresh(iPlane, find(waveL ==1020, 1, 'first'))])
end

subplot(2,3,3)
plot(img1,img2, '.b'); axis square


lims1 = [0, max(makeVec(dsRedImg(:,:,iPlane)))*0.8];
lims2 = [0, max(makeVec(mCherryImg(:,:,iPlane)))*0.8];

lims = [0, max(lims1(2), lims1(2))];

r1 = subplot(2,3,4);
imagesc(dsRedImg(:,:,iPlane)); axis image; caxis(lims1);colorbar;
r2 = subplot(2,3,5);
imagesc(mCherryImg(:,:,iPlane)); axis image;caxis(lims2); colorbar;
linkaxes([r1 r2 r3 r4], 'xy')
subplot(2,3,6)
plot(dsRedImg(:,:,iPlane),mCherryImg(:,:,iPlane), '.b'); axis square
pause;
end

end