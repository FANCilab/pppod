%% create and save ImagingPlane file for each experiment;

clear;
% npRet.FR141_dbNpRF_1;
% npRet.FR141_dbNpRF_2;
% npRet.FR140_dbNpRF_1;
npRet.FR143_dbNpRF_1;

root = 'D:\OneDrive - University College London\Data\2P';

[frameG, frameR, micronsX, micronsY, micronsZ] = s2pUtils.loadMultiPlaneMeanStack(db);

for iDb = 1:numel(db)
    folder = buildExpPath(db(iDb));
    load(fullfile(root, folder, 'mimgSingleExp'));
    
    for iExp = 1:numel(db(iDb).expts)
        
        info(iDb,iExp) = ppbox.infoPopulateTempLFR(db(iDb).mouse_name, db(iDb).date, db(iDb).expts(iExp));
        laserPw{iDb}(:,iExp) = ppbox.getLaserPower(info(iDb,iExp));
        
    end
    %     m = matfile(fullfile(root, folder, 'mimgSingleExp'),'Writable',true);
    waveL{iDb} = db(iDb).waveL ;
    prismR{iDb} = mimgSingleExp.R;
    prismG{iDb} = mimgSingleExp.G;
    
    for iPlane = 1:info(iDb).nPlanes
        planeCount = iPlane + (iDb-1)*info(iDb).nPlanes;
        zPlane(planeCount) = mean(micronsZ{iDb}(:, iPlane));   
    end
      
end

saveTo = 'C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\Data\rvRetinotopy';
save(fullfile(saveTo, ['ImagingPlanes_', db(1).mouse_name, '_', num2str(db(1).starterID)]), 'frameG', 'frameR','micronsX', 'micronsY', 'micronsZ', 'db',...
    'waveL', 'laserPw', 'prismR', 'prismG', 'info','zPlane', '-v7.3');

% save(fullfile(saveTo, ['ImagingPlanes_', db(1).mouse_name, '_', num2str(db(1).starterID)]),  'zPlane', '-append');


%% load ImagingPlanes for one dataset
clear;
cameraLucida.db_INsorting2;

nDb = numel(db);

if ispc
    saveTo = 'C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\Data\rvRetinotopy';
    addpath('C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\FedericoBox\2P\Tools');
else
    saveTo = '/Users/Federico/Google Drive/CarandiniLab/CarandiniLab_MATLAB/Data/rvRetinotopy';
    addpath('/Users/Federico/Google Drive/CarandiniLab/CarandiniLab_MATLAB/FedericoBox/2P/Tools');
end

iDb =3;

file = ['ImagingPlanes_', db(iDb).animal, '_', num2str(db(iDb).starterID)];
load(fullfile(saveTo,  file), 'frameG', 'frameR','micronsX', 'micronsY', 'micronsZ', ...
    'waveL','prismR', 'prismG', 'laserPw', 'zPlane');
imgDb = load(fullfile(saveTo,  file), 'db');


%% Split mCherry and dsRed fluorescence

%FR141_2
rectifySource = [1 1 0 0 0]; 
forceRatio = [0 0 0 0 0];
zLimit = Inf; 
useAllPx = [0 0 0 0 0]; %plane 20 1, plane 26 = 1 

%FR141_1
rectifySource = [0 0 0 0 0];
forceRatio = [0 0 1 1 1];
zLimit = 510;
useAllPx = [0 0 1 1 0];

%FR140
rectifySource = [1 0 0 0 0]; 
forceRatio = [0 0 0 0 0];
zLimit = Inf; 
useAllPx = [0 0 0 1 1];

%FR143
rectifySource = [0 0 0 0 0]; 
forceRatio = [0 0 0 0 0]; 
zLimit = Inf; 
useAllPx = [0 0 0 0 1];


%for supplementary fig use FR140 exp 3
for iExp = 1:numel(prismR)
    nPlanes = size(prismR{iExp},3);
    regrOps(iExp).rectifySource = rectifySource(iExp);
    regrOps(iExp).forceRatio = forceRatio(iExp);
    regrOps(iExp).zLimit = zLimit;
    regrOps(iExp).useAllPx = useAllPx(iExp);
    for iPlane = 1:nPlanes
        planeCount = iPlane + (iExp-1)*nPlanes;
        redCh=prismR{iExp}(:,:,iPlane, :);
        greenCh=prismG{iExp}(:,:,iPlane, :);
        [dsRedImg(:,:,planeCount), mCherryImg(:,:,planeCount), curedRed{planeCount}, mixing{planeCount}, iterErr(:,planeCount)] =...
            prism.splitColors(redCh, greenCh, waveL{iExp}, laserPw{iExp}(iPlane,:), zPlane(planeCount),  regrOps(iExp));
    end
end


% saveastiff(int16(mat2gray(dsRedImg).*(2.^15-1)), 'r' );
% saveastiff(int16(mat2gray(mCherryImg).*(2.^15-1)), 'm' );
file = ['RedCherry_', db(iDb).animal, '_', num2str(db(iDb).starterID)]; 
save(fullfile(saveTo, file),  'dsRedImg', 'mCherryImg', 'regrOps', 'curedRed', 'mixing', 'iterErr', '-v7.3');

%%

file = ['RedCherry_', db(iDb).animal, '_', num2str(db(iDb).starterID)]; 
load(fullfile(saveTo, file),  'dsRedImg', 'mCherryImg', 'regrOps', 'curedRed', 'mixing', 'iterErr');


 flatRedImg = dsRedImg; 
 flatCherryImg = mCherryImg; 
 
for iPlane = 1:numel(mixing)
    totFluo(iPlane) = mean(sum(mixing{iPlane},2));
    
    flatRedImg(:,:,iPlane) = zstack.removeBackground(dsRedImg(:,:,iPlane), 1, [10 10 1]);
    
    flatCherryImg(:,:,iPlane) = zstack.removeBackground(mCherryImg(:,:,iPlane), 1, [10 10 1]);
    
end
% flatRedImg(flatRedImg<0) = 0;
% flatCherryImg(flatCherryImg<0) = 0;

figure;
fitFluoVsZ = fit(zPlane', totFluo', 'exp1');
plot(zPlane,totFluo, 'o'); hold on; plot(1:max(zPlane), fitFluoVsZ(1:max(zPlane)))
fitFluoVsZ = fitFluoVsZ(zPlane);

nflatRedImg = bsxfun(@rdivide, flatRedImg, reshape(fitFluoVsZ, 1,1,[]))*mean(fitFluoVsZ);
nflatCherryImg = bsxfun(@rdivide, flatCherryImg, reshape(fitFluoVsZ, 1,1,[]))*mean(fitFluoVsZ);

file = ['RedCherry_', db(iDb).animal, '_', num2str(db(iDb).starterID)]; 
save(fullfile(saveTo, file),  'nflatRedImg', 'nflatCherryImg', 'fitFluoVsZ', '-append');

% saveastiff(int16(mat2gray(flatRedImg).*(2.^15-1)), 'r' );
% saveastiff(int16(mat2gray(flatCherryImg).*(2.^15-1)), 'm' );

% %% process images by volume
% for iExp = 1:numel(prismR)
%     planeCount = (iExp-1)*nPlanes;
%     nPlanes = size(prismR{iExp},3);
%     redCh=prismR{iExp};
%     greenCh=prismG{iExp};
%     [dsRedImg(:,:,planeCount+1:planeCount+nPlanes), mCherryImg(:,:,planeCount+1:planeCount+nPlanes)] =...
%         prism.splitColorsVolume(redCh, greenCh, waveL{iExp}, [],  forceZero(iExp), rectifySource(iExp));
% end



%%
file = ['RedCherry_', db(iDb).animal, '_', num2str(db(iDb).starterID)]; 
load(fullfile(saveTo, file),  'nflatRedImg', 'nflatCherryImg', 'fitFluoVsZ');

gRedImg = mat2gray(nflatRedImg, [0 max(nflatRedImg(:))]);
gCherryImg = mat2gray(nflatCherryImg, [0 max(nflatCherryImg(:))]);

for img = 1: size(mCherryImg, 3)
[dsRedCells(:,:, img), nbRed(:,:,img),redList{img}, rr{img}, tR(img)] = ...
    zstack.forgrBlobs_dev(gRedImg(:,:, img),  [10 10 1], 4, 0);

[mCherryCells(:,:, img),nbCherry(:,:,img), redList{img}, mc{img}, tM(img)] = ...
    zstack.forgrBlobs_dev(gCherryImg(:,:, img), [8 8 1], 5, 0);
end

nbRed = mat2gray(gRedImg);
nbCherry = mat2gray(gCherryImg);


 
% multiplets = findNeuronMultiplets_dev(xyz, xyT, zT, planePos, resps, corrThr);



%%

% load(fullfile(saveTo, ['ImagingPlanes_', db(1).animal, '_', num2str(db(1).starterID)]), 'prismG');

gCampImg =[];
for iExp = 1:numel(prismG)
    gCampImg = cat(3, gCampImg, mean(prismG{iExp}(:,:,:, waveL{iExp} ~= 780), 4));
end

zPx = repmat(reshape(cat(3, micronsZ{:}), 512, 1, []), [1 512, 1]);

[cropG, cropR, cropC, cropM, choices, neunCc, accepted, multiplets, allMulti]= ...
    cameraLucida.classifyRedCherry(gCampImg,dsRedImg, mCherryImg, rr,zPx, max(gRedImg,[], 3), imgDb.db(1).starterYXPlane(1:2));

file = ['RedCherry_', db(iDb).animal, '_', num2str(db(iDb).starterID)]; 
save(fullfile(saveTo, file), 'cropG', 'cropR', 'cropC', 'cropM', 'choices', 'neunCc', 'accepted', 'multiplets','allMulti','-append');
% load(fullfile(saveTo, file), 'cropR', 'cropC', 'cropM', 'choices', 'neunCc', 'accepted', 'multiplets','allMulti');

% save(fullfile(saveTo, 'INsorting'), 'thStack', 'zcStack', 'G', 'R', 'CC', 'neurons', 'bwBlob', 'pxSz', '-v7.3')
% % 
%  [allMulti] = prism.fixAccepted(dsRedImg, mCherryImg, rr,zPx, max(gRedImg,[], 3), imgDb.db(1).starterYXPlane(1:2));
% for im = 1:numel(multiplets)
%     [~, n] = min(abs(sum(reshape(bsxfun(@minus, allMulti.G{im},cropR(:,:,im)),[], numel(multiplets{im})),1)));
%     accepted2(im, :)= allMulti.pos(multiplets{im}(n),:);    
% end
% save(fullfile(saveTo, file), 'cropR', 'cropC', 'cropM', 'choices', 'neunCc', 'accepted', 'multiplets','allMulti', '-append');


%%
allIN = find(choices =='p');
figure;
for in = 1:numel(allIN)
       clf;
        g=subplot(1,3,1);        
        imagesc(cropG(:,:,allIN(in))); axis image; colormap(g, Green); hold on
         caxis([ 0 prctile(makeVec(cropG(:,:,allIN(in))), 99)])
        contour(any(cropM(:,:,allIN(in)),3),1,'r', 'Linewidth', 0.2)
        
        r=subplot(1,3,2);  
       imagesc(cropR(:,:,allIN(in))); axis image; colormap(r,  Red); hold on
        caxis([ 0 prctile(makeVec(cropR(:,:,allIN(in))), 99)])
        contour(any(cropM(:,:,allIN(in)),3),1,'r', 'Linewidth', 0.2)
          
        c=subplot(1,3,3);      

        imagesc(cropC(:,:,allIN(in))); axis image; colormap(c, Red); hold on
        caxis([ 0 prctile(makeVec(cropC(:,:,allIN(in))), 99)])
        contour(any(cropM(:,:,allIN(in)),3),1,'r', 'Linewidth', 0.2)
        title(sprintf('%d', allIN(in)));
        pause;
end
%%
in = [9 16 27 42 57 80 94 102 105];
pc = [114 116 131 149 160 161 167 169 255 331];
saveastiff(int16(cropG(:,:, pc)*8000), 'FR140_1_pcsG')

%%

% 
% cC = [];
% cR = [];
% rR = [];
% rC = [];
% redCoord =[];
% cherryCoord = [];
% 
% nShuffle = 100;
% for iS = 1:nShuffle
% cherrySh(:, iS) = randperm(numel(mc));
% end
% 
% for iR = 1:numel(rr)
%     redCoord(iR, :) = rr(iR).Centroid;
%     valid = rr(iR).PixelList(:,3) == round(rr(iR).Centroid(3));
%     thisPx = sub2ind(size(dsRedImg), rr(iR).PixelList(valid,2), rr(iR).PixelList(valid,1),rr(iR).PixelList(valid,3));
%     %     thisPx = sub2ind(size(dsRedImg), rr(iR).PixelList(:,2), rr(iR).PixelList(:,1),rr(iR).PixelList(:,3));
%     rC(iR) = mean(nbCherry (thisPx));
%     rR(iR) = mean(nbRed(thisPx));
% end
% for iC = 1:numel(mc)
%     cherryCoord(iC, :) = mc(iC).Centroid;
%     
%     valid = mc(iC).PixelList(:,3) == round(mc(iC).Centroid(3));
%     thisPx = sub2ind(size(dsRedImg), mc(iC).PixelList(valid ,2), mc(iC).PixelList(valid ,1),mc(iC).PixelList(valid ,3));
%     %     thisPx = sub2ind(size(dsRedImg), mc(iC).PixelList(: ,2), mc(iC).PixelList(:,1),mc(iC).PixelList(:,3));
%     
%     cC(iC) = mean(nbCherry (thisPx));
%     cR(iC) = mean(nbRed(thisPx));
%     
%     for iS = 1:nShuffle
%      
%      shZ= repmat(round(mc(cherrySh(iC, iS)).Centroid(3)), numel(thisPx),1);
%      thisPx = sub2ind(size(dsRedImg), mc(iC).PixelList(valid ,2), mc(iC).PixelList(valid ,1),shZ);
%      shC(iC, iS) = mean(nbCherry (thisPx));
%      shR(iC, iS) = mean(nbRed(thisPx));
%     end
% end
% 
% 
% redCoord(:, 3) = cameraLucida.zfromplane(redCoord, micronsZ, 50, 0);
% cherryCoord(:, 3) = cameraLucida.zfromplane(cherryCoord, micronsZ, 50, 0);
% 
% 
% [rf, envelope, envelopeEdges] = quantileReg(shR(:),shC(:), [0 25 50 75 100]);
% %  [rf, envelope, envelopeEdges] =lowEnvelopeReg(shR(:),shC(:), 10);
% figure; hold on; plot(shR(:), shC(:), 'ok');axis image;
% plot([0 1], [0 1]*rf(2) +rf(1), 'r-')
% 
% figure; hold on; plot(rR, rC, 'or');plot(cR, cC, 'ob');axis image; 
% plot([0 1], [0 1]*rf(2) +rf(1), 'k-')
% 
% rIsIN = rC > rR*rf(2) +rf(1);
% 
% % rIsIN = rC>0.03 & rR <0.2;
% 
% figure; plot(redCoord(~rIsIN,1), 513-redCoord(~rIsIN,2), 'or'); hold on
% plot(redCoord(rIsIN,1), 513-redCoord(rIsIN,2), 'ob'); axis image; 
% xlim([0 512]); ylim([0 512]);
% 
% figure; plot3(redCoord(~rIsIN,1), 513-redCoord(~rIsIN,2), -redCoord(~rIsIN,3), 'or'); hold on
% plot3(redCoord(rIsIN,1), 513-redCoord(rIsIN,2),-redCoord(rIsIN,3), 'ob'); axis image; 
% xlim([0 512]); ylim([0 512]);
% 
% % figure; plot3(cherryCoord(:,1), 513-cherryCoord(:,2), -cherryCoord(:,3), 'or'); hold on
% % xlim([0 512]); ylim([0 512]);
% % 
% % figure; hist(cherryCoord(:,3), 10)
% 
% file = ['RedCherry_', db(iDb).animal, '_', num2str(db(iDb).starterID)]; 
% save(fullfile(saveTo, file),  'redCoord', 'cherryCoord', 'rC','cC','rR','cR','rIsIN','shC','shR','rr','cc','-append');

%%

%% pool all data together
% load(fullfile(saveTo, 'INsorting'));

%% load ImagingPlanes for one dataset
clear;
cameraLucida.db_INsorting2;

nDb = numel(db);

if ispc
    saveTo = 'C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\Data\rvRetinotopy';
    addpath('C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\FedericoBox\2P\Tools');
else
    saveTo = '/Users/Federico/Google Drive/CarandiniLab/CarandiniLab_MATLAB/Data/rvRetinotopy';
    addpath('/Users/Federico/Google Drive/CarandiniLab/CarandiniLab_MATLAB/FedericoBox/2P/Tools');
end



for iDb = 1:numel(db)
    INsort(iDb) = load(fullfile(saveTo, ['RedCherry_', db(iDb).animal, '_', num2str(db(iDb).starterID)]), ...
        'cropR', 'cropC', 'cropM', 'choices', 'neunCc', 'accepted', 'multiplets','allMulti');
%     load(fullfile(saveTo, ['ImagingPlanes_', db(iDb).animal, '_', num2str(db(iDb).starterID)]), 'micronsZ', 'frameG');
%     if iDb ==1
%         INsort(iDb).accepted(:, 3)= cameraLucida.zfromplane(INsort(iDb).accepted, micronsZ, size(frameG,3), 0);
%     else
%         INsort(iDb).accepted(:, 3)= cameraLucida.zfromplane(INsort(iDb).accepted, micronsZ, size(frameG,3), 1);
%         
%     end

end

for iDb = 1:numel(db)
file = ['ImagingPlanes_', db(iDb).animal, '_', num2str(db(iDb).starterID)];
imgDb = load(fullfile(saveTo,  file), 'db');

  INsort(iDb).starterYX = imgDb.db(1).starterYXPlane(1:2);
  clear imgDb
end
% 
% for iDb = 1:numel(db)    
%      load(fullfile(saveTo, ['ImagingPlanes_', db(iDb).animal, '_', num2str(db(iDb).starterID)]), 'micronsZ', 'framesG');
%      Z= cameraLucida.zfromplane(micronsZ, size(frameG,3));
% end 

%%

rIns = [];
rPc = [];
cIns = [];
cPc = [];
% figure;
for iDb = 1:numel(db)
    
    for iN = 1:numel(INsort(iDb).multiplets)
        clf;
        thisR = INsort(iDb).cropR(:,:,iN);
        thisC = INsort(iDb).cropC(:,:,iN);
        thisM = INsort(iDb).cropM(:,:,iN);
        [ii,jj] = meshgrid ( -35:35, -35:35 );
        thisRing = (ii.^2+jj.^2) < 10^2 & (ii.^2+jj.^2) >= 2;
        thisRing = thisRing & ~thisM;
        rR(iN) = mean(makeVec(thisR(INsort(iDb).cropM(:,:,iN))))/mean(makeVec(thisR(thisRing)));
        rC(iN) = mean(makeVec(thisC(INsort(iDb).cropM(:,:,iN))))/mean(makeVec(thisC(thisRing)));
        
%         r = subplot(1,2,1); 
%         rm =  cbrewer('seq', 'Reds',100);
%         imagesc(thisR); axis image; colormap gray; caxis([0, 0.7]); colormap(r, rm)
%         set(gca, 'XTickLabel', [], 'YTickLabel', []); formatAxes;title('dsRed');
%         c = subplot(1,2,2);
%         cm =  cbrewer('seq', 'Blues',100);
%         imagesc(thisC);axis image; colormap gray; caxis([0, 0.7]); colormap(c, cm)
%         set(gca, 'XTickLabel', [], 'YTickLabel', []); formatAxes;title('mCherry');
%         pause;

    end
    
    rPc = cat(2, rPc, rR(INsort(iDb).choices == 'p'));
    cPc = cat(2, cPc, rC(INsort(iDb).choices == 'p'));
    rIns = cat(2, rIns, rR(INsort(iDb).choices == 'i'));
    cIns = cat(2, cIns, rC(INsort(iDb).choices == 'i'));
    clear rR rC
end

figure; 
plot(rPc, cPc, 'or'); hold on 
plot(rIns, cIns, 'ob'); hold on; axis square
xlabel('dsRed')
ylabel('mCherry')


%%

% figure;
% r = subplot(1,2,1);
% rm =  cbrewer('seq', 'Reds',100);
% imagesc(thisR); axis image; colormap gray; caxis([0, prctile(thisR(:), 99)]); colormap(r, rm)
% set(gca, 'XTickLabel', [], 'YTickLabel', []); formatAxes;
% c = subplot(1,2,2);
% cm =  cbrewer('seq', 'Blues',100);
% imagesc(thisC);axis image; colormap gray; caxis([0, prctile(thisC(:), 99)]); colormap(c, cm)
% set(gca, 'XTickLabel', [], 'YTickLabel', []); formatAxes;
% pause;

%%
poolChoice = [];
poolNeurons = [];
poolImgG = [];
poolImgR = [];
poolImgM = [];
poolInCoord = [];
poolPcCoord = [];


poolChoice = [INsort.choices];
poolImgG = cat(3,INsort.cropR);
poolImgR = cat(3,INsort.cropC);
poolImgM = cat(3,INsort.cropM);
 
%  goodDb = ~isnan(zoom);
 
for iDb = 1:numel(db)
[pxSzX, pxSzY] = ppbox.zoom2fov(db(iDb).zoom);

% neuron and starter coordinates are in ij system. neuron(:,1) is j,
% neuron(:,2) is i. StarterYX(1) is i, Starter(2) is j.

INsort(iDb).accepted(:, 2) = 513 - INsort(iDb).accepted(:, 2);
INsort(iDb).starterYX(1) = 513 - INsort(iDb).starterYX(1);

relCoord = bsxfun(@minus, INsort(iDb).accepted(:, 1:2),INsort(iDb).starterYX([2 1]));
zCoord = INsort(iDb).accepted(:, 3);

relCoord(:,1) = relCoord(:,1)*pxSzX/512;
relCoord(:,2) = relCoord(:,2)*pxSzY/512;
poolNeurons = cat(1,poolNeurons, relCoord);

acceptedChoices = INsort(iDb).choices(INsort(iDb).choices == 'i' | INsort(iDb).choices == 'p' | INsort(iDb).choices == 'x'| INsort(iDb).choices == 'd');
% coordinate where rows * cols, now transform in X and Y in cortical space
inCoord{iDb} = relCoord(acceptedChoices == 'i', [1 2]);
pcCoord{iDb} = relCoord(acceptedChoices == 'p',  [1 2]);

inCoord{iDb}(:,3) = zCoord(acceptedChoices == 'i');
pcCoord{iDb}(:,3) = zCoord(acceptedChoices== 'p');

poolInCoord = [poolInCoord; inCoord{iDb}];
poolPcCoord = [poolPcCoord; pcCoord{iDb}];

end

% cameraLucida.plotInPcClust(InPcKlust);

% cameraLucida.plotInPcNumber(inCoord, pcCoord);
%% plot individual networks

savePlot = 0;
normFlag = 1;

for iDb = 1:numel(db)
    if savePlot
        latFig = fullfile(saveTo, [db(iDb).animal,'_',num2str(db(iDb).starterID),'_latDist']);
        layerFig = fullfile(saveTo, [db(iDb).animal,'_',num2str(db(iDb).starterID),'_layerDist']);
    else
        latFig = [];
        layerFig = [];
    end
    [hatN(:,:,iDb),hatIN(:,:,iDb), hatPC(:,:,iDb), mexHat(:,:,iDb)]=...
        cameraLucida.plotInPcKlust(inCoord{iDb}(:, 1:2), pcCoord{iDb}(:, 1:2), normFlag, latFig);
    
    [zDistN(:,:,iDb),zDistIN(:,:,iDb), zDistPC(:,:,iDb), zDistDiff(:,:,iDb)] =...
        cameraLucida.plotInPcKlustZ(inCoord{iDb}, pcCoord{iDb}, db(iDb).starterZ, normFlag, layerFig);
    
    
end
%% plot average

savePlot = 0;
normFlag = 1;

if savePlot
    latFig = fullfile(saveTo, 'Pool_latDist');
    layerFig = fullfile(saveTo, 'Pool_layerDist');
    latProfileFig = fullfile(saveTo, 'Profile_latDist');
    layerProfileFig =fullfile(saveTo, 'Profile_layerDist');
else
    latFig = [];
    layerFig = [];
    latProfileFig = [];
    layerProfileFig = [];
end
[poolHatN,poolHatIN, poolHatPC, poolMexHat, latBins] = ...
    cameraLucida.plotInPcKlust(poolInCoord(:, 1:2), poolPcCoord(:, 1:2), normFlag, latFig);

[zPoolDistN,zPoolDistIN, zPoolDistPC, zPoolDistDiff, zBins] =...
    cameraLucida.plotInPcKlustZ(poolInCoord, poolPcCoord,[db.starterZ], normFlag, layerFig);

cameraLucida.plot_singleVSavgNet(hatN, hatIN, hatPC, mexHat, poolHatN, poolHatIN, poolHatPC, poolMexHat, latBins, 0,latProfileFig )

cameraLucida.plot_singleVSavgNet(zDistN,zDistIN, zDistPC, zDistDiff, zPoolDistN,zPoolDistIN, zPoolDistPC, zPoolDistDiff, zBins, 1, layerProfileFig)

%% convert coord to Ret coord
spatConst = [134.0884, 133.7699, 176.9610, 145.3904, 136.0761, 85.9435, 117.5160, 92.1770, 128.9015, 151.1016, 169.6434, 140.6265, 136.9703]; %presynaptics;

allRotIn =  [];
allRotPc = [];
allRotRet = [];
allRetIn =  [];
allRetPc = [];
allRetRet = [];
allRetCoord = [];
done = 0;
for iDb = 1:numel(db)
    
    %     if ~isnan(db(iDb).prefDir) && ~isnan(spatConst(iDb))
    if ~isnan(db(iDb).prefDir)
        
        done = done + 1;
        coord= inCoord{iDb};
        coord(:,2) = -coord(:,2); %invert y coordinate as image microns coordinates are in ij
        
        [retInX{iDb}, retInY{iDb}, retColX{iDb}, retColY{iDb}, dxGrid, dyGrid,~,~,~, sRetX(iDb),sRetY(iDb), sX(iDb),sY(iDb)] ...
            = npRet.getRetDistFromCoords(db(iDb), coord,[], db(iDb).prefDir, 0);

        coord= pcCoord{iDb};
        coord(:,2) = -coord(:,2);
        [retPcX{iDb}, retPcY{iDb}] ...
            = npRet.getRetDistFromCoords(db(iDb), coord, [], db(iDb).prefDir, 0);
        
        theta = db(iDb).prefDir; % to rotate counterclockwise
        R{iDb} = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
        % Rotate ret coords
        
%                 ret{iDb} = [makeVec(retColX{iDb}-sRetX(iDb)), makeVec(-retColY{iDb}+sRetY(iDb))];
        [retColX{iDb}, retColY{iDb}, ~,~,xc{iDb}, yc{iDb}] = npRet.nullRetDist(db(iDb), 100000, 1, spatConst(iDb));
        ret{iDb} = [makeVec(retColX{iDb}), makeVec(-retColY{iDb})];
        retCoord{iDb} = [makeVec(xc{iDb}), makeVec(yc{iDb})];
        %
        retIn{iDb} = [retInX{iDb}, -retInY{iDb}]; % invert Y ret cos mpep maps upside down
        retPc{iDb} = [retPcX{iDb}, -retPcY{iDb}];
        
        rotRet{iDb} = R{iDb}*ret{iDb}';
        rotIn{iDb} = R{iDb}*retIn{iDb}';
        rotPc{iDb} = R{iDb}*retPc{iDb}';
        
        allRotRet = cat(2, allRotRet, rotRet{iDb});
        allRotIn =  cat(2, allRotIn, rotIn{iDb});
        allRotPc = cat(2, allRotPc, rotPc{iDb});
        
        allRetRet = cat(1, allRetRet, ret{iDb});
        allRetIn =  cat(1, allRetIn, retIn{iDb});
        allRetPc = cat(1, allRetPc, retPc{iDb});
        
        allRetCoord = cat(1,allRetCoord, retCoord{iDb});
        
        [rotMapPC(:,:, done), rotMapIN(:,:, done),rotMapDiff(:,:, done)] = cameraLucida.plotInPcKlustRet(rotIn{iDb}', rotPc{iDb}',rotRet{iDb}', 0.5,0);
        %         cameraLucida.plotInPcKlustRet(retIn{iDb}, retPc{iDb},0.5,1)
        %             drawnow; pause;
        
    else
%           done = done + 1;
        coord= inCoord{iDb};
        coord(:,2) = -coord(:,2); %invert y coordinate as image microns coordinates are in ij
        
        [retInX{iDb}, retInY{iDb}, retColX{iDb}, retColY{iDb}, dxGrid, dyGrid,~,~,~, sRetX(iDb),sRetY(iDb), sX(iDb),sY(iDb)] ...
            = npRet.getRetDistFromCoords(db(iDb), coord,[], [], 0);

        coord= pcCoord{iDb};
        coord(:,2) = -coord(:,2);
        [retPcX{iDb}, retPcY{iDb}] ...
            = npRet.getRetDistFromCoords(db(iDb), coord, [], [], 0);
        
      
        [retColX{iDb}, retColY{iDb}, ~,~,xc{iDb}, yc{iDb}] = npRet.nullRetDist(db(iDb), 100000, 1, spatConst(iDb));
        ret{iDb} = [makeVec(retColX{iDb}), makeVec(-retColY{iDb})];
        retCoord{iDb} = [makeVec(xc{iDb}), makeVec(yc{iDb})];
        %
        retIn{iDb} = [retInX{iDb}, -retInY{iDb}]; % invert Y ret cos mpep maps upside down
        retPc{iDb} = [retPcX{iDb}, -retPcY{iDb}];
        

        allRetRet = cat(1, allRetRet, ret{iDb});
        allRetIn =  cat(1, allRetIn, retIn{iDb});
        allRetPc = cat(1, allRetPc, retPc{iDb});
        
        allRetCoord = cat(1,allRetCoord, retCoord{iDb});
                %        
    end
end

for iDb = 1:size(rotMapDiff,3)
    
    %     if ~isnan(db(iDb).prefDir) && ~isnan(spatConst(iDb))
    if ~isnan(db(iDb).prefDir)
        scaledRotMapDiff(:,:,iDb) = rotMapDiff(:,:,iDb)./max(abs(makeVec(rotMapDiff(:,:,iDb))));
    end
end
% figure; imagesc(mean(imgaussfilt(scaledRotMapDiff,0.2), 3)); colormap(BlueWhiteRed); axis image; caxis([-0.01 0.01])

% figure; imagesc(mean(imgaussfilt(rotMapDiff,0.2), 3)); colormap(BlueWhiteRed); axis image; caxis([-0.01 0.01])
% figure; imagesc(mean(imgaussfilt(-rotMapIN,1), 3)); colormap(BlueWhiteRed); axis image; caxis([-0.02 0.02])
% figure; imagesc(mean(imgaussfilt(rotMapPC,1), 3)); colormap(BlueWhiteRed); axis image; caxis([-0.01 0.01])
% figure; imagesc(mean(imgaussfilt(rotMapPC-rotMapIN,1), 3)); colormap(BlueWhiteRed); axis image; caxis([-0.01 0.01])


% retFig = fullfile(saveTo, 'Pool_retDist');
% rotFig = fullfile(saveTo, 'Pool_retRotDist');
% 
% cameraLucida.plotInPcKlustRet(allRotIn', allRotPc',allRotRet',0.5,1, [])
% cameraLucida.plotInPcKlustRet(allRetIn, allRetPc,allRetRet,0.5,1)

% % cameraLucida.angularAnisotropy(pcCoord, retPc);
% % cameraLucida.angularAnisotropy(inCoord, retIn);
% % 
% % cameraLucida.angularAnisotropy([inCoord; pcCoord], [retIn; retPc]);

%%

[nE, nA, nD] = size(rotMapIN);

[gridA, gridE] = meshgrid(1:nE, 1:nA);
% figure;
for iD = 1:nD

coord = cat(3, gridA,gridE);
data = imgaussfilt(double(rotMapIN(:,:,iD)),2);%+rotMapPC(:,:,iD)

options=optimset('Display','off');

%pars = [Amp, x0, y0, w]
pars0 = [0.1, 21, 21, 5];
parsUP = [inf, 30,  30, 15];
parsLW = [0, 10, 10, 1];

fitParsIN(:, iD) = lsqcurvefit(@Ret.gaussian2DRound, pars0, coord, data,parsLW ,parsUP, options);

end
%%


% figure;
center = 21;%9


long = fitParsIN(4,:);

stretch = mean(long);

for iD = 1:nD
    

thisStretch(iD) = stretch/long(iD);

cmA(iD) =fitParsIN(2,iD);
cmE(iD) = fitParsIN(3,iD);


alignedPC(:,:,iD) = interp2( gridA-(cmA(iD)-center), gridE-(cmE(iD)-center),rotMapPC(:,:,iD), gridA,gridE,'linear',0);
alignedIN(:,:,iD) = interp2( gridA-(cmA(iD)-center), gridE-(cmE(iD)-center),rotMapIN(:,:,iD), gridA,gridE,'linear',0);
alignedPC(:,:,iD) = interp2( (gridA-center)*thisStretch(iD) , (gridE-center)*thisStretch(iD) , alignedPC(:,:,iD), (gridA-center),(gridE-center),'linear',0);
alignedIN(:,:,iD) = interp2( (gridA-center)*thisStretch(iD) , (gridE-center)*thisStretch(iD) , alignedIN(:,:,iD), (gridA-center),(gridE-center),'linear',0);

alignedPC(:,:,iD) = interp2( gridA-(-cmA(iD)+center), gridE,alignedPC(:,:,iD), gridA,gridE,'linear',0);
alignedIN(:,:,iD) = interp2( gridA-(-cmA(iD)+center), gridE,alignedIN(:,:,iD), gridA,gridE,'linear',0);

alignedPC(:,:,iD) = alignedPC(:,:,iD)/max(makeVec(alignedPC(:,:,iD)));
alignedIN(:,:,iD) = alignedIN(:,:,iD)/max(makeVec(alignedIN(:,:,iD)));


end
 
angleSector = cart2pol(gridA-center,gridE-center);

aveAlignedIN = mean(imgaussfilt(alignedIN,0.75), 3);
centeredIN = -aveAlignedIN/sum(aveAlignedIN(:));

cmIN = sum(makeVec(-centeredIN.*(gridA-center)));

aveAlignedPC = mean(imgaussfilt(alignedPC,0.75), 3);
centeredPC = aveAlignedPC/sum(aveAlignedPC(:));
cmPC = sum(makeVec(centeredPC.*(gridA-center)));


%%
cameraLucida.swipeInPcRF(centeredPC, centeredIN, gridA, gridE,center, pi/2)

cameraLucida.offsetInPc(alignedPC, alignedIN, gridA, gridE, center, pi/2)

