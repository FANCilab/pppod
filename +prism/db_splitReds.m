function db_splitReds(db, regrOps)

root = 'D:\OneDrive - University College London\Data\2P';

[frameG, frameR, micronsX, micronsY, micronsZ] = s2pUtils.loadMultiPlaneMeanStack(db);


folder = buildExpPath(db);
load(fullfile(root, folder, 'mimgSingleExp'));

for iExp = 1:numel(db.expts)
    
    info(iExp) = ppbox.infoPopulateTempLFR(db.mouse_name, db.date, db.expts(iExp));
    laserPw(:,iExp) = ppbox.getLaserPower(info(iExp));
    
end
%     m = matfile(fullfile(root, folder, 'mimgSingleExp'),'Writable',true);
waveL = db.waveL ;
prismR = mimgSingleExp.R;
prismG = mimgSingleExp.G;

for iPlane = 1:info(1).nPlanes
    
    zPlane(iPlane) = mean(micronsZ{1}(:, iPlane));
    
    redCh=prismR(:,:,iPlane, :);
    greenCh=prismG(:,:,iPlane, :);
    [dsRedImg(:,:,iPlane), mCherryImg(:,:,iPlane), curedRed{iPlane}, mixing{iPlane}, iterErr(:,iPlane)] =...
        prism.splitColors(redCh, greenCh, waveL, laserPw(iPlane,:), zPlane(iPlane),  regrOps);
    
end
    
    
flatRedImg = dsRedImg;
flatCherryImg = mCherryImg;


for iPlane = 1:numel(mixing)
    totFluo(iPlane) = mean(sum(mixing{iPlane},2));
    
    flatRedImg(:,:,iPlane) = zstack.removeBackground(dsRedImg(:,:,iPlane), 1, [10 10 1]);
    
    flatCherryImg(:,:,iPlane) = zstack.removeBackground(mCherryImg(:,:,iPlane), 1, [10 10 1]);
    
end

if numel(unique(waveL))>1
fitFluoVsZ = fit(zPlane', totFluo', 'exp1');
fitFluoVsZ = fitFluoVsZ(zPlane);
nflatRedImg = bsxfun(@rdivide, flatRedImg, reshape(fitFluoVsZ, 1,1,[]))*mean(fitFluoVsZ);
nflatCherryImg = bsxfun(@rdivide, flatCherryImg, reshape(fitFluoVsZ, 1,1,[]))*mean(fitFluoVsZ);

else
 nflatRedImg   =  flatRedImg;
 nflatCherryImg = flatCherryImg;
 fitFluoVsZ = [];
end

save(fullfile(root, folder, 'splitReds'), 'frameG', 'frameR','micronsX', 'micronsY', 'micronsZ', 'db',...
    'waveL', 'laserPw', 'prismR', 'prismG', 'info','zPlane', 'nflatRedImg', 'nflatCherryImg', 'fitFluoVsZ',...
    'dsRedImg', 'mCherryImg', 'regrOps', 'curedRed', 'mixing', 'iterErr','-v7.3');

end