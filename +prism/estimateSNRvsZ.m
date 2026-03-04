function [snrZ, nPhotons, dsRed_exc]  = estimateSNRvsZ(redCh, greenCh, waveL,laserPw, doPlot)

if nargin <5
    
   doPlot = false;
end

nPhotons = bsxfun(@times, laserPw/100, waveL.*prism.lookUpPower(waveL)); % nphotons ~ pw*waveL
nPhotons = nPhotons.^2;

[dsRed_exc] = prism.excitationSpectrum('dsRed2', waveL);

% nPhotons = nPhotons.*dsRed_exc;

nPhotons = nPhotons./mean(nPhotons(:));

redCh = bsxfun(@rdivide, redCh, reshape(nPhotons, [1 1 size(nPhotons)]));
greenCh = bsxfun(@rdivide, greenCh, reshape(nPhotons, [1 1 size(nPhotons)]));

% % nPhotons = ones(size(nPhotons));
% 
% [mCherry_exc] = prism.excitationSpectrum('mCherry', waveL);
% 
% %bscope red channel 607/70
% [mCherry_ems] = prism.emissionSpectrum('mCherry', 572:642);
% [dsRed_ems] = prism.emissionSpectrum('dsRed2', 572:642);
% 
% [G6Ca_ems, pmt] = prism.emissionSpectrum('dsRed2', 572:642);
% [G6free_ems] = prism.emissionSpectrum('dsRed2', 572:642);
% 
% mCherry = mCherry_exc.*nansum(mCherry_ems.*pmt);
% dsRed = dsRed_exc.*nansum(dsRed_ems.*pmt);
% 
% % ch1 and ch2 are nPxX X nPxY X nPlanes X nWaves
% nCh = 2; 

[nX, nY, nPlanes, nWave] = size(redCh);

refWave = find(waveL == 890);
otherWave = find(waveL ~= 890);

waveOrder = [refWave, otherWave];

for iPlane = 1:nPlanes
    refFound = 1;
    for iW = waveOrder
        rr= imgaussfilt(redCh(:,:,iPlane,iW), 0.1);
        gg= imgaussfilt(greenCh(:,:,iPlane,iW), 0.1);
        [img, ~, ~,~, globalThresh(iPlane, iW)] = prism.bleedCure(rr, gg,1);
        %         thStack = zstack.removeBackground(curedR{iPlane}(:,:, iW), 0.1, [8 8 1]);
        if refFound == 1
            [bwBloob, redList{iPlane, iW}, cc] = zstack.forgrBlobs(img, globalThresh(iPlane, iW), [10 10 1], 1);
        end
        signal = img(logical(bwBloob));
        snrZ(iPlane, iW) = mean(signal);
        refFound = 0;
        
    end
    
%     [~, idx] = sort(waveOrder, 'ascend');
%     snrZ = snrZ(:, sort(idx));
%     snrZ = snrZ./(dsRed_exc/mean(dsRed_exc));
    
end
end
