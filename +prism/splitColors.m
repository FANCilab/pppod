function [dsRedImg, mCherryImg, curedR, mixing, iterErr, regrOps] =splitColors(redCh, greenCh, waveL, laserPw, z, regrOps)

if nargin < 6
    regrOps.rectifySource = 0;
    regrOps.forceRatio = 0;
    regrOps.zLimit = Inf;
    regrOps.useAllPx = 1;
end



addpath(genpath('C:\Users\Federico\Documents\MATLAB\FastICA_2.5\FastICA_25'));

if z>regrOps.zLimit
    laserSNR = prism.laserPenetration(waveL, z);
    nPhotons = bsxfun(@times, laserPw/100, waveL.*prism.lookUpPower(waveL)); % nphotons ~ pw*waveL
    nPhotons = nPhotons.^2;
    % % nPhotons = ones(size(laserSNR));
    laserSNR = laserSNR.*nPhotons/mean(laserSNR.*nPhotons);
    %         laserSNR = nPhotons/mean(nPhotons);
    
else
    laserSNR = ones(size(waveL));
    
end

[fpR, fpG] = prism.mixFP({'dsRed2', 'mCherry'}, waveL);

% ch1 and ch2 are nPxX X nPxY X nPlanes X nWaves

nCh = 2;

[nX, nY, nPlanes, nWave] = size(redCh);

for iPlane = 1:nPlanes
    for iW = 1:nWave
        
        rr= imgaussfilt(redCh(:,:,iPlane,iW), .1);
        gg= imgaussfilt(greenCh(:,:,iPlane,iW), .1);
        [curedR{iPlane}(:,:, iW), ~,redList(iPlane, iW),~,globalThresh(iPlane, iW)] = prism.bleedCureNL(rr, gg, [], [], [], 0);
        if z >regrOps.zLimit
            curedR{iPlane}(:,:, iW) = zstack.removeBackground(curedR{iPlane}(:,:, iW), 0.5, [10 10 1]);
        end
        %         [bwBlob{iPlane, iW}, redList{iPlane, iW}] = zstack.forgrBlobs_dev(curedR{iPlane}(:,:, iW), [10 10 1]);
    end
    
    if numel(unique(waveL))>1
        
        redPx{iPlane} = unique(cat(1, redList{:}));
        
        redPxAcross{iPlane} = [];
        
        for iC = 1:nWave
            redPxAcross{iPlane} = cat(1,redPxAcross{iPlane}, redPx{iPlane} + (iC-1)*nX*nY);
        end
        
        colorPx{iPlane} = reshape(curedR{iPlane}, nX*nY, nWave);
        colorPx{iPlane} = bsxfun(@rdivide, colorPx{iPlane}, laserSNR);
        
        if regrOps.useAllPx
            nonZeroPx = find(sum(colorPx{iPlane}>0, 2));
        else
            nonZeroPx = redPx{iPlane};
        end
        if z<regrOps.zLimit
            [~, mixing, iterErr] = prism.learnSources_iterative(colorPx{iPlane}(nonZeroPx, :)', fpR,regrOps.forceRatio, 1);
            sources = pinv(mixing)*colorPx{iPlane}';
            
        else
            %     sources = pinv(fpR)*nbcolorPx{iPlane}';
            [~, mixing, iterErr] = prism.learnSources_iterative(colorPx{iPlane}(nonZeroPx, :)', fpR,regrOps.forceRatio, 1);
            sources = pinv(mixing)*colorPx{iPlane}';
        end
        
        if regrOps.rectifySource
            
            nonZeroCherry = sources(2,:)>0;
            [lowEnvFit, fitY, lowY] = lowEnvelopeReg(sources(2,nonZeroCherry),sources(1,nonZeroCherry),100, 50, 99, 50);
            theta = atan(lowEnvFit(2));
            figure;
            plot(sources(2,nonZeroCherry),sources(1,nonZeroCherry), '.'); hold on
            plot(fitY, lowY, 'r')
            
            %         if abs(theta) > pi/4
            %             source2Px = sources(1,nonZeroPx)< mean(sources(1,nonZeroPx));
            %             [lowEnvFit, fitY, lowY] = lowEnvelopeReg(sources(2,nonZeroPx(source2Px)),sources(1,nonZeroPx(source2Px)));
            %             theta = atan(lowEnvFit(2));
            %         end
            
            rotMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            sources = rotMatrix*sources;
            
        end
        
        dsRedImg = reshape(sources(1,:), nX, nY);
        mCherryImg= reshape(sources(2,:), nX, nY);
   
    %%
    
    figure;
    
    rm =  cat(1, [1 1 1], cbrewer('seq', 'Reds',100));
    bm =  cat(1, [1 1 1], cbrewer('seq', 'Blues',100));
    gm =  cat(1, [1 1 1], cbrewer('seq', 'Greens',100));
    pm = cat(1, [1 1 1], cbrewer('seq', 'RdPu',100));
    
    r3 = subplot(2,3,1);
    
    img1 = curedR{iPlane}(:,:, find(waveL >=890 & waveL <=1000, 1, 'first'));
    imagesc(img1); axis image; colorbar
    caxis(r3,[0 3*globalThresh(iPlane, find(waveL >=890 & waveL <=1000, 1, 'first'))])
    title(sprintf('%d nm', waveL(find(waveL >=890 & waveL <=1000, 1, 'first'))))
    
    r4 = subplot(2,3,2);
    try
        img2   = curedR{iPlane}(:,:, find(waveL <=780, 1, 'first'));
        imagesc(img2); axis image; colorbar;
        caxis(r4,[0 3*globalThresh(iPlane, find(waveL ==1020, 1, 'first'))])
        title(sprintf('%d nm', waveL(find(waveL<=780, 1, 'first'))))
        
    catch
        img2 = curedR{iPlane}(:,:, find(waveL >=1000, 1, 'first'));
        imagesc(img2); axis image; colorbar;
        caxis(r4,[0 3*globalThresh(iPlane, find(waveL >=1000, 1, 'first'))])
        title(sprintf('%04d nm', waveL(find(waveL >=1000, 1, 'first'))))
        
    end
    
    %     subplot(2,3,3)
    % %     plot(img1,img2, '.b'); axis square
    % %
    % %
    lims1 = [0, max(dsRedImg(:))*0.1];
    lims2 = [0.05, max(mCherryImg(:))*0.2];
    %
    
    r1 = subplot(2,3,4);
    imagesc(dsRedImg); axis image; caxis(lims1);colorbar;
    colormap(r1, rm);
    formatAxes
    title('dsRed')
    
    r2 = subplot(2,3,5);
    imagesc(imgaussfilt(mCherryImg)); axis image;caxis(lims2); colorbar;
    colormap(r2, pm);
    linkaxes([r1 r2 r3 r4], 'xy')
    formatAxes
    title('mCherry')
    
    
    subplot(2,3,6)
    hold on
    lims = [min(min(dsRedImg(:)),min(mCherryImg(:))), max(max(dsRedImg(:)), max(mCherryImg(:)))];
    
    [xy, xbins, ybins] = cameraLucida.plot_Density2D (dsRedImg, mCherryImg, 0.1, 1, [lims, lims], 0, 0);
    xy = xy/sum(xy(:));
    zax=imagesc(xbins, ybins, imgaussfilt(xy,1)); axis image; hold on;
    caxis([0 0.00001])
    %     colormap(1-gray);
    colormap(cat(1, [1 1 1], RedWhite(100),flip(Reds(100))));
    %     nanMask = xy;
    %     nanMask(nanMask == 0) = NaN;
    %     set(zax, 'alphadata', ~isnan(nanMask))
    %     xlim([-0.5 5]); ylim([-0.5 5])
    xlabel('dsRed')
    ylabel('mCherry')
    
    xlim(lims)
    ylim(lims)
    formatAxes
    
     else
         dsRedImg = mean(curedR{iPlane}, 3);
         mCherryImg = dsRedImg;
         mixing = mean(dsRedImg(:));
         iterErr = 0;
         
         rm =  cat(1, [1 1 1], cbrewer('seq', 'Reds',100));
    bm =  cat(1, [1 1 1], cbrewer('seq', 'Blues',100));
    gm =  cat(1, [1 1 1], cbrewer('seq', 'Greens',100));
    pm = cat(1, [1 1 1], cbrewer('seq', 'RdPu',100));
%         lims1 = [0, max(dsRedImg(:))*0.1];
    lims2 = [0.05, max(mCherryImg(:))*0.2];
    lims1 = lims2;
         figure;
         
         r1 = subplot(1,2,1);
    imagesc(dsRedImg); axis image; caxis(lims1);colorbar;
    colormap(r1, rm);
    formatAxes
    title('dsRed')
    
    r2 = subplot(1,2,2);
    imagesc(imgaussfilt(mCherryImg)); axis image;caxis(lims2); colorbar;
    colormap(r2, pm);
    linkaxes([r1 r2], 'xy')
    formatAxes
    title('mCherry')
    end
    
    
end
end