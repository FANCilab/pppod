function [imgG, imgR, imgM, bkgG, bkgR, cellG, cellR] = ...
    cutRoiImg(ROI, frameG, frameR, PxSz, cropSize, scaleFlag)

if nargin < 5
    cropSize = 25;
end

if nargin <6
    scaleFlag = 1;
end

if isempty(PxSz)
    frameHW = round(cropSize);
else
    frameHW = round(cropSize/PxSz);
end

nROI = numel(ROI);

[y,x] = size(frameG);

newFrameG = ones(y + 2*frameHW, x + 2*frameHW)*mean(frameG(:));
newFrameG(frameHW+1:y+frameHW,frameHW+1:x+frameHW) =  frameG;

newFrameR = ones(y + 2*frameHW, x + 2*frameHW)*mean(frameR(:));
newFrameR(frameHW+1:y+frameHW,frameHW+1:x+frameHW) =  frameR;

% f=figure;
for ir = 1: nROI
    mask = zeros(y,x);
    
    pix = ROI{ir};
    
    mask(pix) = 1;
    
    newMask = zeros(y + 2*frameHW, x + 2*frameHW);
    newMask(frameHW+1:y+frameHW,frameHW+1:x+frameHW) = mask;
    
    [yc, xc] = find(newMask);
    
    yc = round(mean(yc));
    xc = round(mean(xc));
    
    % cropG = imadjust(mat2gray(newFrameG(yc-frameHW:yc+frameHW, xc-frameHW:xc+frameHW)));
    % cropR = imadjust(mat2gray(newFrameR(yc-frameHW:yc+frameHW, xc-frameHW:xc+frameHW)));
    if scaleFlag
        cropG = mat2gray(newFrameG(yc-frameHW:yc+frameHW, xc-frameHW:xc+frameHW));
        cropR = mat2gray(newFrameR(yc-frameHW:yc+frameHW, xc-frameHW:xc+frameHW));
    else
        cropG = newFrameG(yc-frameHW:yc+frameHW, xc-frameHW:xc+frameHW);
        cropR = newFrameR(yc-frameHW:yc+frameHW, xc-frameHW:xc+frameHW);
    end
    
    imgG(:,:, ir) = cropG;
    imgR(:,:, ir) = cropR;
    imgM(:,:, ir) = newMask(yc-frameHW:yc+frameHW, xc-frameHW:xc+frameHW);
    
    bkgG(ir) = nanmedian(makeVec(cropG(logical(~imgM(:,:, ir)))));
    bkgR(ir) = nanmedian(makeVec(cropR(logical(~imgM(:,:, ir)))));
    
    cellG(ir) =  mean(cropG(logical(imgM(:,:, ir))));
    cellR(ir) =  mean(cropR(logical(imgM(:,:, ir))));
    
%     figure(f)
%     clf;
%     subplot(1,2,1)
%     imagesc(imgG(:,:,ir)); axis image; colormap gray; hold on
%     contour(imgM,'r', 'Linewidth', 0.5)
%     
%     subplot(1,2,2)
%     imagesc(imgR(:,:,ir)); axis image; colormap gray;
%     pause;
    
end

end