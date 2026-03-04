function twoChanICA(ch1, ch2)

% ch1 and ch2 are nPxX X nPxY X nPlanes X nWaves
nCh = 2; 

[nX, nY, nPlanes, nWave] = size(ch1);

for iW = 1:nWave
    
    [curedR, bleedThrough, redList(:, iW), redAll{iW}] = prism.bleedCure(ch2(:,:,:,iW), ch1(:,:,:,iW));
    
end

redPx = unique(cat(1, redAll{:}));

redPxAcross = [];

for iC = 1:nWave*nCh
    
    redPxAcross = cat(1,redPxAcross, redPx + (iC-1)*nX* nY* nPlanes);
    
end

colors = cat(4, ch1, ch2);

colorPx = reshape(colors(redPxAcross), numel(redPx), nWave*nCh);

[A,W] = fastica(colorPx'); 

colors = reshape(colors, nX*nY*nPlanes, nWave);
icaImages = transpose(W*colors');

icaImages = reshape(icaImages,[ nX,nY,nPlanes,nWave]);

colors = reshape(colors,[ nX,nY,nPlanes,nWave]);

figure;
for iPlane = 1:nPlanes
    
    for iM = 1:8
    subplot(2, 4,iM)
    imagesc(icaImages(:,:,iPlane,iM));axis image
    end
    
end

figure;
for iPlane = 1:nPlanes
    
    for iM = 1:8
    subplot(2, 4,iM)
    imagesc(colors(:,:,iPlane,iM));axis image
    end
    
end


end