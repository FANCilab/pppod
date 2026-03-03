  function fixMissedFrame( iPlane, nPlanes, rep, nMissed)
% iPlane is the plane that is missing 1 frame
% rep is the # of rep before dropped frame rep
% assumes 2 channels were recorded
% nPlanes total number of planes
% nMissed number of missed frames

nCh = 2;
addpath('C:\Users\Federico\Google Drive\CarandiniLab\CarandiniLab_MATLAB\FedericoBox\2P\Tools');

[fn, root] = uigetfile('\\zserver.cortexlab.net\Data\Subjects\*');

fixFn = [fn(1:end-4), '_fix.tif'];
    
[mov, ~] = img.loadFrames(fullfile(root, fn));

[ny, nx, nt] = size(mov);

fixMov = zeros(ny, nx, nt + nCh*nMissed, 'like', mov);

lastR = nPlanes*nCh*(rep-1) +iPlane *nCh;
lastG = lastR-1;

missR = nPlanes*nCh*rep+iPlane *nCh;
missG = missR - 1;

fixMov(:,:, 1:missG - 1) = mov(:,:, 1:missG - 1);

fixMov(:,:, missG + nCh*nMissed:end) = mov(:,:, missG:end);

for iMiss = 0: nMissed-1
fixMov(:,:,missG + iMiss*nCh) = (fixMov(:,:, lastG+iMiss*nCh) + fixMov(:,:, lastG+nCh*nPlanes*2+ iMiss*nCh))/2;
fixMov(:,:,missR+ iMiss*nCh) = (fixMov(:,:, lastR+ iMiss*nCh) + fixMov(:,:, lastR+nCh*nPlanes*2+ iMiss*nCh))/2;
end

img.saveFrames(fixMov, fullfile(root, fixFn))


end