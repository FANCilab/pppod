function mov = planes2mov(db)

info = ppbox.infoPopulateTempLFR(db.mouse_name, db.date, db.expts(db.expID));
info.fovX = info.scanPixelsPerLine;
info.fovY = info.scanLinesPerFrame;

allFrameTimes = ppbox.getFrameTimes(info, 'all');

numberOfFrames = 600;


mov = nan(512*3, 512*3,numberOfFrames);

row = [ 1 1 1 2 2 2 3 3 3];
col = [1 2 3 1 2 3 1 2 3];


for iPlane = 2:info.nPlanes % discard flyback
    
    nmov = iPlane -1;

    slot_x{nmov} = 1+(col(nmov)-1)*512:col(nmov)*512;
    slot_y{nmov} = 1+(row(nmov)-1)*512:row(nmov)*512;

    [root, refF, refNeu, refSVD] = starter.getAnalysisRefs(db.mouse_name, db.date, db.expts, iPlane);
    
    svd = load(fullfile(root, refSVD));
    
    load(fullfile(root, [refF(1:end-9), '.mat']), 'ops');
    
    [nY, nX, nBasis] = size(svd.U);
    
    nFrames = size(svd.Vcell{db.expID},2);
               
%     S = imgaussfilt(svd.U, 3*sigma); % smooth spatial basis with gaussian 30 um sigma
    S = svd.U; 
    S = reshape(S, [], nBasis); % nPx * nB
    T = svd.Vcell{db.expID} ; % nB * nT
%     planeFrames = iPlane:info.nPlanes:(size(T, 2)*info.nPlanes);
%     allFrameTimes = allFrameTimes(planeFrames);

    one_mov = reshape(S*T(:, 1:numberOfFrames),nY, nX,numberOfFrames);
    one_mov = mat2gray(imgaussfilt3(one_mov, [0.05 0.05 1]));
%         implay(reshape(mov, nY, nX,300))
%                     Fneu = S(iPx, :)*svd.Vcell{db.expID};

    mov(slot_y{nmov}(ops.yrange), slot_x{nmov}(ops.xrange), :) = mat2gray(one_mov, [0.1, 0.8]);


end

mov(isnan(mov)) = 1;

% implay(mov)

%%
new_file = 'planes_montage';
vw = VideoWriter(fullfile(root, new_file),'MPEG-4');
open(vw); % open for writing to disk

nX = 1024;
nY = 1024;


pan_x = [513:1024];
pan_y = [513:1024];

dPan = round(512/(numberOfFrames/3));
nPan = floor(512/dPan);

step = 0;
edit main
for iF = 1:numberOfFrames % for every frame
    
    
    
    if iF>numberOfFrames/3
        
        step=step+1;
        
      
        
    end
    
    
      start_x = 513+round(513*0.2);
        start_y = 513+round(513*0.2);
        end_x = 1024-round(513*0.2);
        end_y = 1024-round(513*0.2);

        pan_x = [start_x-dPan*step:end_x+dPan*step];
        pan_y =  [start_y-dPan*step:end_y+dPan*step];

        if pan_x(1)<=0
            pan_x = [1:512*3];
            pan_y = [1:512*3];
        end
    
    
        thisFrame = mov(pan_x,pan_y, iF); 
        
        nPx = 750;
  xvec = linspace(1, nPx, size(thisFrame, 1));
    yvec = linspace(1,nPx, size(thisFrame, 2))';
 qx = 1:nPx;
 qy = makeVec(1:nPx);
    newFrame = mat2gray(interp2(xvec, yvec, thisFrame, qy, qx), [0.05 0.9]);
    writeVideo(vw,   newFrame); % write the frame to the mj2 file
    
end

	close(vw); % close videowriter object

%%

end