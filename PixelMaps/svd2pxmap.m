function [map, rr] = svd2pxmap(db, targetPlane, stim_type, rmv_svd1)
%% load a svd compresesed 2P recordinga adn computes stim triggered pixelMaps
%INPUTS
% 

%OUTPUTS
% - map: contain tuning pixel maps. Struct with fields:
%       - dir: pixel map of preferred direction
%       - ori: pixel map of preferred orientation
%       - mimg: average image
%       - p: protocol file 

if nargin <4

   rmv_svd1 = true;

end

%% load SVD compressed recording

root_folder = db.root_storage;
info = getExpInfo(db.subject, db.date, db.exp(db.expID));

switch db.s2p_version
    case 'python'
targetPlane = targetPlane-1;
s2p_folder = fullfile(info.folder2praw,sprintf('plane%d', targetPlane));

s2p_file = sprintf('%s/Fall.mat', s2p_folder);

svd_file = sprintf('%s/SVD_%s_%s_plane%d.mat', s2p_folder, ...
    db.subject, db.date, targetPlane);

    case 'matlab'

s2p_folder = fullfile(root_folder, db.mouse_name, db.date, ...
    sprintf('%d_', db.expts));
s2p_folder = s2p_folder(1:end-1);

svd_file = sprintf('%s/SVD_%s_%s_plane%d.mat', s2p_folder, ...
    db.mouse_name, db.date, targetPlane);

end


svd = load(svd_file);
[nY, nX, nBasis] = size(svd.U);
nFrames = size(svd.Vcell{db.expID},2);

S = svd.U;
S = reshape(S, [], nBasis); % nPx * nB
T = svd.Vcell{db.expID} ; % nB * nT

map.all_st.mimg = reshape(S*mean(T,2), nY,nX)';

if rmv_svd1
    S = S(:, 2:end);
    T = T(2:end,:);
    nSVD = svd.ops.nSVD-1; % minimum hard coded, perhaps make input
else
    nSVD = svd.ops.nSVD; % minimum hard coded, perhaps make input
end

%% load stimulus info and compute STA on temporal components


switch stim_type
    case 'gratings' % sta across stimuli of same direction, averaged across sf and tf if multiple present


        % load stimulus info
        switch db.s2p_version
            case 'python'
        nFrames = svd.ops.frames_per_folder(db.expID);
            case 'matlab'
                nFrames = svd.ops.Nframes(db.expID);

        end
        planeFrames = targetPlane:db.n_planes:(nFrames*info.nPlanes); % check if it works for multiplane recs

        %%

        event = bonsai.load_events(db);

        frameTimes = event.frame.on;

        frameRate = (1/mean(diff(frameTimes)));

        stimTimes = event.grating.stimTimes;

        stimSequence = event.grating.stimSequence;

        stimMatrix = event.grating.stimMatrix(event.grating.stimIdx, targetplane:info.nPlanes:end);
        stimMatrix_blank = event.grating.stimMatrix(event.grating.blankStimIdx, targetplane:info.nPlanes:end);


        p.nDir = sum(event.grating.is_stim);
        p.dirs = 0:30:330;

      
        dirs = p.dirs; %dirs = -dirs;
        dirs = deg2rad(dirs);
        oris = dirs;
        % oris = dirs - pi/2;
        oris(oris >= pi) = oris(oris >=pi) -pi;
        % oris(oris < 0) = oris(oris <0) +pi;
        oris = oris -pi/2;
        oris = oris*2;

        %% now compute average across SfTf

        % compute stimulus triggered responses
       
        % filter in time; HARDCODED change for different indicators
        sT = gaussFilt(T(1:nSVD, :)',2);
        
        % measure sta response timecourse 
        [resp, aveResp, ~, kernelTime] = ...
            ppbox.getStimulusSweepsLFR(sT, stimTimes, stimMatrix,frameRate); % responses is (nroi, nStim, nResp, nT)
        aveResp = aveResp(:,:,kernelTime>-1 & kernelTime<3);

         [resp_blank, aveResp_blank] = ...
            ppbox.getStimulusSweepsLFR(sT, stimTimes, stimMatrix_blank,frameRate); % responses is (nroi, nStim, nResp, nT)
        aveResp_blank = aveResp_blank(:,:,kernelTime>-1 & kernelTime<3);
       
        % measure response in respWindow
        respWin = [0, 2];
        [resPeak, aveResPeak] = ...
            ppbox.gratingOnResp(resp, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

         [resPeak_blank, aveResPeak_blank] = ...
            ppbox.gratingOnResp(resp_blank, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

         % gentle spatial smoothing of the spatial components
        sS =reshape( S(:, 1:nSVD), nY, nX, nSVD);
        sS = imgaussfilt(sS, 0.5);
        sS = reshape(sS, nY*nX, nSVD);

        % save US and V to reconstruct resp for each px, px =
        % sS*sT;
        nStim = size(resPeak,2); nRep = size(resPeak,3);
        resPeak = reshape(resPeak, nSVD, nStim*nRep);
        
        rr.trial_stim_resp = single(sS*resPeak);
        rr.stim_resp = single(sS*aveResPeak);
        rr.stimSequence = stimSequence;


        nStim_blank = size(resPeak_blank,2); nRep_blank = size(resPeak_blank,3);
        resPeak_blank = reshape(resPeak_blank, nSVD, nStim_blank*nRep_blank);
        rr.trial_blank_resp = single(sS*resPeak_blank);
        rr.blank_resp = single(sS*aveResPeak_blank);

        %pixel maxps for direction
        svdTun_dir = aveResPeak.*exp(1i*dirs);
        svdTun_dir = mean(svdTun_dir, 2);
        pxTun_dir = sS*svdTun_dir(1:nSVD);
        map.all_st.dir = permute(reshape(pxTun_dir, nY, nX), [2,1]); % for some reason the SVDs are transposed

        %pixel maxps for orientation
        svdTun_ori = aveResPeak.*exp(1i*oris);
        svdTun_ori = mean(svdTun_ori, 2);
        pxTun_ori = sS*svdTun_ori(1:nSVD);
        map.all_st.ori = permute(reshape(pxTun_ori, nY, nX), [2,1]); % for some reason the SVDs are transposed

        map.all_st.mimg = mat2gray(map.all_st.mimg); % for some reason the SVDs are transposed

        plot_px_map(map.all_st, stim_type, 1, [10 90],0);
        print(fullfile(s2p_folder, sprintf('plane_%d_%s_all_st_px_maps', targetPlane, stim_type)), '-dpdf', '-vector', '-bestfit')

    case 'sparsenoise'

        % need to integrate old code from function 'RF_mapping_mov'
end

map.p = p;
rr.p = p;
rr.nX = nX;
rr.nY = nY;
rr.nRep = nRep;
rr.nStim = nStim;
rr.nRep_blank = nRep_blank;
rr.oris = oris;
rr.dirs = dirs;

end




