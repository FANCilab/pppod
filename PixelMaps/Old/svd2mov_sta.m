function [sta_mov, map] = svd2mov_sta(db, targetPlane, stim_type)

%% load SVD compressed recording

root_folder = db.root_folder;

s2p_folder = fullfile(root_folder, db.mouse_name, db.date, ...
    sprintf('%d', db.expts),'suite2P', sprintf('plane%d', targetPlane-1));

svd_file = sprintf('%s/SVD_%s_%s_plane%d.mat', s2p_folder, ...
    db.mouse_name, db.date, targetPlane-1);

svd = load(svd_file);

[nY, nX, nBasis] = size(svd.U);

nFrames = size(svd.Vcell{db.expID},2);

%     S = imgaussfilt(svd.U, 3*sigma); % smooth spatial basis with gaussian 30 um sigma
S = svd.U;
S = reshape(S, [], nBasis); % nPx * nB
T = svd.Vcell{db.expID} ; % nB * nT



%% load stimulus info and compute STA on temporal components

info = ppbox.infoPopulateTempLFR(db.mouse_name, db.date, db.expts(db.expID));

switch stim_type
    case 'gratings' % sta across stimuli of same direction, averaged across sf and tf if multiple present


        % load stimulus info
        nFrames = svd.ops.frames_per_folder(db.expID);
        planeFrames = targetPlane:info.nPlanes:(nFrames*info.nPlanes); % check if it works for multiplane recs
        %         allFrameTimes = allFrameTimes(planeFrames);

        frameTimes = ppbox.getFrameTimes(info, planeFrames);
        frameRate = (1/mean(diff(frameTimes)));

        stimTimes = ppbox.getStimTimes(info);
        [stimSequence, p] = ppbox.getStimSequence_LFR(info);

        nSfTf = floor(p.nstim/p.nDir);
        stimSequence.seq = ceil(stimSequence.seq/nSfTf);
        stimSequence.labels = stimSequence.labels(1:nSfTf:end);

        stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, frameTimes);

        % remove blank from stimMatrix
        stimMatrix(p.nDir+1,:) = [];

        % compute stimulus triggered responses
        nSVD = min(1000, svd.ops.nSVD); % hard coded, perhaps make input

        sT = gaussFilt(T(1:nSVD, :)',2);
        [resp, aveResp, ~, kernelTime] = ...
            ppbox.getStimulusSweepsLFR(sT, stimTimes, stimMatrix,frameRate); % responses is (nroi, nStim, nResp, nT)
        % measure response in respWindow
        respWin = [0, 2];
        [~, aveResPeak] = ...
            ppbox.gratingOnResp(resp, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

        aveResp = aveResp(:,:,kernelTime>-1 & kernelTime<3);

  
        %%
        dirs = p.dirs; %dirs = -dirs;
        dirs = deg2rad(dirs);

        oris = dirs;
        oris(oris >= pi) = oris(oris >=pi) -pi;
        oris = oris*2 + pi/2;

        svdTun_dir = aveResPeak.*exp(1i*dirs);
        svdTun_dir = mean(svdTun_dir, 2);
        pxTun_dir = S(:, 1:nSVD)*svdTun_dir(1:nSVD);
        map.dir = permute(reshape(pxTun_dir, nY, nX), [2,1]); % for some reason the SVDs are transposed

        svdTun_ori = aveResPeak.*exp(1i*oris);
        svdTun_ori = mean(svdTun_ori, 2);
        pxTun_ori = S(:, 1:nSVD)*svdTun_ori(1:nSVD);
        map.ori = permute(reshape(pxTun_ori, nY, nX), [2,1]); % for some reason the SVDs are transposed


        map.mimg = mat2gray(reshape(S(:,1), nY, nX))'; % for some reason the SVDs are transposed

        plot_px_map(map, stim_type);
        print(fullfile(s2p_folder, sprintf('%s_px_maps', stim_type)), '-dpdf', '-vector', '-bestfit')

    case 'sparsenoise'

        % need to integrate old code from function 'RF_mapping_mov'
end


%% now generate and save the sta movie

newT = reshape(permute(aveResp, [3,2,1]), [], nSVD)'; %nSVD*nT

sta_mov =  S(:, 1:nSVD)*newT(1:nSVD, :);
% sta_mov = zscore(sta_mov')';
sta_mov = reshape(sta_mov, nY, nX, []);
% sta_mov = bsxfun(@minus, sta_mov, min(sta_mov, [], 2));
% sta_mov = bsxfun(@rdivide, sta_mov, max(sta_mov, [], 2));

sta_mov = permute(imgaussfilt3(sta_mov, [0.5, 0.5, 1]), [2 1 3]); % for some reason the SVDs are transposed

saveastiff(uint16(mat2gray(sta_mov)*(2^15-1)), fullfile(s2p_folder, sprintf('%s_sta_mov.tif', stim_type)));

end




