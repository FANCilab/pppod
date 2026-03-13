function [sta_mov, map] = svd2mov_sta_dev(db, targetPlane, stim_type)
%% load a svd compresesed 2P recordinga adn computes stim triggered movies and pixelMaps
%OUTPUTS
% - sta_mov: stimulus triggered movie, concat across stimuli
% - map: contain tuning pixel maps. Struct with fields:
%       - dir: pixel map of preferred direction
%       - dir_Sf_Tf: pixel map of preferred combo of SfTf for orientation
%       - ori: pixel map of preferred orientation
%       - ori_Sf_Tf: pixel map of preferred combo of SfTf for orientation
%       - mimg: average image
%       - st: nSfTf-dim struct with fields dir, ori, mimg (like above) for each sftf
%       combo
%       - all_st: struct with fields dir, ori, mimg computed by averaging
%       across sftf combos
%       - p: protocol file 


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

        dirs = p.dirs; %dirs = -dirs;
        dirs = deg2rad(dirs);
        oris = dirs;
        oris(oris >= pi) = oris(oris >=pi) -pi;
        oris = oris*2 + pi/2;

        %% compute pixel map for each SfTf combo
        if nSfTf >1

            stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, frameTimes);

            for iST = 1:nSfTf
                %select stims
                stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, frameTimes);
                % remove blank from stimMatrix
                stimMatrix(p.blankstim,:) = [];
                % select this SfTf - hardcoded, improve by using protocol info
                stimMatrix = stimMatrix(iST:nSfTf:end, :);
                % sta of temporal components
                nSVD = min(1000, svd.ops.nSVD); % hard coded, perhaps make input
                % filter in time; HARDCODED change for different indicators
                sT = gaussFilt(T(1:nSVD, :)',2);
                [resp, ~, ~, kernelTime] = ppbox.getStimulusSweepsLFR(sT, stimTimes, stimMatrix,frameRate); % responses is (nroi, nStim, nResp, nT)

                % measure response in respWindow
                respWin = [0, 2];
                [~, aveResPeak] = ...
                    ppbox.gratingOnResp(resp, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

                % gentle spatial smoothing of the spatial components
                sS =reshape( S(:, 1:nSVD), nY, nX, nSVD);
                sS = imgaussfilt(sS, 0.5);
                sS = reshape(sS, nY*nX, nSVD);
                
                %pixel maxps for direction
                svdTun_dir = aveResPeak.*exp(1i*dirs);
                svdTun_dir = mean(svdTun_dir, 2);
                pxTun_dir = sS*svdTun_dir(1:nSVD);
                map.st(iST).dir = permute(reshape(pxTun_dir, nY, nX), [2,1]); % for some reason the SVDs are transposed

                %pixel maxps for orientation
                svdTun_ori = aveResPeak.*exp(1i*oris);
                svdTun_ori = mean(svdTun_ori, 2);
                pxTun_ori = sS*svdTun_ori(1:nSVD);
                map.st(iST).ori = permute(reshape(pxTun_ori, nY, nX), [2,1]); % for some reason the SVDs are transposed

                %average image
                map.st(iST).mimg = 1-mat2gray(reshape(sS(:,1), nY, nX))'; % for some reason the SVDs are transposed

                % plot
                plot_px_map(map.st(iST), stim_type, 1, [2 4]);
                print(fullfile(s2p_folder, sprintf('%s_st%d_px_maps', stim_type, iST)), '-dpdf', '-vector', '-bestfit')

            end
            
            % for each pixel, pick the sftf combo gave the strongest resp
            amp = abs(cat(3, map.st(:).ori));
            [~, map.ori_SfTf] = max(amp, [], 3, 'linear');
            amp = reshape(amp(map.ori_SfTf(:)),size(map.ori_SfTf));

            ang = angle(cat(3, map.st(:).ori));
            ang = reshape(ang(map.ori_SfTf(:)),size(map.ori_SfTf));
            map.ori = amp.*exp(1i*ang);

            amp = abs(cat(3, map.st(:).dir));
            [~, map.dir_SfTf] = max(amp, [], 3, 'linear');
            amp = reshape(amp(map.dir_SfTf(:)),size(map.ori_SfTf));

            ang = angle(cat(3, map.st(:).dir));
            ang = reshape(ang(map.dir_SfTf(:)),size(map.ori_SfTf));
            map.dir = amp.*exp(1i*ang);

            map.mimg = mean(cat(3, map.st(:).mimg),3);
            plot_px_map(map, stim_type, 1, [2 4]);
            print(fullfile(s2p_folder, sprintf('%s_px_maps', stim_type)), '-dpdf', '-vector', '-bestfit')

        end

        %% now compute average across SfTf
        stimSequence.seq = ceil(stimSequence.seq/nSfTf); % treat all sftf as the same
        stimSequence.labels = stimSequence.labels(1:nSfTf:end);

        stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, frameTimes);

        % remove blank from stimMatrix
        stimMatrix(p.nDir+1,:) = [];

        % compute stimulus triggered responses
        nSVD = min(1000, svd.ops.nSVD); % hard coded, perhaps make input

        % filter in time; HARDCODED change for different indicators
        sT = gaussFilt(T(1:nSVD, :)',2);
        
        % measure sta response timecourse 
        [resp, aveResp, ~, kernelTime] = ...
            ppbox.getStimulusSweepsLFR(sT, stimTimes, stimMatrix,frameRate); % responses is (nroi, nStim, nResp, nT)
        aveResp = aveResp(:,:,kernelTime>-1 & kernelTime<3);
       
        % measure response in respWindow
        respWin = [0, 2];
        [~, aveResPeak] = ...
            ppbox.gratingOnResp(resp, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

         % gentle spatial smoothing of the spatial components
        sS =reshape( S(:, 1:nSVD), nY, nX, nSVD);
        sS = imgaussfilt(sS, 0.5);
        sS = reshape(sS, nY*nX, nSVD);

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

        map.all_st.mimg = 1-mat2gray(reshape(sS(:,1), nY, nX))'; % for some reason the SVDs are transposed

        plot_px_map(map.all_st, stim_type);
        print(fullfile(s2p_folder, sprintf('%s_all_st_px_maps', stim_type)), '-dpdf', '-vector', '-bestfit')

    case 'sparsenoise'

        % need to integrate old code from function 'RF_mapping_mov'
end

map.p = p;
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




