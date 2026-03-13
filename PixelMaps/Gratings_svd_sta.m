function sta_mov = svd2mov_sta(db, targetPlane, stim_type)

if nargin <2
nTile = 5;
end

info = ppbox.infoPopulateTempLFR(db.mouse_name, db.date, db.expts(db.expID));

allFrameTimes = ppbox.getFrameTimes(info, 'all');



%%

switch stim_type
    case 'gratings'
%% load stimulus info
        info = ppbox.infoPopulateTempLFR(db.mouse_name{iExp}, db.date{iExp}, db.expts{iExp}(db.expID{iExp}));
        nFrames = numel(dF{iExp});
        planeFrames = db.plane(iExp):info.nPlanes:(nFrames*info.nPlanes);
        frameTimes{iExp} = ppbox.getFrameTimes(info, planeFrames);
        frameRate = (1/mean(diff(frameTimes{iExp})));
        stimTimes = ppbox.getStimTimes(info);
        stimSequence = ppbox.getStimSequence_LFR(info);
        stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, frameTimes{iExp});
        % select a subset of stimuli if specified in the db
        if isfield(db, 'stimList')
            stimSet = db.stimList{iExp};
            stimLabels{iExp} = stimSequence.labels(stimSet);
        end
        
        %% compute stimulus triggered responses
        % STA of neuropil
        [~, aveRespNeu] = ...
            getStimulusSweepsLFR(Fneu, stimTimes, stimMatrix,frameRate, stimSet); % responses is (nroi, nStim, nResp, nT)
        % this is the average response of neuropil across stimuli
        F0_neu(iExp) = max(mean(squeeze(aveRespNeu), 1));
        
        %STA of neuron
        [responses{iExp}, aveResp{iExp}, seResp{iExp}, kernelTime{iExp}, stimDur{iExp}] = ...
            getStimulusSweepsLFR(dF{iExp}, stimTimes, stimMatrix,frameRate, stimSet); % responses is (nroi, nStim, nResp, nT)
        % measure response in respWindow
        [resPeak{iExp}, aveResPeak{iExp}, seResPeak{iExp}] = ...
            Ori.gratingOnResp(responses{iExp}, kernelTime{iExp}, respWin);  % resPeak is (nroi, nStim, nResp)
        %     plotSweepResp(responses{iExp}, kernelTime{iExp}, stimDur{iExp});









end


%%



stimTimes = ssLocal.getStimTimes(info);

 [stimFrames, stimFrameTimes, stimPosition] = whiteNoise.getStimulusFrames(info);

% [rangeY, rangeX, ~] = size(stimFrames);

% retX = NaN(info.fovY, info.fovX, info.nPlanes);
% retY = NaN(info.fovY, info.fovX, info.nPlanes);
% 
% interpRetX = retX;
% interpRetY = retY;

% sigma = 10/(max(micronsX(:))/info.fovX);

iPlane = targetPlane;
    
    [root, refF, refNeu, refSVD] = starter.getAnalysisRefs(db.mouse_name, db.date, db.expts, iPlane);
    
    svd = load(fullfile(root, refSVD));
    
    load(fullfile(root, [refF(1:end-9), '.mat']), 'ops');
    
%     if isfield(ops, 'mimgRED')
%         fovR(:,:,iPlane) = ops.mimgRED;
%     end
%     
%     fovG(:,:,iPlane) = svd.ops.mimg;
    
    [nY, nX, nBasis] = size(svd.U);
    
    nFrames = size(svd.Vcell{db.expID},2);
    
    
%     frameTimes = ppbox.getFrameTimes(info, planeFrames);
        
%     S = imgaussfilt(svd.U, 3*sigma); % smooth spatial basis with gaussian 30 um sigma
    S = svd.U; 
    S = reshape(S, [], nBasis); % nPx * nB
    T = svd.Vcell{db.expID} ; % nB * nT
    planeFrames = iPlane:info.nPlanes:(size(T, 2)*info.nPlanes);
    allFrameTimes = allFrameTimes(planeFrames);

    
    
%                 Fneu = S(iPx, :)*svd.Vcell{db.expID};

nRep = numel(stimTimes.onset);
[w, h, nFrames] = size(stimFrames);

% switch respType
%     case 'on'
%         igood = stimFrames(:,:,2:nFrames) == 1 & stimFrames(:,:,1:nFrames-1) == 0;
%     case 'off'
%         igood = stimFrames(:,:,2:nFrames) == -1 & stimFrames(:,:,1:nFrames-1) == 0;
%     case 'onoff'
%         igoodOn = stimFrames(:,:,2:nFrames) == 1 & stimFrames(:,:,1:nFrames-1) == 0;
%         igoodOff = stimFrames(:,:,2:nFrames) == -1 & stimFrames(:,:,1:nFrames-1) == 0;
%         igood = igoodOn | igoodOff;
%     case 'any'
%         igood = stimFrames(:,:,2:nFrames) ~= stimFrames(:,:,1:nFrames-1);
% end

igood = abs(stimFrames(:,:,2:nFrames)) == 1 & abs(stimFrames(:,:,1:nFrames-1)) == 0;
igood = reshape(igood, w*h, []);

onidx = cell(w,h);
for iRep = 1:nRep
    for ix = 1: w*h
    onidx{ix} = [onidx{ix}, stimFrameTimes(find(igood(ix, :))+1) + stimTimes.onset(iRep)];
    end
end

target = [4, 8, 16, 23];
onidx = onidx(target(1):target(2),target(3):target(4));

periT = 0:0.3:2.1;
nt = numel(periT);
ETAmat = [];
for ix = 1:numel(onidx)
    
 [~, ETAmat(:,ix,:)] = magicETA2(allFrameTimes, T(2:end, :)', onidx{ix}, periT);
 
 ix/numel(onidx)
end

newT = reshape(ETAmat, [], nBasis-1)';
% 
%    nTile = 5;
%   tileX = round(linspace(1,nX, nTile +1));
%     tileY = round(linspace(1,nY, nTile+1));
%       for iTlX = 2
%         for iTlY = 3
%             [xPx, yPx] = meshgrid(  tileX(iTlX):tileX(iTlX+1), tileY(iTlY):tileY(iTlY+1));
%             iPx = sub2ind([nY,nX], yPx(:), xPx(:));
%         end
%     end
% 
% sta_mov =  S(iPx, :)*newT;
% 
% sta_mov = reshape(sta_mov, size(xPx,1), size(xPx,2), []); 

sta_mov =  S(:, 2:end)*newT;
sta_mov = zscore(sta_mov')';
sta_mov = reshape(sta_mov, nY, nX, []);
% sta_mov = bsxfun(@minus, sta_mov, min(sta_mov, [], 2));
% sta_mov = bsxfun(@rdivide, sta_mov, max(sta_mov, [], 2));

sta_mov = reshape(sta_mov, nY, nX, [], numel(onidx));
for ix = 1:numel(onidx)
% sta_mov(:,:,:,ix)  = sta_mov(:,:,:,ix)/max(makeVec(sta_mov(:,:,:,ix)));
sta_mov(:,:,:,ix) = imgaussfilt3(sta_mov(:,:,:,ix), [0.1 0.1 1]);
end
sta_mov = reshape(sta_mov, nY, nX, []);

end




