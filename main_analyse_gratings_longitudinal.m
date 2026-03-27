clear;

%% database of the recording
% 
i = 0;

i = i+1;
db(i).subject    = 'HS028'; % animal name
db(i).date          = {'20260211', '20260216', '20260217', '20260218', '20260224'}; % date of the recording
db(i).exp         = {[2],[2],[2],[2],[2]} ; % all the experiments in the recording
db(i).expID         = [1 1 1 1 1]; % the experiment you want to compute pixel map of
db(i).n_planes       = 5;
db(i).fun_channel   = 1;
db(i).n_channels    = 2;
db(i).s2p_version   = 'python';
db(i).root_storage   ='Z:\Data\2P';
db(i).stim_type  = 'gratings';
db(i).nDays = numel(db(i).date);

%% Set path to relevant code

if ispc
    code_repo = 'C:\Users\User\Documents\Code\pppod';
else
    % code_repo = 'D:\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\Code\retinotopy\';
end

cd(code_repo);
addpath(genpath(code_repo));
addpath(genpath('C:\Users\User\Documents\Code\Suite2P_Matlab'))
addpath('C:\Users\User\Documents\Code\FedBox');
addpath(genpath('C:\Users\User\Documents\Code\npy-matlab'));
% addpath('C:\Users\User\Documents\Code\FedBox\rastermap_matlab');


%% load the data
% edit this function to point to your data folders

iDay = 1;

info= getExpInfo(db(i).subject , db(i).date{iDay} , db(i).exp{iDay}(db(i).expID(iDay)), 1);
info.folderDay = fullfile(info.localData, db(i).subject , db(i).date{iDay} , num2str(db(i).exp{iDay}(db(i).expID(iDay))));

targetPlane = 5;
s2p_folder = 'D:\Data\suite2p\HS028\20260211_20260217_20260218_20260224\2_2_2_2\plane4';
s2p_file = sprintf('%s/Fall.mat', s2p_folder);
s2p_iscell = sprintf('%s/iscell.npy', s2p_folder);
s2p_redcell = sprintf('%s/redcell.npy', s2p_folder);


load(s2p_file);
readNPY(s2p_iscell);
readNPY(s2p_redcell);

%% Extract only frames of target Day
%need to write a function that does this properly

if iDay == 1
    thisDayStartFrame = 1;
    thisDayEndFrame = ops.frames_per_folder(iDay);
else
thisDayStartFrame = sum(ops.frames_per_folder(1:iDay-1));
    thisDayEndFrame = thisDayStartFrame + ops.frames_per_folder(iDay)-1;
end

F = F(:, thisDayStartFrame:thisDayEndFrame);
Fneu = Fneu(:, thisDayStartFrame:thisDayEndFrame);

%% Plot neuropil corrected fluorescent traces

bad_cells = sum(F, 2) == 0;

iscell(bad_cells,1) = 0;
redcell(bad_cells,:) = 0;

neurons_raw = F(logical(iscell(:,1)), :);
iscellred = redcell(logical(iscell(:,1)),2);
iscellred = iscellred>0.7;
%% subtract neuropil

neuropil = Fneu(logical(iscell(:,1)), :);

neurons = s2pUtils.estimateNeuropil_LFR(neurons_raw, neuropil);

[nN, nFrames] = size(neurons);

neurons = zscore(neurons, [], 2);

%% load the stimulus data

event = bonsai.load_events(info);

planeFrameTimes = event.frame.on(targetPlane:info.nPlanes:end);

planeRate = info.volumeRate;

stimTimes = event.grating.stimTimes;

stimMatrix = event.grating.stimMatrix(event.grating.stimIdx, targetPlane:info.nPlanes:end);
stimMatrix_blank = event.grating.stimMatrix(event.grating.blankStimIdx, targetPlane:info.nPlanes:end);

%% now compute average across stims

% measure sta response timecourse
[resp, aveResp, ~, kernelTime] = ...
    ppbox.getStimulusSweepsLFR(neurons', stimTimes, stimMatrix,planeRate); % responses is (nroi, nStim, nResp, nT)

[resp_blank, aveResp_blank] = ...
    ppbox.getStimulusSweepsLFR(neurons', stimTimes, stimMatrix_blank,planeRate); % responses is (nroi, nStim, nResp, nT)

% measure response in respWindow
respWin = [0, 2];
[resPeak, aveResPeak] = ...
    ppbox.gratingOnResp(resp, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

[resPeak_blank, aveResPeak_blank] = ...
    ppbox.gratingOnResp(resp_blank, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

%% for other neurons
neurons = ~iscellred;

data.tc =   permute(resp(neurons,:,:,:), [2 3 4 1]); % [dir, size, tf, sf, repeat, time, neuron]
data.amp = permute(resPeak(neurons, :,:), [2 3 1]); % [dir, size, tf, sf, repeat, neuron]
data.blankTc = permute(resp_blank(neurons,:,:,:), [2 3 4 1]);
data.blankAmp = permute(resPeak_blank(neurons, :,:), [2 3 1]);
data.t = kernelTime;
data.directions =   event.grating.uniquePars.direction(1:end-1);
data.sizes   =   event.grating.uniquePars.size;
data.tfs  =   event.grating.uniquePars.tf;
data.sfs  =   event.grating.uniquePars.sf;
data.activeParams = event.grating.activePars;
%
targetFolder  = fullfile(info.folderDay, 'Results', sprintf('plane_%d', targetPlane-1));
%
opts.alpha = 0.05;
opts.saveExt = 'pdf';
opts.visible = 'off';

results = analyze_grating_experiment(data, targetFolder, opts);
    figDist = plot_results_distributions(results, targetFolder, opts);

%% for red neurons
neurons = iscellred;
targetFolder  = fullfile(info.folderDay, 'Results_red', sprintf('plane_%d', targetPlane-1));

data.tc =   permute(resp(neurons,:,:,:), [2 3 4 1]); % [dir, size, tf, sf, repeat, time, neuron]
data.amp = permute(resPeak(neurons, :,:), [2 3 1]); % [dir, size, tf, sf, repeat, neuron]
data.blankTc = permute(resp_blank(neurons,:,:,:), [2 3 4 1]);
data.blankAmp = permute(resPeak_blank(neurons, :,:), [2 3 1]);

results = analyze_grating_experiment(data, targetFolder, opts);

