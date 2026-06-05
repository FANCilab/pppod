clear;

%% database of the recording

i = 0;

i = i+1;
db(i).subject    = 'NM035'; % animal name
db(i).date          = '20260513'; % date of the recording
db(i).exp         = [1]; % all the experiments in the recording
db(i).expID         = 1; % the experiment you want to compute pixel map of
db(i).n_planes       = 1;
db(i).fun_channel   = 1;
db(i).n_channels    = 1;
db(i).s2p_version   = 'python';
db(i).root_storage   ='Z:\Data\2P';
db(i).stim_type  = 'sparse_noise';

%% Set path to relevant code

if ispc
    code_repo = 'C:\Users\User\Documents\Code\pppod';
else
    % code_repo = 'D:\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\Code\retinotopy\';
end

addpath(genpath(code_repo));
addpath(genpath('C:\Users\User\Documents\Code\npy-matlab'));

%% load the data
% edit this function to point to your data folders
info= getExpInfo(db(i).subject , db(i).date , db(i).exp(db(i).expID), 1);

targetPlane = 1;
% targetplane = 'combined';

switch db.s2p_version
    case 'python'

        if isnumeric (targetPlane)
            % python indexes from 0
            s2p_folder = fullfile(info.folder2p,sprintf('plane%d', targetPlane-1));
        else
            s2p_folder = fullfile(info.folder2p,targetPlane-1);
        end
end

try
    s2p_file = sprintf('%s/Fall.mat', s2p_folder);
    load(s2p_file);
catch

    % s2p_ops = fullfile(s2p_folder, 'ops.npy');
    % readNPY(s2p_ops)

    F = fullfile(s2p_folder, 'F.npy');
    readNPY(F);
    Fneu = fullfile(s2p_folder, 'Fneu.npy');
    readNPY(Fneu);
    s2p_iscell = sprintf('%s/iscell.npy', s2p_folder);
    readNPY(s2p_iscell);

end
%% Plot neuropil corrected fluorescent traces

bad_cells = sum(F, 2) == 0;

iscell(bad_cells,1) = 0;

neurons_raw = F(logical(iscell(:,1)), :);

%% subtract neuropil

neuropil = Fneu(logical(iscell(:,1)), :);

neurons = s2pUtils.estimateNeuropil_LFR(neurons_raw, neuropil);

[nN, nFrames] = size(neurons);

neurons = zscore(neurons, [], 2);

%% load the stimulus data

event = bonsai.load_events(db);

planeFrameTimes = event.frame.on(targetPlane:info.nPlanes:end); % time stamps of imaging frames

planeRate = info.volumeRate;

stimFrames = event.sparse_noise.frames; % time stamps of sparse noise frames

stimFrameTimes = event.frame.on'; 

stimTimes.onset = 0;

stimPosition = [-135 135 -40 40]; % boundaris of the stim grid in visual degrees (##HARDCODED, FIX)

