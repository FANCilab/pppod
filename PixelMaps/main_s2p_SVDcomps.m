%% CODE to GENERATE SVD COMPRESSION from S2P REGISTERED MOVIES
clear;
%% database of the recording
i = 0;

% i = i+1;
% db(i).subject    = 'HS019'; % animal name
% db(i).date          = '20251021'; % date of the recording
% db(i).exp         = [3]; % all the experiments in the recording
% db(i).expID         = 1; % the experiment you want to compute pixel map of
% db(i).n_planes       = 1;
% db(i).fun_channel   = 1;
% db(i).n_channels    = 2;
% db(i).s2p_version   = 'python';
% db(i).root_storage   ='Z:\Data\2P';


i = i+1;
db(i).subject    = 'NM023'; % animal name
db(i).date          = '20260226'; % date of the recording
db(i).exp         = [1]; % all the experiments in the recording
db(i).expID         = 1; % the experiment you want to compute pixel map of
db(i).n_planes       = 1;
db(i).fun_channel   = 1;
db(i).n_channels    = 1;
db(i).s2p_version   = 'python';
db(i).root_storage   ='Z:\Data\2P';
db(i).stim_type  = 'gratings';

%% Set path to relevant code
if ispc
    code_repo = 'C:\Users\User\Documents\Code\retinotopy';
else
    % code_repo = 'D:\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\Code\retinotopy\';
end

cd(code_repo);
addpath(genpath(code_repo));
addpath(genpath('C:\Users\User\Documents\Code\Suite2P_Matlab'))
% script setting the paths we need
% set_paths(); % edit the paths pointing to the code, contains which packages to download

%% Use the .bin registered file from suite2p to compute SVD decomposition of each recording
for iExp = 1:numel(db)
    %takes a while to run (few mins per db)
    add_SVDcomps(db(iExp))
end

%% USE SVD COMPRESSED MOVIES TO GENERATE STAs MOVIE and TUNING PIXEL MAPS
clear;

%% Set path to relevant code

if ispc
    code_repo = 'C:\Users\Federico\Documents\GitHub\cameraLucida';
else
    code_repo = '/Users/lfedros/Documents/GitHub/cameraLucida';

end
cd(code_repo);
addpath(genpath(code_repo));
set_dendrite_paths(); % edit the paths pointing to the code, contains which packages to download
addpath('\\znas.cortexlab.net\Code\Spikes');
%% database of the recording
clear;
i = 0;

i = i+1;
db(i).mouse_name    = 'FR213'; % animal name
db(i).date          = '2022-06-15'; % date of the recording
db(i).expts         = [4]; % all the experiments in the recording
db(i).expID         = 1; % the experiment you want to compute pixel map of
db(i).expType  = 'gratings';
db(i).nplanes       = 1;
db(i).s2p_version = 'python';
db(i).root_folder ='D:\OneDrive - University College London\Data\2P\';


%% Use SVD compressed recording to make pixel maps

for iExp = 1:numel(db)
    for iPlane = 1:db(i).nplanes
    % pxmap
    [px_map, px_rr]= svd2pxmap(db, iPlane, db(i).expType);
    % run if you want to generate sta movie, otherwise comment out for
    % speed
%     sta_mov= svd2mov_sta(db, 2, db(i).expType);
    end
end

%manually select the dendrites of interest
px_map.doi = draw.freePolyROI(px_map.mimg);

% save px_map
% add here your preferred way of saving px_map, rr, (sta_mov)



%% some test code

load('/Users/federico/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Data/Dendrites/Dendrites paper/Test/dendrite.mat')

trial_stim_resps = px_rr.all_st.trial_stim_resp;

nStim = 12; nRep = 24; %nTfSf = numel(px_rr.st);

trial_stim_resps = reshape(trial_stim_resps, [], nStim, nRep);

%%
trial_perm = randperm(nRep);
even_resp = mean(trial_stim_resps(:, :, trial_perm(2:2:end)),3);
odd_resp = mean(trial_stim_resps(:, :, trial_perm(1:2:end)),3);
all_resp = mean(trial_stim_resps,3);

max_resp = max(all_resp, [],2);
max_thres = prctile(max_resp, 99.9);
max_sg = max_resp >max_thres;

[even_max, even_sort] = max(even_resp, [],2);

[ ~,even_sort] = sort(even_sort, 'ascend');

even_resp = even_resp(even_sort(max_sg), :); 
odd_resp = odd_resp(even_sort(max_sg), :); 

figure; 
subplot(1,2,1)
imagesc(even_resp)
caxis([0 5])
subplot(1,2,2)
imagesc(odd_resp)
caxis([0 5])

%%
resps = px_rr.all_st.stim_resp;

% resps = mean(trial_stim_resps,3);
[max_resp, dir_sort] = max(resps, [],2);

[ ~,dir_sort] = sort(dir_sort, 'ascend');

max_thres = prctile(max_resp, 95);

max_sg = max_resp >max_thres;

resps = resps(dir_sort(max_sg), :); 

% resps = bsxfun(@minus, resps, min(resps, [], 2));
% resps = bsxfun(@rdivide, resps, max(resps, [], 2));

figure; 
imagesc(resps)
caxis([0 5])

%%
figure;

[max_resp, dir_sort] = max(px_rr.st(2).stim_resp, [],2);

[ ~,dir_sort] = sort(dir_sort, 'ascend');

max_thres = prctile(max_resp, 99.9);

max_sg = max_resp >max_thres;

for iSt = 1:numel(px_rr.st)

resps = px_rr.st(iSt).stim_resp;

resps = resps(dir_sort(max_sg), :); 

% resps = bsxfun(@minus, resps, min(resps, [], 2));
% resps = bsxfun(@rdivide, resps, max(resps, [], 2));

subplot(1,numel(px_rr.st), iSt);
imagesc(resps)
caxis([0 5])
end



