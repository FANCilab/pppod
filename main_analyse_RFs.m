clear;

%% database of the recording

i = 0;

i = i+1;
db(i).subject    = 'HS028'; % animal name
db(i).date          = '20260216'; % date of the recording
db(i).exp         = [2]; % all the experiments in the recording
db(i).expID         = 1; % the experiment you want to compute pixel map of
db(i).n_planes       = 5;
db(i).fun_channel   = 1;
db(i).n_channels    = 1;
db(i).s2p_version   = 'python';
db(i).root_storage   ='Z:\Data\2P';
db(i).stim_type  = 'gratings';

%% Set path to relevant code

if ispc
    code_repo = 'C:\Users\User\Documents\Code';
else
    % code_repo = 'D:\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\Code\retinotopy\';
end

cd(code_repo);
addpath(genpath(code_repo));
addpath(genpath('C:\Users\User\Documents\Code\Suite2P_Matlab'))
addpath('C:\Users\User\Documents\Code\FedBox');
addpath('C:\Users\User\Documents\Code\FedBox\rastermap_matlab');


%% load the data
% edit this function to point to your data folders
info= getExpInfo(db(i).subject , db(i).date , db(i).exp(db(i).expID), 1);

targetPlane = 4;
% targetplane = 'combined';

switch db.s2p_version
    case 'python'

        if isnumeric (targetPlane)
            % python indexes from 0
            s2p_folder = fullfile(info.folder2p,sprintf('plane%d', targetPlane-1));
        else
            s2p_folder = fullfile(info.folder2p,targetPlane-1);
        end

        s2p_file = sprintf('%s/Fall.mat', s2p_folder);

end

load(s2p_file);

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

planeFrameTimes = event.frame.on(targetPlane:info.nPlanes:end);

planeRate = info.volumeRate;

stimTimes = event.grating.stimTimes;

%% Parameters
% for correcting baseline drifts of calcium traces at start of experiments
win_decay = 20; % in s, window to test whether baseline is higher than normal
thresh_decay = 1.5; % in std, threshold for drift
win_correct = 150; % in s, window to fit exponential

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously), and fitting running
% kernels and RFs simultaneously
lambdas = logspace(-4, -1, 4);
rf_timeLimits = [0 0.4];
crossFolds = 10;

% for evaluation of receptive fields (significance/goodness)
minExplainedVariance = 0.01;
maxPVal = 0.05;
% for plotting RFs
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);
titles = {'ON field','OFF field'};

% load data
caData = io.getCalciumData(f);
stimData = io.getVisNoiseInfo(f);
t_stim = stimData.times;
tBin_stim = median(diff(t_stim));
t_rf = (floor(rf_timeLimits(1) / tBin_stim) : ...
    ceil(rf_timeLimits(2) / tBin_stim)) .* tBin_stim;

% prepare calcium traces
% interpolate calcium traces to align all to same time
t_ind = caData.time > t_stim(1) - 10 & ...
    caData.time < t_stim(end) + 10;
caTraces = caData.traces(t_ind,:);
t_ca = caData.time(t_ind);
[caTraces, t_ca] = traces.alignSampling(caTraces, t_ca, ...
    caData.planes, caData.delays);

% remove strong baseline decay at start of experiment in cells that
% show it
caTraces = traces.removeDecay(caTraces, t_ca, win_decay, ...
    win_correct, thresh_decay);

stimFrames = stimData.frames(stimData.stimOrder,:,:);

% map RF
% rFields: [rows x columns x t_rf x ON/OFF x units]
[rFields, ev] = ...
    whiteNoise.getReceptiveField(caTraces, t_ca, ...
    stimFrames, t_stim, t_rf./tBin_stim, ...
    lambdas, crossFolds);

v = squeeze(mean(ev,3)); % [neuron x lamStim]
% average EV across cross-folds
[maxEV, maxStimLam] = max(v,[],2);
bestLambdas = lambdas(maxStimLam)';

% test signficance of each RF (note: resulting ev are not
% cross-validated, while maxEV are)
[ev, ev_shift] = ...
    whiteNoise.receptiveFieldShiftTest( ...
    caTraces, t_ca, stimFrames, t_stim, ...
    t_rf./tBin_stim, rFields, bestLambdas, 500);
pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
pvals(isnan(ev)) = NaN;

% save results
results.maps = permute(rFields, [5 1 2 3 4]);
results.explVars = maxEV;
results.lambdas = bestLambdas;
results.pValues = pvals;
results.timestamps = t_rf;
results.edges = stimData.edges;
io.writeNoiseRFResults(results, f);


for iUnit = 1:length(results.explVars)
    if isnan(results.explVars(iUnit)) || ...
            results.explVars(iUnit) < minExplainedVariance || ...
            results.pValues(iUnit) > maxPVal
        continue
    end
    % rf: [rows x columns x time x ON/OFF];
    rf = squeeze(results.maps(iUnit,:,:,:,:));
    rf(:,:,:,2) = -rf(:,:,:,2);
    [mx,mxTime] = max(max(abs(rf),[],[1 2 4]));
    stimPos = results.edges; % [left, right, top (negative), bottom]
    squW = diff(stimPos(1:2)) / size(rf,1);
    squH = diff(stimPos(3:4)) / size(rf,2);

    figure('Position', [75 195 1470 475])
    for sf = 1:2
        subplot(1,2,sf)
        imagesc([stimPos(1)+squW/2 stimPos(2)-squW/2], ...
            [stimPos(3)+squH/2 stimPos(4)-squH/2], ...
            rf(:,:,mxTime,sf),[-mx mx])
        axis image
        set(gca, 'box', 'off')
        colormap(gca, cms(:,:,sf))
        title(titles{sf})
        colorbar
    end
    sgtitle(sprintf(...
        'ROI %d (lam: %.0e, t: %.2fs, EV: %.3f, pVal: %.3f)', ...
        iUnit, results.lambdas(iUnit), ...
        results.timestamps(mxTime), results.explVars(iUnit), ...
        results.pValues(iUnit)))

    saveas(gcf, fullfile(fPlots, sprintf('Unit%03d_noise.jpg', iUnit)));
    close gcf
end
end