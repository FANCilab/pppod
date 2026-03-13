clear;

%% database of the recording
% i = 0;
% 
% i = i+1;
% db(i).mouse_name    = 'NM023'; % animal name
% db(i).date          = '20251224'; % date of the recording
% db(i).exp_n         = [5]; % all the experiments in the recording
% db(i).expID         = 1; % the experiment you want to compute pixel map of
% db(i).n_planes       = 1;
% db(i).fun_channel   = 1;
% db(i).n_channels    = 1;
% db(i).s2p_version   = 'python';
% db(i).root_storage   ='Z:\Data\2P';
% db(i).expType  = 'gratings';

i = 0;

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
    code_repo = 'C:\Users\User\Documents\Code';
else
    % code_repo = 'D:\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\Code\retinotopy\';
end

cd(code_repo);
addpath(genpath(code_repo));
addpath(genpath('C:\Users\User\Documents\Code\Suite2P_Matlab'))
addpath('C:\Users\User\Documents\Code\FedBox');
addpath('C:\Users\User\Documents\Code\FedBox\rastermap_matlab');


%% Use SVD compressed recording to make pixel maps
info= getExpInfo(db(i).subject , db(i).date , db(i).exp(db(i).expID), 1);

for iExp = 1:numel(db)
    for iPlane = 1:info.nPlanes
    % pxmap
    [px_map, px_rr]= svd2pxmap(db, iPlane, db(i).stim_type);
    % run if you want to generate sta movie, otherwise comment out for
    % speed
%   sta_mov= svd2mov_sta(db, 2, db(i).expType);
    end
end

%% load the data

info= getExpInfo(db(i).subject , db(i).date , db(i).exp(db(i).expID), 1);

targetplane = 4;
% targetplane = 'combined';

switch db.s2p_version
    case 'python'

        if isnumeric (targetplane)
            targetplane = targetplane-1; % python indexes from 0
s2p_folder = fullfile(info.folder2p,sprintf('plane%d', targetPlane));
        else
s2p_folder = fullfile(info.folder2p,targetPlane);
        end

s2p_file = sprintf('%s/Fall.mat', s2p_folder);

    case 'matlab'

s2p_folder = fullfile(root_folder, db.mouse_name, db.date, ...
    sprintf('%d_', db.expts));
s2p_folder = s2p_folder(1:end-1);

end

load(s2p_file);

%% PLot neuropil corrected fluorescent traces

bad_cells = sum(F, 2) == 0;

iscell(bad_cells,1) = 0;

neurons_raw = F(logical(iscell(:,1)), :);

[nN, nFrames] = size(neurons_raw);

frameTime = (1:nFrames)/ops.fs;

%% de-bleach if needed (modify to do after neuropil correction)

[neurons_deBleached, indunits] = s2pUtils.deBleach(neurons_raw', frameTime);

% figure; plot(frameTime, neurons_deBleached(:, 100)); hold on;

%% subtract neuropil

neuropil = Fneu(logical(iscell(:,1)), :);

neurons = s2pUtils.estimateNeuropil_LFR(neurons_raw, neuropil);

[nN, nFrames] = size(neurons);

neurons = zscore(neurons, [], 2);

% interestingN = [18 40 45 53 66 72 81 87 91 102 105 108 217 267];

% figure;
% PlotDisplacedLFRDev(frameTimes, gaussFilt(neurons(interestingN, :)',10),3)

%% raster map

[nN, nT] = size(neurons);

[iclustup, isort, Vout] = activityMap(gaussFilt(neurons',10)');

neurons_sorted = neurons(isort,:);

figure('Color', 'White');
imagesc(frameTime, 1:nN, imgaussfilt(neurons_sorted,1) ); hold on
xlabel('Time(s)');
ylabel('Neurons (zscore)');
%% load the stimulus data

event = bonsai.load_events(db);

planeFrameTimes = event.frame.on(targetplane:info.nPlanes:end);

planeRate = info.volumeRate;

stimTimes = event.grating.stimTimes;

stimMatrix = event.grating.stimMatrix(event.grating.stimIdx, targetplane:info.nPlanes:end);
stimMatrix_blank = event.grating.stimMatrix(event.grating.blankStimIdx, targetplane:info.nPlanes:end);

% dirs = p.dirs; 
% dirs = deg2rad(dirs);
% oris = dirs;
% oris(oris >= pi) = oris(oris >=pi) -pi;
% oris = oris -pi/2;
% oris = oris*2;

%% now compute average across stims

% measure sta response timecourse
[resp, aveResp, ~, kernelTime] = ...
    ppbox.getStimulusSweepsLFR(neurons', stimTimes, stimMatrix,planeRate); % responses is (nroi, nStim, nResp, nT)

[resp_blank, aveResp_blank] = ...
    ppbox.getStimulusSweepsLFR(neurons', stimTimes, stimMatrix_blank,frameRate); % responses is (nroi, nStim, nResp, nT)

% measure response in respWindow
respWin = [0, 2];
[resPeak, aveResPeak] = ...
    ppbox.gratingOnResp(resp, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

[resPeak_blank, aveResPeak_blank] = ...
    ppbox.gratingOnResp(resp_blank, kernelTime, respWin);  % resPeak is (nroi, nStim, nResp)

%%

  data.tc =   permute(resp(:,:,:,:), [2 3 4 1]); % [dir, size, tf, sf, repeat, time, neuron]
  data.amp = permute(resPeak(:, :,:), [2 3 1]); % [dir, size, tf, sf, repeat, neuron]
  data.blankTc = resp_blank;
  data.blankAmp = resPeak_blank;
  data.t = kernelTime;
  data.directions =   event.grating.uniquePars.direction(1:end-1);
  data.sizes   =   event.grating.uniquePars.size;
  data.tfs  =   event.grating.uniquePars.tf;
  data.sfs  =   event.grating.uniquePars.sf;
  data.activeParams = event.grating.activePars;
%
  targetFolder  = fullfile(info.folder2p, 'Results');
%
  opts.alpha = 0.05;
  opts.saveExt = 'png';
  opts.visible = 'on';

results = analyze_grating_experiment(data, targetFolder, opts);


%%

ppbox.plotSweepResp_LFR(permute(resp(interestingN,:,:,:), [2 3 4 1]), kernelTime, 2)

%% 
%% compute average and se of responses
nRep = size(tune.allResp,2);
tune.aveResp = squeeze(nanmean(tune.allResp, 2));
tune.seResp = squeeze(nanstd(tune.allResp, [], 2))/sqrt(nRep);
tune.avePeak = squeeze(nanmean(tune.allPeaks,2));
tune.sePeak = squeeze(nanstd(tune.allPeaks,[], 2))/sqrt(nRep);
tune.aveOriPeak = nanmean(cat(2, tune.allPeaks(1:6, :), tune.allPeaks(7:12, :)),2);
tune.seOriPeak = nanstd(cat(2, tune.allPeaks(1:6, :), tune.allPeaks(7:12, :)),[], 2)/sqrt(nRep*2);

%% fit model tuning curve

toFit= makeVec(tune.allPeaks(1:end-1, :))';
nan_resp = isnan(toFit);
fitDirs = repmat(tune.dirs, 1,nRep);

% % double gaussian (MatteoBox) direction tuning
[tune.dir_pars_g, ~] = fitori(fitDirs (~nan_resp), toFit(~nan_resp));
tune.fit_g = oritune(tune.dir_pars_g, 0:1:359);

%double von Mises direction tuning
tune.dir_pars_vm = mfun.fitTuning(fitDirs, toFit, 'vm2', fixPars.dir);
tune.fit_vm = mfun.vonMises2(tune.dir_pars_vm, 0:1:359);
tune.fit_vm_12 = mfun.vonMises2(tune.dir_pars_vm, 0:30:330);
tune.fit_pt = 0:1:359;
tune.prefDir = tune.dir_pars_vm(1);
tune.prefDir_Ori = tune.prefDir -90;
tune.prefDir_Ori(tune.prefDir_Ori >90) = tune.prefDir_Ori(tune.prefDir_Ori >90) -180;

% von Mises orientation tuning
fitOris = repmat(tune.oris*2, 1,nRep);
tune.ori_pars_vm = mfun.fitTuning(fitOris, toFit, 'vm1', fixPars.ori);
tune.ori_fit_vm = mfun.vonMises(tune.ori_pars_vm, -180:1:179);
tune.ori_fit_pt = (-180:1:179)/2;
tune.prefOri = tune.ori_pars_vm(1)/2;
tune.prefOri = unwrap_angle(tune.prefOri, 1,1);

%% measure parameters from tuning curves fit: prefDir, prefOri, Rp, Rn, Ro, Rb, DS, OS

tune.Ro = min(tune.fit_vm);

% pp = findpeaks(tune.fit_vm, 'SortStr', 'descend');
% if numel(pp) ==2
%     [tune.Rp, tune.Rn] = vecdeal(pp);
% elseif numel(pp)==1
%     [tune.Rp, tune.Rn] = vecdeal([pp, tune.Ro]);
% elseif numel(pp)==0
%     [tune.Rp, tune.Rn] = vecdeal([max(tune.fit_vm), tune.Ro]);
% else
%     warning('Weird tuning curve');
% end

% tune.Ro = tune.dir_pars_vm(4); 

tune.Rp = tune.dir_pars_vm(2) + tune.dir_pars_vm(4); 

tune.Rn = tune.dir_pars_vm(3) + tune.dir_pars_vm(4); 

tune.Rb = tune.avePeak(13); % response to blank

tune.DS = (tune.Rp-tune.Rn)/(tune.Rp+abs(tune.Rn)); %[0 1], works if Rp, Rn >0, or Rp>0, Rn<0, but not Rp,Rn<0

tune.Rp_ori = max(tune.ori_fit_vm); % peak of orientation tuning resp 

tune.Ro_ori = min(tune.ori_fit_vm); % min of orientation tuning resp

% orientation selectivity linear
tune.OS = (tune.Rp_ori-tune.Ro_ori)/(tune.Rp_ori+abs(tune.Ro_ori)); %[0 1], works if Rp, Rn >0, or Rp>0, Rn<0, but not Rp,Rn<0
% tune.OS_circ = circ_var(repmat(tune.oris*2*pi/180, 1,nRep), (toFit-min(toFit)));

% orientation selectivity circular
bsl = min(tune.fit_vm_12);
if bsl>=0
    
%     tune.OS_circ = circ_var(tune.oris*2*pi/180, tune.fit_vm_12');
    tune.OS_circ = circ_r(tune.oris*2*pi/180, tune.fit_vm_12');
    
elseif sum(tune.fit_vm_12<0) <12
%     tune.OS_circ = circ_var(tune.oris*2*pi/180, tune.fit_vm_12' - min(tune.fit_vm_12));
    tune.OS_circ = circ_r(tune.oris*2*pi/180, tune.fit_vm_12' - min(tune.fit_vm_12));

else
%         tune.OS_circ = circ_var(tune.oris*2*pi/180, abs(tune.fit_vm_12'));
        tune.OS_circ = circ_r(tune.oris*2*pi/180, abs(tune.fit_vm_12'));

end

% orientation selectivity circular
bsl = min(tune.fit_vm_12);
if bsl>=0
    
%     tune.OS_circ = circ_var(tune.oris*2*pi/180, tune.fit_vm_12');
    tune.DS_circ = circ_r(tune.dirs*pi/180, tune.fit_vm_12');
    
elseif sum(tune.fit_vm_12<0) <12
%     tune.OS_circ = circ_var(tune.oris*2*pi/180, tune.fit_vm_12' - min(tune.fit_vm_12));
    tune.DS_circ = circ_r(tune.dirs*pi/180, tune.fit_vm_12' - min(tune.fit_vm_12));

else
%         tune.OS_circ = circ_var(tune.oris*2*pi/180, abs(tune.fit_vm_12'));
        tune.DS_circ = circ_r(tune.dirs*pi/180, abs(tune.fit_vm_12'));

end

%%

 %% plot stimulus responses
    
    yval(1) = min([tuning.aveResp(:); tuning.aveResp(:)]);
    yval(2) = max([tuning.aveResp(:); tuning.aveResp(:)]);
    
    [~, peak_ori] = max(tuning.avePeak(1:12)) ;
    stim_vals = 1:12;
    stim_vals = circshift(stim_vals, 7-peak_ori);
    d_dir_vals = -180:30:150;
    for iStim = 1:12
        
        subplot(nCut, nCols, iStim+nCols*(iCut-1))
        
        plot_stim_response(tuning, stim_vals(iStim), [.7 .7 .7], yval);
        if iCut>1
            plot_stim_response(tuning_cut, stim_vals(iStim), line_c, yval);
        end
        if iCut == nCut
            or = (iStim-1)*30;
            
            xlabel(sprintf('%d deg \n(\\Delta pref dir)', d_dir_vals(iStim)));
            set(get(gca,'XLabel'),'Visible','on')
        end
    end

    %%

      %% plot direction tuning
    
    subplot(nCut, nCols, 13 +nCols*(iCut-1))
    
    yval(1) = min([tuning.avePeak(:); tuning_cut.avePeak(:)]);
    yval(2) = max([tuning.avePeak(:); tuning_cut.avePeak(:)]);
    
    plot_dir_tuning(tuning, [.7 .7 .7], yval);
    if iCut>1
        
        plot_dir_tuning(tuning_cut, line_c, yval);
    end
    
    if iCut == nCut
        
        xlabel('Direction')
    end
    ylabel('Norm resp')
    
    % title(sprintf('%d', round(neuron.tuning.prefDir)));
    
    %% plot orientation tuning
    
    subplot(nCut, nCols, 14+nCols*(iCut-1))
    
    yval(1) = min([tuning.aveOriPeak(:); tuning_cut.aveOriPeak(:)]);
    
    yval(2) = max([tuning.aveOriPeak(:); tuning_cut.aveOriPeak(:)]);
    
    plot_ori_tuning(tuning, [.7 .7 .7], yval);
    if iCut>1
        
        plot_ori_tuning(tuning_cut, line_c, yval);
    end
    
    if iCut == nCut
        xlabel('Orientation')
    end
    ylabel('Norm resp')
    % title(sprintf('%d', round(neuron.tuning.prefOri)));
    formatAxes