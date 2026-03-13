% RUN_GRATING_WORKFLOW Example entry-point script.
%
% Replace the placeholder variables below with the arrays from your
% experiment workspace, then run this script.

% -------------------------------------------------------------------------
% Assign your data here
% -------------------------------------------------------------------------
data = struct();

% Required arrays
% Preferred input layout:
%   data.tc  = timeCourseMatrix;        % [nStim, nRep, nTime, nNeuron]
%   data.amp = responseAmplitudeMatrix; % [nStim, nRep, nNeuron]
%
% where nStim = nDir * nSize * nTf * nSf and the stimulus index follows the
% canonical order [dir, size, tf, sf].
%
% Backward-compatible alternative layouts that are also accepted:
%   full canonical tc  = [dir, size, tf, sf, rep, time, neuron]
%   full canonical amp = [dir, size, tf, sf, rep, neuron]
%   legacy squeezed canonical arrays with inactive stimulus dimensions
%   removed before repeat/time/neuron.
%
% data.tc  = timeCourseMatrix;
% data.amp = responseAmplitudeMatrix;

% Optional blank arrays (accepted but not used in the requested analysis)
% data.blankTc = blankTimeCourseMatrix;
% data.blankAmp = blankAmplitudeMatrix;

% Required vectors
% data.t = timeVector;
% data.directions = directions;
% data.sizes = sizes;
% data.tfs = tfs;
% data.sfs = sfs;

% Optional explicit list of active parameters
% Examples:
% data.activeParams = {'direction'};
% data.activeParams = {'direction', 'size'};
% data.activeParams = {'direction', 'size', 'tf', 'sf'};

% -------------------------------------------------------------------------
% Choose output folder and options
% -------------------------------------------------------------------------
targetFolder = fullfile(pwd, 'grating_analysis_output');

opts = struct();
opts.alpha = 0.05;
opts.saveExt = 'png';   % 'png', 'pdf', 'jpg', etc.
opts.visible = 'off';   % use 'on' if you want figures to appear on screen

% -------------------------------------------------------------------------
% Run analysis
% -------------------------------------------------------------------------
results = analyze_grating_experiment(data, targetFolder, opts); %#ok<NASGU>

disp('Analysis complete. Results struct saved to:');
disp(fullfile(targetFolder, 'grating_analysis_results.mat'));
