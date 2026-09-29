%% ALF_RF_pipeline


%clear;
%clc;
%% Path setup
repoRoot = fileparts(mfilename('fullpath'));
addpath(genpath(repoRoot));
opts.npyMatlabPath = 'C:\Users\ichaker\Desktop\RFmapping\npy-matlab\npy-matlab';

ALFroot = 'D:\Pipeline\DataTest'; % ALF export path:                      ALFroot/<subject>/<date>/<session>
resultsRoot = 'D:\Pipeline\ResultsTest'; % RF results path:          resultsRoot/<subject>/<date>/<session>
plotsRoot = 'D:\Pipeline\PlotsTest'; % Grouped plotting path: plotsRoot/<group>/<subject>/<date>/<session>

% s2p output paths to process
s2pPaths = [
    "\\10.233.25.135\FANCiNAS1\Data\2P\NM023\Processed\20260209\1_2"
    %"\\10.233.25.135\FANCiNAS1\Data\2P\NM039\Processed\20260630\1"
    %"\\10.233.25.135\FANCiNAS1\Data\2P\NM039\Processed\20260629\1"
    ];

%% Convert to ALF format and load sessions
% ALFexporter can't yet handle grating stimuli or multi-date s2p runs
[sessions, exportReport] = ALFexporter(s2pPaths, ALFroot);
loadedSessions = ALFloader(ALFroot, sessions);

% >These three field uniquely identify ONE session under ALFroot
% loadedSessions(i).subject 
% loadedSessions(i).date
% loadedSessions(i).session

% >Useful provenance paths
% loadedSessions(i).rawDataPath: path containing raw TIF and stimulus data 
% loadedSessions(i).s2pSessionPath: path of s2p output
% loadedSessions(i).files: list of all loaded ALF dataset file paths

% >Wheel data:
% loadedSessions(i).wheel.timestamps
% loadedSessions(i).wheel.position
% loadedSessions(i).wheel.velocity

% >Sparse noise stimulus data, loader hardcoded for 7x27 stimulus movie shape: 
% loadedSessions(i).sparseNoise.frameTimes
% loadedSessions(i).sparseNoise.frames
% loadedSessions(i).sparseNoise.eventTimes
% loadedSessions(i).sparseNoise.xy
% loadedSessions(i).sparseNoise.edges: hardcoded to [-135 135 -40 40]
% loadedSessions(i).sparseNoise .rawPassiveFolder: raw movie before casting to 7x27 

% >Imaging data:
% loadedSessions(i).nFOVs
% loadedSessions(i).fov(j).nFrames
% loadedSessions(i).fov(j).nROIs
% loadedSessions(i).fov(j).times
% loadedSessions(i).fov(j).F
% loadedSessions(i).fov(j).Fneu
% loadedSessions(i).fov(j).deconvolved
% loadedSessions(i).fov(j).roiTypes: The loaded iscell.npy data 
% loadedSessions(i).fov(j).stackPos: (y,x,zplane) pixel location of ROI centers
% loadedSessions(i).fov(j).masksFile: path to spatial ROI masks



%% RF mapping

opts.runningRegressor = false;
opts.useL1 = false;
opts.rf_timeLimits = [0.2 0.8]; % delay range in seconds
opts.crossFolds = 3;
opts.numShifts = 50;
opts.azimuthCropAngleDeg = -135; % in case we want to crop the design matrix
modality = 'deconv'; % 'traces' use F traces, 'deconv' use s2p deconvolved traces

% Best results use:
%NM041 06/26 deconv
%NM023 02/09 deconv
%NM023 04/28 deconv
%NM034 05/06 deconv
%HS035 05/12 deconv
%NM035 05/13 traces
%NM035 05/15 traces
%NM035 05/19 traces
%NM035 05/20 traces
%NM039 06/10 deconv
%NM044 06/10 traces
%NM040 06/16 deconv
%NM036 05/07 deconv
%NM037 06/07 traces

% right hemisphereFOVs are cropped at -azimuthCropAngleDeg and have their RF azimuth coordinates flipped in plotting
rightHemisphereFOVs = [
    "NM023/2026-02-09/2/FOV_01"
    "NM023/2026-04-28/1/FOV_01"
    "NM034/2026-05-06/1/FOV_01"
    "NM035/2026-05-13/1/FOV_00"
    "NM035/2026-05-15/2/FOV_00"
    "NM035/2026-05-19/2/FOV_00"
    "NM035/2026-05-20/1/FOV_00"
    "NM036/2026-05-07/1/FOV_01"
    "NM037/2026-06-09/1/FOV_00"
    "NM039/2026-06-10/1/FOV_00"
    "NM039/2026-06-10/1/FOV_01"
    "NM044/2026-06-10/1/FOV_00"
    "NM044/2026-06-10/1/FOV_01"
    "NM041/2026-06-26/1/FOV_01"
    "NM036/2026-06-30/1/FOV_01"
    "NM039/2026-06-29/1/FOV_01"
    ];

summaries = ALF_main_analyse_RFs(loadedSessions, modality, resultsRoot, rightHemisphereFOVs,opts);
% Results under resultsRoot/<subject>/<date>/<session>/alf/FOV_XX
% saved in both .npy and .mat
%   _FANCi_rf.withoutRunningRegressor.mat contains results and metadata:

% results.maps
%   [nROIs x nY x nX x nLags x 2]
%   Complete fitted spatiotemporal RF maps.
%   Dimension 5 contains the ON and OFF predictor maps:
%       1 = ON/white subfield
%       2 = OFF/black subfield

% results.explVars
%   [nROIs x 1]
%   RF-only cross-validated explained variance evaluated at the
%   hyperparameter combination selected using the full model.

% results.explVarsFullModel
%   [nROIs x 1]
%   Best cross-validated explained variance from the complete model:
%       RF contribution + running contribution
%   When the running regressor is disabled, this should match explVars.

% results.lambdas
%   [nROIs x 1]
%   Selected stimulus ridge-lambda grid value for each ROI.

% results.runLambdas
%   [nROIs x 1]
%   Selected running-regressor ridge-lambda grid value for each ROI.
%   Contains NaN when the running regressor is disabled.

% results.l1Lambdas
%   [nROIs x 1]
%   Selected L1 penalty for each ROI.
%   Present only when opts.useL1 is true.

% results.l1LambdaGrid
%   [1 x nL1Lambdas]
%   Complete set of tested L1 penalty values.
%   Present only when opts.useL1 is true.

% results.bestL1LambdaIndices
%   [nROIs x 1]
%   Index into l1LambdaGrid selected for each ROI.
%   Present only when opts.useL1 is true.

% results.pValues
%   [nROIs x 1]
%   Circular-shift significance p-value for each ROI.

% results.timestamps
%   [1 x nLags]
%   RF temporal lags in seconds.

% results.gaussPars
%   [nROIs x 6]
%   Two-dimensional Gaussian RF parameters:
%       [amplitude, azimuthCenter, azimuthSigma, ...
%        elevationCenter, elevationSigma, rotationAngle]
%   Spatial parameters are in visual degrees and rotationAngle is in radians.

% results.peakToNoise
%   [nROIs x 1]
%   Gaussian amplitude divided by the standard deviation of the residual
%   between the RF and its Gaussian fit.

% results.bestSubfields
%   [nROIs x 1]
%   RF representation selected for optimal-delay and Gaussian fitting:
%       1 = ON subfield
%       2 = OFF subfield
%       3 = sign-aligned average of ON and OFF
%   ROIs excluded before subfield analysis remain NaN.

% results.subfieldSigns
%   [nROIs x 2]
%   Sign of the strongest response in the ON and OFF subfields:
%       column 1 = ON sign
%       column 2 = OFF sign

% results.optimalDelays
%   [nROIs x 1]
%   One-based index of the selected temporal-lag slice for each ROI.
%   Convert to seconds using:
%       delaySeconds = results.timestamps(results.optimalDelays(iROI));
%   ROIs excluded before delay selection remain NaN.

% results.edges
%   [1 x 4]
%   Stimulus bounds actually used for fitting and interpolation:
%       [left, right, bottom, top]
%   These are the post-cropping edges.

% results.stimulusCrop
%   Scalar struct describing azimuth cropping, (needs update)

% results.useRunningRegressor
%   Logical scalar indicating whether the running regressor was enabled.

% results.running
%   Scalar struct containing running-regressor outputs

% results.predictionsFullModel
%   [nToeplitzTimes x nROIs]
%   Prediction from the complete RF-plus-running model.
%   Rows excluded during model fitting are stored as NaN.

% results.predictionsRFOnly
%   [nToeplitzTimes x nROIs]
%   Prediction generated using only the RF coefficients.
%   Rows excluded during model fitting are stored as NaN.

% results.evRFOnly
%   [nROIs x 1]
%   RF-only explained variance from the final fit using all valid model rows.
%   This is not cross-validated and is the observed statistic used in the
%   circular-shift test.

% results.evShiftRFOnly
%   [nROIs x numShifts]
%   RF-only explained variance from every circularly shifted null fit.

% results.goodRFs_peakToNoise
%   [nGoodRFs x 1]
%   One-based indices into the curated ROI ordering for ROIs satisfying:
%       results.peakToNoise > results.minPeak_shift

% results.minPeak_shift
%   Scalar peak-to-noise threshold used to define goodRFs_peakToNoise.

% results.modality
%   Character vector containing either:
%       'deconv'
%       'traces'

% results.L1
%   Logical scalar indicating whether L1/FISTA fitting was enabled.

% results.subject
%   Subject identifier.

% results.date
%   Session date.

% results.session
%   Session identifier.

% results.fov
%   FOV name, for example 'FOV_00'.

% results.roiIndexMatlab
%   [nROIs x 1]
%   One-based source ROI indices in the original FOV arrays after curation.

% results.roiIndex0
%   [nROIs x 1]
%   Zero-based source ROI indices after curation.

% results.curatedMask
%   [nOriginalROIs x 1]
%   Logical mask identifying which original FOV ROIs were retained.

% results.roiTypes
%   [nROIs x 1]
%   ROI-type values for the retained ROIs.

% results.cellClassifier
%   [nROIs x 1] or []
%   Cell-classifier values for retained ROIs when available.

% results.tracesPreprocessed
%   [nModelTimes x nROIs], stored as single
%   Activity traces actually used for RF modeling after stimulus-window
%   cropping, interpolation, per-ROI z-scoring, missing-value handling, and
%   model-row cleanup.

% results.tracesPreprocessed_times
%   [nModelTimes x 1]
%   Timestamps corresponding to tracesPreprocessed rows.

% results.tracesPreprocessed_roiIndex0
%   [nROIs x 1]
%   Zero-based source ROI indices corresponding to tracesPreprocessed columns.

% results.tracesPreprocessed_roiIndexMatlab
%   [nROIs x 1]
%   One-based source ROI indices corresponding to tracesPreprocessed columns.

% results.preprocessing
%   Scalar struct describing trace preprocessing.
%
%   For modality='deconv':
%       results.preprocessing.sourceDataset
%       results.preprocessing.description
%       results.preprocessing.useRunningRegressor
%
%   For modality='traces':
%       results.preprocessing.sourceDataset
%       results.preprocessing.neuropilCorrection
%       results.preprocessing.description
%       results.preprocessing.useRunningRegressor

% results.neuropilCorrection
%   Neuropil-correction parameters or diagnostics returned by
%   s2pUtils.estimateNeuropil_LFR.
%   Present only when modality='traces'.

% db struct contains further metadata 

%% Grouped Plotting

plotGroups = struct();
plotGroups.groupA = ["NM023", "NM034", "NM036","NM041"]; % MO
plotGroups.groupB = ["HS035","HS036","HS037","HS038"]; % V1
plotGroups.groupC = ["NM035", "NM037", "NM039", "NM044"]; %Axons
plotGroups.groupD = ["NM040","NM042","NM043"]; %V1 to MO
plotGroup = 'groupA'; % which group to plot, empty plots all four
forceReplot = false; % true overwrites existing plot outputs


RF_plots(plotsRoot, resultsRoot, plotGroup, rightHemisphereFOVs, plotGroups, [], forceReplot);






%% Helper functions


function status = RF_plots(plotsRoot, resultsRoot, plotGroup, flipFOVs, plotGroups, pythonCommand, forceReplot)
if nargin < 1 || isempty(plotsRoot)
    plotsRoot = 'D:\Pipeline\Plots';
end
if nargin < 2 || isempty(resultsRoot)
    resultsRoot = 'D:\Pipeline\Results';
end
if nargin < 3 || isempty(plotGroup)
    plotGroup = 'all';
end
if nargin < 4 || isempty(flipFOVs)
    flipFOVs = strings(0, 1);
end
usedOldSignature = false;
if nargin >= 5 && (ischar(plotGroups) || isstring(plotGroups)) && ...
        ~(isstruct(plotGroups) || isa(plotGroups, 'containers.Map'))
    if nargin >= 6
        forceReplot = pythonCommand;
    else
        forceReplot = false;
    end
    pythonCommand = plotGroups;
    plotGroups = localDefaultPlotGroups();
    usedOldSignature = true;
end
if nargin < 5 || isempty(plotGroups)
    plotGroups = localDefaultPlotGroups();
end
if nargin < 6 || isempty(pythonCommand)
    pythonCommand = 'python';
end
if (~usedOldSignature && nargin < 7) || isempty(forceReplot)
    forceReplot = false;
end

repoRoot = fileparts(mfilename('fullpath'));
scriptPath = fullfile(repoRoot, 'RF_plots.py');
if ~isfile(scriptPath)
    error('ALF_RF_pipeline:MissingPythonScript', 'Could not find RF_plots.py at %s', scriptPath);
end

flipArgs = "";
flipFOVs = string(flipFOVs(:));
flipFOVs = flipFOVs(strlength(strtrim(flipFOVs)) > 0);
for iFlip = 1:numel(flipFOVs)
    flipArgs = flipArgs + " --flip-fov " + string(localQuoteCommandArg(flipFOVs(iFlip)));
end
forceArg = "";
if forceReplot
    forceArg = " --force-replot";
end
groupSubjectsArg = " --group-subjects " + string(localQuoteCommandArg(localEncodeGroupSubjectsArg(plotGroups)));

pythonCommandText = sprintf('%s %s --results-root %s --plots-root %s --group %s%s%s%s', ...
    localQuoteCommandArg(pythonCommand), ...
    localQuoteCommandArg(scriptPath), ...
    localQuoteCommandArg(resultsRoot), ...
    localQuoteCommandArg(plotsRoot), ...
    localQuoteCommandArg(plotGroup), ...
    char(groupSubjectsArg), ...
    char(flipArgs), ...
    char(forceArg));
if ispc
    command = sprintf('set "PYTHONUNBUFFERED=1" && %s', pythonCommandText);
else
    command = sprintf('PYTHONUNBUFFERED=1 %s', pythonCommandText);
end
fprintf('[RF_plots] %s\n', command);
[status, output] = system(command, '-echo');
if status ~= 0
    if ~isempty(output)
        fprintf('%s', output);
    end
    error('ALF_RF_pipeline:PythonPlotsFailed', 'RF_plots.py failed with status %d.', status);
end
end

function text = localQuoteCommandArg(value)
text = char(string(value));
text = ['"' strrep(text, '"', '\"') '"'];
end

function plotGroups = localDefaultPlotGroups()
plotGroups = struct();
plotGroups.groupA = ["NM023", "NM034", "NM036"];
plotGroups.groupB = "HS035";
plotGroups.groupC = ["NM035", "NM037", "NM039", "NM044"];
plotGroups.groupD = "NM040";
end

function text = localEncodeGroupSubjectsArg(plotGroups)
if isa(plotGroups, 'containers.Map')
    groupNames = string(keys(plotGroups));
    getSubjects = @(name) plotGroups(char(name));
elseif isstruct(plotGroups)
    groupNames = string(fieldnames(plotGroups));
    getSubjects = @(name) plotGroups.(char(name));
else
    error('ALF_RF_pipeline:InvalidPlotGroups', ...
        'plotGroups must be a struct or containers.Map of group names to subject lists.');
end

parts = strings(numel(groupNames), 1);
for iGroup = 1:numel(groupNames)
    groupName = strtrim(groupNames(iGroup));
    subjects = string(getSubjects(groupName));
    subjects = subjects(:);
    subjects = strtrim(subjects);
    subjects = subjects(strlength(subjects) > 0);
    if strlength(groupName) == 0 || isempty(subjects)
        continue;
    end
    parts(iGroup) = groupName + "=" + strjoin(subjects, ",");
end
parts = parts(strlength(parts) > 0);
if isempty(parts)
    error('ALF_RF_pipeline:EmptyPlotGroups', 'plotGroups must contain at least one non-empty group.');
end
text = char(strjoin(parts, ";"));
end
