function summaries = ALF_main_analyse_RFs(loadedSessions, modality, resultsRoot, runningRegressor, useL1)
%ALF_MAIN_ANALYSE_RFS Run sparse-noise RF mapping from ALFloader output.
%
% summaries = ALF_main_analyse_RFs(loadedSessions, modality, resultsRoot, runningRegressor, useL1)
%
% Inputs
%   loadedSessions    Struct array returned by ALFloader.m.
%   modality          'deconv', 'traces', or 'best'. Default: 'deconv'.
%   resultsRoot       Output root. Default: 'D:\Results'.
%   runningRegressor  Include wheel-speed temporal regressors. Default: false.
%   useL1             Use global L1/FISTA solver after modality is known. Default: false.
%
% Results are written as resultsRoot/subject/date/session/alf/FOV_XX with
% the same _FANCi_rf ALF names used by the deprecated RF scripts. Additional
% JSON ALF metadata files report selected modality and L1 use.

if nargin < 1 || isempty(loadedSessions)
    error('ALF_main_analyse_RFs:MissingSessions', 'loadedSessions from ALFloader.m is required.');
end
if nargin < 2 || isempty(modality)
    modality = 'deconv';
end
if nargin < 3 || isempty(resultsRoot)
    resultsRoot = 'D:\Results';
end
if nargin < 4 || isempty(runningRegressor)
    runningRegressor = false;
end
if nargin < 5 || isempty(useL1)
    useL1 = false;
end

modality = validatestring(char(string(modality)), {'deconv', 'traces', 'best'}, mfilename, 'modality');
resultsRoot = char(string(resultsRoot));
runningRegressor = logical(runningRegressor);
useL1 = logical(useL1);

localEnsurePaths();

global RF_VERBOSE;
RF_VERBOSE = true;

opts = localDefaultOptions();
opts.resultsRoot = resultsRoot;
opts.runningRegressor = runningRegressor;
opts.useL1 = useL1;
opts.requestedModality = modality;
opts.saveResults = true;
opts.verbose = true;

loadedSessions = loadedSessions(:);
summaries = repmat(localEmptySummary(), 0, 1);

fprintf('\n========== ALF_main_analyse_RFs ==========' );
fprintf('\n[PROGRESS] sessions=%d requestedModality=%s resultsRoot=%s runningRegressor=%d L1=%d\n', ...
    numel(loadedSessions), modality, resultsRoot, runningRegressor, useL1);

for iSession = 1:numel(loadedSessions)
    sessionData = loadedSessions(iSession);
    localValidateSession(sessionData);

    subject = char(string(sessionData.subject));
    date = char(string(sessionData.date));
    session = char(string(sessionData.session));
    sessionResultsPath = fullfile(resultsRoot, subject, date, session);
    if ~isfolder(sessionResultsPath)
        mkdir(sessionResultsPath);
    end

    logFile = fullfile(sessionResultsPath, sprintf('%s_%s.log', mfilename, datestr(now, 'yyyymmdd_HHMMSS')));
    diary(logFile);
    diaryCleanup = onCleanup(@() diary('off')); %#ok<NASGU>

    fprintf('\n========== RF session %d/%d: %s / %s / %s ==========' , ...
        iSession, numel(loadedSessions), subject, date, session);
    fprintf('\n[PROGRESS] command-window log: %s\n', logFile);

    if strcmp(modality, 'best')
        [selectedModality, selectionInfo] = localSelectBestModality(sessionData, opts);
    else
        selectedModality = modality;
        selectionInfo = struct('requestedModality', modality, 'selectedModality', modality, ...
            'testedFOVs', string.empty(0, 1), 'deconvGoodRFs', [], 'tracesGoodRFs', [], ...
            'rule', 'caller selected modality');
    end

    fprintf('[PROGRESS] selected modality for %s/%s/%s: %s\n', ...
        subject, date, session, selectedModality);
    if useL1
        fprintf('[PROGRESS] L1 requested;\n');
    end

    fovNames = localFovNames(sessionData);
    for iFOV = 1:numel(sessionData.fov)
        fovOpts = opts;
        fovOpts.selectedModality = selectedModality;
        fovOpts.selectionInfo = selectionInfo;
        fovOpts.logFile = logFile;
        summary = localRunOneFov(sessionData, iFOV, fovNames, selectedModality, fovOpts);
        summaries(end + 1, 1) = summary; %#ok<AGROW>
    end

    clear diaryCleanup;
end

fprintf('\n[PROGRESS] ALF_main_analyse_RFs complete. FOV results=%d\n', numel(summaries));
end

function opts = localDefaultOptions()
opts = struct();
opts.stimPosition = [-135 135 -40 40];
opts.lambdas = logspace(-4, -1, 4);
opts.l1Lambdas = logspace(-5, -2, 4);
opts.rf_timeLimits = [0.2 1];
opts.run_timeLimits = [-5 5];
opts.crossFolds = 3;
opts.minEV = 0.01;
opts.maxPVal = 0.05;
opts.minEV_shift = 0.04;
opts.minPeak_shift = 7.7;
opts.thresh_subfield = 0.7;
opts.cropAzimuthBelowEnabled = false;
opts.cropAzimuthBelowDeg = -40;
opts.numShifts = 10;
opts.resultsRoot = 'D:\Results';
opts.runningRegressor = false;
opts.useL1 = false;
opts.requestedModality = 'deconv';
opts.selectedModality = 'deconv';
opts.saveResults = true;
opts.verbose = true;
opts.selectionInfo = struct();
opts.logFile = '';
end

function summary = localEmptySummary()
summary = struct('subject', string.empty(0, 0), 'date', string.empty(0, 0), ...
    'session', string.empty(0, 0), 'fov', string.empty(0, 0), ...
    'requestedModality', string.empty(0, 0), 'selectedModality', string.empty(0, 0), ...
    'runningRegressor', false, 'L1', false, 'nROIs', NaN, 'validUnits', NaN, ...
    'validRFs', NaN, 'goodRFs', NaN, 'resultFile', string.empty(0, 0));
end

function localEnsurePaths()
codeRepo = fileparts(mfilename('fullpath'));
if isempty(codeRepo)
    codeRepo = pwd;
end
addpath(genpath(codeRepo));
npyMatlabPath = 'C:\Users\ichaker\Desktop\RFmapping\npy-matlab\npy-matlab';
if isfolder(npyMatlabPath)
    addpath(genpath(npyMatlabPath));
end
end

function localValidateSession(sessionData)
required = {'ALFroot', 'subject', 'date', 'session', 'sparseNoise', 'fov'};
for i = 1:numel(required)
    if ~isfield(sessionData, required{i})
        error('ALF_main_analyse_RFs:BadSession', 'loadedSessions entry is missing field: %s', required{i});
    end
end
if isempty(sessionData.fov)
    error('ALF_main_analyse_RFs:NoFOV', 'Session %s/%s/%s has no FOVs.', ...
        string(sessionData.subject), string(sessionData.date), string(sessionData.session));
end
if ~isfield(sessionData.sparseNoise, 'frameTimes') || isempty(sessionData.sparseNoise.frameTimes) || ...
        ~isfield(sessionData.sparseNoise, 'frames') || isempty(sessionData.sparseNoise.frames)
    error('ALF_main_analyse_RFs:MissingStimulus', 'Session %s/%s/%s has no loaded sparse-noise movie.', ...
        string(sessionData.subject), string(sessionData.date), string(sessionData.session));
end
if ~isfield(sessionData, 'wheel')
    sessionData.wheel = struct('timestamps', [], 'velocity', []);
end
end

function fovNames = localFovNames(sessionData)
fovNames = cell(numel(sessionData.fov), 1);
for i = 1:numel(sessionData.fov)
    fovNames{i} = char(string(sessionData.fov(i).name));
end
end

function [selectedModality, info] = localSelectBestModality(sessionData, opts)
selectOpts = opts;
selectOpts.useL1 = false;
selectOpts.saveResults = false;
selectOpts.selectedModality = 'best-selection';

fprintf('[PROGRESS] modality=best: testing deconv and traces without L1.\n');
fovNames = localFovNames(sessionData);
deconvCounts = NaN(numel(fovNames), 1);
tracesCounts = NaN(numel(fovNames), 1);
tested = strings(0, 1);

for iFOV = 1:numel(fovNames)
    tested(end + 1, 1) = string(fovNames{iFOV}); %#ok<AGROW>
    for iMod = 1:2
        candidate = {'deconv', 'traces'};
        candidate = candidate{iMod};
        if ~localCanUseModality(sessionData.fov(iFOV), candidate)
            fprintf('[PROGRESS] best-selection %s %s skipped: required arrays missing.\n', fovNames{iFOV}, candidate);
            continue
        end
        fprintf('[PROGRESS] best-selection testing %s on %s with L1=0.\n', candidate, fovNames{iFOV});
        try
            summary = localRunOneFov(sessionData, iFOV, fovNames, candidate, selectOpts);
            if strcmp(candidate, 'deconv')
                deconvCounts(iFOV) = summary.goodRFs;
            else
                tracesCounts(iFOV) = summary.goodRFs;
            end
        catch exc
            warning('ALF_main_analyse_RFs:BestSelectionFailed', ...
                'Best-selection %s failed for %s: %s', candidate, fovNames{iFOV}, exc.message);
        end
    end

    deconvTotal = sum(deconvCounts, 'omitnan');
    tracesTotal = sum(tracesCounts, 'omitnan');
    fprintf('[PROGRESS] best-selection after %s: deconvGood=%g tracesGood=%g\n', ...
        fovNames{iFOV}, deconvTotal, tracesTotal);
    if isfinite(deconvTotal) && isfinite(tracesTotal) && deconvTotal ~= tracesTotal
        break
    end
end

deconvTotal = sum(deconvCounts, 'omitnan');
tracesTotal = sum(tracesCounts, 'omitnan');
if tracesTotal > deconvTotal
    selectedModality = 'traces';
elseif deconvTotal > tracesTotal
    selectedModality = 'deconv';
elseif localCanUseAny(sessionData, 'deconv')
    selectedModality = 'deconv';
    fprintf('[PROGRESS] best-selection tie; defaulting to deconv.\n');
else
    selectedModality = 'traces';
    fprintf('[PROGRESS] best-selection tie with no deconv arrays; using traces.\n');
end

info = struct();
info.requestedModality = 'best';
info.selectedModality = selectedModality;
info.testedFOVs = tested;
info.deconvGoodRFs = deconvCounts;
info.tracesGoodRFs = tracesCounts;
info.rule = 'first strict good-RF count winner across full-FOV no-L1 tests; ties test more FOVs, final tie prefers deconv';
end

function tf = localCanUseAny(sessionData, modality)
tf = false;
for i = 1:numel(sessionData.fov)
    tf = tf || localCanUseModality(sessionData.fov(i), modality);
end
end

function tf = localCanUseModality(fov, modality)
switch modality
    case 'deconv'
        tf = isfield(fov, 'deconvolved') && ~isempty(fov.deconvolved);
    case 'traces'
        tf = isfield(fov, 'F') && ~isempty(fov.F) && isfield(fov, 'Fneu') && ~isempty(fov.Fneu);
    otherwise
        tf = false;
end
end

function summary = localRunOneFov(sessionData, iFOV, fovNames, selectedModality, opts)
subject = char(string(sessionData.subject));
date = char(string(sessionData.date));
session = char(string(sessionData.session));
alfRoot = char(string(sessionData.ALFroot));
resultsRoot = opts.resultsRoot;
useRunningRegressor = opts.runningRegressor;
useL1 = opts.useL1;
stimPosition = opts.stimPosition;
lambdas = opts.lambdas;
if useL1
    l1Lambdas = opts.l1Lambdas;
else
    l1Lambdas = 0;
end
rf_timeLimits = opts.rf_timeLimits;
run_timeLimits = opts.run_timeLimits;
crossFolds = opts.crossFolds;
minEV = opts.minEV;
maxPVal = opts.maxPVal;
minEV_shift = opts.minEV_shift; %#ok<NASGU>
minPeak_shift = opts.minPeak_shift;
thresh_subfield = opts.thresh_subfield;
cropAzimuthBelowEnabled = opts.cropAzimuthBelowEnabled;
cropAzimuthBelowDeg = opts.cropAzimuthBelowDeg;
numShifts = opts.numShifts;
verbose = opts.verbose;

expNumeric = str2double(session);
if isnan(expNumeric)
    exp = session;
else
    exp = expNumeric;
end
expLabel = char(string(session));

sessionPath = fullfile(alfRoot, subject, date, session);
sessionResultsPath = fullfile(resultsRoot, subject, date, session);

stimData = sessionData.sparseNoise;
if ~isfield(stimData, 'edges') || isempty(stimData.edges)
    stimData.edges = localInferEdgesFromFrames(stimData.frames, stimPosition);
end
if useRunningRegressor
    wheel = sessionData.wheel;
else
    wheel = struct('timestamps', [], 'velocity', []);
end

fov = sessionData.fov(iFOV);
fovName = char(string(fov.name));
db = struct();
db.subject = subject;
db.date = date;
db.session = session;
db.exp = exp;
db.expID = exp;
db.stim_type = 'sparse_noise';
db.stimPosition = stimPosition;
db.alfRoot = alfRoot;
db.resultsRoot = resultsRoot;
db.useRunningRegressor = useRunningRegressor;
db.selectedModality = selectedModality;
db.requestedModality = opts.requestedModality;
db.L1 = useL1;
db.cropAzimuthBelowEnabled = cropAzimuthBelowEnabled;
db.cropAzimuthBelowDeg = cropAzimuthBelowDeg;
db.fov = fovName;

folder = struct();
folder.results = fullfile(sessionResultsPath, 'alf', fovName);
folder.plots = fullfile(sessionResultsPath, 'plots', fovName);
if opts.saveResults
    if ~isfolder(folder.results)
        mkdir(folder.results);
    end
    if ~isfolder(folder.plots)
        mkdir(folder.plots);
    end
end

[calcium, caData, roiTypes, curatedMask, roiIndexMatlab, roiIndex0, preprocessing, neuropilCorrection] = ...
    localPrepareCalcium(fov, selectedModality, verbose);

fprintf('\n=== Processing %s (%d/%d), modality=%s, L1=%d ===\n', ...
    fovName, iFOV, numel(fovNames), selectedModality, useL1);
fprintf('[PROGRESS] Using %d/%d curated ROIs for RF fitting.\n', numel(roiIndexMatlab), numel(roiTypes));
dbgRFMsg('[DEBUG] session path: %s\n', sessionPath);
dbgRFMsg('[DEBUG] output results folder: %s\n', folder.results);
dbgRFPrintTime('_ibl_passiveRFM.times / stimData.frameTimes', stimData.frameTimes);
dbgRFPrintStimMovie('stimData.frames from ALFloader', stimData.frames);
dbgRFPrintTime('mpci.times / calcium.timestamps', calcium.timestamps);
dbgRFPrintMatrix('caData.traces selected modality input', caData.traces);
dbgRFPrintOverlap('stim/calcium', stimData.frameTimes, caData.time);
dbgRFPrintOverlap('stim/wheel', stimData.frameTimes, wheel.timestamps);
dbgRFPrintOverlap('calcium/wheel', caData.time, wheel.timestamps);

[cm_ON, cm_OFF] = colmaps.getRFMaps;
ellipse_x = linspace(-pi, pi, 100);
RFtypes = {'ON', 'OFF', 'ON+OFF'};

% ALF stimulus.frames is already the frame-by-frame sparse-noise movie.
t_stim = stimData.frameTimes;
tBin_stim = median(diff(t_stim));
rfBins = floor(rf_timeLimits(1) / tBin_stim) : ...
    ceil(rf_timeLimits(2) / tBin_stim);
t_rf = rfBins .* tBin_stim;
dbgRFMsg('[DEBUG] tBin_stim=%.12g, rfBins=%s, t_rf=%s\n', tBin_stim, mat2str(rfBins), mat2str(t_rf, 6));
modelLabel = 'without running regressor';
if useRunningRegressor
    modelLabel = 'with running regressor';
end
try
fprintf('ALF GLMRF %s for %s %s exp %d\n', ...
    modelLabel, subject, date, exp);
catch
    fprintf('ALF GLMRF %s for %s %s exp %s\n', ...
        modelLabel, subject, date, expLabel);

end
fprintf('Loaded %d ROIs, %d stimulus frames, %d wheel samples.\n', ...
    size(caData.traces, 2), numel(t_stim), numel(wheel.timestamps));

%% Step 1. Prepare stimulus matrix

% a. Use sparse-noise edges as exported: [left right bottom top].
edges = double(stimData.edges);


% c. Use the ALF frame-by-frame sparse-noise movie directly.
stimMatrix = double(stimData.frames);
stimulusCrop = localEmptyStimulusCrop();
if cropAzimuthBelowEnabled
    [stimMatrix, edges, stimulusCrop] = localCropStimulusAzimuthBelow( ...
        stimMatrix, edges, cropAzimuthBelowDeg);
    fprintf('[PROGRESS] Cropped stimulus azimuth < %.6g deg: kept cols %d/%d, edges=%s\n', ...
        cropAzimuthBelowDeg, stimulusCrop.nColsKept, stimulusCrop.nColsOriginal, mat2str(edges, 6));
else
    stimulusCrop.enabled = false;
    stimulusCrop.thresholdDeg = cropAzimuthBelowDeg;
    stimulusCrop.edgesOriginal = edges;
    stimulusCrop.edgesCropped = edges;
    stimulusCrop.nColsOriginal = size(stimMatrix, 3);
    stimulusCrop.nColsKept = size(stimMatrix, 3);
end
stimSize = [size(stimMatrix, 2), size(stimMatrix, 3)];
dbgRFPrintStimMovie('stimMatrix passed to whiteNoise.makeStimToeplitz', stimMatrix);
dbgRFMsg('[DEBUG] stimSize=[%d %d], edges=%s\n', stimSize(1), stimSize(2), mat2str(edges, 6));

% d. Build the stimulus Toeplitz matrix over rfBins.
[toeplitz, t_toeplitz] = whiteNoise.makeStimToeplitz( ...
    stimMatrix, t_stim, rfBins);
toeplitzRaw = toeplitz;
dbgRFPrintTime('t_toeplitz returned by whiteNoise.makeStimToeplitz', t_toeplitz);
dbgRFPrintMatrix('toeplitz returned by whiteNoise.makeStimToeplitz', toeplitz);
dbgRFMsg('[DEBUG] toeplitz nonzero fraction=%.6g (%d/%d)\n', nnz(toeplitz) / numel(toeplitz), nnz(toeplitz), numel(toeplitz));
ignoreStimTimes = false(size(t_toeplitz));

%% Step 2. Prepare running-speed design matrix on t_toeplitz

if useRunningRegressor
    % The manuscript running regressor is implemented here. Running is included
    % as a temporal filter over speeds from -5 to +5 s around each response
    % sample. Because wheel.velocity is signed, abs(velocity) is used as running
    % speed.
    runBins = floor(run_timeLimits(1) / median(diff(t_toeplitz))) : ...
        ceil(run_timeLimits(2) / median(diff(t_toeplitz)));
    t_run = runBins .* median(diff(t_toeplitz));

    wheelTime = wheel.timestamps(:);
    runSpeed = abs(wheel.velocity(:));
    [wheelTime, uniqueWheelInd] = unique(wheelTime, 'stable');
    runSpeed = runSpeed(uniqueWheelInd);
    runSpeed = interp1(wheelTime, runSpeed, t_toeplitz, 'linear', NaN);
    runMean = mean(runSpeed, 'omitnan');
    runStd = std(runSpeed, 0, 'omitnan');
    if ~isfinite(runStd) || runStd == 0
        runSpeedZ = zeros(size(runSpeed));
    else
        runSpeedZ = (runSpeed - runMean) ./ runStd;
    end

    runDesign = NaN(length(t_toeplitz), length(runBins));
    for r = 1:length(runBins)
        src = (1:length(t_toeplitz))' + runBins(r);
        ok = src >= 1 & src <= length(t_toeplitz);
        runDesign(ok, r) = runSpeedZ(src(ok));
    end
    runDesign = fillmissing(runDesign, 'constant', 0);
    dbgRFPrintVector('runSpeed interpolated/zscored on t_toeplitz', runSpeedZ);
    dbgRFPrintMatrix('runDesign after fillmissing before normalization', runDesign);
else
    runBins = [];
    t_run = [];
    runSpeedZ = [];
    runDesign = zeros(length(t_toeplitz), 0);
    dbgRFMsg('[DEBUG] running regressor disabled: runDesign is [%d x 0]\n', size(runDesign, 1));
end
runDesignRaw = runDesign;

%% Step 3. Prepare deconvolved calcium traces

% a. Crop deconvolved calcium to the sparse-noise stimulus window plus 10 s padding.
t_ind = caData.time > t_stim(1) - 10 & caData.time < t_stim(end) + 10;
caTraces = caData.traces(t_ind, :);
t_ca = caData.time(t_ind);
dbgRFMsg('[DEBUG] calcium crop to stim +/-10s: kept %d/%d samples\n', sum(t_ind), numel(t_ind));
dbgRFPrintTime('t_ca after crop', t_ca);
dbgRFPrintMatrix('deconvolved caTraces after crop', caTraces);

% b. The only trace preprocessing step: high-pass filter the deconvolved traces.
%dbgRFPrintMatrix('deconvolved caTraces after traces.highPassFilter', caTraces);

% c. Put deconvolved traces on the model timebase and
% z-score each neuron independently for regression.  No alignSampling,
% removeDecay, or stimulus-frame smoothing is applied.
tBin_ca = median(diff(t_ca));
tBin_toeplitz = median(diff(t_toeplitz));
dbgRFMsg('[DEBUG] model sampling: tBin_ca=%.12g, tBin_toeplitz=%.12g\n', tBin_ca, tBin_toeplitz);

caTraces = interp1(t_ca, caTraces, t_toeplitz);
zTraces = (caTraces - mean(caTraces, 1, 'omitnan')) ./ ...
    std(caTraces, 0, 1, 'omitnan');
modelTimes = t_toeplitz(:);
dbgRFPrintMatrix('deconvolved caTraces after interp1 onto t_toeplitz', caTraces);
dbgRFPrintMatrix('zTraces before cleanup', zTraces);

% f. Joint time/unit/stimulus/run cleanup. If the running regressor is enabled,
% running rows are removed wherever the corresponding neural/stimulus rows are removed.
validTimes = ~all(isnan(zTraces), 2);
dbgRFMsg('[DEBUG] validTimes initial: %d/%d valid rows before deleting invalid times\n', sum(validTimes), numel(validTimes));
toeplitz(~validTimes, :) = [];
runDesign(~validTimes, :) = [];
zTraces(~validTimes, :) = [];
ignoreStimTimes(~validTimes) = [];
modelTimes(~validTimes) = [];

ind = any(isnan(zTraces), 1) & ...
    sum(isnan(zTraces), 1) / size(zTraces, 1) <= 0.1;
if sum(ind) > 0
    zTraces(:, ind) = fillmissing(zTraces(:, ind), 'constant', 0);
end
validUnits = ~any(isnan(zTraces), 1)';
dbgRFMsg('[DEBUG] validUnits after NaN cleanup: %d/%d valid ROIs; filledUnitColumns=%d\n', sum(validUnits), numel(validUnits), sum(ind));

if any(ignoreStimTimes)
    stimPars = prod(stimSize);
    for b = 1:length(rfBins)
        ignore = [true(rfBins(b), 1); ignoreStimTimes(1:end-rfBins(b))];
        toeplitz(ignore, (b-1)*stimPars + (1:stimPars)) = 0;
    end
end

indTime = all(toeplitz == 0, 2);
indVal = find(validTimes);
validTimes(indVal(indTime)) = false;
zTraces(indTime, :) = [];
toeplitz(indTime, :) = [];
runDesign(indTime, :) = [];
modelTimes(indTime) = [];
dbgRFMsg('[DEBUG] removed all-zero toeplitz rows: %d; final modeling rows=%d\n', sum(indTime), size(zTraces, 1));
dbgRFPrintMatrix('toeplitz after cleanup', toeplitz);
dbgRFPrintMatrix('zTraces after cleanup', zTraces);
dbgRFPrintMatrix('runDesign after cleanup before normalization', runDesign);

% Split the signed Toeplitz matrix into ON and OFF predictors.
s = toeplitz;
s(toeplitz < 0) = 0;
stim = s;
s = toeplitz;
s(toeplitz > 0) = 0;
stim = [stim, s];
validStim = ~all(stim == 0, 1);
dbgRFMsg('[DEBUG] stim ON/OFF design before normalization: size=[%d %d], validStim=%d/%d, nnz=%d\n', size(stim,1), size(stim,2), sum(validStim), numel(validStim), nnz(stim));

% Normalize the stimulus predictors globally, as in the tutorial.
stimMean = mean(stim(:), 'omitnan');
stimStd = std(stim(:), 'omitnan');
stim = (stim - stimMean) ./ stimStd;
dbgRFMsg('[DEBUG] stim normalization: mean=%.12g std=%.12g finite=%d/%d\n', stimMean, stimStd, sum(isfinite(stim(:))), numel(stim));

% Normalize the full running design over all running predictors and samples,
% analogous to the global stimulus normalization after constructing lags.
if useRunningRegressor
    runDesignMean = mean(runDesign(:), 'omitnan');
    runDesignStd = std(runDesign(:), 'omitnan');
    if isfinite(runDesignStd) && runDesignStd > 0
        runDesign = (runDesign - runDesignMean) ./ runDesignStd;
    else
        runDesign(:) = 0;
    end
    dbgRFMsg('[DEBUG] runDesign normalization: mean=%.12g std=%.12g finite=%d/%d\n', runDesignMean, runDesignStd, sum(isfinite(runDesign(:))), numel(runDesign));
else
    runDesignMean = NaN;
    runDesignStd = NaN;
end

%% Step 4. Prepare regularization

nStimPred = size(stim, 2);
nRunPred = size(runDesign, 2);

% Separate lambda scalings are used for stimulus RFs and, optionally, the
% running filter. When running is disabled, lamRun is a dummy single value so
% the existing lambda-loop structure still fits a stimulus-only model.
lamStim = sqrt(lambdas .* sum(validTimes) .* nStimPred);
if useRunningRegressor
    lamRun = sqrt(lambdas .* sum(validTimes) .* nRunPred);
else
    lamRun = 1;
end

% Global L1 regularization is handled by localFitL1, not by augmented rows.
% The augmented design below still contains the existing smoothness-L2 rows.

lamMatrix_stim = krnl.makeLambdaMatrix( ...
    [stimSize, length(rfBins)], ones(1, length(stimSize) + 1));
lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);

% Temporal smoothness regularizer for the running filter.
if useRunningRegressor
    lamMatrix_run = krnl.makeLambdaMatrix(length(runBins), 1);
else
    lamMatrix_run = zeros(0, 0);
end
nSmoothRegRows = size(lamMatrix_stim, 1) + size(lamMatrix_run, 1);
nRegRows = nSmoothRegRows;
dbgRFMsg('[DEBUG] regularization: nStimPred=%d nRunPred=%d nSmoothRegRows=%d nRegRows=%d smoothLambdas=%s l1Lambdas=%s\n', ...
    nStimPred, nRunPred, nSmoothRegRows, nRegRows, mat2str(lambdas, 6), mat2str(l1Lambdas, 6));

%% Step 5. Prepare cross-validation

nPerFold = ceil(size(stim, 1) / crossFolds);
indPerFold = reshape(1:(crossFolds*nPerFold), crossFolds, [])';
indPerFold(indPerFold > size(stim, 1)) = NaN;

% The complete model EV is used to choose the RF lambda and, when enabled,
% the running lambda. The RF-only EV is stored separately and used for RF
% inclusion, matching the manuscript statement that RF significance is judged
% without the running filter contribution.
explainedVarianceFull = NaN(size(caTraces, 2), length(lamStim), ...
    length(lamRun), length(l1Lambdas), crossFolds);
explainedVarianceRF = NaN(size(caTraces, 2), length(lamStim), ...
    length(lamRun), length(l1Lambdas), crossFolds);
fprintf('[PROGRESS] CV setup: rows=%d ROIs=%d validUnits=%d crossFolds=%d nPerFold=%d\n', size(stim,1), size(caTraces,2), sum(validUnits), crossFolds, nPerFold);

%% Step 6. Lambda selection using the complete model

for fold = 1:crossFolds
    fprintf('[PROGRESS] CV fold %d/%d starting\n', fold, crossFolds);
    ind = indPerFold(:, fold);
    ind(isnan(ind)) = [];

    j = true(size(zTraces, 1), 1);
    j(ind) = false;

    if crossFolds > 1
        y_train = padarray(zTraces(j, validUnits), ...
            nRegRows, 'post');
        y_mean = mean(zTraces(j, validUnits), 1);
        xStim_train = stim(j, :);
        xRun_train = runDesign(j, :);
    else
        y_train = padarray(zTraces(~j, validUnits), ...
            nRegRows, 'post');
        y_mean = mean(zTraces(~j, validUnits), 1);
        xStim_train = stim(~j, :);
        xRun_train = runDesign(~j, :);
    end

    y_test = zTraces(~j, validUnits);
    xStim_test = stim(~j, :);
    xRun_test = runDesign(~j, :);

    for lamS = 1:length(lamStim)
        lms = lamMatrix_stim .* lamStim(lamS);
        for lamR = 1:length(lamRun)
            lmr = lamMatrix_run .* lamRun(lamR);
            for lamL1 = 1:length(l1Lambdas)
                A = [xStim_train, xRun_train; ...
                    lms, zeros(size(lms, 1), nRunPred); ...
                    zeros(size(lmr, 1), nStimPred), lmr];

                B = localFitModel(A, y_train, useL1, l1Lambdas(lamL1));
                predFull = [xStim_test, xRun_test] * B;
                predRF = xStim_test * B(1:nStimPred, :);

                explainedVarianceFull(validUnits, lamS, lamR, lamL1, fold) = 1 - ...
                    sum((y_test - predFull) .^ 2, 1) ./ ...
                    sum((y_test - y_mean) .^ 2, 1);
                explainedVarianceRF(validUnits, lamS, lamR, lamL1, fold) = 1 - ...
                    sum((y_test - predRF) .^ 2, 1) ./ ...
                    sum((y_test - y_mean) .^ 2, 1);
            end
        end
    end
end

vFull = mean(explainedVarianceFull, 5);
vFull2 = reshape(vFull, size(caTraces, 2), []);
[maxEVFull, bestCombo] = max(vFull2, [], 2);
[bestStimLams, bestRunLams, bestL1Lams] = ind2sub( ...
    [length(lamStim), length(lamRun), length(l1Lambdas)], bestCombo);

vRF = mean(explainedVarianceRF, 5);
vRF2 = reshape(vRF, size(caTraces, 2), []);
rowInd = (1:size(caTraces, 2))';
maxEV = vRF2(sub2ind(size(vRF2), rowInd, bestCombo));

bestLambdas = lambdas(bestStimLams)';
if useRunningRegressor
    bestRunLambdas = lambdas(bestRunLams)';
else
    bestRunLambdas = NaN(size(bestRunLams));
end
bestL1Lambdas = l1Lambdas(bestL1Lams)';
dbgRFPrintVector('maxEVFull after CV lambda selection', maxEVFull);
dbgRFPrintVector('maxEV RF-only at selected lambdas', maxEV);
dbgRFMsg('[DEBUG] CV thresholds: count maxEV>minEV %.4g is %d/%d; valid maxEV finite=%d\n', minEV, sum(maxEV > minEV), numel(maxEV), sum(isfinite(maxEV)));
dbgRFPrintVector('bestLambdas stimulus', bestLambdas);
dbgRFPrintVector('bestL1Lambdas global L1', bestL1Lambdas);
if useRunningRegressor
    dbgRFPrintVector('bestRunLambdas running', bestRunLambdas);
end

%% Step 7. Final fitting on all valid samples

rFieldsFlat = NaN(nStimPred, size(caTraces, 2));
runKernels = NaN(nRunPred, size(caTraces, 2));
modelPredictionsValid = NaN(size(stim, 1), size(caTraces, 2));
rfPredictionsValid = NaN(size(stim, 1), size(caTraces, 2));

for lamS = 1:length(lamStim)
    lms = lamMatrix_stim .* lamStim(lamS);
    for lamR = 1:length(lamRun)
        lmr = lamMatrix_run .* lamRun(lamR);
        for lamL1 = 1:length(l1Lambdas)
            ind = bestStimLams == lamS & bestRunLams == lamR & ...
                bestL1Lams == lamL1 & validUnits;
            if sum(ind) == 0
                continue
            end

            A = [stim, runDesign; ...
                lms, zeros(size(lms, 1), nRunPred); ...
                zeros(size(lmr, 1), nStimPred), lmr];
            tr = padarray(zTraces(:, ind), nRegRows, 'post');

            B = localFitModel(A, tr, useL1, l1Lambdas(lamL1));
            rFieldsFlat(:, ind) = B(1:nStimPred, :);
            runKernels(:, ind) = B(nStimPred+1:end, :);
            modelPredictionsValid(:, ind) = [stim, runDesign] * B;
            rfPredictionsValid(:, ind) = stim * B(1:nStimPred, :);
        end
    end
end

predictions = NaN(length(t_toeplitz), size(caTraces, 2));
rfPredictions = NaN(length(t_toeplitz), size(caTraces, 2));
predictions(validTimes, :) = modelPredictionsValid;
rfPredictions(validTimes, :) = rfPredictionsValid;

rFieldsFlatForPrediction = rFieldsFlat;
rFieldsFlat(~validStim, :) = NaN;
rFields = reshape(rFieldsFlat, ...
    [stimSize, length(rfBins), 2, size(caTraces, 2)]);
dbgRFMsg('[DEBUG] final fit: rFieldsFlat finite=%d/%d, runKernels finite=%d/%d, modelPredictions finite=%d/%d\n', ...
    sum(isfinite(rFieldsFlat(:))), numel(rFieldsFlat), sum(isfinite(runKernels(:))), numel(runKernels), sum(isfinite(predictions(:))), numel(predictions));

%% Step 8. RF-only shift test

% Real RF-only EV from the RF component of the complete model.
ev = NaN(size(caTraces, 2), 1);
validUnitInds = find(validUnits);
rfOnlyPred = stim * rFieldsFlatForPrediction;
ev(validUnits) = 1 - sum((zTraces(:, validUnits) - ...
    rfOnlyPred(:, validUnits)) .^ 2, 1) ./ ...
    sum((zTraces(:, validUnits) - mean(zTraces(:, validUnits), 1)) .^ 2, 1);
dbgRFPrintVector('RF-only EV before shift test / ev', ev);

% Shifted null distribution. When the running regressor is enabled, each
% shifted fit includes the same running design and the unit's selected RF/run
% lambda/L1 combination; the EV is still computed from the fitted RF component only.
ev_shift = NaN(size(caTraces, 2), numShifts);
shifts = randi(size(zTraces, 1), numShifts, 1);

batches = ceil(size(zTraces, 1) * numShifts * numel(validUnitInds) / 1000000000);
batchSize = ceil(numel(validUnitInds) / batches);
fprintf('[PROGRESS] shift test: numShifts=%d validUnitInds=%d batches=%d batchSize=%d\n', numShifts, numel(validUnitInds), batches, batchSize);

% parfor cannot classify direct writes such as ev_shift(cellID, :) because
% cellID is computed inside the loop.  Keep each batch's null EVs in a
% sliced cell output, then copy them back into ev_shift after the parfor.
ev_shift_batches = cell(batches, 1);
ev_shift_batch_ids = cell(batches, 1);

parfor b = 1:batches
    fprintf('[PROGRESS] shift batch %d/%d starting\n', b, batches);
    batchLocal = (1:batchSize) + (b-1)*batchSize;
    batchLocal(batchLocal > numel(validUnitInds)) = [];
    indBatch = validUnitInds(batchLocal);

    ev_shift_batch = NaN(numel(indBatch), numShifts);
    shiftedTraces = NaN(size(zTraces, 1), numShifts, numel(indBatch));
    for sh = 1:numShifts
        shiftedTraces(:, sh, :) = circshift(zTraces(:, indBatch), shifts(sh), 1);
    end

    for lamS = 1:length(lamStim)
        lms = lamMatrix_stim .* lamStim(lamS);
        for lamR = 1:length(lamRun)
            lmr = lamMatrix_run .* lamRun(lamR);
            for lamL1 = 1:length(l1Lambdas)
                indNeurons = find(bestStimLams(indBatch) == lamS & ...
                    bestRunLams(indBatch) == lamR & ...
                    bestL1Lams(indBatch) == lamL1 & validUnits(indBatch));
                if isempty(indNeurons)
                    continue
                end

                A = [stim, runDesign; ...
                    lms, zeros(size(lms, 1), nRunPred); ...
                    zeros(size(lmr, 1), nStimPred), lmr];

                for iCell = 1:length(indNeurons)
                    localCellID = indNeurons(iCell);
                    tr = shiftedTraces(:, :, localCellID);
                    B = localFitModel(A, padarray(tr, nRegRows, 'post'), useL1, l1Lambdas(lamL1));
                    pred = stim * B(1:nStimPred, :);
                    ev_shift_batch(localCellID, :) = 1 - ...
                        sum((tr - pred) .^ 2, 1) ./ ...
                        sum((tr - mean(tr, 1)) .^ 2, 1);
                end
            end
        end
    end

    ev_shift_batches{b} = ev_shift_batch;
    ev_shift_batch_ids{b} = indBatch;
end

for b = 1:batches
    ev_shift(ev_shift_batch_ids{b}, :) = ev_shift_batches{b};
end

pvals = sum(ev_shift > ev, 2) ./ size(ev_shift, 2);
pvals(isnan(ev)) = NaN;
validRF = find(pvals < maxPVal & maxEV > minEV);
dbgRFPrintVector('pvals from RF-only shift test', pvals);
dbgRFMsg('[DEBUG] RF inclusion: p<%.4g count=%d/%d; maxEV>%.4g count=%d/%d; validRF=%d\n', ...
    maxPVal, sum(pvals < maxPVal), numel(pvals), minEV, sum(maxEV > minEV), numel(maxEV), numel(validRF));

%% Step 9. Subfield selection

rfGaussPars = NaN(length(maxEV), 6);
peakNoiseRatio = NaN(length(maxEV), 1);
bestSubFields = NaN(length(maxEV), 1);
subFieldSigns = NaN(length(maxEV), 2);
optimalDelays = NaN(length(maxEV), 1);
selectedRF = NaN(stimSize(1), stimSize(2), length(rfBins), length(maxEV));

for cellID = validRF'
    rf = squeeze(rFields(:, :, :, :, cellID));
    rf(:, :, :, 2) = -rf(:, :, :, 2);

    rf_tmp = squeeze(mean(rf, 3));
    signs = NaN(1, 2);
    subs = NaN(1, 3);

    for sub = 1:2
        r = rf_tmp(:, :, sub);
        [subs(sub), ind] = max(abs(r), [], 'all');
        signs(sub) = sign(r(ind));
    end

    rf_tmp(:, :, 3) = (rf_tmp(:, :, 1) .* signs(1) + ...
        rf_tmp(:, :, 2) .* signs(2)) ./ 2;
    subs(3) = max(rf_tmp(:, :, 3), [], 'all');

    [m, mxSub] = max(subs);
    if subs(3) > thresh_subfield * m
        mxSub = 3;
    end
    bestSubFields(cellID) = mxSub;
    subFieldSigns(cellID, :) = signs;

    if mxSub < 3
        rf = rf(:, :, :, mxSub) .* signs(mxSub);
    else
        rf = (rf(:, :, :, 1) .* signs(1) + ...
            rf(:, :, :, 2) .* signs(2)) ./ 2;
    end

    selectedRF(:, :, :, cellID) = rf;
end

%% Step 10. Optimal delay selection

optimalRF = NaN(stimSize(1), stimSize(2), length(maxEV));
for cellID = validRF'
    rf = selectedRF(:, :, :, cellID);
    peakPerDelay = squeeze(max(rf, [], [1 2]));
    [~, mxTime] = max(peakPerDelay);
    optimalDelays(cellID) = mxTime;
    optimalRF(:, :, cellID) = rf(:, :, mxTime);
end

%% Step 11. 2D Gaussian fit

for cellID = validRF'
    rf = optimalRF(:, :, cellID);

    [rf_visDeg, xx, yy] = whiteNoise.interpolateRFtoVisualDegrees(rf, edges);
    [fitPars, rf_gauss] = whiteNoise.fit2dGaussRF( ...
        rf_visDeg, false, xx, yy);
    fitPars(6) = -fitPars(6);

    noise = std(rf_visDeg - rf_gauss, 0, 'all');
    peakNoiseRatio(cellID) = fitPars(1) / noise;
    rfGaussPars(cellID, :) = fitPars;
end

goodRFs = find(peakNoiseRatio > minPeak_shift);
fprintf('[PROGRESS] RF counts after Gaussian fit %s: validRF=%d goodRFs=%d\n', ...
    fovName, numel(validRF), numel(goodRFs));
dbgRFPrintVector('peakNoiseRatio after Gaussian fit', peakNoiseRatio);
dbgRFMsg('[DEBUG] Gaussian/good RFs: peakNoiseRatio>%.4g count=%d/%d; validRF entering gaussian=%d\n', ...
    minPeak_shift, numel(goodRFs), numel(peakNoiseRatio), numel(validRF));

%% Step 12. Save results

dims = 1:ndims(rFields);
results.maps = permute(rFields, dims([end 1:end-1]));
results.explVars = maxEV;
results.explVarsFullModel = maxEVFull;
results.lambdas = bestLambdas;
results.runLambdas = bestRunLambdas;
if useL1
    results.l1Lambdas = bestL1Lambdas;
    results.l1LambdaGrid = l1Lambdas;
    results.bestL1LambdaIndices = bestL1Lams(:);
end
results.pValues = pvals;
results.timestamps = t_rf;
results.gaussPars = rfGaussPars;
results.peakToNoise = peakNoiseRatio;
results.bestSubfields = bestSubFields;
results.subfieldSigns = subFieldSigns;
results.optimalDelays = optimalDelays;
results.edges = edges;
results.stimulusCrop = stimulusCrop;
results.useRunningRegressor = useRunningRegressor;
results.running.enabled = useRunningRegressor;
results.running.timestamps = t_toeplitz;
results.running.lags = t_run;
results.running.designRaw = runDesignRaw;
results.running.design = runDesign;
results.running.kernels = runKernels;
results.predictionsFullModel = predictions;
results.predictionsRFOnly = rfPredictions;
results.evRFOnly = ev;
results.evShiftRFOnly = ev_shift;
results.goodRFs_peakToNoise = goodRFs;
results.minPeak_shift = minPeak_shift;
results.requestedModality = opts.requestedModality;
results.selectedModality = selectedModality;
results.L1 = useL1;
results.modalitySelection = opts.selectionInfo;
results.subject = subject;
results.date = date;
results.session = session;
results.fov = fovName;
results.roiIndexMatlab = roiIndexMatlab(:);
results.roiIndex0 = roiIndex0(:);
results.curatedMask = curatedMask(:);
results.roiTypes = roiTypes(curatedMask);
if ~isempty(calcium.cellClassifier)
    results.cellClassifier = calcium.cellClassifier(curatedMask);
else
    results.cellClassifier = [];
end
results.tracesPreprocessed = single(zTraces);
results.tracesPreprocessed_times = double(modelTimes(:));
results.tracesPreprocessed_roiIndex0 = roiIndex0(:);
results.tracesPreprocessed_roiIndexMatlab = roiIndexMatlab(:);
results.preprocessing = preprocessing;
results.preprocessing.useRunningRegressor = useRunningRegressor;
if strcmp(selectedModality, 'traces')
    results.neuropilCorrection = neuropilCorrection;
end
if useRunningRegressor
    rfResultFile = fullfile(folder.results, '_FANCi_rf.withRunningRegressor.mat');
else
    rfResultFile = fullfile(folder.results, '_FANCi_rf.withoutRunningRegressor.mat');
end
metadata = struct();
metadata.object = '_FANCi_rf';
metadata.requestedModality = opts.requestedModality;
metadata.selectedModality = selectedModality;
metadata.runningRegressor = useRunningRegressor;
metadata.L1 = useL1;
metadata.modalitySelectionUsedL1 = false;
metadata.modalitySelectionRule = 'L1 is never used during modality=best dry-run comparisons.';
metadata.modalitySelection = opts.selectionInfo;
metadata.stimulusCrop = stimulusCrop;
metadata.sourceFunction = mfilename;
metadata.resultMatFile = rfResultFile;

if opts.saveResults
    save(rfResultFile, 'results', 'db', '-v7.3');
    localWriteJson(fullfile(folder.results, '_FANCi_rf.analysisMetadata.json'), metadata);
    localWriteJson(fullfile(folder.results, '_FANCi_rf.selectedModality.json'), selectedModality);
    localWriteJson(fullfile(folder.results, '_FANCi_rf.usedL1.json'), useL1);
    fprintf('Saved RF results %s to %s\n', modelLabel, rfResultFile);
else
    rfResultFile = '';
    fprintf('[PROGRESS] Dry-run result not saved for %s modality=%s.\n', fovName, selectedModality);
end
fprintf('[PROGRESS] Summary %s: ROIs=%d validUnits=%d validRF=%d goodRFs=%d maxEVmax=%.6g peakNoiseMax=%.6g\n', ...
    fovName, size(caData.traces, 2), sum(validUnits), numel(validRF), numel(goodRFs), max(maxEV, [], 'omitnan'), max(peakNoiseRatio, [], 'omitnan'));

if opts.saveResults && exist('writeNPY', 'file') == 2
    writeNPY(results.maps, fullfile(folder.results, '_FANCi_rf.maps.npy'));
    writeNPY(results.explVars, fullfile(folder.results, '_FANCi_rf.explainedVariance.npy'));
    writeNPY(results.explVarsFullModel, fullfile(folder.results, '_FANCi_rf.explainedVarianceFullModel.npy'));
    writeNPY(results.lambdas, fullfile(folder.results, '_FANCi_rf.lambdas.npy'));
    writeNPY(results.runLambdas, fullfile(folder.results, '_FANCi_rf.runLambdas.npy'));
    if useL1
        writeNPY(results.l1Lambdas, fullfile(folder.results, '_FANCi_rf.l1Lambdas.npy'));
    end
    writeNPY(results.pValues, fullfile(folder.results, '_FANCi_rf.pValues.npy'));
    writeNPY(results.timestamps, fullfile(folder.results, '_FANCi_rf.temporalLags.npy'));
    writeNPY(results.gaussPars, fullfile(folder.results, '_FANCi_rf.gaussPars.npy'));
    writeNPY(results.peakToNoise, fullfile(folder.results, '_FANCi_rf.peakToNoise.npy'));
    if useRunningRegressor
        writeNPY(results.running.lags, fullfile(folder.results, '_FANCi_rf.runningLags.npy'));
        writeNPY(results.running.kernels, fullfile(folder.results, '_FANCi_rf.runningKernels.npy'));
    end
    writeNPY(int64(results.roiIndex0), fullfile(folder.results, '_FANCi_rf.roiIndex.npy'));
    writeNPY(int64(results.roiIndexMatlab), fullfile(folder.results, '_FANCi_rf.roiIndexMatlab.npy'));
    writeNPY(double(results.cellClassifier), fullfile(folder.results, '_FANCi_rf.cellClassifier.npy'));
    writeNPY(single(zTraces), fullfile(folder.results, '_FANCi_rf.tracesPreprocessed.npy'));
    writeNPY(double(modelTimes(:)), fullfile(folder.results, '_FANCi_rf.tracesPreprocessed_times.npy'));
    writeNPY(int64(roiIndex0(:)), fullfile(folder.results, '_FANCi_rf.tracesPreprocessed_roiIndex.npy'));
    writeNPY(int64(roiIndexMatlab(:)), fullfile(folder.results, '_FANCi_rf.tracesPreprocessed_roiIndexMatlab.npy'));
end

summary = localEmptySummary();
summary.subject = string(subject);
summary.date = string(date);
summary.session = string(session);
summary.fov = string(fovName);
summary.requestedModality = string(opts.requestedModality);
summary.selectedModality = string(selectedModality);
summary.runningRegressor = useRunningRegressor;
summary.L1 = useL1;
summary.nROIs = size(caData.traces, 2);
summary.validUnits = sum(validUnits);
summary.validRFs = numel(validRF);
summary.goodRFs = numel(goodRFs);
summary.resultFile = string(rfResultFile);
end


function [calcium, caData, roiTypes, curatedMask, roiIndexMatlab, roiIndex0, preprocessing, neuropilCorrection] = localPrepareCalcium(fov, modality, verbose)
calcium = struct();
calcium.alfFolder = char(string(fov.path));
calcium.timestamps = fov.times(:);
calcium.badFrames = fov.badFrames;
calcium.cellClassifier = fov.cellClassifier;
calcium.roiTypes = fov.roiTypes(:);
roiTypes = calcium.roiTypes(:);
neuropilCorrection = [];

switch modality
    case 'deconv'
        if isempty(fov.deconvolved)
            error('ALF_main_analyse_RFs:MissingDeconv', 'FOV %s has no deconvolved activity.', string(fov.name));
        end
        calcium.deconvolved = fov.deconvolved;
        curatedMask = roiTypes == 1;
        roiIndexMatlab = find(curatedMask);
        roiIndex0 = roiIndexMatlab - 1;
        tracesCurated = calcium.deconvolved(:, curatedMask);
        caData.time = calcium.timestamps;
        caData.traces = tracesCurated;
        preprocessing = struct();
        preprocessing.sourceDataset = 'mpci.ROIActivityDeconvolved.npy';
        preprocessing.description = ['tracesPreprocessed contains curated mpci.ROIActivityDeconvolved traces after ', ...
            'crop, interp1 onto modelTimes, z-scoring, and model-row cleanup. ', ...
            'No raw F loading, neuropil correction, alignSampling, removeDecay, or stimulus-frame smoothing was applied.'];

    case 'traces'
        if isempty(fov.F) || isempty(fov.Fneu)
            error('ALF_main_analyse_RFs:MissingTraces', 'FOV %s has no F/Fneu activity.', string(fov.name));
        end
        calcium.raw = fov.F;
        calcium.neuropil = fov.Fneu;
        badRois = sum(calcium.raw, 1).' == 0;
        curatedMask = roiTypes == 1 & ~badRois;
        roiIndexMatlab = find(curatedMask);
        roiIndex0 = roiIndexMatlab - 1;
        rawCurated = calcium.raw(:, curatedMask);
        neuropilCurated = calcium.neuropil(:, curatedMask);
        sigmaSamples = 2;
        rawCurated = smoothGaussianSamples(rawCurated, sigmaSamples);
        neuropilCurated = smoothGaussianSamples(neuropilCurated, sigmaSamples);
        fprintf('[PROGRESS] Running LFR neuropil correction for %d curated ROIs after Gaussian smoothing sigma=%g samples.\n', ...
            numel(roiIndexMatlab), sigmaSamples);
        if verbose
            [neuropilCorrectedT, neuropilCorrection] = s2pUtils.estimateNeuropil_LFR(rawCurated.', neuropilCurated.');
        else
            evalc('[neuropilCorrectedT, neuropilCorrection] = s2pUtils.estimateNeuropil_LFR(rawCurated.'', neuropilCurated.'');');
        end
        caData.time = calcium.timestamps;
        caData.traces = neuropilCorrectedT.';
        preprocessing = struct();
        preprocessing.sourceDataset = 'mpci.ROIActivityF.npy + mpci.ROINeuropilActivityF.npy';
        preprocessing.neuropilCorrection = 's2pUtils.estimateNeuropil_LFR';
        preprocessing.description = ['tracesPreprocessed contains curated LFR neuropil-corrected traces after ', ...
            'Gaussian smoothing sigma=2 samples, crop, interp1 onto modelTimes, z-scoring, and model-row cleanup. ', ...
            'No alignSampling, removeDecay, or stimulus-frame smoothing was applied.'];

    otherwise
        error('ALF_main_analyse_RFs:BadModality', 'Unsupported modality: %s', modality);
end
end

function edges = localInferEdgesFromFrames(frames, stimPosition)
if ~isempty(stimPosition) && numel(stimPosition) == 4
    edges = double(stimPosition);
else
    edges = [1 size(frames, 3) 1 size(frames, 2)];
end
end

function stimulusCrop = localEmptyStimulusCrop()
stimulusCrop = struct();
stimulusCrop.enabled = false;
stimulusCrop.thresholdDeg = NaN;
stimulusCrop.edgesOriginal = [];
stimulusCrop.edgesCropped = [];
stimulusCrop.nColsOriginal = NaN;
stimulusCrop.nColsKept = NaN;
stimulusCrop.keptColumns = [];
stimulusCrop.keptAzimuthCentersDeg = [];
end

function [stimMatrix, edges, stimulusCrop] = localCropStimulusAzimuthBelow(stimMatrix, edges, thresholdDeg)
stimulusCrop = localEmptyStimulusCrop();
stimulusCrop.enabled = true;
stimulusCrop.thresholdDeg = thresholdDeg;
stimulusCrop.edgesOriginal = edges;
stimulusCrop.nColsOriginal = size(stimMatrix, 3);

if numel(edges) ~= 4
    error('ALF_main_analyse_RFs:BadStimulusEdges', ...
        'Sparse-noise edges must be [left right bottom top] before azimuth cropping.');
end
if size(stimMatrix, 3) < 1
    error('ALF_main_analyse_RFs:EmptyStimulusColumns', ...
        'Stimulus movie has no azimuth columns to crop.');
end

azimuthEdges = linspace(edges(1), edges(2), size(stimMatrix, 3) + 1);
azimuthCenters = azimuthEdges(1:end-1) + diff(azimuthEdges) ./ 2;
keepCols = azimuthCenters >= thresholdDeg;
if ~any(keepCols)
    error('ALF_main_analyse_RFs:AzimuthCropRemovedAllColumns', ...
        'Azimuth crop threshold %.6g deg removed all %d stimulus columns.', ...
        thresholdDeg, size(stimMatrix, 3));
end

keptColIdx = find(keepCols);
stimMatrix = stimMatrix(:, :, keepCols);
edges(1) = azimuthEdges(keptColIdx(1));
edges(2) = azimuthEdges(keptColIdx(end) + 1);

stimulusCrop.edgesCropped = edges;
stimulusCrop.nColsKept = size(stimMatrix, 3);
stimulusCrop.keptColumns = keptColIdx(:).';
stimulusCrop.keptAzimuthCentersDeg = azimuthCenters(keepCols);
end

function B = localFitModel(A, Y, useL1, lambdaL1)
if useL1
    B = localFitL1(A, Y, lambdaL1);
else
    B = A \ Y;
end
end

function localWriteJson(path, value)
fid = fopen(path, 'w');
if fid < 0
    error('ALF_main_analyse_RFs:JsonWriteFailed', 'Could not write %s', path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(value), 'char');
end

function tracesOut = smoothGaussianSamples(tracesIn, sigmaSamples)
% smoothGaussianSamples NaN-safe Gaussian smoothing along time.
%
% tracesIn is [time x ROI].
% sigmaSamples is in samples, not seconds.

if nargin < 2 || isempty(sigmaSamples)
    sigmaSamples = 2;
end

if isempty(tracesIn)
    tracesOut = tracesIn;
    return
end

halfWidth = ceil(4 * sigmaSamples);
kernelX = -halfWidth:halfWidth;
gaussKernel = exp(-0.5 * (kernelX ./ sigmaSamples) .^ 2);
gaussKernel = gaussKernel ./ sum(gaussKernel);

tracesIn = double(tracesIn);
finiteMask = isfinite(tracesIn);

tracesZero = tracesIn;
tracesZero(~finiteMask) = 0;

smoothNumerator = conv2(tracesZero, gaussKernel(:), 'same');
smoothDenominator = conv2(double(finiteMask), gaussKernel(:), 'same');

tracesOut = smoothNumerator ./ smoothDenominator;
tracesOut(smoothDenominator == 0) = NaN;
end

function B = localFitL1(A, Y, lambdaL1)
% localFitL1 Fit L1-regularized coefficients with the existing smoothness rows.
%
% A is already augmented with the real design rows plus smoothness-L2 rows.
% Y is already padded with matching zero targets for those smoothness rows.
% This solves, without an intercept:
%   min_B (1/(2*N)) * ||A*B - Y||_F^2 + lambdaL1 * sum(abs(B(:)))
%
% The optimizer is proximal-gradient/FISTA with soft thresholding, so no
% backslash/least-squares optimizer is used for model fitting in v3.

try
    A = double(gather(A));
catch
    A = double(A);
end
try
    Y = double(gather(Y));
catch
    Y = double(Y);
end
lambdaL1 = double(lambdaL1);

if ~isscalar(lambdaL1) || ~isfinite(lambdaL1) || lambdaL1 < 0
    error('lambdaL1 must be a finite nonnegative scalar');
end
if size(A, 1) ~= size(Y, 1)
    error('A and Y must have the same number of rows');
end

nPred = size(A, 2);
nUnits = size(Y, 2);
B = NaN(nPred, nUnits);
if nPred == 0 || nUnits == 0
    return
end

goodA = all(isfinite(A), 2);
if all(goodA) && all(isfinite(Y(:)))
    B = localFistaL1NoIntercept(A, Y, lambdaL1);
    return
end

for iUnit = 1:nUnits
    goodRows = goodA & isfinite(Y(:, iUnit));
    if sum(goodRows) == 0
        continue
    end
    B(:, iUnit) = localFistaL1NoIntercept(A(goodRows, :), Y(goodRows, iUnit), lambdaL1);
end
end


function B = localFistaL1NoIntercept(A, Y, lambdaL1)
% localFistaL1NoIntercept Proximal-gradient solver for squared loss + L1.

nRows = size(A, 1);
nPred = size(A, 2);
nUnits = size(Y, 2);
B = zeros(nPred, nUnits);
if nRows == 0 || nPred == 0 || nUnits == 0
    return
end

AtA = (A' * A) ./ nRows;
AtY = (A' * Y) ./ nRows;
L = norm(AtA, 2);
if ~isfinite(L) || L <= 0
    return
end

step = 1 ./ L;
Z = B;
t = 1;
maxIter = 500;
relTol = 1e-4;

for iter = 1:maxIter %#ok<NASGU>
    Bprev = B;
    grad = AtA * Z - AtY;
    B = localSoftThreshold(Z - step .* grad, lambdaL1 .* step);

    tNew = (1 + sqrt(1 + 4 * t^2)) / 2;
    Z = B + ((t - 1) / tNew) .* (B - Bprev);

    if norm(B(:) - Bprev(:)) <= relTol * max(1, norm(Bprev(:)))
        break
    end
    t = tNew;
end
end


function X = localSoftThreshold(X, threshold)
X = sign(X) .* max(abs(X) - threshold, 0);
end


%% Debug-only stdout helpers.
function tf = dbgRFVerbose()
global RF_VERBOSE;
if isempty(RF_VERBOSE)
    tf = true;
else
    tf = logical(RF_VERBOSE);
end
end

function dbgRFMsg(varargin)
if ~dbgRFVerbose()
    return
end
fprintf(varargin{:});
end

function dbgRFPrintTime(name, t)
if ~dbgRFVerbose()
    return
end
try
    t = double(t(:));
    finite = isfinite(t);
    nFinite = sum(finite);
    if numel(t) < 2 || nFinite == 0
        dbgRFMsg('[DEBUG TIME] %s: n=%d finite=%d/%d\n', name, numel(t), nFinite, numel(t));
        return
    end
    dt = diff(t);
    dbgRFMsg('[DEBUG TIME] %s: n=%d range=[%.12g %.12g] duration=%.12g finite=%d/%d monotonicViolations=%d\n', ...
        name, numel(t), min(t(finite)), max(t(finite)), max(t(finite))-min(t(finite)), nFinite, numel(t), sum(dt <= 0));
    dbgRFMsg('[DEBUG TIME] %s dt: median=%.12g mean=%.12g min=%.12g max=%.12g std=%.12g\n', ...
        name, median(dt, 'omitnan'), mean(dt, 'omitnan'), min(dt, [], 'omitnan'), max(dt, [], 'omitnan'), std(dt, 0, 'omitnan'));
catch ME
    dbgRFMsg('[DEBUG ERROR] dbgRFPrintTime failed for %s: %s\n', name, ME.message);
end
end

function dbgRFPrintVector(name, x)
if ~dbgRFVerbose()
    return
end
try
    if isempty(x)
        dbgRFMsg('[DEBUG VECTOR] %s: empty\n', name);
        return
    end
    x = double(x(:));
    finite = isfinite(x);
    if any(finite)
        dbgRFMsg('[DEBUG VECTOR] %s: n=%d finite=%d nan=%d min=%.12g max=%.12g mean=%.12g std=%.12g median=%.12g\n', ...
            name, numel(x), sum(finite), sum(isnan(x)), min(x(finite)), max(x(finite)), mean(x(finite)), std(x(finite)), median(x(finite)));
    else
        dbgRFMsg('[DEBUG VECTOR] %s: n=%d finite=0 nan=%d\n', name, numel(x), sum(isnan(x)));
    end
catch ME
    dbgRFMsg('[DEBUG ERROR] dbgRFPrintVector failed for %s: %s\n', name, ME.message);
end
end

function dbgRFPrintLogical(name, x)
if ~dbgRFVerbose()
    return
end
try
    if isempty(x)
        dbgRFMsg('[DEBUG LOGICAL] %s: empty\n', name);
        return
    end
    x = logical(x(:));
    dbgRFMsg('[DEBUG LOGICAL] %s: n=%d true=%d false=%d\n', name, numel(x), sum(x), sum(~x));
catch ME
    dbgRFMsg('[DEBUG ERROR] dbgRFPrintLogical failed for %s: %s\n', name, ME.message);
end
end

function dbgRFPrintMatrix(name, x)
if ~dbgRFVerbose()
    return
end
try
    if isempty(x)
        dbgRFMsg('[DEBUG MATRIX] %s: empty\n', name);
        return
    end
    sz = size(x);
    x = double(x(:));
    finite = isfinite(x);
    if any(finite)
        dbgRFMsg('[DEBUG MATRIX] %s: size=%s finite=%d/%d nan=%d min=%.12g max=%.12g mean=%.12g std=%.12g median=%.12g\n', ...
            name, mat2str(sz), sum(finite), numel(x), sum(isnan(x)), min(x(finite)), max(x(finite)), mean(x(finite)), std(x(finite)), median(x(finite)));
    else
        dbgRFMsg('[DEBUG MATRIX] %s: size=%s finite=0/%d nan=%d\n', name, mat2str(sz), numel(x), sum(isnan(x)));
    end
catch ME
    dbgRFMsg('[DEBUG ERROR] dbgRFPrintMatrix failed for %s: %s\n', name, ME.message);
end
end

function dbgRFPrintStimMovie(name, movie)
if ~dbgRFVerbose()
    return
end
try
    sz = size(movie);
    vals = unique(movie(:)).';
    dbgRFMsg('[DEBUG STIM] %s: size=%s unique=%s\n', name, mat2str(sz), mat2str(double(vals)));
    n = numel(movie);
    for v = vals
        c = sum(movie(:) == v);
        dbgRFMsg('[DEBUG STIM]   value %.12g: %d/%d (%.6f)\n', double(v), c, n, c / n);
    end
    if ndims(movie) == 3 && size(movie, 1) > 1
        frameFlat = reshape(movie, size(movie, 1), []);
        changed = any(diff(frameFlat, 1, 1) ~= 0, 2);
        nnzPerFrame = sum(frameFlat ~= 0, 2);
        dbgRFMsg('[DEBUG STIM]   frames with any change from previous: %d/%d\n', sum(changed), numel(changed));
        dbgRFMsg('[DEBUG STIM]   nonzero squares/frame: median=%.12g min=%.12g max=%.12g mean=%.12g\n', ...
            median(nnzPerFrame), min(nnzPerFrame), max(nnzPerFrame), mean(nnzPerFrame));
    end
catch ME
    dbgRFMsg('[DEBUG ERROR] dbgRFPrintStimMovie failed for %s: %s\n', name, ME.message);
end
end

function dbgRFPrintSparseEvents(stimData)
if ~dbgRFVerbose()
    return
end
try
    if isfield(stimData, 'eventTimes')
        dbgRFPrintTime('_ibl_sparseNoise.times / stimData.eventTimes', stimData.eventTimes);
    end
    if isfield(stimData, 'xy') && ~isempty(stimData.xy)
        xy = double(stimData.xy);
        dbgRFMsg('[DEBUG SPARSE EVENTS] stimData.xy size=%s xRange=[%.12g %.12g] yRange=[%.12g %.12g] uniqueX=%d uniqueY=%d\n', ...
            mat2str(size(xy)), min(xy(:,1)), max(xy(:,1)), min(xy(:,2)), max(xy(:,2)), numel(unique(xy(:,1))), numel(unique(xy(:,2))));
    end
    if isfield(stimData, 'frames')
        dbgRFMsg('[DEBUG SPARSE EVENTS] nnz(stimData.frames)=%d\n', nnz(stimData.frames));
    end
catch ME
    dbgRFMsg('[DEBUG ERROR] dbgRFPrintSparseEvents failed: %s\n', ME.message);
end
end

function dbgRFPrintOverlap(name, a, b)
if ~dbgRFVerbose()
    return
end
try
    a = double(a(:)); b = double(b(:));
    a = a(isfinite(a)); b = b(isfinite(b));
    if isempty(a) || isempty(b)
        dbgRFMsg('[DEBUG OVERLAP] %s: empty finite input\n', name);
        return
    end
    lo = max(min(a), min(b));
    hi = min(max(a), max(b));
    dbgRFMsg('[DEBUG OVERLAP] %s: overlap=[%.12g %.12g] duration=%.12g; aRange=[%.12g %.12g]; bRange=[%.12g %.12g]\n', ...
        name, lo, hi, hi-lo, min(a), max(a), min(b), max(b));
catch ME
    dbgRFMsg('[DEBUG ERROR] dbgRFPrintOverlap failed for %s: %s\n', name, ME.message);
end
end

function dbgRFPrintRawOrientation(sessionPath, frameTimes, currentMovie)
if ~dbgRFVerbose()
    return
end
try
    rawFile = fullfile(sessionPath, 'raw_passive_data', '_iblrig_RFMapStim.raw.bin');
    if ~isfile(rawFile)
        dbgRFMsg('[DEBUG RAW] raw sparse movie not found: %s\n', rawFile);
        return
    end
    fid = fopen(rawFile, 'rb');
    if fid < 0
        dbgRFMsg('[DEBUG RAW] could not open raw sparse movie: %s\n', rawFile);
        return
    end
    raw = fread(fid, inf, 'uint8=>uint8');
    fclose(fid);
    nFrames = numel(frameTimes);
    shapes = [10 30; 7 27];
    nRows = NaN; nCols = NaN;
    for iShape = 1:size(shapes, 1)
        if numel(raw) == nFrames * shapes(iShape, 1) * shapes(iShape, 2)
            nRows = shapes(iShape, 1);
            nCols = shapes(iShape, 2);
            break
        end
    end
    rawSigned = double(raw);
    rawSigned(rawSigned >= 128) = rawSigned(rawSigned >= 128) - 256;
    dbgRFMsg('[DEBUG RAW] %s\n', rawFile);
    dbgRFMsg('[DEBUG RAW] nBytes=%d inferredFrames=%d inferredRows=%g inferredCols=%g uniqueUint8=%s uniqueInt8=%s\n', ...
        numel(raw), nFrames, nRows, nCols, mat2str(double(unique(raw(:)).')), mat2str(unique(rawSigned(:)).'));
    if isnan(nRows)
        return
    end
    decoded = dbgRFDecodeSparseNoise(rawSigned);
    candidateMatlab = reshape(decoded, [nRows, nCols, nFrames]);
    candidateMatlab = permute(candidateMatlab, [3 1 2]);
    candidatePython = reshape(decoded, [nFrames, nCols, nRows]);
    candidatePython = permute(candidatePython, [1 3 2]);
    currentMovie = int8(currentMovie);
    if isequal(size(currentMovie), size(candidateMatlab))
        matchMatlab = mean(currentMovie(:) == int8(candidateMatlab(:)));
        matchPython = mean(currentMovie(:) == int8(candidatePython(:)));
        candidateMatch = mean(int8(candidateMatlab(:)) == int8(candidatePython(:)));
        dbgRFMsg('[DEBUG RAW ORIENTATION] current stimData.frames match old MATLAB reshape candidate: %.6f\n', matchMatlab);
        dbgRFMsg('[DEBUG RAW ORIENTATION] current stimData.frames match Python/C-order candidate: %.6f\n', matchPython);
        dbgRFMsg('[DEBUG RAW ORIENTATION] old MATLAB candidate vs Python/C-order candidate equality: %.6f\n', candidateMatch);
    else
        dbgRFMsg('[DEBUG RAW ORIENTATION] currentMovie size=%s candidate size=%s; cannot compare\n', mat2str(size(currentMovie)), mat2str(size(candidateMatlab)));
    end
catch ME
    dbgRFMsg('[DEBUG ERROR] dbgRFPrintRawOrientation failed: %s\n', ME.message);
end
end

function decoded = dbgRFDecodeSparseNoise(rawSigned)
values = unique(rawSigned(:)).';
if all(ismember(values, [-1 0 1])) && ~ismember(-128, values)
    decoded = int8(rawSigned);
    return
end
if all(ismember(values, [-128 -1 0]))
    decoded = zeros(size(rawSigned), 'int8');
    decoded(rawSigned == -128) = 0;
    decoded(rawSigned == -1) = -1;
    decoded(rawSigned == 0) = 1;
    return
end
if all(ismember(values, [-128 0 127]))
    decoded = zeros(size(rawSigned), 'int8');
    decoded(rawSigned == -128) = -1;
    decoded(rawSigned == 0) = 0;
    decoded(rawSigned == 127) = 1;
    return
end
error('Unexpected sparse-noise raw values: %s', mat2str(values));
end
