function main_mapRF_ALF_deconv_single(subject, date, session, alfRoot, resultsRoot, rfParams, verbose)
% main_mapRF_ALF_deconv_single Run deconvolved sparse-noise RF mapping for one session.
%
% Required inputs:
%   subject     Subject ID, e.g. 'NM023'
%   date        Session date, e.g. '2026-04-28'
%   session     Session number/string, e.g. '1'
%   alfRoot     Root containing subject/date/session ALF folders
%   resultsRoot Root where RF results are written
%
%   rfParams    RF/preprocessing configuration supplied by masterPipeline_single
%
% Optional input:
%   verbose     Logical flag for [DEBUG ...] stdout. Default is true.

if nargin < 6
    error(['Usage: main_mapRF_ALF_deconv_single(subject, date, session, ', ...
        'alfRoot, resultsRoot, rfParams[, verbose])']);
end
if nargin < 7 || isempty(verbose)
    verbose = true;
end

requiredParams = {'lambdas', 'rf_timeLimits', 'run_timeLimits', 'crossFolds', ...
    'minEV', 'maxPVal', 'minEV_shift', 'minPeak_shift', 'thresh_subfield', ...
    'numShifts', 'stimPosition', 'useRunningRegressor', ...
    'stimulusPaddingSeconds', 'maxNaNFraction', 'shiftBatchElementLimit', ...
    'npyMatlabPath', 'fovs'};
for iParam = 1:numel(requiredParams)
    if ~isfield(rfParams, requiredParams{iParam})
        error('rfParams.%s is required.', requiredParams{iParam});
    end
end

subject = char(string(subject));
date = char(string(date));
session = char(string(session));
alfRoot = char(string(alfRoot));
resultsRoot = char(string(resultsRoot));
verbose = logical(verbose);

lambdas = double(rfParams.lambdas);
rf_timeLimits = double(rfParams.rf_timeLimits);
run_timeLimits = double(rfParams.run_timeLimits);
crossFolds = double(rfParams.crossFolds);
minEV = double(rfParams.minEV);
maxPVal = double(rfParams.maxPVal);
minEV_shift = double(rfParams.minEV_shift); %#ok<NASGU> retained tutorial parameter
minPeak_shift = double(rfParams.minPeak_shift);
thresh_subfield = double(rfParams.thresh_subfield);
numShifts = double(rfParams.numShifts);
stimPosition = double(rfParams.stimPosition);
useRunningRegressor = logical(rfParams.useRunningRegressor);
stimulusPaddingSeconds = double(rfParams.stimulusPaddingSeconds);
maxNaNFraction = double(rfParams.maxNaNFraction);
shiftBatchElementLimit = double(rfParams.shiftBatchElementLimit);

validateattributes(lambdas, {'numeric'}, {'vector', 'real', 'finite', 'positive', 'nonempty'});
validateattributes(rf_timeLimits, {'numeric'}, {'vector', 'numel', 2, 'real', 'finite'});
validateattributes(run_timeLimits, {'numeric'}, {'vector', 'numel', 2, 'real', 'finite'});
validateattributes(crossFolds, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(numShifts, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(stimPosition, {'numeric'}, {'vector', 'numel', 4, 'real', 'finite'});
validateattributes(stimulusPaddingSeconds, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(maxNaNFraction, {'numeric'}, {'scalar', 'real', 'finite', '>=', 0, '<=', 1});
validateattributes(shiftBatchElementLimit, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});

global RF_VERBOSE;
RF_VERBOSE = logical(verbose);

expNumeric = str2double(session);
if isnan(expNumeric)
    exp = session;
else
    exp = expNumeric;
end
expLabel = char(string(session));

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
db.rfParams = rfParams;

% Start a session-level command-window log before any verbose/progress output.
sessionResultsPath = fullfile(resultsRoot, subject, date, session);
if ~isfolder(sessionResultsPath)
    mkdir(sessionResultsPath);
end
logFile = fullfile(sessionResultsPath, sprintf('%s_%s.log', mfilename, datestr(now, 'yyyymmdd_HHMMSS')));
diary(logFile);
diaryCleanup = onCleanup(@() diary('off')); %#ok<NASGU>

%% Set path to relevant code

code_repo = fileparts(mfilename('fullpath')); % automatically gets the folder where this script is saved
if isempty(code_repo)
    code_repo = pwd;
end
addpath(genpath(code_repo));

npyMatlabPath = char(string(rfParams.npyMatlabPath));
if isfolder(npyMatlabPath)
    addpath(genpath(npyMatlabPath));
end

fprintf('\n========== main_mapRF_ALF_deconv_single ==========' );
fprintf('\n[PROGRESS] verbose debug stdout: %d\n', verbose);
fprintf('[PROGRESS] command-window log: %s\n', logFile);
fprintf('[PROGRESS] RF parameters supplied by masterPipeline_single:\n');
disp(rfParams);
dbgRFMsg('[DEBUG] Deconvolved RF pipeline: loads mpci.ROIActivityDeconvolved.npy through alfLoaders_deconv.\n');
dbgRFMsg('[DEBUG] script path: %s\n', mfilename('fullpath'));
dbgRFMsg('[DEBUG] code_repo: %s\n', code_repo);
dbgRFMsg('[DEBUG] alfLoaders_deconv resolves to: %s\n', which('alfLoaders_deconv'));
dbgRFMsg('[DEBUG] readNPY resolves to: %s\n', which('readNPY'));
dbgRFMsg('[DEBUG] writeNPY resolves to: %s\n', which('writeNPY'));

loaders = alfLoaders_deconv();
sessionPath = loaders.sessionPath(alfRoot, subject, date, session);
sessionResultsPath = loaders.sessionPath(resultsRoot, subject, date, session);

%% Step 0. Setup, constants, and ALF loading

db.useRunningRegressor = useRunningRegressor;

% Session-level ALF data are shared by all FOV collections.
stimData = loaders.loadStimulusData(sessionPath, stimPosition);
if useRunningRegressor
    wheel = loaders.loadWheelData(sessionPath);
else
    wheel = struct('timestamps', [], 'velocity', []);
end
fovNames = loaders.listFOVs(sessionPath);
requestedFovs = cellstr(string(rfParams.fovs));
requestedFovs = requestedFovs(~cellfun(@isempty, requestedFovs));
if ~isempty(requestedFovs)
    missingFovs = setdiff(requestedFovs, fovNames);
    if ~isempty(missingFovs)
        warning('Requested FOV(s) not found in ALF session: %s', strjoin(missingFovs, ', '));
    end
    fovNames = fovNames(ismember(fovNames, requestedFovs));
end
if isempty(fovNames)
    error('No requested FOV collections were found under %s.', sessionPath);
end

fprintf('ALF session path: %s\n', sessionPath);
fprintf('Found %d FOV collection(s).\n', numel(fovNames));
dbgRFMsg('[DEBUG] results session path: %s\n', sessionResultsPath);
dbgRFMsg('[DEBUG] FOV collections: %s\n', strjoin(fovNames, ', '));
dbgRFPrintTime('_ibl_passiveRFM.times / stimData.frameTimes', stimData.frameTimes);
dbgRFPrintStimMovie('stimData.frames returned by alfLoaders.loadStimulusData', stimData.frames);
dbgRFPrintSparseEvents(stimData);
if useRunningRegressor
    dbgRFPrintTime('_ibl_wheel.timestamps / wheel.timestamps', wheel.timestamps);
    dbgRFPrintVector('_ibl_wheel.velocity / wheel.velocity', wheel.velocity);
else
    fprintf('[PROGRESS] Running regressor disabled; wheel data will not be used.\n');
end
dbgRFPrintRawOrientation(sessionPath, stimData.frameTimes, stimData.frames);

for iFOV = 1:numel(fovNames)
    fovName = fovNames{iFOV};
    db.fov = fovName;

    folder = struct();
    folder.results = fullfile(sessionResultsPath, 'alf', fovName);
    folder.plots = fullfile(sessionResultsPath, 'plots', fovName);
    if ~isfolder(folder.results)
        mkdir(folder.results);
    end
    if ~isfolder(folder.plots)
        mkdir(folder.plots);
    end

    % Load ALF calcium from one FOV collection. Do not concatenate FOVs.
    calcium = loaders.loadCalciumTraces(sessionPath, fovName);

    dbgRFMsg('\n========== DEBUG FOV %s: loaded through alfLoaders.loadCalciumTraces ==========\n', fovName);
    dbgRFMsg('[DEBUG] calcium ALF folder: %s\n', calcium.alfFolder);
    dbgRFMsg('[DEBUG] output results folder: %s\n', folder.results);
    dbgRFMsg('[DEBUG] output plots folder: %s\n', folder.plots);
    dbgRFPrintTime('mpci.times / calcium.timestamps', calcium.timestamps);
    dbgRFPrintMatrix('mpci.ROIActivityDeconvolved / calcium.deconvolved', calcium.deconvolved);
    if isfield(calcium, 'badFrames')
        dbgRFPrintLogical('mpci.badFrames / calcium.badFrames', calcium.badFrames);
    end
    if isfield(calcium, 'cellClassifier')
        dbgRFPrintVector('mpciROIs.cellClassifier / calcium.cellClassifier', calcium.cellClassifier);
    end
    if isfield(calcium, 'roiTypes')
        dbgRFPrintVector('mpciROIs.mpciROITypes / calcium.roiTypes', calcium.roiTypes);
    end

    % Keep only curated Suite2p/ALF ROIs and use the exported deconvolved
    % activity directly.  Do not load raw F, do not load neuropil F, and do
    % not run any neuropil-correction step here.
    roiTypes = calcium.roiTypes(:);
    curatedMask = roiTypes == 1;
    roiIndexMatlab = find(curatedMask);
    roiIndex0 = roiIndexMatlab - 1;

    dbgRFMsg('[DEBUG ROI] curatedMask = roiTypes == 1: %d/%d kept\n', ...
        numel(roiIndexMatlab), numel(roiTypes));
    dbgRFPrintVector('curated ROI indices MATLAB 1-based', roiIndexMatlab);
    dbgRFPrintVector('curated ROI indices Python/Suite2p 0-based', roiIndex0);

    deconvolvedCurated = calcium.deconvolved(:, curatedMask);
    dbgRFPrintMatrix('curated mpci.ROIActivityDeconvolved before high-pass', deconvolvedCurated);

    caData.time = calcium.timestamps;
    caData.traces = deconvolvedCurated;
    dbgRFPrintMatrix('caData.traces = curated deconvolved ROIs', caData.traces);
    dbgRFPrintOverlap('stim/calcium', stimData.frameTimes, caData.time);
    dbgRFPrintOverlap('stim/wheel', stimData.frameTimes, wheel.timestamps);
    dbgRFPrintOverlap('calcium/wheel', caData.time, wheel.timestamps);

    fprintf('\n=== Processing %s (%d/%d) ===\n', fovName, iFOV, numel(fovNames));
    fprintf('[PROGRESS] Using %d/%d curated ROIs for RF fitting.\n', numel(roiIndexMatlab), numel(roiTypes));

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

% a. Crop deconvolved calcium to the sparse-noise stimulus window plus configured padding.
t_ind = caData.time > t_stim(1) - stimulusPaddingSeconds & ...
    caData.time < t_stim(end) + stimulusPaddingSeconds;
caTraces = caData.traces(t_ind, :);
t_ca = caData.time(t_ind);
dbgRFMsg('[DEBUG] calcium crop to stim +/-%.6gs: kept %d/%d samples\n', stimulusPaddingSeconds, sum(t_ind), numel(t_ind));
dbgRFPrintTime('t_ca after crop', t_ca);
dbgRFPrintMatrix('deconvolved caTraces after crop', caTraces);

% b. The only trace preprocessing step: high-pass filter the deconvolved traces.
%dbgRFPrintMatrix('deconvolved caTraces after traces.highPassFilter', caTraces);

% c. Put high-pass-filtered deconvolved traces on the model timebase and
% z-score each neuron independently for regression.  No alignSampling,
% removeDecay, or stimulus-frame smoothing is applied.
tBin_ca = median(diff(t_ca));
tBin_toeplitz = median(diff(t_toeplitz));
dbgRFMsg('[DEBUG] model sampling: tBin_ca=%.12g, tBin_toeplitz=%.12g\n', tBin_ca, tBin_toeplitz);

caTraces = interp1(t_ca, caTraces, t_toeplitz);
zTraces = (caTraces - mean(caTraces, 1, 'omitnan')) ./ ...
    std(caTraces, 0, 1, 'omitnan');
modelTimes = t_toeplitz(:);
dbgRFPrintMatrix('high-pass deconvolved caTraces after interp1 onto t_toeplitz', caTraces);
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
    sum(isnan(zTraces), 1) / size(zTraces, 1) <= maxNaNFraction;
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

lamMatrix_stim = krnl.makeLambdaMatrix( ...
    [stimSize, length(rfBins)], ones(1, length(stimSize) + 1));
lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);

% Temporal smoothness regularizer for the running filter.
if useRunningRegressor
    lamMatrix_run = krnl.makeLambdaMatrix(length(runBins), 1);
else
    lamMatrix_run = zeros(0, 0);
end
nRegRows = size(lamMatrix_stim, 1) + size(lamMatrix_run, 1);
dbgRFMsg('[DEBUG] regularization: nStimPred=%d nRunPred=%d nRegRows=%d lambdas=%s\n', nStimPred, nRunPred, nRegRows, mat2str(lambdas, 6));

%% Step 5. Prepare cross-validation

nPerFold = ceil(size(stim, 1) / crossFolds);
indPerFold = reshape(1:(crossFolds*nPerFold), crossFolds, [])';
indPerFold(indPerFold > size(stim, 1)) = NaN;

% The complete model EV is used to choose the RF lambda and, when enabled,
% the running lambda. The RF-only EV is stored separately and used for RF
% inclusion, matching the manuscript statement that RF significance is judged
% without the running filter contribution.
explainedVarianceFull = NaN(size(caTraces, 2), length(lamStim), ...
    length(lamRun), crossFolds);
explainedVarianceRF = NaN(size(caTraces, 2), length(lamStim), ...
    length(lamRun), crossFolds);
fprintf('[PROGRESS] CV setup: rows=%d ROIs=%d validUnits=%d crossFolds=%d nPerFold=%d\n', size(stim,1), size(caTraces,2), sum(validUnits), crossFolds, nPerFold);

%% Step 6. Lambda selection using the complete model

for fold = 1:crossFolds
    fprintf('[PROGRESS] CV fold %d/%d starting\n', fold, crossFolds);
    ind = indPerFold(:, fold);
    ind(isnan(ind)) = [];

    j = true(size(zTraces, 1), 1);
    j(ind) = false;

    if crossFolds > 1
        y_train = gpuArray(padarray(zTraces(j, validUnits), ...
            nRegRows, 'post'));
        y_mean = mean(zTraces(j, validUnits), 1);
        xStim_train = stim(j, :);
        xRun_train = runDesign(j, :);
    else
        y_train = gpuArray(padarray(zTraces(~j, validUnits), ...
            nRegRows, 'post'));
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
            A = gpuArray([xStim_train, xRun_train; ...
                lms, zeros(size(lms, 1), nRunPred); ...
                zeros(size(lmr, 1), nStimPred), lmr]);

            B = gather(A \ y_train);
            predFull = [xStim_test, xRun_test] * B;
            predRF = xStim_test * B(1:nStimPred, :);

            explainedVarianceFull(validUnits, lamS, lamR, fold) = 1 - ...
                sum((y_test - predFull) .^ 2, 1) ./ ...
                sum((y_test - y_mean) .^ 2, 1);
            explainedVarianceRF(validUnits, lamS, lamR, fold) = 1 - ...
                sum((y_test - predRF) .^ 2, 1) ./ ...
                sum((y_test - y_mean) .^ 2, 1);
        end
    end
end

vFull = mean(explainedVarianceFull, 4);
vFull2 = reshape(vFull, size(caTraces, 2), []);
[maxEVFull, bestPair] = max(vFull2, [], 2);
[bestStimLams, bestRunLams] = ind2sub([length(lamStim), length(lamRun)], bestPair);

vRF = mean(explainedVarianceRF, 4);
vRF2 = reshape(vRF, size(caTraces, 2), []);
rowInd = (1:size(caTraces, 2))';
maxEV = vRF2(sub2ind(size(vRF2), rowInd, bestPair));

bestLambdas = lambdas(bestStimLams)';
if useRunningRegressor
    bestRunLambdas = lambdas(bestRunLams)';
else
    bestRunLambdas = NaN(size(bestRunLams));
end
dbgRFPrintVector('maxEVFull after CV lambda selection', maxEVFull);
dbgRFPrintVector('maxEV RF-only at selected lambdas', maxEV);
dbgRFMsg('[DEBUG] CV thresholds: count maxEV>minEV %.4g is %d/%d; valid maxEV finite=%d\n', minEV, sum(maxEV > minEV), numel(maxEV), sum(isfinite(maxEV)));
dbgRFPrintVector('bestLambdas stimulus', bestLambdas);
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
        ind = bestStimLams == lamS & bestRunLams == lamR & validUnits;
        if sum(ind) == 0
            continue
        end

        lmr = lamMatrix_run .* lamRun(lamR);
        A = [stim, runDesign; ...
            lms, zeros(size(lms, 1), nRunPred); ...
            zeros(size(lmr, 1), nStimPred), lmr];
        tr = padarray(zTraces(:, ind), nRegRows, 'post');

        B = gather(gpuArray(A) \ gpuArray(tr));
        rFieldsFlat(:, ind) = B(1:nStimPred, :);
        runKernels(:, ind) = B(nStimPred+1:end, :);
        modelPredictionsValid(:, ind) = [stim, runDesign] * B;
        rfPredictionsValid(:, ind) = stim * B(1:nStimPred, :);
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
% lambda pair; the EV is still computed from the fitted RF component only.
ev_shift = NaN(size(caTraces, 2), numShifts);
shifts = randi(size(zTraces, 1), numShifts, 1);

batches = max(1, ceil(size(zTraces, 1) * numShifts * numel(validUnitInds) / shiftBatchElementLimit));
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
            indNeurons = find(bestStimLams(indBatch) == lamS & ...
                bestRunLams(indBatch) == lamR & validUnits(indBatch));
            if isempty(indNeurons)
                continue
            end

            lmr = lamMatrix_run .* lamRun(lamR);
            A = gpuArray([stim, runDesign; ...
                lms, zeros(size(lms, 1), nRunPred); ...
                zeros(size(lmr, 1), nStimPred), lmr]);

            for iCell = 1:length(indNeurons)
                localCellID = indNeurons(iCell);
                tr = shiftedTraces(:, :, localCellID);
                B = gather(A \ gpuArray(padarray(tr, nRegRows, 'post')));
                pred = stim * B(1:nStimPred, :);
                ev_shift_batch(localCellID, :) = 1 - ...
                    sum((tr - pred) .^ 2, 1) ./ ...
                    sum((tr - mean(tr, 1)) .^ 2, 1);
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
results.pValues = pvals;
results.timestamps = t_rf;
results.gaussPars = rfGaussPars;
results.peakToNoise = peakNoiseRatio;
results.bestSubfields = bestSubFields;
results.subfieldSigns = subFieldSigns;
results.optimalDelays = optimalDelays;
results.edges = edges;
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
results.rfParams = rfParams;
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
results.preprocessing.sourceDataset = 'mpci.ROIActivityDeconvolved.npy';
results.preprocessing.useRunningRegressor = useRunningRegressor;
results.preprocessing.description = ['tracesPreprocessed contains curated mpci.ROIActivityDeconvolved traces after ', ...
    'crop, highPassFilter, interp1 onto modelTimes, z-scoring, and model-row cleanup. ', ...
    'No raw F loading, neuropil correction, alignSampling, removeDecay, or stimulus-frame smoothing was applied.'];
if useRunningRegressor
    rfResultFile = fullfile(folder.results, '_FANCi_rf.withRunningRegressor.mat');
else
    rfResultFile = fullfile(folder.results, '_FANCi_rf.withoutRunningRegressor.mat');
end
save(rfResultFile, 'results', 'db', '-v7.3');



fprintf('Saved RF results %s to %s\n', modelLabel, rfResultFile);
fprintf('[PROGRESS] Summary %s: ROIs=%d validUnits=%d validRF=%d goodRFs=%d maxEVmax=%.6g peakNoiseMax=%.6g\n', ...
    fovName, size(caData.traces, 2), sum(validUnits), numel(validRF), numel(goodRFs), max(maxEV, [], 'omitnan'), max(peakNoiseRatio, [], 'omitnan'));

if exist('writeNPY', 'file') == 2
    writeNPY(results.maps, fullfile(folder.results, '_FANCi_rf.maps.npy'));
    writeNPY(results.explVars, fullfile(folder.results, '_FANCi_rf.explainedVariance.npy'));
    writeNPY(results.explVarsFullModel, fullfile(folder.results, '_FANCi_rf.explainedVarianceFullModel.npy'));
    writeNPY(results.lambdas, fullfile(folder.results, '_FANCi_rf.lambdas.npy'));
    writeNPY(results.runLambdas, fullfile(folder.results, '_FANCi_rf.runLambdas.npy'));
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

%% Step 13. Plotting

fprintf(['Per-neuron RF plots disabled. Use ' ...
    'plot_RFsALFRunning_FOVOverlay.py for FOV overlays.\n']);

end

end


%% Debug-only stdout helpers. These do not load analysis data and do not replace alfLoaders_deconv.
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
