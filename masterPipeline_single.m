%% masterPipeline_single.m
% Run the RF pipeline for exactly one configured subject/date/session.
% Multiple FOVs inside that session are supported. All user-editable behavior
% is in the configuration sections immediately below.

clear;
clc;

%% ========================================================================
% A. SESSION
% ========================================================================

subject = 'NM039';
date = '2026-06-30';
session = '1';

%% ========================================================================
% B. PATHS
% ========================================================================

outputRoot = 'D:\SinglePipeline';
alfRoot = fullfile(outputRoot, 'Data', 'ALF');
deconvResultsRoot = fullfile(outputRoot, 'Results', 'Deconv');
tracesResultsRoot = fullfile(outputRoot, 'Results', 'Traces');
deconvPlotsRoot = fullfile(outputRoot, 'Plots', 'Deconv');
tracesPlotsRoot = fullfile(outputRoot, 'Plots', 'Traces');
logsRoot = fullfile(outputRoot, 'Logs');

% Raw data are read-only. The exporter also rejects output/log paths under it.
raw2pRoot = '\\10.233.25.135\FANCiNAS1\Data\2P';
suite2pRoot = '\\10.233.25.135\FANCiNAS1\Processed\2P';

% MATLAB's Python interface uses this interpreter. An executable name or an
% absolute executable path can be supplied. pythonModuleRoot=[] means the
% directory containing this script.
pythonExecutable = 'C:\Users\ichaker\AppData\Local\Python\pythoncore-3.14-64\python.exe';
pythonExecutionMode = 'OutOfProcess';
pythonModuleRoot = '';

% MATLAB NPY reader/writer dependency used by both mapper copies.
npyMatlabPath = 'C:\Users\ichaker\Desktop\RFmapping\npy-matlab\npy-matlab';

%% ========================================================================
% C. RF MAPPING PARAMETERS
% ========================================================================

rfParams = struct();
rfParams.lambdas = logspace(-4, -1, 4);
rfParams.rf_timeLimits = [0.2 1];
rfParams.run_timeLimits = [-5 5];       % seconds; running-speed filter window
rfParams.crossFolds = 3;
rfParams.minEV = 0.01;
rfParams.maxPVal = 0.05;
rfParams.minEV_shift = 0.04;            % retained tutorial parameter
rfParams.minPeak_shift = 7.7;
rfParams.thresh_subfield = 0.7;
rfParams.numShifts = 5;

% Existing scientific/preprocessing behavior, now controlled here.
rfParams.stimPosition = [-135 135 -40 40]; % [left right bottom top], degrees
rfParams.useRunningRegressor = false;
rfParams.stimulusPaddingSeconds = 10;
rfParams.maxNaNFraction = 0.1;
rfParams.shiftBatchElementLimit = 1e9;
rfParams.traceGaussianSigmaSamples = 2;

%% ========================================================================
% D. EXPORT CONFIGURATION
% ========================================================================

encoderCountsPerRevolution = 65536;
forceReexport = false;
exportDryRun = false;

%% ========================================================================
% E. PLOTTING CONFIGURATION (what is generated, not styling)
% ========================================================================

makeGaussianCoverage = true;
makePerNeuronRFMaps = true;
makeFOVOverlays = true;
makeNeutralROIOverlays = true;
makeShiftEVHistogram = true;
makePredictedVsActual = true;
makePeakTimeCourses = true;
makeRFSizeHistogram = true;
makeCombinedSessionPlots = true;
makeCombinedFOVOverlays = true;

fovs = {};                              % empty = all FOVs in this session
flipFovs = {};                          % empty = no FOV azimuth flipping
maxNeuronRFPlots = [];                  % empty = all qualifying neurons
tracePlotStride = 1;                    % 1 = keep every trace sample

%% ========================================================================
% F. EXECUTION / VERBOSITY
% ========================================================================

rfVerbose = true;

%% ========================================================================
% Resolve paths and validate configuration
% ========================================================================

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
if isempty(pythonModuleRoot)
    pythonModuleRoot = scriptDir;
end

exportHelperPy = fullfile(scriptDir, 'ONEpy', 'rf_master_export_single.py');
plotHelperPy = fullfile(scriptDir, 'RFPlot_single.py');
mapperDeconv = 'main_mapRF_ALF_deconv_single';
mapperTraces = 'main_mapRF_ALF_traces_single';

rfParams.npyMatlabPath = npyMatlabPath;
rfParams.fovs = fovs;

exportOptions = struct();
exportOptions.encoder_counts_per_revolution = encoderCountsPerRevolution;
exportOptions.sparse_noise_edges_deg = rfParams.stimPosition;
exportOptions.force_reexport = forceReexport;
exportOptions.dry_run = exportDryRun;

basePlotOptions = struct();
basePlotOptions.fovs = fovs;
basePlotOptions.flip_fovs = flipFovs;
basePlotOptions.min_ev = rfParams.minEV;
basePlotOptions.max_p_value = rfParams.maxPVal;
basePlotOptions.min_peak_to_noise = rfParams.minPeak_shift;
basePlotOptions.max_neuron_plots = maxNeuronRFPlots;
basePlotOptions.trace_plot_stride = tracePlotStride;
basePlotOptions.make_gaussian_coverage = makeGaussianCoverage;
basePlotOptions.make_per_neuron_rf_maps = makePerNeuronRFMaps;
basePlotOptions.make_fov_overlays = makeFOVOverlays;
basePlotOptions.make_neutral_roi_overlays = makeNeutralROIOverlays;
basePlotOptions.make_shift_ev_histogram = makeShiftEVHistogram;
basePlotOptions.make_predicted_vs_actual = makePredictedVsActual;
basePlotOptions.make_peak_time_courses = makePeakTimeCourses;
basePlotOptions.make_rf_size_histogram = makeRFSizeHistogram;
basePlotOptions.make_combined_session_plots = makeCombinedSessionPlots;
basePlotOptions.make_combined_fov_overlays = makeCombinedFOVOverlays;

validateSingleConfiguration(subject, date, session, raw2pRoot, suite2pRoot, outputRoot, ...
    exportHelperPy, plotHelperPy, mapperDeconv, mapperTraces, rfParams, ...
    flipFovs, tracePlotStride, encoderCountsPerRevolution);

addpath(genpath(scriptDir));
ensureDir(alfRoot);
ensureDir(deconvResultsRoot);
ensureDir(tracesResultsRoot);
ensureDir(deconvPlotsRoot);
ensureDir(tracesPlotsRoot);
ensureDir(logsRoot);

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
runLogDir = fullfile(logsRoot, ['masterPipeline_single_' timestamp]);
ensureDir(runLogDir);
pipelineLogFile = fullfile(runLogDir, ['masterPipeline_single_' timestamp '.log']);
diary(pipelineLogFile);
diaryCleanup = onCleanup(@() diary('off')); %#ok<NASGU>

try
    fprintf('\n========== masterPipeline_single ==========%s', newline);
    fprintf('Started:                  %s%s', char(datetime('now')), newline);
    fprintf('Selected session:         %s / %s / %s%s', subject, date, session, newline);
    fprintf('Output root:              %s%s', outputRoot, newline);
    fprintf('ALF root:                 %s%s', alfRoot, newline);
    fprintf('Deconv results / plots:   %s | %s%s', deconvResultsRoot, deconvPlotsRoot, newline);
    fprintf('Traces results / plots:   %s | %s%s', tracesResultsRoot, tracesPlotsRoot, newline);
    fprintf('Logs:                     %s%s', runLogDir, newline);
    fprintf('Raw 2P root (READ ONLY):  %s%s', raw2pRoot, newline);
    fprintf('Suite2p root (READ ONLY): %s%s', suite2pRoot, newline);
    fprintf('Python executable:        %s%s', pythonExecutable, newline);
    fprintf('Python module root:       %s%s', pythonModuleRoot, newline);
    fprintf('Requested FOVs:           %s%s', displayCellList(fovs, 'all'), newline);
    fprintf('Azimuth-flipped FOVs:     %s%s', displayCellList(flipFovs, 'none'), newline);
    fprintf('Mapper verbose:           %d%s', rfVerbose, newline);
    fprintf('RF parameters:%s', newline);
    disp(rfParams);

    [exporterModule, plotterModule] = configurePython( ...
        pythonExecutable, pythonExecutionMode, pythonModuleRoot);

    %% 1. Export exactly one selected session
    fprintf('\n========== STAGE 1/5: EXPORT SINGLE SESSION ==========%s', newline);
    exportResult = exporterModule.export_rf_session_single( ...
        suite2pRoot, raw2pRoot, alfRoot, subject, date, session, runLogDir, ...
        matlabStructToPyDict(exportOptions));
    fprintf('[EXPORT RETURN]%s%s%s', newline, pyObjectToJson(exportResult), newline);
    if exportDryRun
        error('exportDryRun=true: stopping before RF mapping because no export was requested.');
    end

    [alfOK, alfMessage, detectedFovs] = validateAlfSession( ...
        alfRoot, subject, date, session, fovs);
    fprintf('[VALIDATE ALF] %s%s', alfMessage, newline);
    if ~alfOK
        error('ALF validation failed: %s', alfMessage);
    end
    fprintf('[DETECTED FOVs] %s%s', strjoin(detectedFovs, ', '), newline);

    unknownFlips = setdiff(cellstr(string(flipFovs)), detectedFovs);
    if ~isempty(unknownFlips)
        warning('Configured flipFovs not found in selected session: %s', strjoin(unknownFlips, ', '));
    end

    %% 2. Deconvolved RF mapping
    fprintf('\n========== STAGE 2/5: RF MAPPING - DECONV ==========%s', newline);
    runMapperVisible(mapperDeconv, subject, date, session, alfRoot, ...
        deconvResultsRoot, rfParams, rfVerbose, pipelineLogFile);
    [deconvOK, deconvMessage] = validateRfResults( ...
        deconvResultsRoot, subject, date, session, detectedFovs);
    fprintf('[VALIDATE DECONV RESULTS] %s%s', deconvMessage, newline);
    if ~deconvOK
        error('Deconv result validation failed: %s', deconvMessage);
    end

    %% 3. Deconvolved plotting through MATLAB's Python interface
    fprintf('\n========== STAGE 3/5: PLOTTING - DECONV ==========%s', newline);
    deconvPlotOptions = basePlotOptions;
    deconvPlotOptions.analysis_label = 'Deconv';
    deconvPlotResult = plotterModule.plot_rf_session_single( ...
        deconvResultsRoot, subject, date, session, deconvPlotsRoot, ...
        alfRoot, suite2pRoot, matlabStructToPyDict(deconvPlotOptions));
    fprintf('[DECONV PLOT RETURN]%s%s%s', newline, pyObjectToJson(deconvPlotResult), newline);
    [deconvPlotsOK, deconvPlotsMessage] = validatePlotOutputs( ...
        deconvPlotsRoot, subject, date, session, detectedFovs, basePlotOptions);
    fprintf('[VALIDATE DECONV PLOTS] %s%s', deconvPlotsMessage, newline);
    if ~deconvPlotsOK
        error('Deconv plot validation failed: %s', deconvPlotsMessage);
    end

    %% 4. Raw-trace RF mapping
    fprintf('\n========== STAGE 4/5: RF MAPPING - TRACES ==========%s', newline);
    runMapperVisible(mapperTraces, subject, date, session, alfRoot, ...
        tracesResultsRoot, rfParams, rfVerbose, pipelineLogFile);
    [tracesOK, tracesMessage] = validateRfResults( ...
        tracesResultsRoot, subject, date, session, detectedFovs);
    fprintf('[VALIDATE TRACES RESULTS] %s%s', tracesMessage, newline);
    if ~tracesOK
        error('Traces result validation failed: %s', tracesMessage);
    end

    %% 5. Raw-trace plotting through MATLAB's Python interface
    fprintf('\n========== STAGE 5/5: PLOTTING - TRACES ==========%s', newline);
    tracesPlotOptions = basePlotOptions;
    tracesPlotOptions.analysis_label = 'Traces';
    tracesPlotResult = plotterModule.plot_rf_session_single( ...
        tracesResultsRoot, subject, date, session, tracesPlotsRoot, ...
        alfRoot, suite2pRoot, matlabStructToPyDict(tracesPlotOptions));
    fprintf('[TRACES PLOT RETURN]%s%s%s', newline, pyObjectToJson(tracesPlotResult), newline);
    [tracesPlotsOK, tracesPlotsMessage] = validatePlotOutputs( ...
        tracesPlotsRoot, subject, date, session, detectedFovs, basePlotOptions);
    fprintf('[VALIDATE TRACES PLOTS] %s%s', tracesPlotsMessage, newline);
    if ~tracesPlotsOK
        error('Traces plot validation failed: %s', tracesPlotsMessage);
    end

    summary = struct();
    summary.subject = subject;
    summary.date = date;
    summary.session = session;
    summary.fovs = detectedFovs;
    summary.alfRoot = alfRoot;
    summary.deconvResultsRoot = deconvResultsRoot;
    summary.deconvPlotsRoot = deconvPlotsRoot;
    summary.tracesResultsRoot = tracesResultsRoot;
    summary.tracesPlotsRoot = tracesPlotsRoot;
    summary.pipelineLog = pipelineLogFile;
    summary.completed = true;
    summary.finished = char(datetime('now'));
    summaryPath = fullfile(runLogDir, 'pipeline_summary.json');
    writeTextFile(summaryPath, jsonencode(summary, 'PrettyPrint', true));

    fprintf('\n========== SINGLE-SESSION PIPELINE COMPLETE ==========%s', newline);
    fprintf('Session:          %s / %s / %s%s', subject, date, session, newline);
    fprintf('FOVs:             %s%s', strjoin(detectedFovs, ', '), newline);
    fprintf('ALF:              %s%s', fullfile(alfRoot, subject, date, session, 'alf'), newline);
    fprintf('Deconv results:   %s%s', fullfile(deconvResultsRoot, subject, date, session), newline);
    fprintf('Deconv plots:     %s%s', fullfile(deconvPlotsRoot, subject, date, session), newline);
    fprintf('Traces results:   %s%s', fullfile(tracesResultsRoot, subject, date, session), newline);
    fprintf('Traces plots:     %s%s', fullfile(tracesPlotsRoot, subject, date, session), newline);
    fprintf('Run log:          %s%s', pipelineLogFile, newline);
    fprintf('Summary:          %s%s', summaryPath, newline);

catch ME
    try
        diary(pipelineLogFile);
    catch
    end
    fprintf('\n[PIPELINE FAILED] %s / %s / %s%s', subject, date, session, newline);
    fprintf('%s%s', exceptionText(ME), newline);
    rethrow(ME);
end

%% ========================================================================
% Local helpers (no user-editable configuration below this line)
% ========================================================================

function validateSingleConfiguration(subject, date, session, raw2pRoot, suite2pRoot, outputRoot, ...
        exportHelperPy, plotHelperPy, mapperDeconv, mapperTraces, rfParams, ...
        flipFovs, tracePlotStride, encoderCountsPerRevolution)
    if isempty(strtrim(char(subject))) || isempty(strtrim(char(session)))
        error('subject and session must be non-empty.');
    end
    if isempty(regexp(char(date), '^\d{4}-\d{2}-\d{2}$', 'once'))
        error('date must use YYYY-MM-DD format; received %s.', char(date));
    end
    if normalizedPathStartsWith(outputRoot, raw2pRoot)
        error('Output root must not be inside the read-only raw 2P root.');
    end
    if normalizedPathStartsWith(outputRoot, suite2pRoot)
        error('Output root must not be inside the read-only Suite2p root.');
    end
    if ~isfile(exportHelperPy)
        error('Missing single-session exporter: %s', exportHelperPy);
    end
    if ~isfile(plotHelperPy)
        error('Missing single-session plotter: %s', plotHelperPy);
    end
    if isempty(which(mapperDeconv)) && ~isfile(fullfile(fileparts(mfilename('fullpath')), [mapperDeconv '.m']))
        error('Missing deconv mapper: %s.m', mapperDeconv);
    end
    if isempty(which(mapperTraces)) && ~isfile(fullfile(fileparts(mfilename('fullpath')), [mapperTraces '.m']))
        error('Missing traces mapper: %s.m', mapperTraces);
    end
    if ~isfolder(raw2pRoot)
        warning('Raw 2P root is not currently reachable: %s', raw2pRoot);
    end
    if ~isfolder(suite2pRoot)
        warning('Suite2p root is not currently reachable: %s', suite2pRoot);
    end
    if ~isfolder(rfParams.npyMatlabPath)
        warning('npy-matlab path is not currently reachable: %s', rfParams.npyMatlabPath);
    end
    validateattributes(tracePlotStride, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(encoderCountsPerRevolution, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(rfParams.numShifts, {'numeric'}, {'scalar', 'integer', 'positive'});
    invalidFlips = cellstr(string(flipFovs));
    invalidFlips = invalidFlips(~cellfun(@isempty, invalidFlips) & ...
        ~startsWith(upper(string(invalidFlips)), 'FOV_'));
    if ~isempty(invalidFlips)
        error('flipFovs entries must start with FOV_: %s', strjoin(invalidFlips, ', '));
    end
end


function [exporterModule, plotterModule] = configurePython(executable, executionMode, moduleRoot)
    env = pyenv;
    if env.Status == "NotLoaded"
        fprintf('[PYTHON] Loading %s (%s).%s', executable, executionMode, newline);
        pyenv('Version', executable, 'ExecutionMode', executionMode);
        env = pyenv;
    elseif isfile(executable) && ~strcmpi(char(env.Executable), char(executable))
        error(['MATLAB Python is already loaded from %s, but pythonExecutable is %s. ', ...
            'Restart MATLAB or make the top-level setting match.'], char(env.Executable), executable);
    end
    fprintf('[PYTHON] Version: %s%s', char(env.Version), newline);
    fprintf('[PYTHON] Executable: %s%s', char(env.Executable), newline);
    fprintf('[PYTHON] Execution mode: %s%s', char(env.ExecutionMode), newline);

    pythonPath = py.sys.path();
    if int64(pythonPath.count(moduleRoot)) == 0
        insert(pythonPath, int32(0), moduleRoot);
    end
    py.importlib.invalidate_caches();
    exporterModule = py.importlib.import_module('ONEpy.rf_master_export_single');
    exporterModule = py.importlib.reload(exporterModule);
    plotterModule = py.importlib.import_module('RFPlot_single');
    plotterModule = py.importlib.reload(plotterModule);
    fprintf('[PYTHON] Imported ONEpy.rf_master_export_single and RFPlot_single.%s', newline);
end


function pyDict = matlabStructToPyDict(value)
    pyDict = py.json.loads(jsonencode(value));
end


function textValue = pyObjectToJson(value)
    try
        textValue = char(py.json.dumps(value, pyargs('indent', int32(2), 'sort_keys', true)));
    catch
        textValue = char(string(value));
    end
end


function runMapperVisible(mapperName, subject, date, session, alfRoot, ...
        resultsRoot, rfParams, verbose, pipelineLogFile)
    fprintf('[MAPPER] %s%s', mapperName, newline);
    fprintf('[MAPPER] ALF root: %s%s', alfRoot, newline);
    fprintf('[MAPPER] Results root: %s%s', resultsRoot, newline);
    resumeCleanup = onCleanup(@() resumeDiary(pipelineLogFile)); %#ok<NASGU>
    feval(mapperName, subject, date, session, alfRoot, resultsRoot, rfParams, verbose);
    fprintf('[MAPPER] %s completed.%s', mapperName, newline);
end


function resumeDiary(pathText)
    try
        diary(pathText);
    catch ME
        fprintf('[WARN] Could not resume master diary %s: %s%s', pathText, ME.message, newline);
    end
end


function [tf, message, fovNames] = validateAlfSession(alfRoot, subject, date, session, requestedFovs)
    alfSession = fullfile(alfRoot, subject, date, session, 'alf');
    sessionDir = fileparts(alfSession);
    requiredSession = {
        fullfile(alfSession, '_ibl_passiveRFM.times.npy'), ...
        fullfile(alfSession, '_ibl_sparseNoise.times.npy'), ...
        fullfile(alfSession, '_ibl_sparseNoise.xy.npy'), ...
        fullfile(alfSession, '_ibl_wheel.position.npy'), ...
        fullfile(alfSession, '_ibl_wheel.timestamps.npy'), ...
        fullfile(alfSession, '_ibl_wheel.velocity.npy'), ...
        fullfile(sessionDir, 'raw_passive_data', '_iblrig_RFMapStim.raw.bin')};
    missing = requiredSession(~cellfun(@isfile, requiredSession));

    entries = dir(fullfile(alfSession, 'FOV_*'));
    fovNames = sort({entries([entries.isdir]).name});
    requested = cellstr(string(requestedFovs));
    requested = requested(~cellfun(@isempty, requested));
    if ~isempty(requested)
        missingRequested = setdiff(requested, fovNames);
        if ~isempty(missingRequested)
            missing{end+1} = ['requested FOVs missing: ' strjoin(missingRequested, ', ')];
        end
        fovNames = fovNames(ismember(fovNames, requested));
    end

    requiredFov = {'mpci.times.npy', 'mpci.ROIActivityF.npy', ...
        'mpci.ROINeuropilActivityF.npy', 'mpci.ROIActivityDeconvolved.npy'};
    for iFov = 1:numel(fovNames)
        for iFile = 1:numel(requiredFov)
            pathText = fullfile(alfSession, fovNames{iFov}, requiredFov{iFile});
            if ~isfile(pathText)
                missing{end+1} = pathText; %#ok<AGROW>
            end
        end
    end
    if isempty(fovNames)
        missing{end+1} = ['no FOV_* collections under ' alfSession];
    end

    tf = isempty(missing);
    if tf
        message = sprintf('Complete ALF inputs for %d FOV(s): %s', numel(fovNames), strjoin(fovNames, ', '));
    else
        message = sprintf('Missing ALF input(s): %s', strjoin(missing, ' | '));
    end
end


function [tf, message] = validateRfResults(resultsRoot, subject, date, session, fovNames)
    resultSession = fullfile(resultsRoot, subject, date, session, 'alf');
    required = {'_FANCi_rf.maps.npy', '_FANCi_rf.explainedVariance.npy', ...
        '_FANCi_rf.pValues.npy', '_FANCi_rf.gaussPars.npy', ...
        '_FANCi_rf.peakToNoise.npy'};
    missing = {};
    for iFov = 1:numel(fovNames)
        for iFile = 1:numel(required)
            pathText = fullfile(resultSession, fovNames{iFov}, required{iFile});
            if ~isfile(pathText)
                missing{end+1} = pathText; %#ok<AGROW>
            end
        end
    end
    tf = isempty(missing);
    if tf
        message = sprintf('Required RF outputs exist for %d/%d FOV(s).', numel(fovNames), numel(fovNames));
    else
        message = sprintf('Missing RF result(s): %s', strjoin(missing, ' | '));
    end
end


function [tf, message] = validatePlotOutputs(plotsRoot, subject, date, session, fovNames, options)
    plotSession = fullfile(plotsRoot, subject, date, session, 'alf');
    missing = {};
    if ~isfile(fullfile(plotSession, 'plot_summary.json'))
        missing{end+1} = fullfile(plotSession, 'plot_summary.json');
    end
    anyPerFovFamily = options.make_gaussian_coverage || ...
        options.make_per_neuron_rf_maps || options.make_fov_overlays || ...
        options.make_neutral_roi_overlays || options.make_shift_ev_histogram || ...
        options.make_predicted_vs_actual || options.make_peak_time_courses || ...
        options.make_rf_size_histogram;
    for iFov = 1:numel(fovNames)
        fovDir = fullfile(plotSession, fovNames{iFov});
        if ~isfile(fullfile(fovDir, 'plot_input_summary.json'))
            missing{end+1} = fullfile(fovDir, 'plot_input_summary.json'); %#ok<AGROW>
        end
        if anyPerFovFamily && ~hasFigureRecursive(fovDir)
            missing{end+1} = [fovDir ' (no PNG/SVG figures)']; %#ok<AGROW>
        end
    end
    if options.make_combined_session_plots
        combinedDir = fullfile(plotSession, 'combined');
        if ~isfile(fullfile(combinedDir, 'combined_plot_summary.json'))
            missing{end+1} = fullfile(combinedDir, 'combined_plot_summary.json');
        end
        if ~hasFigureRecursive(combinedDir)
            missing{end+1} = [combinedDir ' (no PNG/SVG figures)'];
        end
    end
    tf = isempty(missing);
    if tf
        message = sprintf('Plot summaries and requested figures exist for %d FOV(s).', numel(fovNames));
    else
        message = sprintf('Missing plot output(s): %s', strjoin(missing, ' | '));
    end
end


function tf = hasFigureRecursive(folderPath)
    tf = false;
    if ~isfolder(folderPath)
        return;
    end
    if ~isempty(dir(fullfile(folderPath, '*.png'))) || ~isempty(dir(fullfile(folderPath, '*.svg')))
        tf = true;
        return;
    end
    entries = dir(folderPath);
    for iEntry = 1:numel(entries)
        if entries(iEntry).isdir && ~any(strcmp(entries(iEntry).name, {'.', '..'}))
            if hasFigureRecursive(fullfile(entries(iEntry).folder, entries(iEntry).name))
                tf = true;
                return;
            end
        end
    end
end


function ensureDir(pathText)
    if ~isfolder(pathText)
        mkdir(pathText);
    end
end


function writeTextFile(pathText, textValue)
    fid = fopen(pathText, 'w');
    if fid < 0
        error('Could not open file for writing: %s', pathText);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s', textValue);
end


function textValue = exceptionText(ME)
    try
        textValue = getReport(ME, 'extended', 'hyperlinks', 'off');
    catch
        textValue = sprintf('%s: %s', class(ME), ME.message);
    end
end


function textValue = displayCellList(values, emptyLabel)
    list = cellstr(string(values));
    list = list(~cellfun(@isempty, list));
    if isempty(list)
        textValue = emptyLabel;
    else
        textValue = strjoin(list, ', ');
    end
end


function tf = normalizedPathStartsWith(pathText, rootText)
    pathText = lower(strrep(char(pathText), '/', '\'));
    rootText = lower(strrep(char(rootText), '/', '\'));
    rootText = regexprep(rootText, '\\+$', '');
    tf = strcmp(pathText, rootText) || startsWith(pathText, [rootText '\']);
end
