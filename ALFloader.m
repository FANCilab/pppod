function ALF = ALFloader(ALFroot, sessions)
%ALFLOADER Load ALF/ONE data as a simple session struct array.
%
% Usage
%   ALF = ALFloader(ALFroot)
%   ALF = ALFloader(ALFroot, sessions)
%
% Inputs
%   ALFroot  Root folder organized as subject/date/session/alf/...
%   sessions Optional struct array with string scalar fields subject, date,
%            and session. Legacy tables with variables subject, date, and
%            session are also accepted. When omitted, all valid sessions under
%            ALFroot are loaded.
%
% Output
%   ALF(i) is one session. If sessions is provided, ALF(i) matches
%   sessions(i). No global ID tables and no calcium block table are created.
%
% Common access
%   ALF(i).ALFroot
%   ALF(i).rawDataPath
%   ALF(i).s2pSessionPath
%   ALF(i).wheel.timestamps
%   ALF(i).sparseNoise.frameTimes
%   ALF(i).fov(j).times
%   ALF(i).fov(j).F
%   ALF(i).fov(j).Fneu
%   ALF(i).fov(j).deconvolved

if nargin < 1 || isempty(ALFroot)
    error('ALFloader:MissingRoot', 'ALFroot is required.');
end

ALFroot = char(ALFroot);
if ~isfolder(ALFroot)
    error('ALFloader:BadRoot', 'ALFroot does not exist or is not a folder: %s', ALFroot);
end

readNPYPath = localEnsureReadNPYOnPath();
if readNPYPath ~= ""
    fprintf('[ALFloader] Added readNPY path: %s\n', readNPYPath);
end

if nargin < 2 || (isnumeric(sessions) && isempty(sessions))
    sessions = localDiscoverSessions(ALFroot);
    fprintf('[ALFloader] Discovered %d session(s) under %s\n', numel(sessions), ALFroot);
else
    sessions = localNormalizeSessions(sessions);
    fprintf('[ALFloader] Loading %d provided session(s) from %s\n', numel(sessions), ALFroot);
end

nSessions = numel(sessions);
ALF = repmat(localEmptySession(), nSessions, 1);

for iSession = 1:nSessions
    subject = sessions(iSession).subject;
    dateValue = sessions(iSession).date;
    sessionValue = sessions(iSession).session;

    sessionPath = fullfile(ALFroot, char(subject), char(dateValue), char(sessionValue));
    alfPath = fullfile(sessionPath, 'alf');

    ALF(iSession).ALFroot = string(ALFroot);
    ALF(iSession).subject = subject;
    ALF(iSession).date = dateValue;
    ALF(iSession).session = sessionValue;
    ALF(iSession).exists = isfolder(sessionPath);
    ALF(iSession).hasAlf = isfolder(alfPath);

    fprintf('[ALFloader] Session %d/%d: %s / %s / %s\n', ...
        iSession, nSessions, subject, dateValue, sessionValue);

    if ~isfolder(sessionPath)
        msg = sprintf('Session folder not found: %s', sessionPath);
        ALF(iSession).errorMessage = string(msg);
        ALF(iSession).errors = localAddError(ALF(iSession).errors, 'session', sessionPath, msg);
        fprintf('  [missing] %s\n', msg);
        continue
    end

    if ~isfolder(alfPath)
        msg = sprintf('ALF folder not found: %s', alfPath);
        ALF(iSession).errorMessage = string(msg);
        ALF(iSession).errors = localAddError(ALF(iSession).errors, 'session', alfPath, msg);
        fprintf('  [missing] %s\n', msg);
        continue
    end

    ALF(iSession).files = localListFiles(alfPath);
    sourcePaths = localLoadFanciSourcePaths(alfPath);
    ALF(iSession).rawDataPath = sourcePaths.rawDataPath;
    ALF(iSession).s2pSessionPath = sourcePaths.s2pSessionPath;

    [ALF(iSession).sparseNoise, stimErrors] = localLoadSparseNoise(sessionPath);
    ALF(iSession).errors = localConcatErrors(ALF(iSession).errors, stimErrors);
    ALF(iSession).hasSparseNoiseStimulus = ALF(iSession).sparseNoise.loaded;

    [ALF(iSession).wheel, wheelErrors] = localLoadWheel(sessionPath);
    ALF(iSession).errors = localConcatErrors(ALF(iSession).errors, wheelErrors);
    ALF(iSession).hasWheel = ALF(iSession).wheel.loaded;

    fovNames = localListFOVs(sessionPath);
    ALF(iSession).fov = repmat(localEmptyFov(), numel(fovNames), 1);

    if isempty(fovNames)
        fprintf('  [calcium] no FOV_* folders found\n');
    end

    for iFov = 1:numel(fovNames)
        fovName = fovNames{iFov};
        fovFolder = fullfile(alfPath, fovName);
        fprintf('  [FOV %d/%d] %s\n', iFov, numel(fovNames), fovName);

        [ALF(iSession).fov(iFov), fovErrors] = localLoadFov(fovName, fovFolder);
        ALF(iSession).errors = localConcatErrors(ALF(iSession).errors, fovErrors);
        ALF(iSession).files = [ALF(iSession).files; localListFiles(fovFolder)]; %#ok<AGROW>

        fprintf('    loaded=%d, frames=%d, rois=%d\n', ...
            ALF(iSession).fov(iFov).loaded, ...
            ALF(iSession).fov(iFov).nFrames, ...
            ALF(iSession).fov(iFov).nROIs);
    end

    ALF(iSession).nFOVs = numel(ALF(iSession).fov);
    ALF(iSession).nROIs = sum([ALF(iSession).fov.nROIs]);
    ALF(iSession).hasCalcium = any([ALF(iSession).fov.loaded]);
    ALF(iSession).loaded = true;
end

fprintf('[ALFloader] Done. Sessions=%d\n', numel(ALF));
end

function session = localEmptySession()
session = struct( ...
    'ALFroot', string.empty(0, 0), ...
    'rawDataPath', string.empty(0, 0), ...
    's2pSessionPath', string.empty(0, 0), ...
    'subject', string.empty(0, 0), ...
    'date', string.empty(0, 0), ...
    'session', string.empty(0, 0), ...
    'exists', false, ...
    'hasAlf', false, ...
    'loaded', false, ...
    'hasSparseNoiseStimulus', false, ...
    'hasWheel', false, ...
    'hasCalcium', false, ...
    'nFOVs', 0, ...
    'nROIs', 0, ...
    'wheel', localEmptyWheel(), ...
    'sparseNoise', localEmptySparseNoise(), ...
    'fov', repmat(localEmptyFov(), 0, 1), ...
    'files', strings(0, 1), ...
    'errors', localEmptyErrors(), ...
    'errorMessage', string.empty(0, 0));
end

function fov = localEmptyFov()
fov = struct( ...
    'name', string.empty(0, 0), ...
    'path', string.empty(0, 0), ...
    'loaded', false, ...
    'nFrames', 0, ...
    'nROIs', 0, ...
    'times', [], ...
    'badFrames', [], ...
    'F', [], ...
    'Fneu', [], ...
    'deconvolved', [], ...
    'cellClassifier', [], ...
    'roiTypes', [], ...
    'stackPos', [], ...
    'meanImage', [], ...
    'masksFile', string.empty(0, 0), ...
    'files', strings(0, 1), ...
    'missing', strings(0, 1), ...
    'errorMessage', string.empty(0, 0));
end

function sparseNoise = localEmptySparseNoise()
sparseNoise = struct( ...
    'loaded', false, ...
    'frameTimes', [], ...
    'frames', [], ...
    'eventTimes', [], ...
    'xy', [], ...
    'edges', [], ...
    'rawPassiveFolder', string.empty(0, 0), ...
    'files', strings(0, 1), ...
    'missing', strings(0, 1), ...
    'errorMessage', string.empty(0, 0));
end

function wheel = localEmptyWheel()
wheel = struct( ...
    'loaded', false, ...
    'timestamps', [], ...
    'position', [], ...
    'velocity', [], ...
    'files', strings(0, 1), ...
    'missing', strings(0, 1), ...
    'errorMessage', string.empty(0, 0));
end

function errors = localEmptyErrors()
errors = struct('dataset', {}, 'path', {}, 'message', {});
end

function addedPath = localEnsureReadNPYOnPath()
addedPath = "";
if exist('readNPY', 'file') == 2
    return
end

thisFolder = fileparts(mfilename('fullpath'));
repoParent = fileparts(thisFolder);
candidates = unique([
    string(fullfile(repoParent, 'npy-matlab', 'npy-matlab'))
    "C:\Users\ichaker\Desktop\FANCiRepos\npy-matlab\npy-matlab"
    ], 'stable');

for i = 1:numel(candidates)
    candidate = char(candidates(i));
    if isfolder(candidate) && isfile(fullfile(candidate, 'readNPY.m'))
        addpath(candidate);
        addedPath = string(candidate);
        return
    end
end

error('ALFloader:MissingReadNPY', ...
    ['readNPY.m is required to load ALF .npy files. Add npy-matlab to the MATLAB path, ', ...
    'for example: addpath(''C:\Users\ichaker\Desktop\FANCiRepos\npy-matlab\npy-matlab'').']);
end

function sourcePaths = localLoadFanciSourcePaths(alfPath)
sourcePaths = struct( ...
    'rawDataPath', string.empty(0, 0), ...
    's2pSessionPath', string.empty(0, 0));

sourcePaths.rawDataPath = localReadJsonString( ...
    fullfile(alfPath, '_FANCi_source.rawDataPath.json'), ...
    'rawDataPath');
sourcePaths.s2pSessionPath = localReadJsonString( ...
    fullfile(alfPath, '_FANCi_source.s2pSessionPath.json'), ...
    's2pSessionPath');
end

function value = localReadJsonString(path, fieldName)
value = string.empty(0, 0);
if ~isfile(path)
    return
end

decoded = jsondecode(fileread(path));
if isstruct(decoded) && isfield(decoded, fieldName)
    decoded = decoded.(fieldName);
end

if ischar(decoded) || (isstring(decoded) && isscalar(decoded))
    value = string(decoded);
else
    error('ALFloader:BadFanciSourceMetadata', ...
        'Expected %s to contain a JSON string or object field "%s".', path, fieldName);
end
end

function sessions = localNormalizeSessions(sessions)
required = {'subject', 'date', 'session'};

if istable(sessions)
    sessions = localSessionsFromTable(sessions, required);
elseif isstruct(sessions)
    sessions = localSessionsFromStruct(sessions, required);
else
    error('ALFloader:BadSessions', ...
        'sessions must be a struct array or legacy table with subject, date, and session.');
end
end

function sessions = localSessionsFromTable(sessionTable, required)
for i = 1:numel(required)
    if ~ismember(required{i}, sessionTable.Properties.VariableNames)
        error('ALFloader:BadSessions', 'sessions is missing required variable: %s', required{i});
    end
end

nSessions = height(sessionTable);
sessions = localEmptySessionStructArray(nSessions);
for iSession = 1:nSessions
    for iField = 1:numel(required)
        fieldName = required{iField};
        sessions(iSession).(fieldName) = localSessionStringScalar( ...
            sessionTable.(fieldName)(iSession), fieldName, iSession);
    end
end
end

function sessions = localSessionsFromStruct(sessionStruct, required)
for i = 1:numel(required)
    if ~isfield(sessionStruct, required{i})
        error('ALFloader:BadSessions', 'sessions is missing required field: %s', required{i});
    end
end

sessionStruct = sessionStruct(:);
nSessions = numel(sessionStruct);
sessions = localEmptySessionStructArray(nSessions);
for iSession = 1:nSessions
    for iField = 1:numel(required)
        fieldName = required{iField};
        sessions(iSession).(fieldName) = localSessionStringScalar( ...
            sessionStruct(iSession).(fieldName), fieldName, iSession);
    end
end
end

function value = localSessionStringScalar(value, fieldName, iSession)
value = string(value);
if ~isscalar(value)
    error('ALFloader:BadSessions', ...
        'sessions(%d).%s must be a string scalar.', iSession, fieldName);
end
end

function sessions = localDiscoverSessions(ALFroot)
subjects = strings(0, 1);
dates = strings(0, 1);
sessionValues = strings(0, 1);

subjectDirs = localSubdirs(ALFroot);
for iSubject = 1:numel(subjectDirs)
    subjectName = subjectDirs{iSubject};
    dateDirs = localSubdirs(fullfile(ALFroot, subjectName));
    for iDate = 1:numel(dateDirs)
        dateName = dateDirs{iDate};
        sessionDirs = localSubdirs(fullfile(ALFroot, subjectName, dateName));
        for iSession = 1:numel(sessionDirs)
            sessionName = sessionDirs{iSession};
            alfPath = fullfile(ALFroot, subjectName, dateName, sessionName, 'alf');
            if isfolder(alfPath)
                subjects(end + 1, 1) = string(subjectName); %#ok<AGROW>
                dates(end + 1, 1) = string(dateName); %#ok<AGROW>
                sessionValues(end + 1, 1) = string(sessionName); %#ok<AGROW>
            end
        end
    end
end

sessions = localSessionStructArray(subjects, dates, sessionValues);
end

function sessions = localSessionStructArray(subjects, dates, sessionValues)
if isempty(subjects)
    sessions = localEmptySessionStructArray();
    return
end

sessionTable = table(subjects(:), dates(:), sessionValues(:), ...
    'VariableNames', {'subject', 'date', 'session'});
sessions = localSessionsFromTable(sessionTable, {'subject', 'date', 'session'});
end

function sessions = localEmptySessionStructArray(nSessions)
if nargin < 1
    nSessions = 0;
end

sessions = repmat(struct( ...
    'subject', string.empty(0, 0), ...
    'date', string.empty(0, 0), ...
    'session', string.empty(0, 0)), nSessions, 1);
end

function names = localSubdirs(folder)
d = dir(folder);
d = d([d.isdir]);
names = {d.name};
names = names(~ismember(names, {'.', '..'}));
names = sort(names);
end

function fovNames = localListFOVs(sessionPath)
d = dir(fullfile(sessionPath, 'alf', 'FOV_*'));
d = d([d.isdir]);
fovNames = sort({d.name});
end

function files = localListFiles(folder)
files = strings(0, 1);
if ~isfolder(folder)
    return
end
d = dir(fullfile(folder, '*'));
d = d(~[d.isdir]);
for i = 1:numel(d)
    files(end + 1, 1) = string(fullfile(folder, d(i).name)); %#ok<AGROW>
end
end

function [fov, errors] = localLoadFov(fovName, fovFolder)
fov = localEmptyFov();
errors = localEmptyErrors();
fov.name = string(fovName);
fov.path = string(fovFolder);

knownFiles = { ...
    'mpci.times.npy', ...
    'mpci.badFrames.npy', ...
    'mpci.ROIActivityF.npy', ...
    'mpci.ROINeuropilActivityF.npy', ...
    'mpci.ROIActivityDeconvolved.npy', ...
    'mpciROIs.cellClassifier.npy', ...
    'mpciROIs.mpciROITypes.npy', ...
    'mpciROIs.stackPos.npy', ...
    'mpciROIs.masks.sparse_npz', ...
    'mpciMeanImage.images.npy'};

for i = 1:numel(knownFiles)
    path = fullfile(fovFolder, knownFiles{i});
    if isfile(path)
        fov.files(end + 1, 1) = string(path);
    else
        fov.missing(end + 1, 1) = string(path);
    end
end

[pathOk, value, err] = localTryReadColumn(fullfile(fovFolder, 'mpci.times.npy'), true);
if pathOk
    fov.times = value;
    fov.nFrames = numel(value);
elseif isfile(fullfile(fovFolder, 'mpci.times.npy'))
    errors = localAddError(errors, 'mpci.times', fullfile(fovFolder, 'mpci.times.npy'), err);
end

nSamples = fov.nFrames;
[pathOk, value, err] = localTryReadMatrix(fullfile(fovFolder, 'mpci.ROIActivityF.npy'), nSamples);
if pathOk
    fov.F = value;
else
    errors = localAddErrorIfPresent(errors, 'mpci.ROIActivityF', fullfile(fovFolder, 'mpci.ROIActivityF.npy'), err);
end

[pathOk, value, err] = localTryReadMatrix(fullfile(fovFolder, 'mpci.ROINeuropilActivityF.npy'), nSamples);
if pathOk
    fov.Fneu = value;
else
    errors = localAddErrorIfPresent(errors, 'mpci.ROINeuropilActivityF', fullfile(fovFolder, 'mpci.ROINeuropilActivityF.npy'), err);
end

[pathOk, value, err] = localTryReadMatrix(fullfile(fovFolder, 'mpci.ROIActivityDeconvolved.npy'), nSamples);
if pathOk
    fov.deconvolved = value;
else
    errors = localAddErrorIfPresent(errors, 'mpci.ROIActivityDeconvolved', fullfile(fovFolder, 'mpci.ROIActivityDeconvolved.npy'), err);
end

[pathOk, value, err] = localTryReadColumn(fullfile(fovFolder, 'mpci.badFrames.npy'), false);
if pathOk
    fov.badFrames = logical(value);
else
    errors = localAddErrorIfPresent(errors, 'mpci.badFrames', fullfile(fovFolder, 'mpci.badFrames.npy'), err);
end

[pathOk, value, err] = localTryReadArray(fullfile(fovFolder, 'mpciROIs.cellClassifier.npy'));
if pathOk
    fov.cellClassifier = double(value);
else
    errors = localAddErrorIfPresent(errors, 'mpciROIs.cellClassifier', fullfile(fovFolder, 'mpciROIs.cellClassifier.npy'), err);
end

[pathOk, value, err] = localTryReadArray(fullfile(fovFolder, 'mpciROIs.mpciROITypes.npy'));
if pathOk
    fov.roiTypes = double(value(:));
else
    errors = localAddErrorIfPresent(errors, 'mpciROIs.mpciROITypes', fullfile(fovFolder, 'mpciROIs.mpciROITypes.npy'), err);
end

[pathOk, value, err] = localTryReadArray(fullfile(fovFolder, 'mpciROIs.stackPos.npy'));
if pathOk
    fov.stackPos = double(value);
else
    errors = localAddErrorIfPresent(errors, 'mpciROIs.stackPos', fullfile(fovFolder, 'mpciROIs.stackPos.npy'), err);
end

[pathOk, value, err] = localTryReadArray(fullfile(fovFolder, 'mpciMeanImage.images.npy'));
if pathOk
    fov.meanImage = double(value);
else
    errors = localAddErrorIfPresent(errors, 'mpciMeanImage.images', fullfile(fovFolder, 'mpciMeanImage.images.npy'), err);
end

masksFile = fullfile(fovFolder, 'mpciROIs.masks.sparse_npz');
if isfile(masksFile)
    fov.masksFile = string(masksFile);
end

fov.nROIs = localInferNRois(fov);
fov.loaded = ~isempty(fov.times) || ~isempty(fov.F) || ~isempty(fov.Fneu) || ~isempty(fov.deconvolved);
if ~fov.loaded
    fov.errorMessage = "No supported calcium arrays loaded.";
end
end

function nRois = localInferNRois(fov)
candidates = [];
if ~isempty(fov.F) && ndims(fov.F) == 2
    candidates(end + 1) = size(fov.F, 2); %#ok<AGROW>
end
if ~isempty(fov.Fneu) && ndims(fov.Fneu) == 2
    candidates(end + 1) = size(fov.Fneu, 2); %#ok<AGROW>
end
if ~isempty(fov.deconvolved) && ndims(fov.deconvolved) == 2
    candidates(end + 1) = size(fov.deconvolved, 2); %#ok<AGROW>
end
if ~isempty(fov.cellClassifier)
    candidates(end + 1) = size(fov.cellClassifier, 1); %#ok<AGROW>
end
if ~isempty(fov.roiTypes)
    candidates(end + 1) = numel(fov.roiTypes); %#ok<AGROW>
end
if ~isempty(fov.stackPos)
    candidates(end + 1) = size(fov.stackPos, 1); %#ok<AGROW>
end
if isempty(candidates)
    nRois = 0;
else
    nRois = max(candidates);
end
end

function [sparseNoise, errors] = localLoadSparseNoise(sessionPath)
sparseNoise = localEmptySparseNoise();
errors = localEmptyErrors();
alfFolder = fullfile(sessionPath, 'alf');
rawPassiveFolder = fullfile(sessionPath, 'raw_passive_data');
sparseNoise.rawPassiveFolder = string(rawPassiveFolder);

frameTimesPath = fullfile(alfFolder, '_ibl_passiveRFM.times.npy');
rawMoviePath = fullfile(rawPassiveFolder, '_iblrig_RFMapStim.raw.bin');
eventTimesPath = fullfile(alfFolder, '_ibl_sparseNoise.times.npy');
xyPath = fullfile(alfFolder, '_ibl_sparseNoise.xy.npy');
knownFiles = {frameTimesPath, rawMoviePath, eventTimesPath, xyPath};

for i = 1:numel(knownFiles)
    if isfile(knownFiles{i})
        sparseNoise.files(end + 1, 1) = string(knownFiles{i});
    else
        sparseNoise.missing(end + 1, 1) = string(knownFiles{i});
    end
end

[pathOk, value, err] = localTryReadColumn(frameTimesPath, true);
if pathOk
    sparseNoise.frameTimes = value;
else
    errors = localAddErrorIfPresent(errors, '_ibl_passiveRFM.times', frameTimesPath, err);
end

if isfile(rawMoviePath) && ~isempty(sparseNoise.frameTimes)
    try
        sparseNoise.frames = localReadSparseNoiseRawBin(rawMoviePath, numel(sparseNoise.frameTimes));
    catch exc
        errors = localAddError(errors, '_iblrig_RFMapStim.raw', rawMoviePath, exc.message);
    end
end

[pathOk, value, err] = localTryReadColumn(eventTimesPath, false);
if pathOk
    sparseNoise.eventTimes = value;
else
    errors = localAddErrorIfPresent(errors, '_ibl_sparseNoise.times', eventTimesPath, err);
end

[pathOk, value, err] = localTryReadArray(xyPath);
if pathOk
    sparseNoise.xy = double(value);
    if size(sparseNoise.xy, 2) ~= 2 && size(sparseNoise.xy, 1) == 2
        sparseNoise.xy = sparseNoise.xy.';
    end
    sparseNoise.edges = localInferEdges(sparseNoise.xy);
else
    errors = localAddErrorIfPresent(errors, '_ibl_sparseNoise.xy', xyPath, err);
end

sparseNoise.loaded = ~isempty(sparseNoise.frameTimes) || ~isempty(sparseNoise.frames) || ...
    ~isempty(sparseNoise.eventTimes) || ~isempty(sparseNoise.xy);
if ~sparseNoise.loaded && any(cellfun(@isfile, knownFiles))
    sparseNoise.errorMessage = "Sparse-noise files were present but no supported arrays loaded.";
end

if sparseNoise.loaded
    fprintf('  [sparseNoise] loaded frameTimes=%d, events=%d\n', ...
        numel(sparseNoise.frameTimes), numel(sparseNoise.eventTimes));
else
    fprintf('  [sparseNoise] no supported sparse-noise files loaded\n');
end
end

function [wheel, errors] = localLoadWheel(sessionPath)
wheel = localEmptyWheel();
errors = localEmptyErrors();
alfFolder = fullfile(sessionPath, 'alf');

positionPath = fullfile(alfFolder, '_ibl_wheel.position.npy');
timestampsPath = fullfile(alfFolder, '_ibl_wheel.timestamps.npy');
velocityPath = fullfile(alfFolder, '_ibl_wheel.velocity.npy');
knownFiles = {positionPath, timestampsPath, velocityPath};

for i = 1:numel(knownFiles)
    if isfile(knownFiles{i})
        wheel.files(end + 1, 1) = string(knownFiles{i});
    else
        wheel.missing(end + 1, 1) = string(knownFiles{i});
    end
end

[pathOk, value, err] = localTryReadColumn(positionPath, false);
if pathOk
    wheel.position = value;
else
    errors = localAddErrorIfPresent(errors, '_ibl_wheel.position', positionPath, err);
end

[pathOk, value, err] = localTryReadColumn(timestampsPath, true);
if pathOk
    wheel.timestamps = value;
else
    errors = localAddErrorIfPresent(errors, '_ibl_wheel.timestamps', timestampsPath, err);
end

[pathOk, value, err] = localTryReadColumn(velocityPath, false);
if pathOk
    wheel.velocity = value;
else
    errors = localAddErrorIfPresent(errors, '_ibl_wheel.velocity', velocityPath, err);
end

wheel.loaded = ~isempty(wheel.position) || ~isempty(wheel.timestamps) || ~isempty(wheel.velocity);
if ~wheel.loaded && any(cellfun(@isfile, knownFiles))
    wheel.errorMessage = "Wheel files were present but no supported arrays loaded.";
end

if wheel.loaded
    fprintf('  [wheel] loaded samples=%d\n', ...
        max([numel(wheel.position), numel(wheel.timestamps), numel(wheel.velocity)]));
else
    fprintf('  [wheel] no supported wheel files loaded\n');
end
end

function [ok, value, message] = localTryReadArray(path)
ok = false;
value = [];
message = '';
if ~isfile(path)
    message = sprintf('File not found: %s', path);
    return
end
try
    value = readNPY(path);
    ok = true;
catch exc
    message = exc.message;
end
end

function [ok, value, message] = localTryReadColumn(path, requireIncreasing)
[ok, value, message] = localTryReadArray(path);
if ~ok
    return
end
try
    value = double(value);
    value = value(:);
    if ~all(isfinite(value))
        error('Vector contains non-finite values.');
    end
    if requireIncreasing && numel(value) > 1 && any(diff(value) <= 0)
        error('Vector is not strictly increasing.');
    end
catch exc
    ok = false;
    message = exc.message;
    value = [];
end
end

function [ok, value, message] = localTryReadMatrix(path, nSamples)
[ok, value, message] = localTryReadArray(path);
if ~ok
    return
end
try
    value = double(value);
    if ndims(value) ~= 2
        error('Expected a 2D matrix, got ndims=%d.', ndims(value));
    end
    if nargin >= 2 && ~isempty(nSamples) && nSamples > 0
        if size(value, 1) ~= nSamples && size(value, 2) == nSamples
            value = value.';
        end
        if size(value, 1) ~= nSamples
            error('Matrix first dimension does not match time vector: %d vs %d.', size(value, 1), nSamples);
        end
    end
catch exc
    ok = false;
    message = exc.message;
    value = [];
end
end

function stimMap = localReadSparseNoiseRawBin(filename, nFrames)
fid = fopen(filename, 'rb');
if fid < 0
    error('Could not open sparse-noise raw movie: %s', filename);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

raw = fread(fid, inf, 'int8=>int8');

shapes = [10 30; 7 27];
for iShape = 1:size(shapes, 1)
    nRows = shapes(iShape, 1);
    nCols = shapes(iShape, 2);
    if numel(raw) == nFrames * nRows * nCols
        stimMap = reshape(raw, [nRows, nCols, nFrames]);
        stimMap = permute(stimMap, [3 1 2]);
        stimMap = localDecodeSparseNoiseValues(stimMap);
        stimMap = int8(stimMap);
        return
    end
end

error('Could not infer sparse-noise movie shape from %s with %d frames.', filename, nFrames);
end

function decoded = localDecodeSparseNoiseValues(rawSigned)
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

function edges = localInferEdges(xy)
edges = [];
if isempty(xy) || size(xy, 2) ~= 2
    return
end
xVals = sort(unique(xy(:, 1)));
yVals = sort(unique(xy(:, 2)));
if isempty(xVals) || isempty(yVals)
    return
end

dx = 1;
dy = 1;
if numel(xVals) > 1
    dx = median(diff(xVals), 'omitnan');
end
if numel(yVals) > 1
    dy = median(diff(yVals), 'omitnan');
end
edges = [min(xVals) - dx / 2, max(xVals) + dx / 2, min(yVals) - dy / 2, max(yVals) + dy / 2];
end

function errors = localAddErrorIfPresent(errors, dataset, path, message)
if isfile(path)
    errors = localAddError(errors, dataset, path, message);
end
end

function errors = localAddError(errors, dataset, path, message)
errors(end + 1, 1) = struct( ... %#ok<AGROW>
    'dataset', string(dataset), ...
    'path', string(path), ...
    'message', string(message));
end

function errors = localConcatErrors(errors, moreErrors)
if isempty(moreErrors)
    return
end
if isempty(errors)
    errors = moreErrors;
else
    errors = [errors; moreErrors];
end
end
