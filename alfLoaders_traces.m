function loaders = alfLoaders_traces()
%ALFLOADERS_TRACES Lightweight loaders for ALF/ONE raw trace 2P datasets.
%
% Usage:
%   loaders = alfLoaders_traces();
%   sessionPath = loaders.sessionPath('D:\Data\ALF2P', 'NM023', '2026-04-28', '1');
%   fovNames = loaders.listFOVs(sessionPath);
%   calcium = loaders.loadCalciumTraces(sessionPath, fovNames{1});
%   stimulus = loaders.loadStimulusData(sessionPath, [-135 135 -40 40]);
%   wheel = loaders.loadWheelData(sessionPath);
%
% The functions assume the ALF files already exist and follow the Python
% exporter layout:
%   session/alf/FOV_XX/mpci.ROIActivityF.npy
%   session/alf/FOV_XX/mpci.ROINeuropilActivityF.npy
%   session/alf/_ibl_wheel.*.npy
%   session/alf/_ibl_passiveRFM.times.npy
%   session/raw_passive_data/_iblrig_RFMapStim.raw.bin

loaders.sessionPath = @localSessionPath;
loaders.alfPath = @localAlfPath;
loaders.listFOVs = @localListFOVs;
loaders.loadCalciumTraces = @localLoadCalciumTraces;
loaders.loadStimulusData = @localLoadStimulusData;
loaders.loadWheelData = @localLoadWheelData;
end

function folder = localSessionPath(root, subject, date, session)
folder = fullfile(char(root), char(subject), char(date), char(string(session)));
end

function folder = localAlfPath(sessionPath)
folder = fullfile(sessionPath, 'alf');
end

function fovNames = localListFOVs(sessionPath)
d = dir(fullfile(localAlfPath(sessionPath), 'FOV_*'));
d = d([d.isdir]);
fovNames = sort({d.name});
end

function calcium = localLoadCalciumTraces(sessionPath, fovName)
alfFolder = fullfile(localAlfPath(sessionPath), char(fovName));

calcium = struct();
calcium.alfFolder = alfFolder;
calcium.fov = char(fovName);
calcium.timestamps = localReadColumn(fullfile(alfFolder, 'mpci.times.npy'));
calcium.raw = localReadMatrix(fullfile(alfFolder, 'mpci.ROIActivityF.npy'), numel(calcium.timestamps));
calcium.neuropil = localReadMatrix(fullfile(alfFolder, 'mpci.ROINeuropilActivityF.npy'), numel(calcium.timestamps));

file = fullfile(alfFolder, 'mpci.ROIActivityDeconvolved.npy');
if isfile(file)
    calcium.deconvolved = localReadMatrix(file, numel(calcium.timestamps));
else
    calcium.deconvolved = [];
end

file = fullfile(alfFolder, 'mpci.badFrames.npy');
if isfile(file)
    calcium.badFrames = logical(localReadColumn(file));
else
    calcium.badFrames = [];
end

file = fullfile(alfFolder, 'mpciROIs.cellClassifier.npy');
if isfile(file)
    calcium.cellClassifier = localReadColumn(file);
else
    calcium.cellClassifier = [];
end

file = fullfile(alfFolder, 'mpciROIs.mpciROITypes.npy');
if isfile(file)
    calcium.roiTypes = localReadColumn(file);
else
    calcium.roiTypes = [];
end

file = fullfile(alfFolder, 'mpciROIs.stackPos.npy');
if isfile(file)
    calcium.stackPos = double(readNPY(file));
else
    calcium.stackPos = [];
end
end

function stimulus = localLoadStimulusData(sessionPath, stimulusEdgesDeg)
alfFolder = localAlfPath(sessionPath);
rawPassiveFolder = fullfile(sessionPath, 'raw_passive_data');

stimulus = struct();
stimulus.alfFolder = alfFolder;
stimulus.rawPassiveFolder = rawPassiveFolder;
stimulus.frameTimes = localReadColumn(fullfile(alfFolder, '_ibl_passiveRFM.times.npy'));
stimulus.frames = localReadSparseNoiseRawBin( ...
    fullfile(rawPassiveFolder, '_iblrig_RFMapStim.raw.bin'), ...
    numel(stimulus.frameTimes));
stimulus.eventTimes = localReadColumn(fullfile(alfFolder, '_ibl_sparseNoise.times.npy'));
stimulus.xy = double(readNPY(fullfile(alfFolder, '_ibl_sparseNoise.xy.npy')));

if nargin >= 2 && ~isempty(stimulusEdgesDeg) && numel(stimulusEdgesDeg) == 4
    stimulus.edges = double(stimulusEdgesDeg(:).');
else
    stimulus.edges = localInferEdges(stimulus.xy);
end
end

function wheel = localLoadWheelData(sessionPath)
alfFolder = localAlfPath(sessionPath);

wheel = struct();
wheel.alfFolder = alfFolder;
wheel.position = localReadColumn(fullfile(alfFolder, '_ibl_wheel.position.npy'));
wheel.timestamps = localReadColumn(fullfile(alfFolder, '_ibl_wheel.timestamps.npy'));
wheel.velocity = localReadColumn(fullfile(alfFolder, '_ibl_wheel.velocity.npy'));
end

function x = localReadColumn(filename)
x = double(readNPY(filename));
x = x(:);
end

function x = localReadMatrix(filename, nSamples)
x = double(readNPY(filename));
if size(x, 1) ~= nSamples && size(x, 2) == nSamples
    x = x.';
end
end

function stimMap = localReadSparseNoiseRawBin(filename, nFrames)
fid = fopen(filename, 'rb');
if fid < 0
    error('Could not open sparse-noise raw movie: %s', filename);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

% Match the legacy MATLAB load_events.m sparse-noise interpretation exactly:
%   stimulus = fread(fid, 'int8');
%   stimulus = reshape(stimulus, rows, cols, []);
%
% The RF analysis expects [frame x row x col], so after the legacy reshape
% we permute [row x col x frame] -> [frame x row x col].
%
% Do not use the Python/C-order reshape here.  The raw bytes are the same,
% but C-order reshape puts sparse-noise squares in different screen pixels.
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
if size(xy, 2) ~= 2 && size(xy, 1) == 2
    xy = xy.';
end
xVals = sort(unique(xy(:, 1)));
yVals = sort(unique(xy(:, 2)));

dx = 1;
dy = 1;
if numel(xVals) > 1
    dx = median(diff(xVals), 'omitnan');
end
if numel(yVals) > 1
    dy = median(diff(yVals), 'omitnan');
end
edges = [min(xVals)-dx/2, max(xVals)+dx/2, min(yVals)-dy/2, max(yVals)+dy/2];
end
