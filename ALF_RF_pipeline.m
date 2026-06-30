%% ALF_RF_pipeline

clear;
clc;

repoRoot = fileparts(mfilename('fullpath'));
addpath(genpath(repoRoot));

ALFroot = 'D:\Pipeline\Data'; % ALF export path. This basically rewrites s2p outputs into standerdized format under ALFroot/<subject>/<date>/<session>
% the implementation so far can split multi session s2p outputs, but not yet multi date s2p outputs

resultsRoot = 'D:\Pipeline\Results'; % RF mapping results, these also follow ALF paths resultsRoot/<subject>/<date>/<session>
plotsRoot = 'D:\Pipeline\Plots'; % Plotting is divided into groups  plotsRoot/<group>/<subject>/<date>/<session>
plotGroup = 'groupA'; % groupA=MO groupB=V1  groupC=Axons  groupD=V1->MO. can be 'all'
modality = 'deconv'; % this can be 'deconv' or 'traces' depending on what format is best for the recording. 'best' automatically selects based on the first FOV
runningRegressor = false;
useL1 = false; % this is very slow but may gain you few good RFs


s2pPaths = ["\\10.233.25.135\FANCiNAS1\Data\2P\NM035\Processed\20260513\1"];

[sessions, exportReport] = ALFexporter(s2pPaths, ALFroot);
loadedSessions = ALFloader(ALFroot, sessions);
summaries = ALF_main_analyse_RFs(loadedSessions, modality, resultsRoot, runningRegressor, useL1);

flipFOVs = [
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
    ]; % subject/date/session/FOV_XX entries whose plotted RF azimuth should be flipped

RF_plots(plotsRoot, resultsRoot, plotGroup, flipFOVs);

















function status = RF_plots(plotsRoot, resultsRoot, plotGroup, flipFOVs, pythonCommand, forceReplot)
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
if nargin < 5 || isempty(pythonCommand)
    pythonCommand = 'python';
end
if nargin < 6 || isempty(forceReplot)
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

pythonCommandText = sprintf('%s %s --results-root %s --plots-root %s --group %s%s%s', ...
    localQuoteCommandArg(pythonCommand), ...
    localQuoteCommandArg(scriptPath), ...
    localQuoteCommandArg(resultsRoot), ...
    localQuoteCommandArg(plotsRoot), ...
    localQuoteCommandArg(plotGroup), ...
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
