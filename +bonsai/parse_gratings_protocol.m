function out = parse_gratings_protocol(filename)
%PARSE_EXPERIMENTAL_PROTOCOL Parse a stimulus protocol spreadsheet.
%
% out = parse_experimental_protocol(filename)
%
% Expected parameters are:
%   direction, tf, sf, size, contrast
%
% Special blank-trial rule:
%   A row with direction == 360 and contrast == 0 is treated as a blank
%   trial. Blank trials are treated as stimuli with a different direction,
%   so blank-only changes in contrast do NOT make contrast an active
%   parameter. Contrast is considered active only if it varies outside
%   those blank trials.
%
% The function:
%   1) finds which parameters vary across stimuli (activePars)
%   2) identifies the unique base set of stimuli
%   3) sorts that base set to match the requested example ordering:
%        direction changes fastest, then size, then tf, then sf
%      (therefore the lexicographic sort keys are reversed for sortrows)
%   4) maps each row of the original experiment back to the sorted stimulus idx
%   5) builds one text label per unique stimulus in the base set
%
% Returned fields:
%   out.protocol                    standardized input table
%   out.activePars                  cell array of varying parameter names
%   out.uniquePars                  struct with effective unique values used
%                                   to determine activePars
%   out.nValuesPerPar               struct with number of effective values
%   out.stimuli                     sorted unique stimulus combinations
%   out.stimLabels                  cell array of labels, one per row of
%                                   out.baseStimuli
%   out.repeatCount                 repetitions per unique stimulus
%   out.stimulusIdxSequence         1-based idx for each presented stimulus row
%   out.blankTrials                 logical column vector, true only for
%                                   blank trials in presentation order
%   out.trialStimIdx                [nTrials x nBaseStimuli] if contiguous
%                                   trials exist
%   out.nRepeats                    number of contiguous trials, when detected
%   out.nStim                       number of unique stimuli
%   out.nActivePars                 number of active pars

%   out.hasCompleteContiguousTrials logical flag
%
% Example:
%   out = parse_experimental_protocol('metadata.csv');
%   baseStimuli = out.baseStimuli;
%   labels      = out.Stimulus_labels;
%   seq         = out.stimulusIdxSequence;
%   blank       = out.blankTrials;
%   trialIdx    = out.trialStimIdx;

% 2026 chatGPT created, edited by LFR

    parNames = {'direction','sf','tf','size','contrast'};
    protocol = readProtocolTable(filename, parNames);
    out = struct();
    out.protocol = protocol;

    %% Blank trials
    blankTrials = identifyBlankTrials(protocol);
    out.blankTrials = blankTrials;

    %% 1) Identify active parameters
    out.activePars = {};
    out.uniquePars = struct();
    out.nValuesPerPar = struct();

    rawContrastVals = sort(unique(protocol.contrast));
    effectiveContrastVals = getEffectiveContrastValues(protocol, blankTrials);
    contrastIgnoredBecauseBlank = ...
        (numel(rawContrastVals) > numel(effectiveContrastVals));

    fprintf('Blank-trial handling:\n');
    fprintf('  identified %d blank trials (direction == 360 and contrast == 0)\n', ...
        nnz(blankTrials));
    if contrastIgnoredBecauseBlank
        fprintf(['  contrast changes caused only by blank trials are ignored ' ...
                 'when determining activePars\n']);
    end

    fprintf('\nChanged parameters (activePars):\n');
    for i = 1:numel(parNames)
        p = parNames{i};

        if strcmp(p, 'contrast')
            vals = effectiveContrastVals;
        else
            vals = sort(unique(protocol.(p)));
        end

        out.uniquePars.(p) = vals;
        out.nValuesPerPar.(p) = numel(vals);

        if numel(vals) > 1
            out.activePars{end+1} = p; %#ok<AGROW>
            fprintf('  %s: %d values -> %s\n', p, numel(vals), mat2str(vals(:).'));
        end
    end
    if isempty(out.activePars)
        fprintf('  none (all parameters are constant)\n');
    end

    fprintf('\nConstant parameters:\n');
    for i = 1:numel(parNames)
        p = parNames{i};
        vals = out.uniquePars.(p);
        if numel(vals) == 1
            if strcmp(p, 'contrast') && contrastIgnoredBecauseBlank
                fprintf('  %s: %s (outside blank trials)\n', p, mat2str(vals(:).'));
            else
                fprintf('  %s: %s\n', p, mat2str(vals(:).'));
            end
        end
    end

    %% 2) Unique base stimulus set
    % Keep the raw table values here so blank stimuli remain identifiable
    % as direction==360, contrast==0 rows in the base set.
    stimuli = unique(protocol(:, parNames), 'rows');

    % To match the user's example, direction varies fastest, then size,
    % then tf, then sf. sortrows() uses the first key as the slowest/
    % outermost grouping, so we reverse those keys here.
    sortKeys = {'sf','tf','size','direction'};

    % Only include contrast in the sort if it is truly active outside blank
    % trials. Otherwise blank rows are treated as stimuli that differ only
    % in direction.
    if out.nValuesPerPar.contrast > 1
        sortKeys{end+1} = 'contrast';
    end

    stimuli = sortrows(stimuli, sortKeys);

    out.stimuli = stimuli;
    out.stimLabels = buildStimulusLabels(stimuli);

%    % % Alternative building the labels (according to the active parameters)
% nActivePars=length(p.activepars);
% for iStim=1:p.nstim
%     for iActivePar=1:nActivePars
%         ind=p.activepars{iActivePar};
%         ind=(ind(min(length(ind), 2)));
%         if iActivePar==1
%             stimSequence.labels{iStim}=sprintf('%s = %d', p.parnames{ind}, p.pars(ind, iStim));
%         else
%             stimSequence.labels{iStim}=sprintf('%s, %s = %d', stimSequence.labels{iStim}, p.parnames{ind}, p.pars(ind, iStim));
%         end
%     end
%     if ismember(iStim, p.blankstims)
%         stimSequence.labels{iStim}='blank';
%     end
% end
    %% 3) Map each presented stimulus row to the sorted base stimulus idx
    [isMember, stimSequence] = ismember(protocol(:, parNames), stimuli, 'rows');
    if ~all(isMember)
        error('Could not map all rows back to the sorted base stimulus set.');
    end

    repeatCount = accumarray(stimSequence, 1, [height(stimuli) 1]);

    out.repeatCount = repeatCount;
    out.stimSequence = stimSequence;

    fprintf('\nBasic stimulus set:\n');
    fprintf('  %d unique stimulus combinations\n', height(stimuli));
    if all(repeatCount == repeatCount(1))
        fprintf('  each unique stimulus is repeated %d times\n', repeatCount(1));
    else
        fprintf('  repeat counts vary across stimuli: min = %d, max = %d\n', ...
            min(repeatCount), max(repeatCount));
    end

    %% 4) Detect contiguous trials, if the rows split cleanly into repeats
    out.trialStimIdx = [];
    out.nRepeats = NaN;
    out.hasCompleteContiguousTrials = false;

    nBase = height(stimuli);
    nRows = numel(stimSequence);

    if mod(nRows, nBase) == 0
        candidateTrialStimIdx = reshape(stimSequence, nBase, []).';
        expected = repmat(1:nBase, size(candidateTrialStimIdx, 1), 1);
        isCompleteTrial = all(sort(candidateTrialStimIdx, 2) == expected, 2);

        if all(isCompleteTrial)
            out.trialStimIdx = candidateTrialStimIdx;
            out.nRepeats = size(candidateTrialStimIdx, 1);
            out.hasCompleteContiguousTrials = true;
            fprintf('  detected %d contiguous trials of %d stimuli each\n', ...
                out.nRepeats, nBase);
        else
            fprintf('  rows do not split into complete contiguous trials\n');
        end
    else
        fprintf('  total row count is not an integer multiple of the base set size\n');
    end

    nShow = min(20, numel(stimSequence));
    fprintf('\nFirst %d stimulus idx in presentation order:\n', nShow);
    fprintf('  %s\n', mat2str(stimSequence(1:nShow).'));

    nBlankShow = min(20, numel(blankTrials));
    fprintf('First %d blank-trial flags in presentation order:\n', nBlankShow);
    fprintf('  %s\n', mat2str(double(blankTrials(1:nBlankShow)).'));

    nLabelShow = min(5, numel(out.stimLabels));
    fprintf('First %d stimulus labels:\n', nLabelShow);
    for i = 1:nLabelShow
        fprintf('  %d -> %s\n', i, out.stimLabels{i});
    end

    out.nStim = numel(out.stimLabels);
    out.nActivePars = numel(out.activePars);
    out.blankStimIdx = unique(out.stimSequence(out.blankTrials));
    out.stimIdx = unique(out.stimSequence(~out.blankTrials));

end

function T = readProtocolTable(filename, parNames)
% Read either a spreadsheet with headers or a plain 5-column numeric file.

    aliases = struct();
    aliases.direction = {'direction','dir'};
    aliases.tf        = {'tf','temporalfrequency','temporalfreq'};
    aliases.sf        = {'sf','spatialfrequency','spatialfreq'};
    aliases.size      = {'size','stimsize','diameter'};
    aliases.contrast  = {'contrast','ctr'};

    T = table();
    useHeader = false;

    try
        T0 = readtable(filename);
        idx = matchColumns(T0.Properties.VariableNames, parNames, aliases);
        if all(~isnan(idx))
            T = T0(:, idx);
            T.Properties.VariableNames = parNames;
            useHeader = true;
        end
    catch
        % Fall back to headerless read below.
    end

    if ~useHeader
        T = readtable(filename, 'ReadVariableNames', false);
        if width(T) < numel(parNames)
            error('File must contain at least %d columns.', numel(parNames));
        end
        T = T(:, 1:numel(parNames));
        T.Properties.VariableNames = parNames;
    end

    % Remove completely empty rows before numeric conversion.
    emptyRow = all(ismissing(T), 2);
    T(emptyRow, :) = [];

    % Convert to numeric.
    for i = 1:numel(parNames)
        p = parNames{i};
        T.(p) = forceNumeric(T.(p), p);
    end
end

function idx = matchColumns(varNames, parNames, aliases)
    idx = nan(1, numel(parNames));
    normVars = cellfun(@normalizeName, varNames, 'UniformOutput', false);

    for i = 1:numel(parNames)
        thisPar = parNames{i};
        aliasList = aliases.(thisPar);
        normAlias = cellfun(@normalizeName, aliasList, 'UniformOutput', false);
        hit = find(ismember(normVars, normAlias), 1, 'first');
        if ~isempty(hit)
            idx(i) = hit;
        end
    end
end

function x = forceNumeric(x, colName)
    if isnumeric(x)
        x = double(x);
        return;
    end

    if islogical(x)
        x = double(x);
        return;
    end

    x = str2double(string(x));
    if any(isnan(x))
        error('Column "%s" contains non-numeric values that could not be parsed.', colName);
    end
end

function s = normalizeName(s)
    s = char(string(s));
    s = lower(s);
    s = regexprep(s, '[^a-z0-9]', '');
end

function blankTrials = identifyBlankTrials(T)
% Blank trials are defined as direction == 360 and contrast == 0.

    tol = 1e-12;
    blankTrials = abs(T.direction - 360) < tol & abs(T.contrast) < tol;
end

function vals = getEffectiveContrastValues(T, blankTrials)
% Contrast is active only if it varies outside blank trials.

    nonBlankMask = ~blankTrials;
    if any(nonBlankMask)
        vals = sort(unique(T.contrast(nonBlankMask)));
    else
        vals = sort(unique(T.contrast));
    end
end

function labels = buildStimulusLabels(baseStimuli)
% One label per unique stimulus in the sorted base set.

    nStim = height(baseStimuli);
    labels = cell(nStim, 1);

    for i = 1:nStim
        labels{i} = sprintf(['Dir = %s, tf = %s, sf = %s, ' ...
                             'size = %s, contrast = %s'], ...
            formatScalar(baseStimuli.direction(i)), ...
            formatScalar(baseStimuli.tf(i)), ...
            formatScalar(baseStimuli.sf(i)), ...
            formatScalar(baseStimuli.size(i)), ...
            formatScalar(baseStimuli.contrast(i)));



    end
end

function s = formatScalar(x)
% Compact numeric formatting for labels.

    s = num2str(x, '%.15g');
end
