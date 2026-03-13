function out = parse_experimental_protocol(filename)
%PARSE_EXPERIMENTAL_PROTOCOL Parse a stimulus protocol spreadsheet.
%
% out = parse_experimental_protocol(filename)
%
% Expected parameters are:
%   direction, tf, sf, size, contrast
%
% The function:
%   1) finds which parameters vary across stimuli (activePars)
%   2) identifies the unique base set of stimuli
%   3) sorts that base set to match the requested example ordering:
%        direction changes fastest, then size, then tf, then sf
%      (therefore the lexicographic sort keys are reversed for sortrows)
%   4) maps each row of the original experiment back to the sorted stimulus idx
%
% Returned fields:
%   out.T                         standardized input table
%   out.activePars                cell array of varying parameter names
%   out.uniqueValues              struct with unique values for each parameter
%   out.nValuesPerPar             struct with number of values for each parameter
%   out.baseStimuli               sorted unique stimulus combinations
%   out.repeatCount               repetitions per unique stimulus
%   out.stimulusIdxSequence       1-based idx for each presented stimulus row
%   out.trialStimIdx              [nTrials x nBaseStimuli] if contiguous trials exist
%   out.nTrials                   number of contiguous trials, when detected
%   out.hasCompleteContiguousTrials logical flag
%
% Example:
%   out = parse_experimental_protocol('metadata.csv');
%   baseStimuli = out.baseStimuli;
%   seq        = out.stimulusIdxSequence;
%   trialIdx   = out.trialStimIdx;

    parNames = {'direction','tf','sf','size','contrast'};
    T = readProtocolTable(filename, parNames);
    out = struct();
    out.T = T;

    %% 1) Identify active parameters
    out.activePars = {};
    out.uniqueValues = struct();
    out.nValuesPerPar = struct();

    fprintf('Changed parameters (activePars):\n');
    for i = 1:numel(parNames)
        p = parNames{i};
        vals = sort(unique(T.(p)));
        out.uniqueValues.(p) = vals;
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
        vals = out.uniqueValues.(p);
        if numel(vals) == 1
            fprintf('  %s: %s\n', p, mat2str(vals(:).'));
        end
    end

    %% 2) Unique base stimulus set
    baseStimuli = unique(T(:, parNames), 'rows');

    % To match the user's example, direction varies fastest, then size,
    % then tf, then sf. sortrows() uses the first key as the slowest/
    % outermost grouping, so we reverse those keys here.
    sortKeys = {'sf','tf','size','direction'};

    % contrast was not specified in the requested ordering; use it only as
    % a final tie-breaker if it varies.
    if out.nValuesPerPar.contrast > 1
        sortKeys{end+1} = 'contrast';
    end

    baseStimuli = sortrows(baseStimuli, sortKeys);

    %% 3) Map each presented stimulus row to the sorted base stimulus idx
    [isMember, stimulusIdxSequence] = ismember(T(:, parNames), baseStimuli, 'rows');
    if ~all(isMember)
        error('Could not map all rows back to the sorted base stimulus set.');
    end

    repeatCount = accumarray(stimulusIdxSequence, 1, [height(baseStimuli) 1]);

    out.baseStimuli = baseStimuli;
    out.repeatCount = repeatCount;
    out.stimulusIdxSequence = stimulusIdxSequence;

    fprintf('\nBasic stimulus set:\n');
    fprintf('  %d unique stimulus combinations\n', height(baseStimuli));
    if all(repeatCount == repeatCount(1))
        fprintf('  each unique stimulus is repeated %d times\n', repeatCount(1));
    else
        fprintf('  repeat counts vary across stimuli: min = %d, max = %d\n', ...
            min(repeatCount), max(repeatCount));
    end

    %% 4) Detect contiguous trials, if the rows split cleanly into repeats
    out.trialStimIdx = [];
    out.nTrials = NaN;
    out.hasCompleteContiguousTrials = false;

    nBase = height(baseStimuli);
    nRows = numel(stimulusIdxSequence);

    if mod(nRows, nBase) == 0
        candidateTrialStimIdx = reshape(stimulusIdxSequence, nBase, []).';
        expected = repmat(1:nBase, size(candidateTrialStimIdx, 1), 1);
        isCompleteTrial = all(sort(candidateTrialStimIdx, 2) == expected, 2);

        if all(isCompleteTrial)
            out.trialStimIdx = candidateTrialStimIdx;
            out.nTrials = size(candidateTrialStimIdx, 1);
            out.hasCompleteContiguousTrials = true;
            fprintf('  detected %d contiguous trials of %d stimuli each\n', ...
                out.nTrials, nBase);
        else
            fprintf('  rows do not split into complete contiguous trials\n');
        end
    else
        fprintf('  total row count is not an integer multiple of the base set size\n');
    end

    nShow = min(20, numel(stimulusIdxSequence));
    fprintf('\nFirst %d stimulus idx in presentation order:\n', nShow);
    fprintf('  %s\n', mat2str(stimulusIdxSequence(1:nShow).'));
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
