function validate_grating_inputs(data)
% VALIDATE_GRATING_INPUTS Basic validation of required fields and sizes.

    requiredFields = {'tc', 'amp', 't', 'directions', 'sizes', 'tfs', 'sfs'};
    for iField = 1:numel(requiredFields)
        if ~isfield(data, requiredFields{iField})
            error('Missing required field data.%s', requiredFields{iField});
        end
    end

    if isempty(data.tc) || isempty(data.amp)
        error('data.tc and data.amp must be non-empty.');
    end

    if ~isvector(data.t) || isempty(data.t)
        error('data.t must be a non-empty vector.');
    end

    paramFields = {'directions', 'sizes', 'tfs', 'sfs'};
    paramNames = {'direction', 'size', 'tf', 'sf'};

    paramCounts = zeros(1, numel(paramFields));
    for iField = 1:numel(paramFields)
        values = data.(paramFields{iField});
        if ~isvector(values) || isempty(values)
            error('data.%s must be a non-empty vector.', paramFields{iField});
        end
        paramCounts(iField) = numel(values);
    end

    if isfield(data, 'activeParams') && ~isempty(data.activeParams)
        validNames = paramNames;

        if isstring(data.activeParams)
            activeParams = cellstr(data.activeParams);
        elseif iscell(data.activeParams)
            activeParams = data.activeParams;
        else
            error('data.activeParams must be a string array or a cell array of character vectors.');
        end

        areTextEntries = cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), activeParams);
        if ~all(areTextEntries)
            error('All entries in data.activeParams must be text scalars.');
        end

        activeParams = cellfun(@char, activeParams, 'UniformOutput', false);
        activeParams = unique(activeParams, 'stable');

        for iName = 1:numel(activeParams)
            thisName = activeParams{iName};
            if ~ismember(thisName, validNames)
                error('Invalid entry in data.activeParams: %s', thisName);
            end
        end

        activeMask = ismember(paramNames, activeParams);
        for iParam = 1:numel(paramNames)
            if activeMask(iParam) && paramCounts(iParam) < 2
                error(['data.activeParams marks %s as active, but data.%s contains only ', ...
                       'one value. Active parameters should vary across stimuli.'], ...
                       paramNames{iParam}, paramFields{iParam});
            end
            if ~activeMask(iParam) && paramCounts(iParam) > 1
                error(['data.activeParams does not include %s, but data.%s contains ', ...
                       '%d values. Inactive parameters should have length 1.'], ...
                       paramNames{iParam}, paramFields{iParam}, paramCounts(iParam));
            end
        end
    end

    info = infer_grating_input_layout(data);

    if info.nTime ~= numel(data.t)
        error('Length of data.t (%d) does not match time dimension in data.tc (%d).', ...
            numel(data.t), info.nTime);
    end

    if info.nRep < 2
        error('At least 2 repeats are required for odd/even responsiveness analysis.');
    end

    expectedAmpNumel = info.nStim * info.nRep * info.nNeurons;
    if numel(data.amp) ~= expectedAmpNumel
        error(['data.amp has %d elements, but expected %d elements based on ', ...
               'parameter vectors, repeats, and neurons.'], numel(data.amp), expectedAmpNumel);
    end

    expectedTcNumel = info.nStim * info.nRep * info.nTime * info.nNeurons;
    if numel(data.tc) ~= expectedTcNumel
        error(['data.tc has %d elements, but expected %d elements based on ', ...
               'parameter vectors, repeats, times, and neurons.'], numel(data.tc), expectedTcNumel);
    end
end
