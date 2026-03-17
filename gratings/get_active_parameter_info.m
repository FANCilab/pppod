function paramInfos = get_active_parameter_info(canon)
% GET_ACTIVE_PARAMETER_INFO Return metadata for all active stimulus parameters.
%
% Orientation is added automatically whenever direction is an active
% parameter. It is a derived parameter with 180-degree periodicity and is
% computed from direction as mod(direction + 90, 180).

    catalog = local_catalog(canon);
    requestedNames = canon.activeParams;

    if canon.hasOrientation
        insertAfterDirection = find(strcmp(requestedNames, 'direction'), 1, 'first');
        if isempty(insertAfterDirection)
            requestedNames{end+1} = 'orientation';
        else
            requestedNames = [requestedNames(1:insertAfterDirection), {'orientation'}, requestedNames(insertAfterDirection+1:end)]; %#ok<AGROW>
        end
    end

    if isempty(requestedNames)
        paramInfos = struct('name', {}, 'prettyName', {}, 'values', {}, ...
            'dim', {}, 'resultField', {}, 'timecourseFolderField', {}, 'isDerived', {});
        return
    end

    paramInfos = repmat(struct( ...
        'name', '', ...
        'prettyName', '', ...
        'values', [], ...
        'dim', [], ...
        'resultField', '', ...
        'timecourseFolderField', '', ...
        'isDerived', false), 1, numel(requestedNames));

    for iParam = 1:numel(requestedNames)
        idx = find(strcmp(requestedNames{iParam}, {catalog.name}), 1, 'first');
        if isempty(idx)
            error('Unknown active parameter: %s', requestedNames{iParam});
        end

        paramInfos(iParam).name = catalog(idx).name;
        paramInfos(iParam).prettyName = catalog(idx).prettyName;
        paramInfos(iParam).values = catalog(idx).values;
        paramInfos(iParam).dim = catalog(idx).dim;
        paramInfos(iParam).resultField = catalog(idx).resultField;
        paramInfos(iParam).timecourseFolderField = catalog(idx).timecourseFolderField;
        paramInfos(iParam).isDerived = catalog(idx).isDerived;
    end
end

function catalog = local_catalog(canon)
% LOCAL_CATALOG Canonical metadata for supported stimulus parameters.

    catalog = struct( ...
        'name', {'direction', 'orientation', 'size', 'tf', 'sf'}, ...
        'prettyName', {'Direction', 'Orientation', 'Size', 'TF', 'SF'}, ...
        'values', {canon.directions, canon.orientations, canon.sizes, canon.tfs, canon.sfs}, ...
        'dim', {1, 0, 2, 3, 4}, ...
        'resultField', {'direction', 'orientation', 'size', 'tf', 'sf'}, ...
        'timecourseFolderField', {'directionTimecourses', 'orientationTimecourses', 'sizeTimecourses', 'tfTimecourses', 'sfTimecourses'}, ...
        'isDerived', {false, true, false, false, false});
end
