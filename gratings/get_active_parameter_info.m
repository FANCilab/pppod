function paramInfos = get_active_parameter_info(canon)
% GET_ACTIVE_PARAMETER_INFO Return metadata for all active stimulus parameters.

    catalog = local_catalog();
    activeNames = canon.activeParams;

    if isempty(activeNames)
        paramInfos = struct('name', {}, 'prettyName', {}, 'values', {}, ...
            'dim', {}, 'resultField', {}, 'timecourseFolderField', {});
        return
    end

    paramInfos = repmat(struct( ...
        'name', '', ...
        'prettyName', '', ...
        'values', [], ...
        'dim', [], ...
        'resultField', '', ...
        'timecourseFolderField', ''), 1, numel(activeNames));

    for iParam = 1:numel(activeNames)
        idx = find(strcmp(activeNames{iParam}, {catalog.name}), 1, 'first');
        if isempty(idx)
            error('Unknown active parameter: %s', activeNames{iParam});
        end

        paramInfos(iParam).name = catalog(idx).name;
        paramInfos(iParam).prettyName = catalog(idx).prettyName;
        paramInfos(iParam).values = canon.(catalog(idx).valueField);
        paramInfos(iParam).dim = catalog(idx).dim;
        paramInfos(iParam).resultField = catalog(idx).resultField;
        paramInfos(iParam).timecourseFolderField = catalog(idx).timecourseFolderField;
    end
end

function catalog = local_catalog()
% LOCAL_CATALOG Canonical metadata for supported stimulus parameters.

    catalog = struct( ...
        'name', {'direction', 'size', 'tf', 'sf'}, ...
        'prettyName', {'Direction', 'Size', 'TF', 'SF'}, ...
        'valueField', {'directions', 'sizes', 'tfs', 'sfs'}, ...
        'dim', num2cell(1:4), ...
        'resultField', {'direction', 'size', 'tf', 'sf'}, ...
        'timecourseFolderField', {'directionTimecourses', 'sizeTimecourses', 'tfTimecourses', 'sfTimecourses'});
end
