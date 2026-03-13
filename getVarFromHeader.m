function values = getVarFromHeader(str, tabs, fields)
% str is the header
% tabs is a cell array of strings with tabs names
% fields is a cell array of strings with variable names. Must be as long as
% tabs
% values is a cell array of corresponding values, they will be strings
if ~iscell(fields)
    fields = cell(fields);
end
values = cell(size(fields));

for iField = 1:numel(fields)

    ff = strsplit(str, {' = ', tabs{iField}});
    ind = find(ismember(ff, fields{iField}));
    values{iField} = ff{ind+1};

end
end