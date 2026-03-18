function db = build_db(varargin)

nDb = numel(varargin{1});
nFields = numel(varargin);

for iDb = 1:nDb
    for iv = 1:nFields
        
        fieldname = inputname(iv);
        db(iDb).(fieldname) = varargin{iv}{iDb};
        
    end
end
end