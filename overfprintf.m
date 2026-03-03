function n = overfprintf(n, varargin)
%overfprintf Summary of this function goes here
%   Detailed explanation goes here

if isfloat(varargin{1}) && isscalar(varargin{1})
  fileID = varargin{1};
  varargin = varargin(2:end);
else
  fileID = 1; %default file ID is the standard output
end

n = fprintf(fileID, [repmat('\b', 1, n) varargin{1}], varargin{2:end}) - n;


end

