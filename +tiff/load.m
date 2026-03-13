function stack = load

%load the .tiff stack

[file, path] = uigetfile;

cd(path);

headers = imfinfo(fullfile(path, file), 'tif');
nFrames = numel(headers);

stack = [];

for iFrame = 1: nFrames
stack(:,:, iFrame) = imread(fullfile(path, file), 'tif', 'Index', iFrame);
end

end

