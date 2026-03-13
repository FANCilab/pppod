% Find the number of frames in a tiff file
%
function n = nFramesLFR(tiff_path)
    
headers = imfinfo(tiff_path, 'tif');
n = numel(headers);

end

