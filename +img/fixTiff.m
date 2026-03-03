function fixTiff(badFile, goodFile, newFile)

% Tested for ScanImgae 4.2 tiff files only, might need adjustments for
% other ScanImage versions
%
% This function will fix a corrupted tiff file. It will do so by replacing 
% the metadata of a corrupted file by the metadata of a 'good' file, 
% while keeping the image data of the corrupted file
%
% fixTiff(badFile, goodFile, newFile)
% badFile - full filename of the corrupted file
% goodFile - full filename of the not corrupted file
%            It is a good idea to use a file from the same experiment, so
%            that the format and the size are exactly the same
% newFile - full filename of the new fixed file
%
% TODO: will not work on the last file of the acquisition, where the number
% of saved frames is different
%
% Written by Michael Krumin, 03-Apr-2018


if nargin<3
    error('Please provide all three filenames')
end

if isequal(newFile, goodFile) || isequal(newFile, badFile)
    error('Let''s not overwirte existing files automatically')
end

if isequal(badFile, goodFile)
    error('Corrupted file and not corrupted file have the same name?!')
end

% get metadata from a not corrupted file
fInfo = imfinfo(goodFile);
nPixels = fInfo(1).Width * fInfo(1).Height;
bytesPerImage = nPixels * fInfo(1).BitsPerSample/8;
firstOffset = fInfo(1).Offset;
interFrameOffset = diff([fInfo(1:2).Offset])-bytesPerImage;
headerSamples = interFrameOffset/fInfo(1).BitsPerSample*8;

% and read the data
fid = fopen(goodFile, 'r');
preDataGoodFile = fread(fid, firstOffset, '*uint8');
dataGoodFile = fread(fid, '*int16');
fclose(fid);

% read the corrupted tiff file in a 'binary mode'
fid = fopen(badFile, 'r');
fseek(fid, firstOffset, 'bof');
dataBadFile = fread(fid, '*int16');
fclose(fid);

dataGoodFile = reshape(dataGoodFile, headerSamples+nPixels, []);
dataBadFile = reshape(dataBadFile, headerSamples+nPixels, []);

% replace the corrupted metadata in the bad file by the metadata from the good file
dataBadFile(1:headerSamples, :) = dataGoodFile(1:headerSamples, :);

% save the new tiff file
fid = fopen(newFile, 'w');
fwrite(fid, preDataGoodFile, 'uint8');
fwrite(fid, dataBadFile(:), 'int16');
fclose(fid);
