function frames = readBinFrames(fname, Lx, Ly, startFrame, endFrame, machinefmt, headerBytes)
% readBinFrames  Read a range of frames from a raw int16 binary image stack.
%   frames = readBinFrames(fname, Lx, Ly, startFrame, endFrame)
%   frames = readBinFrames(..., machinefmt)  % e.g., 'ieee-le' (default) or 'ieee-be'
%   frames = readBinFrames(..., machinefmt, headerBytes)  % headerBytes default = 0
% If the file was written big-endian, pass 'ieee-be' as machinefmt.
% If there’s a fixed header before the pixel data, pass its size in bytes as headerBytes
% Returns frames with size [Ly x Lx x N], where N = endFrame - startFrame + 1.

    if nargin < 6 || isempty(machinefmt), machinefmt = 'ieee-le'; end
    if nargin < 7 || isempty(headerBytes), headerBytes = 0; end

    assert(startFrame >= 1 && endFrame >= startFrame, 'Invalid frame range.');
    bytesPerPixel = 2; % int16
    frameBytes    = int64(Lx) * int64(Ly) * bytesPerPixel;

    fid = fopen(fname, 'r', machinefmt);
    if fid < 0, error('Failed to open %s', fname); end
    c = onCleanup(@() fclose(fid));

    % Check file size vs. frame size
    fseek(fid, 0, 'eof');
    fileBytes = int64(ftell(fid));
    if fileBytes < headerBytes
        error('headerBytes exceeds file size.');
    end
    dataBytes  = fileBytes - int64(headerBytes);
    totalFrames = double(floor(double(dataBytes) / double(frameBytes)));
    if mod(double(dataBytes), double(frameBytes)) ~= 0
        warning('File size is not a multiple of frame size. Check Lx/Ly or headerBytes.');
    end
    if endFrame > totalFrames
        error('Requested endFrame (%d) exceeds frames in file (%d).', endFrame, totalFrames);
    end

    % Seek to the start frame
    offset = int64(headerBytes) + frameBytes * int64(startFrame - 1);
    if fseek(fid, offset, 'bof') ~= 0
        error('Seek failed. Offset may be too large or file is truncated.');
    end

    % Read the requested frames
    nFrames = endFrame - startFrame + 1;
    nElems  = double(Lx) * double(Ly) * double(nFrames);
    [raw, count] = fread(fid, nElems, 'int16=>int16');
    if count < nElems
        warning('Requested %d frames but only read %.3f frames.', nFrames, count/(Lx*Ly));
    end
    nFramesRead = floor(count / (double(Lx*Ly)));

    raw = raw(1 : double(Lx)*double(Ly)*double(nFramesRead));

    % Reshape: file is X-fastest then Y, then frames -> reshape to [Lx Ly N], then permute
    tmp = reshape(raw, [Lx, Ly, nFramesRead]);
    frames = permute(tmp, [2 1 3]);  % -> [Ly x Lx x N]
end