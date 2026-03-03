function writeAviMovie(frames, outFile, fps, scale, quality)
% writeAviFromFrames  Save a 3-D stack [Ly x Lx x N] to an AVI.
%   writeAviFromFrames(frames, outFile, fps)
%   writeAviFromFrames(frames, outFile, fps, scale)
%   writeAviFromFrames(frames, outFile, fps, scale, quality)
%
% Inputs:
%   frames   : [Ly x Lx x N] numeric, typically int16.
%   outFile  : e.g., 'movie.avi'
%   fps      : desired frame rate, e.g., 20
%   scale    : (optional) how to map intensities -> uint8
%              - 'prctile' (default): 1–99% per whole stack
%              - 'minmax' : min/max of whole stack
%              - [low high] : explicit intensity limits (in native units)
%   quality  : (optional) 0–100 (default 95) for Motion JPEG AVI
%
% Notes:
%   - Output is truecolor (RGB) for broad AVI compatibility.
%   - For enormous stacks, consider chunking to limit memory.

    if nargin < 4 || isempty(scale),   scale = 'prctile'; end
    if nargin < 5 || isempty(quality), quality = 95; end

    frames = double(frames);  % for scaling math

    % Determine scaling limits
    switch lower(string(scale))
        case "prctile"
            lo = prctile(frames(:), 1);
            hi = prctile(frames(:), 99.5);
        case "minmax"
            lo = min(frames(:));
            hi = max(frames(:));
        otherwise
            if isnumeric(scale) && numel(scale) == 2
                lo = scale(1); hi = scale(2);
            else
                error('Unrecognized scale option. Use ''prctile'', ''minmax'', or [low high].');
            end
    end
    if hi <= lo
        hi = lo + eps;
    end

    % Create writer
    v = VideoWriter(outFile, 'Motion JPEG AVI');
    v.FrameRate = fps;
    v.Quality   = quality;
    open(v);
    c = onCleanup(@() close(v));

    [Ly, Lx, N] = size(frames);

    % Normalize to [0,1] -> uint8 -> RGB and write
    for k = 1:N
        g = (frames(:,:,k) - lo) / (hi - lo);
        g = min(max(g, 0), 1);            % clip
        g8 = uint8(round(255 * g));       % grayscale 8-bit
        rgb = repmat(g8, [1, 1, 3]);      % truecolor
        writeVideo(v, rgb);
    end
end
