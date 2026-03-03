function writeAviRGBMovie(frames, outFile, fps, scale, quality)
% writeAviFromRGBFrames  Save an RGB stack to an AVI.
%   writeAviFromRGBFrames(frames, outFile, fps)
%   writeAviFromRGBFrames(frames, outFile, fps, scale)
%   writeAviFromRGBFrames(frames, outFile, fps, scale, quality)
%
% Inputs:
%   frames  : RGB frames as [H x W x 3 x N]. Also accepts:
%             - [H x W x 3] (single frame)
%             - [H x W x N] (grayscale) -> replicated to RGB
%             - [H x W x N x 3] -> auto permuted to [H x W x 3 x N]
%             Types: uint8 (preferred), uint16/int16, single/double.
%   outFile : e.g., 'movie_rgb.avi'
%   fps     : frames per second (e.g., 20)
%   scale   : (optional) how to map non-uint8 data -> uint8
%             - 'auto'  (default): uint8 pass-through; floats in [0,1]→*255;
%                                integers → per-channel min/max
%             - 'minmax': per-channel min/max across the stack
%             - [low high]: explicit intensity window (applied to all 3 channels)
%             - 'none' : only clip to [0,255] then cast (use if already 0–255 range)
%   quality : (optional) 0–100 for Motion JPEG (default 95)

    if nargin < 4 || isempty(scale),   scale = 'auto'; end
    if nargin < 5 || isempty(quality), quality = 95;   end

    % Normalize dimensions to [H W 3 N]
    switch ndims(frames)
        case 2
            error('Frames must be at least 3-D (H x W x 3 or H x W x N).');
        case 3
            if size(frames,3) == 3
                frames = reshape(frames, size(frames,1), size(frames,2), 3, 1);
            else
                % Grayscale [H W N] -> replicate to RGB
                frames = repmat(frames, 1, 1, 3, 1);
                frames = permute(frames, [1 2 3 4]); % already correct shape now
            end
        case 4
            if size(frames,3) ~= 3 && size(frames,4) == 3
                frames = permute(frames, [1 2 4 3]); % [H W N 3] -> [H W 3 N]
            elseif size(frames,3) ~= 3
                error('Expected 3 color channels on dim 3.');
            end
        otherwise
            error('Unsupported input dimensions.');
    end

    [H, W, C, N] = size(frames);
    if C ~= 3, error('Expected RGB with 3 channels.'); end

    % Convert to uint8 according to 'scale'
    frames = convertToUint8RGB(frames, scale);

    % Set up writer (widely compatible)
    v = VideoWriter(outFile, 'Motion JPEG AVI');
    v.FrameRate = fps;
    v.Quality   = quality;
    open(v);
    c = onCleanup(@() close(v));

    % Stream frames
    for k = 1:N
        writeVideo(v, frames(:,:,:,k));
    end
end

% -------- helper --------
function out = convertToUint8RGB(in, scale)
    cls = class(in);
    switch cls
        case 'uint8'
            if strcmpi(scale,'none') || strcmpi(scale,'auto') || isempty(scale)
                out = in;
                return;
            end
        otherwise
            % continue to handle scaling/casting
    end

    in = double(in);  % do math in double
    [H,W,~,N] = size(in);

    if isnumeric(scale) && numel(scale) == 2
        lo = scale(1); hi = scale(2);
        if hi <= lo, hi = lo + eps; end
        out = (in - lo) / (hi - lo);
    else
        switch lower(string(scale))
            case "none"
                % Only clip to [0,255] then cast
                out = in / 255;
            case "auto"
                if isfloat(in)
                    % If looks like 0..1, keep; else minmax per channel
                    mn = min(in,[],'all'); mx = max(in,[],'all');
                    if mn >= 0 && mx <= 1
                        out = in;
                    else
                        out = perChannelMinMax(in);
                    end
                else
                    out = perChannelMinMax(in);
                end
            case "minmax"
                out = perChannelMinMax(in);
            otherwise
                error('Unrecognized scale option. Use ''auto'', ''minmax'', ''none'', or [low high].');
        end
    end

    % Clip to [0,1], scale to [0,255] and cast
    out = uint8(round(255 * min(max(out,0),1)));

    % Nested: per-channel min/max scaling over whole stack
    function y = perChannelMinMax(x)
        y = zeros(size(x));
        for ch = 1:3
            xs = x(:,:,ch,:);
            lo = min(xs,[],'all');
            hi = max(xs,[],'all');
            if hi <= lo, hi = lo + eps; end
            y(:,:,ch,:) = (x(:,:,ch,:) - lo) / (hi - lo);
        end
    end
end
