function stack = loadAviRGB(fname, startFrame, endFrame)
% loadAviFramesRGB  Load RGB frames from an AVI into a 4-D uint8 array.
%   stack = loadAviFramesRGB(fname)
%   stack = loadAviFramesRGB(fname, startFrame, endFrame)
%
% Inputs:
%   fname       : path to AVI file.
%   startFrame  : (optional) 1-based index of first frame to load.
%   endFrame    : (optional) 1-based index of last frame to load (inclusive).
%                  If omitted or [], loads to the end of the video.
%
% Output:
%   stack       : uint8 array of size [Height x Width x 3 x N].

    % Open video
    v = VideoReader(fname);
    H = v.Height; W = v.Width;

    % Handle optional inputs
    if nargin < 2 || isempty(startFrame), startFrame = 1; end
    if nargin < 3 || isempty(endFrame),   endFrame   = inf; end

    % Basic checks
    startFrame = max(1, floor(startFrame));
    if ~(isfinite(startFrame) && startFrame >= 1)
        error('startFrame must be a positive integer.');
    end
    if ~(isinf(endFrame) || (isfinite(endFrame) && endFrame >= startFrame))
        error('endFrame must be Inf or an integer >= startFrame.');
    end

    % Seek to start
    v.CurrentTime = (startFrame - 1) / v.FrameRate;

    if isfinite(endFrame)
        % Known count
        N = max(0, floor(endFrame - startFrame + 1));
        stack = zeros(H, W, 3, N, 'uint8');
        n = 0;
        while hasFrame(v) && n < N
            n = n + 1;
            frame = readFrame(v);          % uint8 [H W 3]
            stack(:,:,:,n) = frame;
        end
        stack = stack(:,:,:,1:n);          % trim if short
    else
        % Unknown count: preallocate generously using duration * fps
        estTotal = ceil(v.Duration * v.FrameRate) + 2 - (startFrame - 1);
        estTotal = max(1, estTotal);
        stack = zeros(H, W, 3, estTotal, 'uint8');
        n = 0;
        while hasFrame(v)
            n = n + 1;
            if n > size(stack,4)           % rare, but handle if estimate was low
                stack(:,:,:,end+ceil(n/2)) = uint8(0); %#ok<AGROW>
            end
            stack(:,:,:,n) = readFrame(v);
        end
        stack = stack(:,:,:,1:n);          % trim to actual size
    end
end
