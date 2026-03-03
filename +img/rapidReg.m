function [regMovie, x, y, dx, dy, target] = rapidReg(movie, target, varargin)
% rapidReg Movie frame registration to a target frame using parallelisation
% and array vectorisation for efficiency
%
%   [REG, X, Y, DX, DY, TARGET] = rapidReg(MOVIE, TARGET,...) register 
%   MOVIE (an (X,Y,T) array) to TARGET (either the target image frame (X,Y)
%   or 'auto', to find one automatically). Returns the registered movie, REG,
%   the x and y coordinates of the registered movies pixels relative to the
%   target (useful for when clipping is used), the translations required to
%   register each frame, DX and DY, and TAGRET, the target frame. 
%   Optionally takes, 'noparallel', meaning use single-threaded codepath 
%   instead of parallel, and 'clip', meaning clip the output image to keep
%   only pixels that have valid info at all times (i.e. haven't translated
%   outside target).

% 2013-07 CB created (heavily plagiarised from Mario Dipoppa's code)

[h, w, nFrames] = size(movie);

fprintf('RapidReg..');

%% Setup
%convert movie data to an appropriate floating point type if necessary
dataType = class(movie);
switch dataType
  case {'int8' 'int16' 'uint8' 'uint16' 'int32' 'uint32'}
    %convert all integer types up to 32-bits to single
    movie = single(movie);
    origTypeFun = str2func(dataType);
  case {'int64'  'uint64'}
    %convert 64-bit integer types to double
    movie = double(movie);
    origTypeFun = str2func(dataType);
  case {'single' 'double'}
    %no conversion
    origTypeFun = str2func(dataType);%DS on 13.12.14
  otherwise
    error('''%s'' is not a recognised data type', dataType);
end

%create a Gaussian filter for filtering frames
hGauss = fspecial('gaussian', [5 5], 1);

%look for flag on whether to use parallel codepath
if any(cell2mat(strfind(varargin, 'nopar')) == 1)
  parallel = false;
else
  parallel = true;
end

if strcmpi(target, 'auto')
  %% Compute the best target frame
  fprintf('finding target..');
%   %first compute a smoothed mean of each frame
%   meanF = smooth(mean(reshape(movie, h*w, nFrames)));
%   %now look in the middle third of the image frames for the minimum
%   fromFrame = round(nFrames*1/3);
%   toFrame = round(nFrames*2/3);
%   [~, idx] = min(meanF(fromFrame:toFrame));
%   minFrame = fromFrame + idx;
%   %Gaussian filter the target image
%   target = imfilter(movie(:,:,minFrame), hGauss, 'same', 'replicate');
%   % AR added display output of which frame number is being used as target
%   fprintf(['target frame is ',num2str(minFrame),'..']);
 target = img.selectTarget(2, movie, min(100, nFrames), 20);
end

%% Fourier transform the movie frames, unfiltered and filtered
fprintf('filtering..');
ftMovie = fft2(movie);
ftFilteredMovie = fft2(imfilter(movie, hGauss, 'same', 'replicate'));

%% Compute required displacement and register each frame
dx = zeros(1, nFrames);
dy = zeros(1, nFrames);
nr = ifftshift((-fix(h/2):ceil(h/2) - 1));
nc = ifftshift((-fix(w/2):ceil(w/2) - 1));
[nc, nr] = meshgrid(nc, nr);
regMovie = zeros(h, w, nFrames, 'like', movie);
fprintf('registering..');

if parallel
  %% Register in parallel

  ppool = gcp('nocreate');
 
  if isempty(ppool)
      createPool = true;
      temporaryPool = gcp;
  else
      createPool = false;
  end
  
  
  try
    %do parallel loops in chunks of data to prevent matlab choking
    % currently hardcoded, which is bad
    chunkSize = 500; %frames
    nChunks = ceil(nFrames/chunkSize);
    for i = 0:(nChunks - 1)
      sidx = i*chunkSize + 1;
      eidx = min((i + 1)*chunkSize, nFrames);
      parfor t = sidx:eidx
        %find the best registration translation
        output = dftregistration(fft2(target), ftFilteredMovie(:,:,t), 20);
        dx(t) = output(4);
        dy(t) = output(3);
        %translate the original (i.e. unfiltered) frame
        ftRegFrame = ftMovie(:,:,t).*exp(sqrt(-1)*2*pi*(-dy(t)*nr/h - dx(t)*nc/w));
        regMovie(:,:,t) = abs(ifft2(ftRegFrame));
      end
    end
    if createPool
      delete(temporaryPool); %close worker pool
    end
  catch ex
    if createPool
      %in case of error, ensure temporary worker pool is closed
      delete(temporaryPool);    
    end
    rethrow(ex)
  end
else
  %% Register sequentially
  for t = 1:nFrames
    %find the best registration translation
    output = dftregistration(fft2(target), ftFilteredMovie(:,:,t), 20);
    dx(t) = output(4);
    dy(t) = output(3);
    %translate the original (i.e. unfiltered) frame
    ftRegFrame = ftMovie(:,:,t).*exp(sqrt(-1)*2*pi*(-dy(t)*nr/h - dx(t)*nc/w));
    regMovie(:,:,t) = abs(ifft2(ftRegFrame));
  end
end

%% If requested, clip the frames to the maximum fully valid region
if any(cell2mat(strfind(varargin, 'clip')) == 1)
  fprintf('clipping..');
  dxMax = max(0, ceil(max(dx)));
  dxMin = min(0, floor(min(dx)));
  dyMax = max(0, ceil(max(dy)));
  dyMin = min(0, floor(min(dy)));
  x = (1 + dxMax):(w + dxMin);
  y = (1 + dyMax):(h + dyMin);
  regMovie = regMovie(y,x,:);
else
  x = 1:w;
  y = 1:h;
end
%% Convert the registered movie to its original type
regMovie = origTypeFun(regMovie);

fprintf('.done\n');

end

