clear;

%% set paths
addpath(genpath('D:\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\Code\retinotopy'));
addpath(genpath('D:\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\Code\FedBox'));

addpath(genpath('Z:\Data\2P\'));

ops.root_storage = 'Z:\Data\2P';

ops.mouse_name = 'NM023';
ops.date = '20251224';
ops.exp_n = 4;

ops.n_channels = 1;
ops.fun_channel = 1;
ops.n_planes = 1;
ops.n_tiles = 10;
ops.i_plane = 1;

ops.root_dir = fullfile(ops.root_storage, ops.mouse_name, ops.date, num2str(ops.exp_n));

%% load imaging data

ops.fsroot = [];
ffile = dir(fullfile(ops.root_dir, '*.tif'));
fname = struct2cell(ffile);
fname = fname(1,:)';
[~,index] = sort_nat(fname);
ops.fsroot = ffile(index);
% ops.fsroot = cat(1, ops.fsroot, ...
%     dir(fullfile(ops.root_dir, '*.tiff')));
for k = 1:length(ops.fsroot)
    ops.fsroot(k).name = fullfile(ops.root_dir, ops.fsroot(k).name);
end

fs = ops.fsroot;

% check if there are tiffs in directory
try
    IMG = loadFramesBuff(fs(1).name, 1, 1, 1);
catch
    error('could not find any tif or tiff, check your path');
end
[Ly, Lx, ~, ~] = size(IMG);
ops.Ly = Ly;
ops.Lx = Lx;

tic;

iplane0 = 1:1:ops.n_planes; % identity of planes for first frames in tiff file
nbytes = 0;
frames_done = 0;
mov = [];
mimg = zeros(Ly, Lx);
for j = 1:length(fs)
    % only compute number of frames if size of tiff is different
    % from previous tiff
    if abs(nbytes - fs(j).bytes)>1e3
        nbytes = fs(j).bytes;
        % nFr= nFrames(fs(j).name);
       nFr= nFramesLFR(fs(j).name);
    end

    if mod(nFr, ops.n_channels) ~= 0
        fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
    end

    % identity of planes for first frames in tiff file
    iplane0 = mod(iplane0-1, ops.n_planes) + 1;
    % only load frames of registration channel
    start_frame = (iplane0(ops.i_plane)-1)*ops.n_channels+ops.fun_channel;
    %load the data
    data = loadFramesBuff(fs(j).name, start_frame, nFr, ops.n_planes*ops.n_channels, []);
    n_frames_current = size(data,3);

    mimg = mean(data, 3)*nFr +mimg;
    % mov(:,:, (frames_done+1):(frames_done+n_frames_current)) = data;

    data = reshape(data, Ly*Lx, n_frames_current);
    
    tileX = round(linspace(1,Lx, ops.n_tiles +1));
    tileY = round(linspace(1,Ly, ops.n_tiles+1));
      for iTlX = 1: ops.n_tiles
        for iTlY = 1: ops.n_tiles
            [xPx, yPx] = meshgrid(tileX(iTlX):tileX(iTlX+1), tileY(iTlY):tileY(iTlY+1));
            iPx = sub2ind([Ly,Lx], yPx(:), xPx(:));

            F = single(data(iPx, :));
            tileF(iTlY, iTlX,  (frames_done+1):(frames_done+n_frames_current)) = mean(F, 1);

        end
      end


    fprintf('Tiff %d of %d done in time %2.2f \n', j, numel(fs), toc)

    %keyboard;
    iplane0 = iplane0 - nFr/ops.n_channels;
    frames_done = frames_done+n_frames_current;
end
mimg = mimg/frames_done;
%% load metadata

event = bonsai.load_events(ops);

allFrameTimes = event.frame.on(ops.i_plane:ops.n_planes:end);
stimFrames = event.sparse_noise.frames;
stimFrameTimes = event.frame.on';
stimTimes.onset = 0;

tileF = tileF(:,:, 1:numel(allFrameTimes));

stimPosition = [0 120 -40 40]; % boundaris of the stim grid in visual degrees (##HARDCODED, FIX)

[rangeY, rangeX, ~] = size(stimFrames); % number of squares in the grid

%%

for iTlX = 1: ops.n_tiles
    for iTlY = 1: ops.n_tiles
                
        [ R{iTlY, iTlX}, tKernel(iTlY, iTlX, :), tileFit(iTlX, iTlY, :), ...
            rfContours{iTlY, iTlX}, periT]  = medianSTRF(squeeze(tileF(iTlY, iTlX, :)), ...
            allFrameTimes, stimFrames, stimFrameTimes, stimTimes, 'any', [], 0);
%         
        retC(iTlY, iTlX, :) = tileFit(iTlX, iTlY,[2, 1]);%nTile * 2
    end
end

figure('Color', 'White', 'Position', [152 223 1080 1001]);
for iTlX = 1: ops.n_tiles
    for iTlY = 1: ops.n_tiles
        
        subplot(ops.n_tiles,ops.n_tiles, (iTlY-1)*ops.n_tiles + iTlX)
        cmap = BlueWhiteRed(100, 0.2);
        maxval = max(cellfun(@(x) max(x(:)), R));
        imagesc(R{iTlY, iTlX}/max(maxval)); axis image; colormap(cmap); hold on
        caxis([-1 1])
        axis off
        plot(tileFit(iTlX, iTlY, 2), tileFit(iTlX, iTlY, 1), '*')

    end
end
    

% convert maps to visual degrees
mapX= retC(:,:, 1);
mapX = stimPosition(1) + (stimPosition(2) - stimPosition(1))*(mapX - 1)/rangeX;
mapY= retC(:,:,2);
mapY = stimPosition(3) + (stimPosition(4) - stimPosition(3))*(mapY - 1)/rangeY;
% smooth maps
smapX = imgaussfilt(mapX, 0.2);
% smapY = rescaleVec(imgaussfilt(mapY, 1), max(mapY(:)), min(mapY(:)));
smapY = imgaussfilt(mapY, 0.2); 

% interpolate to px size
Cx = round(linspace(1,Lx, ops.n_tiles));
Cy = round(linspace(1,Ly, ops.n_tiles));
retX = interp2( Cx,  Cy,  smapX, (1:Lx)', 1:Ly);
retY = interp2( Cx ,  Cy, smapY, (1:Lx)', 1:Ly );


% summary plot 
figure('Color', 'White', 'Position', [1307 386 1092 655]);
ax1 = subplot(2,3,1);
imagesc(mimg, 'Parent', ax1); axis(ax1, 'image'); colormap(ax1, 'gray'); colorbar(ax1);

ax2 = subplot(2,3,2);
imagesc(mapX, 'Parent', ax2); axis(ax2, 'image'); colormap(ax2, 'jet'); colorbar(ax2);
caxis([0 120]);

ax3 = subplot(2,3,3);
imagesc(mapY, 'Parent', ax3); axis(ax3, 'image'); colormap(ax3, 'jet'); colorbar(ax3);
caxis([-40 40]);

ax4 = subplot(2,3,4);
imagesc(mimg, 'Parent', ax4); axis(ax4, 'image'); colormap(ax4, 'gray'); colorbar(ax4);

ax5 = subplot(2,3,5);
imagesc(retX, 'Parent', ax5); axis(ax5, 'image'); colormap(ax5, 'jet'); colorbar(ax5);
caxis([0 120]);

ax6 = subplot(2,3,6);
imagesc(retY, 'Parent', ax6); axis(ax6, 'image'); colormap(ax6, 'jet'); colorbar(ax6);
caxis([-40 40]);
