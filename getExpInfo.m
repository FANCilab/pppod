function info= getExpInfo(animal, expDate, exp, get_2p_info)

% This function populates the initial info structure, which is needed to run many
% other 2p analysis

% input arguments' format:
% animal - string with the subject (animal) name
% expDate - string with the date (as in an example below)
% exp - experiment number
%
% Example:
% info= getExpInfo('M140108_MK005', '2014-03-05', 1)


if nargin < 4
    get_2p_info = 1;
end


info.subject=animal;
info.date=expDate;
info.exp=exp;
info.expRef= sprintf('%s_%s_%d', info.subject,info.date, info.exp);

dataFolders = { 'D:\Data\suite2p\','Z:\Data\2P','\\10.233.25.135\Data\2P'}; %edit this so that the first one is your local data folder
for k = 1:length(dataFolders)
    folder = fullfile(dataFolders{k}, info.subject, info.date, num2str(info.exp));
    if exist(folder, 'dir') ~= 0
        info.folder2p = folder;
        thisServer = dataFolders{k};
        break
    end
end
info.folderZstack = fullfile(info.folder2p, 'zStack');
info.basename2p=sprintf('%s_%s_%d_2P', info.subject, info.date, info.exp);

rawDataFolders = {'Z:\Data\2P','\\10.233.25.135\Data\2P'};
for k = 1:length(dataFolders)
    folder = fullfile(rawDataFolders{k}, info.subject, info.date, num2str(info.exp));
    if exist(folder, 'dir') ~= 0
        info.root_storage = rawDataFolders{k};
        info.folder2praw = folder;
        thisServer = dataFolders{k};
        break
    end
end

% end

% get the number of planes from the Timeline data (piezo trace and frame trigger).
% If this fails (usually it does), then the tiff header is used, but this
% number also might be wrong if something was wrong with the acquisition (should not happen).
% !!! Always check your data visually after doing the preliminary analysis !!!

if get_2p_info

    try
        % try
            % allTiffInfo = dir([info.folder2pLocal, filesep, info.basename2p, '*.tif']);
            % tiffName = allTiffInfo(1).name;
            % filename=fullfile(info.folder2pLocal, tiffName);
            % [~, header]=loadFramesBuff(filename, 1, 1, 1);
        % catch
            fprintf('Getting the tiff from the server (local tiffs do not exist)...\n');
            allTiffInfo = dir([info.folder2praw, filesep, info.basename2p,'*.tif']);
            tiffName = allTiffInfo(1).name;
            filename=fullfile(info.folder2praw, tiffName);
            [~, header]= tiff.loadFramesBuff(filename, 1, 1, 1);
        % end
        % getting some parameters from the header
        hh=header{1};

        values = getVarFromHeader(hh, ...
            {'\nSI.hRoiManager.','\nSI.hRoiManager.','\nSI.hRoiManager.','\nSI.hRoiManager.', '\nSI.hRoiManager.','\nSI.hChannels.', '\nSI.hStackManager.', '\nSI.hStackManager.'},...
            {'scanZoomFactor',  'linesPerFrame',    'pixelsPerLine',    'scanFrameRate',   'scanVolumeRate',  'channelSave',      'numSlices', 'zs'});
  

        info.scanZoomFactor = str2double(values{1});
        info.zoomFactor = str2double(values{1});
        info.scanLinesPerFrame = str2double(values{2});
        info.scanPixelsPerLine = str2double(values{3});
        info.frameRate = str2double(values{4});
        info.volumeRate = str2double(values{5});
        info.nChannels = numel(str2num(values{6}));

        info.chData(1).color = 'green';
        if info.nChannels ==2
            info.chData(2).color = 'red';
        end

        info.nPlanes = str2double(values{7});
        info.zs = str2double(values{8});

        if info.nPlanes ~=numel(info.zs)
            info.nPlanes = numel(info.zs);
        end

    catch
        warning('NO IMAGING DATA FOUND, returning basic exp info')
    end
end
end





