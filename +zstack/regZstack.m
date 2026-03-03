function [info, zStack] = regZstack(info, options)

%% setting the default options

if nargin<2 || ~isfield(options, 'noTiff')
    % do not save the tiff file (not really needed for unregistered data)
    options.noTiff=true;
end

if nargin<2 || ~isfield(options, 'fastSave')
    % set to true if you want to use the fast array saving instead of
    % matlab save() function. Will be much faster, but the data is not
    % compressed (you loose approx 20-25% of space for typical BScope data)
    % make sure you have saveArr and loadArr (writeen by Chris) if you use
    % this fastSave option
    options.fastSave=true;
end

if nargin<2 || ~isfield(options, 'iPlane')
    % by default extract the first plane
    options.iPlane=1;
end

if nargin<2 || ~isfield(options, 'channels')
    % by default exctract all channels
    options.channels = 1:info.nChannels;
end

    regCh = options.registrationChannel;

if nargin<2 || ~isfield(options, 'doBidi')
    doBidi = 0;
else
    doBidi = options.doBidi;
end

%% start the processing
% get all the relevant tif-file names in that folder

if isfield(info, 'folder2pLocal')
    localFiles = true;
    allTiffInfo = dir([info.folder2pLocal, filesep, info.basename2p, '*.tif']);
    if isempty(allTiffInfo)
        localFiles = false;
        fprintf('There were no locally saved tiff files at %s\n', info.folder2pLocal);
        fprintf('Will now load the files from the server (slower)\n');
        allTiffInfo = dir([info.folder2p, filesep, info.basename2p, '*.tif']);
    end
else
    localFiles = false;
    allTiffInfo = dir([info.folder2p, filesep, info.basename2p, '*.tif']);
end

nFiles = length(allTiffInfo);
tiffNames = cell(nFiles, 1);
acqNumber = nan(nFiles, 1);
partNumber = nan(nFiles, 1);
for iFile = 1:nFiles
    tiffNames{iFile} = allTiffInfo(iFile).name;
    try
        acqNumber(iFile) = str2num(tiffNames{iFile}(end-10:end-8));
    catch
        acqNumber(iFile) = str2num(tiffNames{iFile}(end-12:end-10));
        
    end
    partNumber(iFile) = str2num(tiffNames{iFile}(end-6:end-4));
end

acquisitions=unique(acqNumber);
nAcqs=length(acquisitions);

if nAcqs>1
    warning('There is more than one acquisition in this folder, something might be wrong');
end

for iAcq=1:nAcqs
    [partsSorted{iAcq} fileIdxSorted{iAcq}]=sort(partNumber(acqNumber==acquisitions(iAcq)));
end


%%
nParts=length(fileIdxSorted{iAcq});
nPlanes=info.nPlanes;
nFrames=cell(nParts, 1);
nFrameSets = info.nPlanes * info.nChannels;

if nAcqs>1
    acqString = sprintf('_acq%03d', iAcq);
else
    acqString = '';
end

% Creating the '\Processed' sub folder, if it doesn't exist
if ~exist(info.folderZstack, 'dir')
    mkdir(info.folderZstack);
end

%     fids = cell(1, length(options.channels));
for iCh = 1:length(options.channels)
    chID = options.channels(iCh);
    if info.nChannels>1
        chString = sprintf('_channel%01d', chID);
    else
        chString = '';
    end
    frameSetFilename = sprintf('%s%s_%s', ...
        info.basename2p, acqString, chString);
%     fids{iCh} = fopen([fullfile(info.folderProcessed, frameSetFilename), '.bin'], 'w');
end

zStack = cell(1, length(options.channels));
for iPlane = 1: info.nPlanes
    
    initChars = overfprintf(0, sprintf('Extracting plane %d/%d ... ', iPlane, nPlanes));
    
    planeData = cell(1, length(options.channels));
     
    for iAcq=1:nAcqs
         
        nFramesAccum=0;
        
            nMsgChars = 0;

        for iPart=1:nParts
            if localFiles
                filename=fullfile(info.folder2pLocal, tiffNames{fileIdxSorted{iAcq}(iPart)});
            else
                filename=fullfile(info.folder2p, tiffNames{fileIdxSorted{iAcq}(iPart)});
            end
            nMsgChars = overfprintf(nMsgChars, sprintf('Loading part %d/%d ', iPart, nParts));
            if isempty(nFrames{iPart})
                if localFiles
                    nFrames{iPart}=img.nFrames(fullfile(info.folder2pLocal, tiffNames{fileIdxSorted{iAcq}(iPart)}));
                else
                    nFrames{iPart}=img.nFrames(fullfile(info.folder2p, tiffNames{fileIdxSorted{iAcq}(iPart)}));
                end
            end
            
            for iCh = 1:length(options.channels)
                
                chID = options.channels(iCh);
                
                % first frame in the current tiff file which belongs to the
                % current plane and current channel
                iFrameSet = (iPlane-1) * info.nChannels + chID;
                firstFrame=mod(iFrameSet-mod(nFramesAccum, nFrameSets), nFrameSets);
                if firstFrame==0
                    firstFrame=nFrameSets;
                end
                
                if firstFrame < nFrames{iPart}
                    
                    frames2load=(firstFrame:nFrameSets:nFrames{iPart})';
                    
                    lastFrame=frames2load(end);
                    
                    [data, headers]=img.loadFrames(filename, firstFrame, lastFrame, nFrameSets);
                    [h, w, nf] = size(data);
                    
                    planeData{iCh} =cat(3,  planeData{iCh}, single(data));
                    
%                     dataPrecision = class(data);
                 
                end
            end
            nFramesAccum=nFramesAccum+nFrames{iPart};

        end
                         overfprintf(0, ' ');

    end

    if doBidi
%     BiDiPhase = img.BiDiPhaseOffsets(planeData{regCh});
        BiDiPhase = img.BiDiPhaseOffsets(planeData{1});
    end
    targetFrame = mean(planeData{regCh}, 3);
    
    [~, ~, ~, dx, dy] = img.rapidReg(planeData{regCh}, ...
        targetFrame);
    
    [~,~, nFrCh] = cellfun(@size, planeData, 'UniformOutput', false);
    
    validFrames = 1:min(nFrCh{:});
    
    
    for iCh = 1:length(options.channels)
        
        if doBidi
        planeData{iCh} = img.ShiftBiDi(BiDiPhase, planeData{iCh}(:,:,validFrames), h, w);
        else
            planeData{iCh} = planeData{iCh}(:,:,validFrames);
        end
        [regPlaneData{iCh}, ~, ~]=img.translate(planeData{iCh}, dx(validFrames), dy(validFrames));
        
        zStack{iCh} = cat(3, zStack{iCh}, mean(regPlaneData{iCh},3));
        
        %                 fwrite(fids{iCh}, zStack{regCh}, dataPrecision);
    end
    
    
    
    clear targetFrame planeData regPlaneData

            
        fprintf('Plane %d complete.\n', iPlane);
end

nCh = length(options.channels);



if nCh ==1
    
%     zStack{iCh} = zstack.registerStackColumn(zStack{iCh});
    
    zStack{iCh} = uint16(mat2gray(zStack{iCh})*(2^16-1));
    
   saveastiff(zStack{iCh}, fullfile(info.folderZstack , [info.expRef, '_zStackMean.tif']));
else
    
%     [zStack{2}, zStack{1}] = zstack.registerStackColumn(zStack{2}, zStack{1}, regCh);

    for iCh = 1:nCh
        

        zStack{iCh} = uint16(mat2gray(zStack{iCh})*(2^16-1));
        switch iCh
            case 1
        saveastiff(zStack{iCh}, fullfile(info.folderZstack , [info.expRef, '_zStackMean_G.tif']));
            case 2
       saveastiff(zStack{iCh}, fullfile(info.folderZstack , [info.expRef, '_zStackMean_R.tif']));

        end
    end
end



end

