function [avgTarget,choosenTargets] = selectTarget(OptionTg,inpMovie,nRandTarg,nRandom2,subpix, gaussianFilter)

% chooseTarget choose the optimal average target for image registration.
% First it extract N random frames from the movie(s) (or receive them as an
% input). Second it computes the average distance of each frame of this
% ensemble with M random frame of the same enseble. It selects the frame
% with the minimal distance from the others and then finds all the frames
% that are aligned within 1/3 of a pixel from this one. The optimal target
% is computed as the average of all these aligned frames.
%
% OptionTg: 1 = input a selection of target frames as a (Y,X,#targets)array
%           2 = input the whole movie (Y,X,#frames) from which randomly
%               select possible targets
%           3 = input a cell of several movies (Y,X,#frames) to be
%               registered together
% inpMovie: targets or movie(s) from which the algoritm extracts the
%             optimal target
% nRandTarg: number of frames selected randomly as possible targets
% nRandom2: number of random targets to be compared
%           with each target of the ensemble
%
% 2014-05-28 - MD Created
% 2014-06-12 - MK modified
% 2014-06-12 - MD added the output choosenTargets corresponding to the
% targets choosen for the average
% 2014-11-27 - MK added the gaussianFilter input argument

if nargin<6 || isempty(gaussianFilter)
    % if not provided, generate a default gaussian filter
    filtWindow=5;
    filtWidth=1; 
    gaussianFilter = fspecial('gaussian',[filtWindow filtWindow],filtWidth);
end

if nargin <5
    subpix = 3;
end
%% first build a targets matrix

if OptionTg==1
    %targets are provided as input
    Targets=imfilter(inpMovie, gaussianFilter, 'replicate');
    
elseif OptionTg==2
    %n targets are randomly choosen from the movie
    [~,~,nFrames]=size(inpMovie);
    randomFrames=randperm(nFrames, nRandTarg);
    Targets = imfilter(inpMovie(:,:,randomFrames), gaussianFilter, 'replicate');
    
elseif OptionTg==3
    %n targets are randomly choosen from each of m movies
    nExps=length(inpMovie);
    [nY,nX,~]=size(inpMovie{1});
    Targets=zeros(nY,nX,nRandTarg*nExps);
    randomFrames=cell(1,nExps);
    for iExp=1:nExps
        [~,~,nFrames]=size(inpMovie{iExp});
        idx = (iExp-1)*nRandTarg+1:iExp*nRandTarg;
        randomFrames{iExp}=randperm(nFrames, nRandTarg);
        Targets(:,:,idx) = ...
        imfilter(inpMovie{iExp}(:,:,randomFrames{iExp}), gaussianFilter, 'replicate');
    end
    
end

%% now do the target selection

nTargets=size(Targets,3);

Drtot=zeros(1,nTargets);

for iTarget=1:nTargets
    nChars = fprintf('Checking target %d/%d', iTarget, nTargets);
    frameInd = randperm(nTargets-1, nRandom2);
    frameInd(frameInd >= nTargets) = frameInd(frameInd >= nTargets) + 1;
    for k=1:nRandom2   
        jTarget = frameInd(k);
        output = dftregistration(fft2(single(squeeze(Targets(:,:,iTarget)))), fft2(single(squeeze(Targets(:,:,jTarget)))),subpix);
        Dr=norm(output(3:4));
        Drtot(iTarget)=Drtot(iTarget)+Dr;
    end
    if (iTarget<nTargets)
        fprintf(repmat('\b', 1, nChars));
    end
end
fprintf('\n');

Drtot=Drtot./nRandom2;


%choose the target that has the minimal distance from all the other targets
%choiceTarget=squeeze(Targets(:,:,minTg));

fprintf('Calculating the average target frame...\n');
[~,minTg]=min(Drtot);

dy=nan(1,nTargets);
dx=nan(1,nTargets);
fftChosen = fft2(single(squeeze(Targets(:,:,minTg))));
for iTarget=1:nTargets
    output = dftregistration(fft2(single(squeeze(Targets(:,:,iTarget)))), fftChosen, subpix);
    dy(iTarget)=output(3);
    dx(iTarget)=output(4);
end

%find frames that are distant less than 1/subpix of a pixel in x and y from the
%optimal frame and average them
alignedTargets=find((abs(dx)<1/subpix).*(abs(dy)<1/subpix));



avgTarget=mean(Targets(:,:,alignedTargets),3);
fprintf('Average target computed from %d/%d frames\n', length(alignedTargets),  nTargets)

%% 

% save which frames were choosen
if OptionTg==1
    choosenTargets=alignedTargets;
elseif OptionTg==2
    choosenTargets=randomFrames(alignedTargets);
elseif OptionTg==3
    choosenTargets=cell(1,nExps);
    idExp=floor((alignedTargets-1)/nRandTarg)+1;
    iRandPerm=alignedTargets-(idExp-1)*nRandTarg;
    for iExp=1:nExps
        choosenTargets{iExp}=randomFrames{iExp}(iRandPerm(idExp==iExp));
    end
end



