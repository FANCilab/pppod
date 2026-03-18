function [responses, reorderedType, blankResp, bStd] = getStimulusTrials(F, frameTimes, stimInfo, Types)

if nargin <4
    Types = {'oriMultiSfTf', 'oris', 'plaids',  'SF', 'TF', 'nat movies', 'size', 'spontaneous', 'sparseNoise'};
end

nStimsTypes = numel({stimInfo.stimType});
[nT, nN] =  size(F);
framesXexp = [0, cumsum(cellfun(@numel, frameTimes))];

for iType = 1:nStimsTypes
    thisF = F(framesXexp(iType)+1:framesXexp(iType+1), :); %
    thisFrameTimes = frameTimes{iType};
    
    %     thisF = npRet.my_conv(thisF',1)';
    
    
    switch stimInfo(iType).stimType
        
        case {'oriMultiSfTf','plaids', 'oris', 'SF', 'TF'}
            stimTimes = stimInfo(iType).stimTimes;
            stimSeq = stimInfo(iType).stimSequence;
            nStim = numel(unique(stimSeq.seq));
            if nStim >13
                if nStim> 40
                    nBlanks = 3;
                else
                    nBlanks = 2;
                    
                end
            else
                nBlanks = 1;
            end
            nStim = nStim-nBlanks;
            
            nRep = sum(stimSeq.seq == 1); % assuming all the stims come with the same nRep
            
            %         case {'nat movies'} % there are always 2 movies
            %             stimTimes = stimInfo(iType).stimTimes;
            %             stimSeq = stimInfo(iType).stimSequence;
            %             nStim = numel(unique(stimSeq.seq));
            %             nRep = sum(stimSeq.seq == 1);
            
        case  {'sparseNoise', 'spontaneous','size', 'nat movies'}
    end
    
    switch stimInfo(iType).stimType
        
        case 'plaids'
            stimSet = 13:24;
            blank = 25:26;
        case 'oris'
            stimSet = 1:nStim;
            blank = 13;
            %         case 'size'
            %             stimSet = 1:10;
            %             blank = 11;
        case 'oriMultiSfTf'
            stimSet = 1:nStim;
            blank = nStim+1:nStim+nBlanks;
            
        case {'SF', 'TF'}
            stimSet = 1:60;
            blank = 61:64;
            %         case 'nat movies'
            %             if nStim == 5
            %                 stimSet2 = 1:2;
            %                 stimSet1 = 3:5;
            %                 blank = [];
            %             else
            %                 stimSet2 = 1;
            %                 stimSet1 = 2:4;
            %                 blank = [];
            %             end
            % %                 stimSet = 1:2;
            % %                 blank = [];
            % %             else
            % %                 stimSet = 1;
            % %                 blank = [];
            % %             end
        case  {'sparseNoise', 'spontaneous','size', 'nat movies'}
    end
    
    % expID = find(strcmp({stimInfo.stimType}, 'plaids'));
    
    
    %% define stimulus parameters
    switch stimInfo(iType).stimType
        
        case {'oriMultiSfTf','plaids', 'oris', 'SF', 'TF'}
            
            
            nStim = numel(stimSet);
            nBlanks = numel(blank);
            
            stimON = zeros(nStim, nRep);
            stimOFF = zeros(nStim, nRep);
            
            blankON = zeros(nBlanks, nRep);
            blankOFF = zeros(nBlanks, nRep);
            
            for iStim = 1:nStim
                stimON(iStim,:) = stimTimes.onset(stimSeq.seq == stimSet(iStim));
                stimOFF(iStim,:) = stimTimes.offset(stimSeq.seq == stimSet(iStim));
            end
            
            for iBl = 1:nBlanks
                blankON(iBl,:) = stimTimes.onset(stimSeq.seq == blank(iBl));
                blankOFF(iBl,:) = stimTimes.offset(stimSeq.seq == blank(iBl));
            end
            
            iSi = min(stimTimes.onset(2:end)-stimTimes.offset(1:(end-1)));
            stimDur = round(mean(stimOFF(:)- stimON(:)));
            
            % chunk up the F in sweeps
            
            signalON =  -1;
%             signalOFF = stimDur+2;
            signalOFF = 2+2;
            signalTimes = linspace(signalON, signalOFF, range([signalON, signalOFF])*3); % interpolate at 3 Hz
            
            if size(thisF, 2)>1
                [signal{iType}, aveSignal, seSignal] = magicETA2(thisFrameTimes, thisF, stimON, signalTimes); % nTrial*nT*nStim*nN, or nT*nStim*nN
            else
                [signal{iType}, aveSignal, seSignal] = magicETA(thisFrameTimes, thisF, stimON, signalTimes); % nTrial*nT*nStim*nN, or nT*nStim*nN
                
            end
            %         stimDur = mean(blankOFF(:)- blankON(:));
            %         signalON =  -1;
            %         signalOFF = stimDur+2;
            signalTimes = linspace(signalON, signalOFF, range([signalON, signalOFF])*3); % interpolate at 3 Hz
            if size(thisF, 2)>1
                [blankSignal{iType}, ~, ~] = magicETA2(thisFrameTimes, thisF, blankON, signalTimes); % nTrial*nT*nStim*nN, or nT*nStim*nN
            else
                [blankSignal{iType}, ~, ~] = magicETA(thisFrameTimes, thisF, blankON, signalTimes); % nTrial*nT*nStim*nN, or nT*nStim*nN
                
            end
            
            blankStd{iType}= std(reshape(blankSignal{iType}, [], nN), [], 1);
            %% #####
            %     case 'nat movies'
            %
            %         stimSet = stimSet1;
            %         nStim = numel(stimSet);
            %         nBlanks = numel(blank);
            %
            %         stimON = zeros(nStim, nRep);
            %         stimOFF = zeros(nStim, nRep);
            %
            %         blankON = zeros(nBlanks, nRep);
            %         blankOFF = zeros(nBlanks, nRep);
            %
            %         for iStim = 1:nStim
            %             stimON(iStim,:) = stimTimes.onset(stimSeq.seq == stimSet(iStim));
            %             stimOFF(iStim,:) = stimTimes.offset(stimSeq.seq == stimSet(iStim));
            %         end
            %
            %         for iBl = 1:nBlanks
            %             blankON(iBl,:) = stimTimes.onset(stimSeq.seq == blank(iBl));
            %             blankOFF(iBl,:) = stimTimes.offset(stimSeq.seq == blank(iBl));
            %         end
            %
            %         iSi = min(stimTimes.onset(2:end)-stimTimes.offset(1:(end-1)));
            %         stimDur = mean(stimOFF(:)- stimON(:));
            %         signalON =  -1;
            %         signalOFF = stimDur+2;
            %         signalTimes = linspace(signalON, signalOFF, range([signalON, signalOFF])*3); % interpolate at 3 Hz
            %         stimON = reshape(stimON, 1, nStim*nRep);
            %         [signal{iType}, ~, ~] = magicETA2(thisFrameTimes, thisF, stimON, signalTimes); % nTrial*nT*nStim*nN, or nT*nStim*nN
            %
            %         stimSet = stimSet2;
            %         nStim = numel(stimSet);
            %         nBlanks = numel(blank);
            %
            %         stimON = zeros(nStim, nRep);
            %         stimOFF = zeros(nStim, nRep);
            %
            %         blankON = zeros(nBlanks, nRep);
            %         blankOFF = zeros(nBlanks, nRep);
            %
            %         for iStim = 1:nStim
            %             stimON(iStim,:) = stimTimes.onset(stimSeq.seq == stimSet(iStim));
            %             stimOFF(iStim,:) = stimTimes.offset(stimSeq.seq == stimSet(iStim));
            %         end
            %
            %         for iBl = 1:nBlanks
            %             blankON(iBl,:) = stimTimes.onset(stimSeq.seq == blank(iBl));
            %             blankOFF(iBl,:) = stimTimes.offset(stimSeq.seq == blank(iBl));
            %         end
            %
            %         iSi = min(stimTimes.onset(2:end)-stimTimes.offset(1:(end-1)));
            %         stimDur = mean(stimOFF(:)- stimON(:));
            %         signalON =  -1;
            %         signalOFF = stimDur+2;
            %         signalTimes = linspace(signalON, signalOFF, range([signalON, signalOFF])*3); % interpolate at 3 Hz
            %         stimON = reshape(stimON, 1, nStim*nRep);
            %         [signal{nStimsTypes+1}, ~, ~] = magicETA2(thisFrameTimes, thisF, stimON, signalTimes); % nTrial*nT*nStim*nN, or nT*nStim*nN
            %         stimInfo(nStimsTypes+1).stimType = 'nat movies';
            
            % #####
            
            
        case {'sparseNoise', 'spontaneous','size', 'nat movies'}
            spontCorrMat{iType} = NaN(nN,nN);
            noiseCorrMat{iType} = NaN(nN,nN);
            signalCorrMat{iType}= NaN(nN,nN);
            signal{iType} = [];
            blankSignal{iType} = [];
            blankStd{iType} = [];
            
    end
    % spontaneousON = +1;
    % spontaneousOFF = iSi;
    % spontaneousTimes = linspace(spontaneousON, spontaneousOFF, range([spontaneousON , spontaneousOFF])*3); % interpolate at 3 Hz
    
    % [spontaneous, aveSpont, seSpont] = magicETA2(thisFrameTimes, thisF, stimOFF, signalTimes);
    
end
% reorder matrices so that all the exp have the same stimulus type sequence
testedTypes = {stimInfo.stimType};
stimTypeOrder = arrayfun(@(k) find(strcmp(testedTypes{k},Types)), 1:numel(testedTypes));
[current, reorder] = sort(stimTypeOrder, 'ascend');

% for iType = 1:numel(Types)
% responses{iType} = [];
% end

% responses(standard) =  signal(reorder);
%
try
    switch intersect(setdiff(1:numel(Types), current), [1,2])
        case 1
            responses = cell(1, numel(signal)+1);
            responses(2:end) = signal(reorder);
            blankResp = cell(1, numel(blankSignal)+1);
            blankResp(2:end) = blankSignal(reorder);
            bStd = cell(1, numel(blankStd)+1);
            bStd(2:end) = blankStd(reorder);
        case 2
            responses = cell(1, numel(signal)+1);
            responses([1, 3:end])=  signal(reorder);
            blankResp = cell(1, numel(blankSignal)+1);
            blankResp([1, 3:end])=  blankSignal(reorder);
            bStd = cell(1, numel(blankStd)+1);
            bStd([1, 3:end])=  blankStd(reorder);
    end
catch
    responses =  signal(reorder);
    blankResp =  blankSignal(reorder);
    bStd =  blankStd(reorder);
    
end

reorderedType = testedTypes(reorder);

try
    if strcmp(reorderedType(end), 'sparseNoise') % if the last stim is sparseNoise drop it.
        responses=  responses(1:end-1);
        blankResp=  blankResp(1:end-1);
        bStd =  bStd(1:end-1);
        
        reorderedType = reorderedType(1:end-1);
    end
    
    if strcmp(reorderedType(end), 'spontaneous') % if the last stim is spontaneous drop it.
        responses=  responses(1:end-1);
        blankResp=  blankResp(1:end-1);
        bStd =  bStd(1:end-1);
        
        reorderedType = reorderedType(1:end-1);
    end
    
    if strcmp(reorderedType(end), 'size') % if the last stim is size drop it.
        responses=  responses(1:end-1);
        blankResp=  blankResp(1:end-1);
        bStd =  bStd(1:end-1);
        
        reorderedType = reorderedType(1:end-1);
    end
    
    if strcmp(reorderedType(end), 'nat movies') % if the last stim is size drop it.
        responses=  responses(1:end-1);
        blankResp=  blankResp(1:end-1);
        bStd =  bStd(1:end-1);
        
        reorderedType = reorderedType(1:end-1);
    end
    
end
end