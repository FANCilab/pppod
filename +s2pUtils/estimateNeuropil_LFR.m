function [signalTrace, neuropCorrPars]=estimateNeuropil_LFR(cellRoiTrace,neuropRoiTrace,opt)
%
% estimateNeuropil estimates the correction factor r to be applied to the
% neuropil subtraction assuming that being C the proper signal of the cell,
% N the contamination from neuropil, and S the measured signal of the cell
% then C = S - rN
%
% INPUTS:
% cellRoiTrace: (S) time traces of fluorescence in the cell ROI
% neuropRoiTrace: (N) time traces of fluorescence surrounding the cell ROI
% opt.numN: is the number of windows in whoch compute the discretized
% values of Neuropil (default is 20)
% opt.minNp: is the lowest percentile of N from which
% compute the windows (default is 10)
% opt.maxNp: is the highest percentile of N from which
% compute the windows (default is 90)
% opt.pCell: is the lowest percentile of S to be applied to each discrete
% window of N (default is 5)
%
% OUTPUTS:
% cellRoiTrace: (C) proper signal of the cell computed by C=S-rN
% neuropCorrPars.fitNeuro: are the discretized values of neuropil used to
% compute the correction factor r
% neuropCorrPars.corrFactor: nCells x 2 matrix, second column (r) is the
% correction factor. It is obtained by a linear fit of fitNeuro and lowCell
% neuropCorrPars.lowCell: is the lowest percentile of S for each discrete
% neuropCorrPars.F0: is the F0 estimated from the model
% window of N (default is 5)
%
% 2015.06.09 Mario Dipoppa - Created
% 2024.01.24 Federico Rossi - added 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<3 || ~isfield(opt, 'numN')
    numN=20; 
else
    numN=opt.numN;
end
if nargin<3 || ~isfield(opt, 'minNp')
    minNp=10;
else
    minNp=opt.minNp;
end
if nargin<3 || ~isfield(opt, 'maxNp')
    maxNp=90;
else
    maxNp=opt.maxNp;
end
if nargin<3 || ~isfield(opt, 'pN')
    pCell=5;
else
    pCell=opt.pCell;
end


[nCells,nT]=size(neuropRoiTrace);

fitNeuro=nan(nCells,numN);
lowCell=nan(nCells,numN);
corrFactor=nan(nCells,2);

signalTrace=nan(nCells,nT);
for iCell=1:nCells
    
    if mod(iCell,5)>0
        fprintf([num2str(iCell) ' '])
    else
        fprintf([num2str(iCell) '\n'])
    end
    
    if all(isnan(neuropRoiTrace(iCell,:))) || all(isnan(cellRoiTrace(iCell,:)))
        continue
    end
    
    minN=prctile(neuropRoiTrace(iCell,:),minNp);
    maxN=prctile(neuropRoiTrace(iCell,:),maxNp);
    discrNeuro=round(numN*(neuropRoiTrace(iCell,:)-minN)/(maxN-minN));
    %discrNeuro are the discretized values of neuropil between minN and
    % maxN, with numN elements
    
    for iN=1:numN
        lowCell(iCell,iN)= prctile(cellRoiTrace(iCell,discrNeuro==iN),pCell);
    end
    
    fitNeuro(iCell,:)=(1:numN).*(maxN-minN)/numN+minN;
    corrFactor(iCell,:) = robustfit(fitNeuro(iCell,:),lowCell(iCell,:));
    %fit between discretized Neuropil and lowest percentile of signal in
    %the cell ROI
    
    signalTrace(iCell,:)=cellRoiTrace(iCell,:)-corrFactor(iCell,2)*neuropRoiTrace(iCell,:)-corrFactor(iCell,1);
%     signalTrace(iCell,:) = signalTrace(iCell,:) - prctile(signalTrace(iCell,:),5);
end
fprintf('\n')

neuropCorrPars.fitNeuro=fitNeuro;
neuropCorrPars.corrFactor=corrFactor;
neuropCorrPars.lowCell=lowCell;
neuropCorrPars.F0 = corrFactor(:, 2).*fitNeuro(:, 1) + corrFactor(:, 1);