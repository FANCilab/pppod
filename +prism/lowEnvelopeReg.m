function [lowEnvFit, fitY, lowY, madY] = lowEnvelopeReg(X,Y, numN, pLow, maxNp, minNp, type)

if nargin <3
    numN = 1000;
end
if nargin < 4 
    pLow = 50;
end

X = X(:); Y = Y(:);

if nargin <7
    type = 'tile';
end

switch type
    case 'tile'
    
if nargin<5
    maxNp =100; %95
end

if nargin <6
minN = max(prctile(X, 5), 0);
% minN =prctile(X,minNp);
else
minN = prctile(X, minNp);
end
maxN=prctile(X,maxNp);

    case 'numeric'
       minN = minNp;
       maxN = maxNp;
end

discrX=round(numN*(X-minN)/(maxN-minN));
%discrX are the discretized values of X between minN and
% maxN, with numN elements

for iN=1:numN
    lowY(iN)= prctile(Y(discrX==iN),pLow);
    madY(iN) = mad(Y(discrX==iN));
end

fitY=(1:numN).*(maxN-minN)/numN+minN;

fitY(isnan(lowY)) = [];
lowY(isnan(lowY)) = [];

lowEnvFit = robustfit(fitY,lowY);

% nQuant = numel(quantiles);
% 
% quantileEdges = prctile(X, quantiles);
% 
% [~,~,idx] = histcounts(X, quantileEdges); 
% 
% 
% for iQ = 1:nQuant-1
%     
%     chosenY = Y(idx == iQ);
%     chosenX = Y(idx == iQ);
%     
%     envelope(iQ) = mean(chosenY(chosenY < prctile(chosenY, 25)));
%     
%     envelopeEdges(iQ)= mean(chosenX);
% 
% end
% 
% 
% rf = robustfit(envelopeEdges, envelope);
% % 

end