function  [sources, mixing, finalErr] = unmixSources( data, mixGuess,waveL, constrainMix, doPlot)

if nargin<3
    waveL = [];
end


if nargin < 4
    constrainMix = 'trueMix';
end

if nargin < 5
    doPlot = false;
end

addpath(genpath('C:\Users\Federico\Documents\MATLAB\nnlslab'));

%  estimate sources with initial mixing guess
sourceGuess = pinv(mixGuess)*data; % data is nCh*nPx, mixing is nCh*nF, sources is nF*nPx
sourceGuess = max(sourceGuess,0);

% sourceGuess = prism.nnlsSeparation_dev(mixGuess, data');

iter = 1;
err(1) = sum(sum((data - mixGuess*sourceGuess).^2, 1),2)./sum(data(:).^2); 
tol = err(1);

if doPlot
e = figure;
end

while  tol > 10^(-10) &&iter <100
    
iter = iter + 1;

% improve guess of mixing matrix
mixing = pinv(sourceGuess')*data'; % data is nPx *nCh, sources is nPx*nF, mixing is nF*nCh
mixing = mixing';
mixing(mixing<0) = 0.000001; % negative weights does not make sense 

% apply constraints to mixing matrix
switch constrainMix
    case 'trueMix'
        mixRatio = bsxfun(@rdivide, mixGuess, mixGuess(:,1));
        mixing(:, 2:end) =  bsxfun(@times, mixing(:, 1), mixRatio(:, 2:end));
%         mixing(waveL <= 970 & waveL >=820 , 2) = 0;
    case 'forceZero'
            mixing(waveL <= 920 & waveL >=820 , 2) = 0;
%             mixing(waveL == 970, 2) = 0;
    case 'free'

end

% recompute sources
sources = pinv(mixing)*data;
sources = max(sources,0);

% calculate reconstruction error

err(iter) = sum(sum((data - mixing*sources).^2, 1),2)./sum(data(:).^2); 
tol =  err(iter-1) - err(iter); 


if doPlot
cla;plot(1:iter, log10(err), '-o'); hold on; drawnow; 
xlabel('iteration #')
ylabel('log10(variance explained)')
end
% pause;
if tol >0
mixGuess = mixing;
sourceGuess = sources;
finalErr = err(iter);
else
    mixing = mixGuess;
    sources = sourceGuess;
    finalErr = err(iter-1);
end

end


end