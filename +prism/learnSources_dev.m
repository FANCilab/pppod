function  [sources, mixing] = learnSources_dev(data, mixGuess, waveL)

addpath(genpath('C:\Users\Federico\Documents\MATLAB\nnlslab'));

%  estimate sources with initial mixing guess

% sourceGuess = pinv(mixGuess)*data; % data is nCh*nPx, mixing is nCh*nF, sources is nF*nPx
sourceGuess = prism.nnlsSeparation_dev(mixGuess, data');
% sourceGuess = max(sourceGuess,0);
% sourceGuess = sourceGuess.*(sourceGuess >0.1);

e = figure;

iter = 1;

err(1) = sum(sum((data - mixGuess*sourceGuess).^2, 1),2)./sum(data(:).^2); 
tol = err(1);

while tol > 10^(-10) && iter <100
    
iter = iter + 1;

% improve guess of mixing matrix

% mixing = pinv(sourceGuess')*data'; % data is nPx *nCh, sources is nPx*nF, mixing is nF*nCh
mixing = prism.nnlsSeparation_dev(sourceGuess', data);


mixing = mixing';
mixing(mixing<0) = 0;
mixing(:,2) =  mixing(:,1).*mixGuess(:,2)./mixGuess(:,1);
mixing(waveL == 890, 2) = 0;
% mixing(waveL == 970, 2) = 0;

% recompute sources

% sources = pinv(mixing)*data;
sources = prism.nnlsSeparation_dev(mixing, data');

% sources = max(sources,0);
% sources = sources.*(sources >0.1);

% calculate reconstruction error

err(iter) = sum(sum((data - mixing*sources).^2, 1),2)./sum(data(:).^2); 
tol =  err(iter-1) - err(iter); 

plot(1:iter, log10(err), '-o'); hold on; drawnow; 
% pause;
if tol >0
mixGuess = mixing;
sourceGuess = sources;
else
    mixing = mixGuess;
    sources = sourceGuess;
end

end


end