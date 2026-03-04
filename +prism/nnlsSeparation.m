function sources = nnlsSeparation(mixing, data, waveL, wRep, nRep)

addpath(genpath('C:\Users\Federico\Documents\MATLAB\nnlslab'));

wRep = find(waveL == wRep);

replicaMix = mean(mixing(wRep, :),1);
replicaImg = mean(data(:, wRep), 2);

if nRep >0
mixing = cat(1, mixing, repmat(replicaMix, nRep, 1));
data = cat(2, data, repmat(replicaImg, 1, nRep));
end

mixing = bsxfun(@rdivide, mixing, mean(mixing(:)));

opt_gen1 = initopt_general('perfmeas', {1,2}, 'maxcpu', 5, 'stopcrit', 2, 'tol', 1E-12);
pg_1 = projgradient(mixing, data', opt_gen1, opt_projgradient());
sources = pg_1.xopt;

% sources = pinv(mixing)*data';

% mixing = bsxfun(@rdivide, mixing, mean(mixing(:)));
% 
% allPx = bsxfun(@rdivide,allPx, max(allPx, [], 1));
 
% % sources = fastica(allPx'); 
end