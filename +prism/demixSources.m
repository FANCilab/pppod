function errFun = demixSources(mixGuess, data)

% mixGuess = cat(2, mixGuess, ones(size(mixGuess,1),1));

sourceGuess = pinv(mixGuess)*data; % data is nCh*nPx, mixing is nCh*nF, sources is nF*nPx
% sourceGuess = max(sourceGuess,0);

err = data - mixGuess*sourceGuess;
err = err(:);

posErr = err(err>0).^2;
negErr = err(err<0).^2;

errReconstruct = (sum(posErr) + sum(negErr))/sum(data(:).^2);

% uncorrPx = setdiff(1:size(sourceGuess,2), corrPx);
% decorrErr = abs(corr(sourceGuess(1,uncorrPx)',sourceGuess(2,uncorrPx)'));
% decorrErr = abs(corr(sourceGuess(1,:)',sourceGuess(2,:)'));

nzs = sourceGuess(1,:)>0 | sourceGuess(2,:);

mutInfo = mutualInformation(sourceGuess(1,nzs), sourceGuess(2,nzs));

negWeigth = sum(sourceGuess(sourceGuess(2,:)<0).^2)/sum(sourceGuess(2,:).^2);

errFun = errReconstruct*mutInfo*negWeigth;
% errFun = abs(mean(sourceGuess(1,:).*sourceGuess(2,:)));

% iter = 1;
% err(1) = sum(sum((data - mixGuess*sourceGuess).^2, 1),2)./sum(data(:).^2); 
% tol = err(1);
% 
% while  tol > 10^(-10) &&iter <100
%     
% iter = iter + 1;
% 
% % improve guess of mixing matrix
% mixing = pinv(sourceGuess')*data'; % data is nPx *nCh, sources is nPx*nF, mixing is nF*nCh
% mixing = mixing';
% mixing(mixing<0) = 0.000001; % negative weights does not make sense 
% 
% % apply constraints to mixing matrix
% 
% mixRatio = bsxfun(@rdivide, mixGuess, mixGuess(:,1));
% mixing(:, 2:end) =  bsxfun(@times, mixing(:, 1), mixRatio(:, 2:end));
% 
% % recompute sources
% sources = pinv(mixing)*data;
% sources = max(sources,0);
% 
% % calculate reconstruction error
% 
% err(iter) = sum(sum((data - mixing*sources).^2, 1),2)./sum(data(:).^2); 
% tol =  err(iter-1) - err(iter); 
% 
% if tol >0
% mixGuess = mixing;
% sourceGuess = sources;
% errFun = err(iter);
% else
%     mixing = mixGuess;
%     sources = sourceGuess;
%     errFun = err(iter-1);
% end
% 
% end
% 

end