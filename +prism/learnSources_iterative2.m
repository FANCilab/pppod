function  [sources, mixing, finalErr] = learnSources_iterative2(data, mixGuess, forceZero, doPlot)


if nargin < 4
    forceZero = false;
end

if nargin < 3
    doPlot = false;
end

% add intercept? 
% mixGuess = cat(2, mixGuess, ones(size(mixGuess, 1),1));
% mixGuess = cat(2, mixGuess, mean(mixGuess, 2));

r0 = mixGuess(:,2)./mixGuess(:,1);
r = r0;


%  estimate sources with initial mixing guess
sourceGuess = pinv(mixGuess)*data; % data is nCh*nPx, mixing is nCh*nF, sources is nF*nPx

if forceZero % if requested, force sources to be non-negative
    sourceGuess = max(sourceGuess,0);
end

mixing = pinv(sourceGuess')*data'; % data is nPx *nCh, sources is nPx*nF, mixing is nF*nCh
mixing = mixing';
mixing(mixing<0) = 0.000001; % negative weights do not make sense

% allow mixing coefficient to scale, but constrain relative ratio to be within range from initial guess
% r = mixing(:,2)./mixing(:,1);
% for iW = 1:size(mixing,1)
%     if r0(iW)/r(iW) <0.1|| r0(iW)/r(iW) >10 % if relative ration changes too much, reset
%         %             mixing(iW, 2:end) =  mixing(iW, 1)*mixRatio(iW, 2:end);
%         mixing(iW, 2) =  mixing(iW, 1)*r0(iW);
%     end
%     %        mixing(:, 2) = bsxfun(@times, mixing(:,1),mixRatio(:,2));
% end

% recompute sources
sourceGuess = pinv(mixing)*data;
if forceZero
    sourceGuess = max(sourceGuess,0);
end

% initialise loss function
errReconstruct = data - mixing*sourceGuess;
errReconstruct = sum(errReconstruct(:).^2)/sum(data(:).^2);
negWeigth = sum(sourceGuess(sourceGuess(2,:)<0).^2)/sum(sourceGuess(2, :).^2);
mutInfo = mutualInformation(sourceGuess(1,:), sourceGuess(2,:));
% backGround = sum(sourceGuess(3,:).^2)/sum(sourceGuess(:).^2);
sourceBlow = sum(sourceGuess(2,:).^2)/sum(sourceGuess(:).^2);

iter = 1;
if forceZero
    err(iter)  = errReconstruct*mutInfo ;
else
    err(iter)  = errReconstruct*mutInfo*negWeigth;
end

tol = 1;

if doPlot
    e = figure;
end

mixGuess = mixing;
while  tol > 10^(-10) && iter <100
    
    iter = iter + 1;
    
    % improve guess of mixing matrix
    mixing = pinv(sourceGuess')*data'; % data is nPx *nCh, sources is nPx*nF, mixing is nF*nCh
    mixing = mixing';
    mixing(mixing<0) = 0.000001; % negative weights do not make sense
    
    % allow mixing coefficient to scale, but constrain relative ratio to be within range from initial guess
%     r = mixing(:,2)./mixing(:,1);
%     for iW = 1:size(mixing,1)
%         
%         if r0(iW)/r(iW) <0.1  % r is much bigger, we can tolerate only s2 getting to 0
%             %             mixing(iW, 2:end) =  mixing(iW, 1)*mixRatio(iW, 2:end);            
%             mixing(iW, 1) =  mixing(iW, 2)/r0(iW);
% 
%         elseif r0(iW)/r(iW) >10        
%             mixing(iW, 2) =  mixing(iW, 1)*r0(iW);
%         end
%         %        mixing(:, 2) = bsxfun(@times, mixing(:,1),mixRatio(:,2));
%     end
    
    % recompute sources
    sources = pinv(mixing)*data;
    if forceZero
        sources = max(sources,0);
    end
    
    % calculate loss function
    
    errReconstruct = data - mixing*sources;
    errReconstruct = sum(errReconstruct(:).^2)/sum(data(:).^2);
    negWeigth = sum(sources(sources(2,:)<0).^2)/sum(sources(2, :).^2);
    mutInfo = mutualInformation(sources(1,:), sources(2,:));
%     backGround = sum(sources(3,:).^2)/sum(sources(:).^2);
    sourceBlow = sum(sources(2,:).^2)/sum(sources(:).^2);

    if forceZero
        err(iter)  = errReconstruct*mutInfo ;
    else
        err(iter)  = errReconstruct*mutInfo*negWeigth*sourceBlow;
    end
    
    % compute improvement from previous iteration
    tol =  err(iter-1) - err(iter);
    
    % plot loss function if requested
    if doPlot
        cla;plot(1:iter, log10(err), '-o'); hold on; drawnow;
        xlabel('iteration #')
        ylabel('log10(variance explained)')
    end
    % pause;
    
    % break iterations when minimum is found or tolerance criterion is met
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