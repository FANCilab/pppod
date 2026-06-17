function [ETAmat, ETA, ETAse, window] = magicETA(t, signal, etaT, periT, bsl)
%%INPUTS

% etaT is a nStim * nReps matrix of onset times
% signal is an nT*1 vector of activity
% t is a nT*1 vector of time stamps
% periT is a 2 element vector that specifies the window around etaT to be considered

%%OUTPUTS
%ETAmat
%ETA
%ETAse
%window

%% one line magic peri-stimulus triggered average, courtesy of CB wisdom
if numel(periT) ==2
    window = linspace(periT(1), periT(2), range(periT)*10);
else
    window = periT;
end

if nargin < 5
    bsl = window(window <0);
else
    if numel(bsl) ==2
        bsl = window(window > bsl(1) & window<bsl(2));
    end
end

% ETAmat = interp1(t,signal,bsxfun(@plus,window,etaT'));
% [nStim, nRep] = size(etaT);
newT = permute(etaT, [2, 3, 1]);

ETAmat = interp1(t,signal,bsxfun(@plus,window,newT));
ETAbsl = prctile(interp1(t,signal,bsxfun(@plus,bsl,newT)), 50,2);
ETAmat = bsxfun(@minus,ETAmat, ETAbsl);

%% compute median resp and se
ETA = squeeze(nanmean(ETAmat,1));
ETAse = squeeze(nanstd(ETAmat,1,1)/sqrt(size(etaT,2)));

% %% plotting
% figure
% shadePlot(window, ETA(:,5), ETAse(:,5),'b');
end