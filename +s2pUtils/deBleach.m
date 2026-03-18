function [newData, indUnits] = deBleach(traces, time, opts)

% for correcting baseline drifts of calcium traces at start of experiments
if nargin <3
    driftWin = 20; % in s, window to test whether baseline is higher than normal
    driftThresh = 1.5; % in std, threshold for drift
    correctWin = 150; % in s, window to fit exponential

else

driftWin =  opts.driftWin; % in s, window to test whether baseline is higher than normal
driftThresh = opts.driftThresh; % in std, threshold for drift
correctWin = opts.correctWin; % in s, window to fit exponential

end

t_ca = time;
timeBin = median(diff(t_ca));

% remove strong baseline decay at start of experiment in cells that
% show it
indUnits = find(mean(traces(1:round(driftWin / timeBin),:), 1, 'omitnan') > ...
    mean(traces, 1, 'omitnan') + driftThresh .* std(traces,0,1, 'omitnan'));
ind = round(correctWin / timeBin);
for iUnit = 1:length(indUnits)
    y = traces(1:ind, indUnits(iUnit));
    y = fillmissing(y, 'linear');
    % fit double exponential to start of trace
    f = fit((1:length(y))', y, ...
        @(a,b,c,d,e,x) a + b .* exp(-x ./ c) + d .* exp(-x ./ e), ...
        'Lower', [0 0 0 0 0], ...
        'Upper', [max(y) max(y) 500 max(y) 500], ...
        'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
    % remove fit
    traces(:, indUnits(iUnit)) = traces(:, indUnits(iUnit)) - ...
        f(1 : size(traces,1)) + f.a;
end

newData = traces; %(traces-mean(traces,1,  'omitnan'))./std(traces, 0, 1, 'omitnan');
end