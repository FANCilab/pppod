function traces = highPassFilter(traces, time, smoothWin)

if nargin < 3 || isempty(smoothWin) || smoothWin <= 0
    return
end

dt = median(diff(time));
winSamples = round(smoothWin / dt);

if winSamples < 1
    return
end

smoothed = smoothdata(traces,1,"movmedian",winSamples,"omitnan");
traces = traces - smoothed;