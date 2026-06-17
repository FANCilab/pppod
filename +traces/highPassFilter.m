function traces = highPassFilter(traces, time, smoothWin)

dt = median(diff(time));
winSamples = round(smoothWin / dt);
smoothed = smoothdata(traces,1,"movmedian",winSamples,"omitnan");
traces = traces - smoothed;