function lpw = lookUpPower(waveL)

measuredWL = [760 780 800 820 840 860 880 890 910 930 950 970 980 1000 1020 1040];
measuredPw = [3670 4100 4140 3920 3590 3280 3020 2740 2340 2010 1600 1380 1280 970 830 740];

WL = 760:1:1040;
Pw = interp1(measuredWL, measuredPw, WL);

% measuredWLObj = [760 780 800 820 840 860 880 890 900 920 940 960 970 990 1000 1020 1040];
% measuredPwObj = [127.9 146.4 149.9 139.6 128 116.8 112.6 98.7 86.4 70.2 62.1 48.8 45.3 37.6 33.8 30.2 25.6];

% PwObj = interp1(measuredWLObj, measuredPwObj, WL);

[~,  query] = ismember(waveL, WL);
lpw = Pw(query);

end