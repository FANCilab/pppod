function [xPref, yPref, rfContours] = getRFContour(RR, c, doPlot)

if nargin < 3
    doPlot = 0;
end

R = RR;

[ny, nx]  = size(R);

if ~isempty(c)
    
    upy = min(round(c(1)) + 5, ny); % hardcoded constraint, should vary with square size
    dwy = max(round(c(1)) - 5, 1);
    upx = min(round(c(2)) + 9, nx);
    dwx = max(round(c(2)) - 9, 1);
    
    R = RR(dwy:upy, dwx:upx);
    
    c = c - [dwy, dwx] +1;
else
    
    upy = ny;
    dwy = 1;
    upx = nx;
    dwx = 1;
    
end

maxR = max(makeVec(R));

[yMax, xMax] = ind2sub(size(R), find(R == maxR));

rfContours = struct('xx', [], 'yy', []);

thR = prctile(makeVec(R), 95);
% thR = prctile(makeVec(R), 98);% hack for Mika


C = contourc(R, [thR, thR]);

CC = getContours(C);

if length(CC)>1
    
    dist = Inf;
    for iC = 1:length(CC)
        dd = [(mean(CC(iC).xx)-xMax)^2, (mean(CC(iC).yy)-yMax)^2];
        if norm(dd) < dist
            dist = norm(dd);
            iClosest = iC;
        end
    end
    CC = CC(iClosest);
end
% xPref = mean(CC.xx)+ dwx-1;
% yPref = mean(CC.yy)+ dwy-1;
rfContours.xx = CC.xx;
rfContours.yy = CC.yy;
rfContours.xx= rfContours.xx + dwx -1;
rfContours.yy= rfContours.yy + dwy -1;

mask = poly2mask(rfContours.xx, rfContours.yy, ny,nx);

[y,x,val] = find(mask.*RR);

xPref = (x'*val)/(sum(val));
yPref = (y'*val)/(sum(val));

if doPlot
    figure;
    imagesc(RR); colormap hot; hold on
    plot(rfContours.xx, rfContours.yy, '-b', 'LineWidth', 1);
    plot(xPref, yPref, '*k')
end

end



function CC = getContours(C)

CC = struct([]);
[~, N] = size(C);
iContour = 0;
startPoint = 1;
while startPoint<N
    iContour = iContour + 1;
    CC(iContour).value = C(1, startPoint);
    nPoints = C(2, startPoint);
    CC(iContour).xx = C(1, startPoint+(1:nPoints))';
    CC(iContour).yy = C(2, startPoint+(1:nPoints))';
    startPoint = startPoint + nPoints + 1;
end
end

