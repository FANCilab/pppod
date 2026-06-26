function [parameters, gaussMap] = fit2dGaussRF(rfield, doPlot, ...
    xcoords, ycoords)
% resp = your measure of the responsiveness of the neuron at each point in
% space
% x = [amplitude, xCenter, xStdev, yCenter, yStdev, rotation, baseline]

if nargin < 3
    % create mesh with x- and y-coordinates of profile
    xcoords = 1:size(rfield,2);
    ycoords = 1:size(rfield,1);
end
[x, y] = meshgrid(xcoords, ycoords);
xdata = cat(3, x, y);

% x0 is the initial guess
[indY, indX] = find(rfield == max(rfield(:)),1);
maxX = xcoords(indX);
maxY = ycoords(indY);
x0 = [1, maxX, mean(diff(xcoords))*5, maxY, mean(diff(ycoords))*5, 0];

% lower and upper bounds of parameters
lb = [0, min(xcoords), 0, min(ycoords), 0, -pi/4];
mx_std = max(range(xcoords(:)), range(ycoords(:))) / 4;
ub = [2 * max(rfield(:)), max(xcoords), mx_std, ...
    max(ycoords), mx_std, pi/4];
options = optimoptions('lsqcurvefit', 'Display', 'off');
parameters = lsqcurvefit(@whiteNoise.D2GaussFunctionRot, x0, xdata, ...
    rfield, lb, ub, options);

gaussMap = whiteNoise.D2GaussFunctionRot(parameters, xdata);

if doPlot
    figure;
    subplot(3,1,1);
    imagesc(xcoords, ycoords, rfield);
    colormap hot;

    subplot(3,1,2);
    imagesc(xcoords, ycoords, gaussMap);

    subplot(3,1,3);
    contour(gaussMap,1,'Color', 'k', 'LineWidth', 2.0);
    set(gca, 'YDir', 'reverse');
end