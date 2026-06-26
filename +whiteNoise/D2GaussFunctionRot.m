function F = D2GaussFunctionRot(x,xdata)
%% x = [Amp, x0, wx, y0, wy, fi]

% x and y vetors (spanning the space) are rotated into opposite direction,
% which results in the final ellipse being rotated in the given direction
rotation = x(6);

xdatarot(:,:,1) = xdata(:,:,1)*cos(rotation) - xdata(:,:,2)*sin(rotation);
xdatarot(:,:,2) = xdata(:,:,1)*sin(rotation) + xdata(:,:,2)*cos(rotation);

x0rot = x(2)*cos(rotation) - x(4)*sin(rotation);
y0rot = x(2)*sin(rotation) + x(4)*cos(rotation);

F = x(1) * exp( -((xdatarot(:,:,1) - x0rot).^2 / (2 * x(3)^2) + ...
    (xdatarot(:,:,2) - y0rot).^2 / (2 * x(5)^2) ) );
