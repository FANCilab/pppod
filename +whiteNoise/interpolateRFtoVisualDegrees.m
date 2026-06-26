function [rf_visDeg, x1, y1] = interpolateRFtoVisualDegrees(rf, stimEdges)

% stimEdges is [left right bottom top]. The y edge order is preserved so RF
% rows map to the exported elevation coordinates without flipping signs.
gridW = diff(stimEdges(1:2)) / size(rf,2);
gridH = diff(stimEdges(3:4)) / size(rf,1);

% gridpoints of rf in visual degrees
gridX = linspace(stimEdges(1)+0.5*gridW, stimEdges(2)-0.5*gridW, size(rf,2));
gridY = linspace(stimEdges(3)+0.5*gridH, stimEdges(4)-0.5*gridH, size(rf,1));
[x0, y0] = meshgrid(gridX, gridY);

% vectors x1 and y1 specify gridlines with distance of 1 degree
x1 = ceil(gridX(1)) : floor(gridX(end));
if gridY(1) <= gridY(end)
    y1 = ceil(gridY(1)) : floor(gridY(end));
else
    y1 = floor(gridY(1)) : -1 : ceil(gridY(end));
end

% interp2 needs query coordinates with matching sizes, so build the full
% visual-degree query grid from the x/y vectors before interpolating.
[xq, yq] = meshgrid(x1, y1);
rf_visDeg = NaN(numel(y1), numel(x1), size(rf,3));
for k = 1:size(rf,3)
    rf_visDeg(:,:,k) = interp2(x0, y0, rf(:,:,k), xq, yq);
end


% % gridpoints of rf
% [x0, y0] = meshgrid((1:size(rf,2))-0.5, (1:size(rf,1))-0.5);
% % vectors x1 and y1 specify gridlines with distance of 1
% % degree (diff(stimPos(...))); values match position of
% % gridlines (edges) in pixels of stimulus row/column
% x1 = linspace(0, size(rf,2), diff(stimEdges(1:2)));
% y1 = linspace(0, size(rf,1), -diff(stimEdges(3:4)));
% 
% % need to delete pixel values outside  given stimulus pixels,
% % so we can use interpolation (rather than extrapolation) when
% % mapping the RF from stimulus pixels to visual degrees
% x2 = x1;
% x2(x1<x0(1) | x1>x0(end)) = [];
% y2 = y1;
% y2(y1<y0(1) | y1>y0(end)) = [];
% [x2, y2] = meshgrid(x2, y2);
% 
% rf_visDeg = NaN(size(x2,1), size(x2,2), size(rf,3));
% for k = 1:size(rf,3)
%     rf_visDeg(:,:,k) = interp2(x0, y0, rf(:,:,k), x2, y2);
% end
