clear;

%% user input to select the files --> needs changing to be automated
root = 'Z:\Data\2P\HS034\20260210\1\zStack';
[Rfile, root] = uigetfile(fullfile(root, '*.tif'));
Gfile = uigetfile(fullfile(root, '*.tif'));

%% Load the data
R = tiff.load(fullfile(root, Rfile));
G = tiff.load(fullfile(root, Gfile));

%% remove the flyback frames;
R = R(:,:, 3:end);
G = G(:,:, 3:end);

%% remove the bleedthrough of GCaMP in R channel
%need to improve, relationship is not purely linear, but nl fits have
%problems with saturating pixels. Piecewise linear maybe better

[curedR, bleedThrough, g2r] = ...
    prism.bleedCureNL3D(R, G, 1, [], 'linear', 1);

%% enhance soma-scale features and remove background
curedR = mat2gray(R);
neuSize = 10; 
sig = 1;
thR = prism.removeBackground(curedR, sig, neuSize);

%% detect somas with cellpose, segmenting 2D max projection

max_img = max(thR, [], 3);
cp = cellpose(Model="cyto2");
labels2D = segmentCells2D(cp,max_img,ImageCellDiameter=10, CellThreshold=0, FlowErrorThreshold= 0.4);
B = labeloverlay(max_img,labels2D);
figure;
imshow(B)

% %% the 3D version sucks a bit
% cp = cellpose(Model='cyto2');
% labels3D = segmentCells3D(cp, thStack, ImageCellDiameter = 12, CellThreshold = -5);
% volshow(thStack,RenderingStyle="SlicePlanes", ...
%     OverlayData=labels, ...
%     OverlayAlphamap=0.5);

%% TODO find the z position of each mask