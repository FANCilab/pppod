function [rgb_image, mimgG, mimgR] = plotRois(db, iPlane, color_opt)

if nargin <3
    color_opt = 1;
end

root = 'D:\OneDrive - University College London\Data\2P\';

if iscell(db.mouse_name)
    info = ppbox.infoPopulateTempLFR(db.mouse_name{1}, db.date{1}, db.expts{1}(1));
else
    info = ppbox.infoPopulateTempLFR(db.mouse_name, db.date, db.expts(1));
    
end

[folder, refF, ~] = starter.getAnalysisRefs(db.mouse_name, db.date, db.expts, iPlane);

load(fullfile(folder, refF)); % load ROIs data

data = dat;

mimgG = data.mimg(:,:,2);

mimgR = data.mimg(:,:,3);

mimgR = reshape(prism.bleedCure(mimgR, mimgG), size(mimgR));

Ly = data.cl.Ly;

Lx = data.cl.Lx;

map = zeros(Ly, Lx);

goodROI = find([data.stat.iscell]);
redROI = find([data.stat.redcell]);

nROI = numel(goodROI);

%     expVar = data{iPlane}.res.M;
expVar = data.res.lambda(:);
varR = zeros(size(expVar));
varRRed = zeros(size(expVar));
h = randperm(nROI);
for ir = 1:nROI
    ipix = data.stat(goodROI(ir)).ipix;
    map(ipix) = h(ir);
    vM = expVar(ipix);
    vM = vM/sum(vM.^2)^.5;
    varR(ipix) = vM;
    
end

for ir = 1:numel(redROI)
    ipix = data.stat(redROI(ir)).ipix;
    mapR(ipix) = h(ir);
    vM = expVar(ipix);
    vM = vM/sum(vM.^2)^.5;
    varRRed(ipix) = vM;
    
end

if color_opt
%    V = max(0, min(.5 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));
V = max(0, min(.2 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));
H = (map-1)/max(map(:));
Sat = ones(size(V));
rgb_image= hsv2rgb(cat(3, H, Sat, V));
rgb_image(:,:,1) = 0;
else
rr = max(0, min(.5 * reshape(varRRed, Ly, Lx)/mean(varR(:)), 2));

R = 1-max(0, min(.15 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));
G = 1-max(0, min(.15 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));
B = 1-max(0, min(.15 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));

R(rr>0) = rr(rr>0);
% R(rr>0) = c;
G(rr>0) = 0;
B(rr>0) = 0;
% R = zeros(size(R));
rgb_image = cat(3, R,G,B);
end
% figure;
% image(rgb_image); axis image;
% set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
% formatAxes

figure;
g = subplot(1,3,1);
imagesc(mat2gray(mimgG)); axis image;
set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
caxis([0.05 0.3]);colormap(g, 'gray');

r = subplot(1,3,2);

imagesc(mat2gray(mimgR)); axis image;
set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
caxis([0.05 0.3]); colormap(r, Red(100));

ror = subplot(1,3,3);
ro = imagesc(rgb_image); axis image;
set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
set(ro, 'alphadata', map ~= 0)
formatAxes

end
