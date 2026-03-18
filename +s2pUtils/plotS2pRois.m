function plotS2pRois(data, c)

if nargin <2
    c = 1;
end

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

% %    V = max(0, min(.5 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));
% V = max(0, min(.2 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));
% H = (map-1)/max(map(:));
% Sat = ones(size(V));
% rgb_image= hsv2rgb(cat(3, H, Sat, V));

rr = max(0, min(.5 * reshape(varRRed, Ly, Lx)/mean(varR(:)), 2));

R = 1-max(0, min(.15 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));
G = 1-max(0, min(.15 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));
B = 1-max(0, min(.15 * reshape(varR, Ly, Lx)/mean(varR(:)), 2));

R(rr>0) = rr(rr>0);
% R(rr>0) = c;
G(rr>0) = 0;
B(rr>0) = 0;

rgb = cat(3, R,G,B);

% figure;
% image(rgb_image); axis image;
% set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
% formatAxes

imagesc(rgb); axis image;
set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [])
formatAxes

end
