function coords = getROIxyz(map, xvec, yvec ,zvec)
%map is a cell array containing px indeces
%xvec, yvec, zvec are marginals coordinates of full FOV

nRoi = numel(map);

x = NaN(nRoi,1);
y = NaN(nRoi,1);
z = NaN(nRoi,1);

for iRoi = 1:nRoi
    
    [r, c] = ind2sub([numel(yvec),numel(xvec)], map{iRoi});
    
    x(iRoi) = mean(xvec(c));
    y(iRoi) = mean(yvec(r));
    z(iRoi) = mean(zvec(r));
    
    
end

 coords = [x(:), y(:), z(:)];
 
end