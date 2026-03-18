function [frameG, frameR, micronsX, micronsY, micronsZ, starterXYZ] = loadMultiPlaneMeanStack(db, doPlot)

if nargin <2
    doPlot = 0;
end


nExp = numel(db);

planeCount = 0;

if doPlot
figure;
end

for iExp = 1:nExp
    
    try
    info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(db(iExp).expID(1)));
    catch
            info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(1));

    end
    [micronsX{iExp}, micronsY{iExp}, micronsZ{iExp}] = ppbox.getPxXYZ(info);
    
%     if iExp == 1
%         if db(1).starterYXPlane(3) > 0
%             starterXYZ = [micronsX{iExp}(db(1).starterYXPlane(2)), ...
%                 micronsY{iExp}(db(1).starterYXPlane(1)),...
%                 micronsZ{iExp}(db(1).starterYXPlane(1), db(1).starterYXPlane(3)) + db(1).ObjZ];
%         elseif isfield(db, 'starterZ')
%             starterXYZ = [micronsX{iExp}(db(1).starterYXPlane(2)), ...
%                 micronsY{iExp}(db(1).starterYXPlane(1)),...
%                 db(1).starterZ];
%             
%         end
%     end
%     
    micronsZ{iExp} = micronsZ{iExp} + db(iExp).ObjZ;

    ny = numel(micronsY{iExp});
    nx = numel(micronsX{iExp});

    for iPlane = 1: info.nPlanes
        
        [root, refF, ~, ~, refRaw] = starter.getAnalysisRefs(db(iExp).mouse_name, db(iExp).date, db(iExp).expts, iPlane);
        
        try
            load(fullfile(root, refF));
            planeCount = planeCount +1;    
            
            dummyG = zeros(ny, nx);
            dummyG(dat.ops.yrange,dat.ops.xrange) = dat.mimg(:,:,2);
            dummyR = zeros(ny, nx);
            dummyR(dat.ops.yrange,dat.ops.xrange) = dat.mimg(:,:,3);
            frameG(:,:,planeCount) = dummyG;
            frameR(:,:,planeCount) = dummyR;
            %             [ fovx, fovy ] = zoom2fov(info.zoomFactor);
            %             PxSz = (fovx+fovy)/(2*512);
            if doPlot
%                 color = [1 1 0.7 0.4];
                r = subplot(nExp, info.nPlanes, (iExp-1)*info.nPlanes + iPlane);
                s2pUtils.plotS2pRois(dat);
            end
        catch
            load(fullfile(root, refRaw));
            planeCount = planeCount +1;    
            dummyG = zeros(ny, nx);
            dummyG(ops.yrange,ops.xrange) = ops.mimg1(ops.yrange, ops.xrange);
            dummyR = zeros(ny, nx);
            frameG(:,:,planeCount) = dummyG;
            if isfield(ops, 'mimgRED')
            dummyR(ops.yrange,ops.xrange) = ops.mimgRED(ops.yrange, ops.xrange);
            frameR(:,:,planeCount) = dummyR;
            else
            frameR(:,:,planeCount) = NaN(size(frameG(:,:,1)));
    
            end
%             [ fovx, fovy ] = zoom2fov(info.zoomFactor);
%             PxSz = (fovx+fovy)/(2*512);
        end
    end
    
  
end


end