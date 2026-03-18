function [imgG, imgR, imgM, bkgG, bkgR, cellG, cellR] = loadROIcrops(db, cropSize)

if nargin <2
    cropSize = 25;
end

nExp = numel(db);



for iExp = 1:nExp
    
    info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(db(iExp).expID(1)));
    
    if isfield(db(iExp).starterID)
        
        [file, folder] = starter.planeRef(db);
        root = fullfile(root,folder);
        load(fullfile(root, file));
        
        
        
    else
    
    [micronsX, micronsY, micronsZ] = ppbox.getPxXYZ(info);
    
    
    for iPlane = 1: info.nPlanes
        
        [root, refF, refNeu] = starter.getAnalysisRefs(db(iExp).mouse_name, db(iExp).date, db(iExp).expts, iPlane);
        
        if exist(fullfile(root, refF))
            load(fullfile(root, refF));
            
            
            [~, nanN] = find(isnan(dat.Fcell{1}'));
            nanN = unique(nanN);

            redCells = setdiff(find([dat.stat.redcell] & [dat.stat.iscell]), nanN);
            
            if sum(redCells>0)
            map = {dat.stat(redCells).ipix};
            
            
            frameG = dat.mimg(:,:,2);
            frameR = dat.mimg(:,:,3);
            try
            frameR = prism.bleedCure(frameR, frameG);
            end
            [ fovx, fovy ] = ppbox.zoom2fov(info.zoomFactor);
            PxSz = (fovx+fovy)/(2*512);
            
            [imgG{iExp, iPlane}, imgR{iExp, iPlane}, imgM{iExp, iPlane},...
                bkgG{iExp, iPlane}, bkgR{iExp, iPlane}, cellG{iExp, iPlane}, cellR{iExp, iPlane}] = ...
                s2pUtils.cutRoiImg(map, frameG, frameR, PxSz, cropSize);
%             coords{iExp, iPlane} = starter.getROIxyz(map, micronsX(dat.ops.xrange), micronsY(dat.ops.yrange) , micronsZ(dat.ops.yrange, iPlane));
%             coords{iExp, iPlane}(:,3) = coords{iExp, iPlane}(:,3) + 125*(iExp-1);
            %     figure; scatter3(coords(:,2), coords(:,1),-coords(:,3));axis xy
            
            end
        
            
        end
    end
    end
    
  
end


end