function [choices] = measureROIRed(db, cropSize, redOnly)

if nargin <2
    cropSize = 25; % microns
end

if nargin <3
    redOnly = 0;
end
nExp = numel(db);

if isfield(db, 'root')
root = db.root;
else
% root = '\\zserver.cortexlab.net\Data\Subjects\';
root = 'D:\OneDrive - University College London\Data\2P\';
end

for iExp = 1:nExp
    
 try
    info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(db(iExp).expID(1)));
    catch
     info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(1));

    end    
    folder = buildExpPath(db(iExp));
        load(fullfile(root, folder, 'splitReds'), 'nflatRedImg', 'nflatCherryImg', 'frameG');
    
    [micronsX, micronsY, micronsZ] = ppbox.getPxXYZ(info);
    
    for iPlane = 1: info.nPlanes
        iPlane
        [~, refF, ~] = s2pUtils.getAnalysisRefs(db(iExp).mouse_name, db(iExp).date, db(iExp).expts, iPlane);
        
        
        
        if exist(fullfile(root, folder, refF))
            load(fullfile(root, folder, refF)); % load ROIs data
            yrange = dat.ops.yrange;
            xrange = dat.ops.xrange;
            
            C = nflatCherryImg(:,:,iPlane);
            R = nflatRedImg(:,:,iPlane);
            G = frameG(:,:,iPlane);
             
            C = C(yrange, xrange);
            R = R(yrange, xrange);
                        G = G(yrange, xrange);

            
            [ fovx, fovy ] = ppbox.zoom2fov(info.zoomFactor);
            PxSz = (fovx+fovy)/(2*512);
            
            if redOnly
                allCells{iExp, iPlane} = find([dat.stat.redcell] & [dat.stat.iscell]);
                redCells{iExp, iPlane} = [dat.stat.redcell];
                redCells{iExp, iPlane} = redCells{iExp, iPlane}(allCells{iExp, iPlane});
                map = {dat.stat(allCells).ipix};
                
                [imgR{iExp, iPlane}, imgC{iExp, iPlane}, imgG{iExp, iPlane}, imgM{iExp, iPlane},...
                    bkgR{iExp, iPlane}, bkgC{iExp, iPlane}, redVal{iExp, iPlane}, cherryVal{iExp, iPlane}] = ...
                    s2pUtils.cutRoiImg3Ch(map, R, C, G,  PxSz, cropSize);
                
                coords{iExp, iPlane} = s2pUtils.getROIxyz(map, micronsX(dat.ops.xrange), micronsY(dat.ops.yrange) , micronsZ(dat.ops.yrange, iPlane));
                coords{iExp, iPlane}(:,3) = coords{iExp, iPlane}(:,3) + db(iExp).ObjZ;
                
            else
                
                allCells{iExp, iPlane} = find([dat.stat.iscell]);
                redCells{iExp, iPlane} = [dat.stat.redcell];
                redCells{iExp, iPlane} = redCells{iExp, iPlane}(allCells{iExp, iPlane});
                planeID{iExp, iPlane} = iPlane*ones(size(allCells{iExp, iPlane})); 
                
                nCells = numel(allCells{iExp, iPlane});
                map = {dat.stat(allCells{iExp, iPlane}).ipix};
                
                [imgR{iExp, iPlane}, imgC{iExp, iPlane}, imgG{iExp, iPlane}, imgM{iExp, iPlane},...
                    bkgR{iExp, iPlane}, bkgC{iExp, iPlane}, redVal{iExp, iPlane}, cherryVal{iExp, iPlane}] = ...
                    s2pUtils.cutRoiImg3Ch(map, R, C, G,  PxSz, cropSize);
                
                 coords{iExp, iPlane} = s2pUtils.getROIxyz(map, micronsX(dat.ops.xrange), micronsY(dat.ops.yrange) , micronsZ(dat.ops.yrange, iPlane));
                coords{iExp, iPlane}(:,3) = coords{iExp, iPlane}(:,3) + db(iExp).ObjZ;
                
                
                %                 choices{iExp, iPlane} = s2pUtils.sortROI(imgR{iExp, iPlane}, imgC{iExp, iPlane}, imgM{iExp, iPlane});
                multiplets = mat2cell(makeVec(1:size(imgR{iExp, iPlane},3)),ones(1,size(imgR{iExp, iPlane},3)));
                [~, choices{iExp, iPlane}, ~] = s2pUtils.sortROIred( multiplets, imgR{iExp, iPlane}, imgC{iExp, iPlane}, imgG{iExp, iPlane},imgM{iExp, iPlane});

                %             for iC = 1:nCells
                %                 cherryVal{iExp, iPlane}(iC) = mean(C(dat.stat(redCells(iC)).ipix));
                %                 redVal{iExp, iPlane}(iC) =  mean(R(dat.stat(redCells(iC)).ipix));
                %             end
                
            end
            
            
        end
        
        
    end
    
    
    thisCherryVal = cat(2,cherryVal{iExp, :});
    thisRedVal = cat(2,redVal{iExp, :});
    thisCoords =cat(1,coords{iExp, :});
    thisImgR = cat(3, imgR{iExp, :});
    thisImgC = cat(3,imgC{iExp, :});
    thisImgM = cat(3,imgM{iExp, :});
    thisImgG = cat(3,imgG{iExp, :});

    thisBkgR =cat(2, bkgR{iExp,:});
    thisBkgC =cat(2, bkgC{iExp, :});
    thisAllCells  = cat(2, allCells{iExp, :});
    thisRedCells =cat(2, redCells{iExp, :});
    thisChoices = cat(2, choices{iExp,:});
    
    thisPlaneID = cat(2, planeID{iExp,:});
    
    if redOnly
        save(fullfile(root, folder, 'presynROI_RedOrCherry'), 'thisCherryVal', 'thisRedVal', 'thisCoords', 'thisImgR', ...
            'thisImgC', 'thisImgG', 'thisImgM','thisBkgR', 'thisBkgC', 'thisAllCells', 'thisRedCells', 'thisChoices', '-v7.3');
        
    else
        
        save(fullfile(root, folder, 'ROI_RedOrCherry'), 'thisCherryVal', 'thisRedVal', 'thisCoords', 'thisImgR', ...
            'thisImgC', 'thisImgG', 'thisImgM','thisBkgR', 'thisBkgC', 'thisAllCells', 'thisRedCells', 'thisChoices', 'thisPlaneID',  '-v7.3');
    end
end




% if redOnly
%     save(fullfile(root, db(1).mouse_name, sprintf('presynROI_RedOrCherry_%d', db(1).starterID)), 'cherryVal', 'redVal', 'coords', 'imgR', 'imgC', 'imgM',...
%         'bkgR', 'bkgC', 'allCells', 'redCells', 'choices', '-v7.3');
%     
% else
%     
%     save(fullfile(root, db(1).mouse_name, sprintf('ROI_RedOrCherry_%d', db(1).starterID)), 'cherryVal', 'redVal', 'coords', 'imgR', 'imgC', 'imgM',...
%         'bkgR', 'bkgC', 'allCells', 'redCells', 'choices', '-v7.3');
% end

end