function [choices] = sort_measureROIRed(db, cropSize)

if nargin <2
    cropSize = 25;
end

if nargin <3
    redOnly = 0;
end
nExp = numel(db);

root = 'D:\OneDrive - University College London\Data\2P\';

for iExp = 1:nExp
    
    try
    info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(db(iExp).expID(1)));
    catch
     info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(1));

    end
    
    folder = buildExpPath(db(iExp));
    try
        load(fullfile(root, folder, 'splitReds'), 'nflatRedImg', 'nflatCherryImg', 'frameG');
    catch

        [frameG, frameR] = s2pUtils.loadMultiPlaneMeanStack(db(iExp));
        [nflatRedImg] = prism.bleedCure(frameR, frameG);
        nflatCherryImg = zeros(size(nflatRedImg));

    end
    [micronsX, micronsY, micronsZ] = ppbox.getPxXYZ(info);
    
    for iPlane = 1: info.nPlanes
        iPlane
        [~, refF, ~] = starter.getAnalysisRefs(db(iExp).mouse_name, db(iExp).date, db(iExp).expts, iPlane);
        
        
        
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
            
            allCells{iExp, iPlane} = find([dat.stat.redcell] & [dat.stat.iscell]);
            redCells{iExp, iPlane} = [dat.stat.redcell];
            redCells{iExp, iPlane} = redCells{iExp, iPlane}(allCells{iExp, iPlane});
            
            map = {dat.stat(allCells{iExp, iPlane}).ipix};
                if ~isempty(map)
                    
              
            [imgR{iExp, iPlane}, imgC{iExp, iPlane}, imgG{iExp, iPlane}, imgM{iExp, iPlane},...
                    bkgR{iExp, iPlane}, bkgC{iExp, iPlane}, bkgG{iExp, iPlane}, ...
                    redVal{iExp, iPlane}, cherryVal{iExp, iPlane}, gcampVal{iExp, iPlane}] = ...
                    s2pUtils.cutRoiImg3Ch(map, R, C, G, PxSz, cropSize);
                else
                    imgR{iExp, iPlane} =[];
                    imgC{iExp, iPlane} =[];
                    imgM{iExp, iPlane} =[];
                    imgG{iExp, iPlane} =[];

                    bkgR{iExp, iPlane} = [];
                    bkgC{iExp, iPlane} = [];
                    bkgG{iExp, iPlane} = [];

                    bkgM{iExp, iPlane} = [];
                    redVal{iExp, iPlane} = [];
                    cherryVal {iExp, iPlane} = [];
                    gcampVal {iExp, iPlane} = [];

                end
            coords{iExp, iPlane} = starter.getROIxyz(map, micronsX(dat.ops.xrange), micronsY(dat.ops.yrange) , micronsZ(dat.ops.yrange, iPlane));
            coords{iExp, iPlane}(:,3) = coords{iExp, iPlane}(:,3) + db(iExp).ObjZ;
        else
            
            coords{iExp, iPlane} =[];
            allCells{iExp, iPlane} = [];
            redCells{iExp, iPlane} = [];
            imgR{iExp, iPlane} =[];
            imgC{iExp, iPlane} =[];
            imgM{iExp, iPlane} =[];
            imgG{iExp, iPlane} =[];         
            bkgR{iExp, iPlane} = [];
            bkgC{iExp, iPlane} = [];
            bkgG{iExp, iPlane} = [];     
            bkgM{iExp, iPlane} = [];
            redVal{iExp, iPlane} = [];
            cherryVal {iExp, iPlane} = [];
            gcampVal {iExp, iPlane} = [];

        end
    end
    neurons = cat(1, coords{iExp, :});

    R = cat(3, imgR{iExp, :});
    C = cat(3, imgC{iExp, :});
    M = cat(3, imgM{iExp, :});
    G =  cat(3, imgG{iExp, :});
    multiplets{iExp} = findNeuronMultiplets_dev(neurons, 5, 20, [],[],[]);  
    
%     if exist(fullfile(root, folder, 'presynROI_RedOrCherry.mat'))
%         
%      thisImgG = imgG{iExp, :};
%      save(fullfile(root, folder, 'presynROI_RedOrCherry'),'thisImgG', '-append');
%   
%     else
    [n2keep{iExp}, choices{iExp}, keeper{iExp}] = s2pUtils.sortROIred(multiplets{iExp}, R, C,G, M);
    
    accum = 0;
    for iPlane = 1:info.nPlanes
        thisPlaneN = size(coords{iExp,iPlane},1);
        bestDuplic{iExp, iPlane} = zeros(thisPlaneN,1);
        bestDuplic{iExp, iPlane} = n2keep{iExp}(accum+1:accum+thisPlaneN);
        accum = accum + thisPlaneN;
    end
            
        
    thisCherryVal = cat(2, cherryVal{iExp, :});
    thisRedVal = cat(2, redVal{iExp, :});
    thisCoords =cat(1,coords{iExp, :});
    thisImgR = cat(3, imgR{iExp, :});
    thisImgC = cat(3,imgC{iExp, :});
    thisImgM = cat(3,imgM{iExp, :});
    thisImgG = cat(3,imgG{iExp, :});

    thisBkgR =cat(2, bkgR{iExp,:});
    thisBkgC =cat(2, bkgC{iExp, :});
    thisAllCells  = cat(2, allCells{iExp, :});
    thisRedCells =cat(2, redCells{iExp, :});
    
    thisChoices = cat(2, choices{iExp});
    thisKeeper = cat(2,keeper{iExp});
    thisBestDuplic = cat(1, bestDuplic{iExp,:});
    thisMultiplets = multiplets{iExp};
    thisN2Keep = n2keep{iExp};
           
   save(fullfile(root, folder, 'presynROI_RedOrCherry'), 'thisCherryVal', 'thisRedVal', 'thisCoords', 'thisImgG', 'thisImgR', ...
            'thisImgC', 'thisImgM','thisBkgR', 'thisBkgC', 'thisAllCells', 'thisRedCells', 'thisChoices', ...
            'thisKeeper','thisBestDuplic','thisMultiplets','thisN2Keep', 'thisPlaneN', '-v7.3');
        
%     end
end


% if exist(fullfile(root, db(1).mouse_name, sprintf('presynROI_RedOrCherry_%d.mat', db(1).starterID)))
%        save(fullfile(root, db(1).mouse_name, sprintf('presynROI_RedOrCherry_%d', db(1).starterID)),'imgG',  '-append');
% 
% else
% 
%     save(fullfile(root, db(1).mouse_name, sprintf('presynROI_RedOrCherry_%d', db(1).starterID)), 'cherryVal', 'redVal', 'coords','imgG',  'imgR', 'imgC', 'imgM',...
%         'bkgR', 'bkgC', 'allCells', 'redCells', 'choices', 'keeper','bestDuplic','multiplets','n2keep','-v7.3');
%     
% end
    
end