function [F, fTimes, redCells, refFrames, xyzPos, nanROIs] = loadF(db, loadFlag)

dataType = inputname(1);

if contains(dataType, 'starter')
    dataType = 'starterRec';
else
    dataType = 'presColumn';
end

if nargin <2
    loadFlag = 0;
end

saveDir= s2pUtils.getAnalysisRefs(db(1).mouse_name);
try
    saveTo = fullfile(saveDir, sprintf('totaldF_%s_%s_%d.mat', dataType, db(1).mouse_name, db(1).starterID));
catch
    saveTo = fullfile(saveDir, sprintf('totaldF_%s_%s.mat', dataType, db(1).mouse_name));
    
end

if exist(saveTo) && loadFlag
    
    load(saveTo,'F', 'fTimes', 'redCells', 'refFrames', 'xyzPos');
%         load(saveTo,'F', 'fTimes', 'redCells', 'refFrames', 'xyzPos', 'nanROIs');

else
    
    nExp = numel(db);
    
    xyzPos = cell(nExp,1);
    
    for iExp = 1:nExp
        
        info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(db(iExp).expID(1)));
        [micronsX, micronsY, micronsZ, xMicronsPerPixel, yMicronsPerPixel] = ppbox.getPxXYZ(info);
        
        
        for iPlane = 1: info.nPlanes
            
            [root, refF, ~] = s2pUtils.getAnalysisRefs(db(iExp).mouse_name, db(iExp).date, db(iExp).expts, iPlane);
            
            if exist(fullfile(root, refF), 'file')
                
                load(fullfile(root, refF), 'dat');
                isCell = logical([dat.stat.iscell]);
                red = [dat.stat.redcell];
                isRed{iExp, iPlane} = red(isCell);
                
                map = {dat.stat(isCell).ipix};
                
                pos = s2pUtils.getROIxyz(map, micronsX(dat.ops.xrange), micronsY(dat.ops.yrange) , micronsZ(dat.ops.yrange, iPlane));
                
                if isfield(db(iExp), 'ObjZ')
                pos(:,3) = pos(:, 3) + db(iExp).ObjZ;
                end
               
                if isfield(db(iExp), 'starterYX')

                pos(:,1) = pos(:, 1) - db(iExp).starterYX(2)*xMicronsPerPixel;
                pos(:,2) = pos(:, 2) - db(iExp).starterYX(1)*yMicronsPerPixel;
                end
                
                xyzPos{iExp} = cat(1, xyzPos{iExp}, pos);
                
                for expID = 1: numel(db(iExp).expID)
                    
                    info = ppbox.infoPopulateTempLFR(db(iExp).mouse_name, db(iExp).date, db(iExp).expts(db(iExp).expID(expID)));
                    Fcell = dat.Fcell{db(iExp).expID(expID)}(isCell, :); % nN*nT
                    Fneu = dat.FcellNeu{db(iExp).expID(expID)}(isCell, :);% nN*nT
                    dF = estimateNeuropil_LFR(Fcell,Fneu); % correct neuropil and then subtracts 5th prctile of trace to equalize eventual drifts across exps
                    
                    allF{iExp, iPlane, expID} = npRet.my_conv(dF, 1)'; % nT*nN
%                     allF{iExp, iPlane, expID} = dF'; % nT*nN

                    [nFrames(iExp, iPlane, expID), ~] = size(allF{iExp, iPlane, expID}); % nT * nNeurons
                    if nFrames(iExp, iPlane, expID) >0
                        planeFrames = iPlane:info.nPlanes:(nFrames(iExp, iPlane, expID)*info.nPlanes);
                        frameTimes{iExp, iPlane, expID} = makeVec(ppbox.getFrameTimes(info, planeFrames));
                        allF{iExp, iPlane, expID} =  allF{iExp, iPlane, expID}(1:numel(frameTimes{iExp, iPlane, expID}) , :);
                    end
                    
                    
                end
            end
        end
        
        
    end
    
    %% concatenate neurons recorded across planes across expID
    F = cell(nExp,1);
    fTimes = cell(nExp,1);
    
    for iExp = 1:nExp
        allPlanes = cell(numel(db(iExp).expID),1);
        for expID = 1: numel(db(iExp).expID)
            
            emptyPlanes = cellfun(@numel, frameTimes(iExp, :, expID));
            fullPlanes = find(emptyPlanes);
            refPlane = fullPlanes(1);
            endPlane = fullPlanes(end);
            dummyT = frameTimes{iExp, refPlane, expID}(1:numel(frameTimes{iExp, endPlane, expID}));
            refFrames{iExp, expID} = dummyT;
            if expID >1
                init = fTimes{iExp}(end);
            else
                init = 0;
            end
            
            %             for iPlane = 1: min(info.nPlanes, size(allF,2))
            for iPlane = 1:size(allF,2)
                
                if ~isempty(allF{iExp, iPlane, expID})
                    dummy = allF{iExp, iPlane, expID};
                    interpF{iExp, iPlane, expID} = interp1(frameTimes{iExp, iPlane, expID}, dummy, dummyT, 'linear', 0);
                else
                    interpF{iExp, iPlane, expID} = [];
                end
            end
            
            dummyT = dummyT - dummyT(1) + mean(diff(dummyT)) + init;
            allPlanes{expID} = cat(2,interpF{iExp, :, expID}); % nT* nNeurons, neurons catenated across planes
            fTimes{iExp} = cat(1, fTimes{iExp}, dummyT);
            
            
        end
        
        redCells{iExp} = cat(2,isRed{iExp, :});
        F{iExp} = cat(1, allPlanes{:}); % nT* nNeurons, neurons catenated across planes, and time points catenated across expID
    end
    
    %remove multiplets and accidentally bad ROIs
    
    for iExp = 1:nExp
        
%         multiplets = findNeuronMultiplets_dev(xyzPos{iExp}, 5, 20,[], F{iExp}, 0.5);
        multiplets = findNeuronMultiplets_dev(xyzPos{iExp}, 1, 1,[], F{iExp}, 1);
        
        score = max(F{iExp},[],1);
        nN = size(F{iExp},2);
        isUnique = zeros(nN,1);
        for im = 1: numel(multiplets)
            [~,chosen] = max(score(multiplets{im})) ;
            isUnique(multiplets{im}(chosen)) = 1;
        end
        
        F{iExp}(:, ~isUnique) = [];
        xyzPos{iExp}(~isUnique, :) = [];
        redCells{iExp}(~isUnique) = [];
                
        
        [~, nanN] = find(isnan(F{iExp}));
        nanN = unique(nanN);
        
        F{iExp}(:, nanN) = [];
        redCells{iExp}(nanN) = [];
        xyzPos{iExp}(nanN, :) = [];
        nanROIs{iExp} = nanN;
    end
    
    save(saveTo, 'F', 'fTimes', 'redCells', 'refFrames', 'xyzPos', 'nanROIs', '-v7.3');
    
end


end