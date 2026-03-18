function plotTraces(db, planes)


if iscell(db.mouse_name)
    info = ppbox.infoPopulateTempLFR(db.mouse_name{1}, db.date{1}, db.expts{1}(1));
else
    info = ppbox.infoPopulateTempLFR(db.mouse_name, db.date, db.expts(1));
    
end
if nargin <2
    planes = 1:info.nPlanes;
end

for iPlane = planes
    
    [root, refF, ~] = starter.getAnalysisRefs(db.mouse_name, db.date, db.expts, iPlane);
    
    if exist(fullfile(root, refF), 'file')
        
        load(fullfile(root, refF));
        
        goodCellsP = find([dat.stat.iscell]);

        figure; hold on;
        tinit = 0;
        for iStim = 1: numel(dat.Fcell)
            
            F = dat.Fcell{iStim}(goodCellsP, :);
            
            Fneu = dat.FcellNeu{iStim}(goodCellsP, :);
            
            dF{iStim} = estimateNeuropil(F,Fneu);
            
            dF{iStim} =  gaussFilt(dF{iStim}',1);
             
            stdi(:, iStim) = std(dF{iStim}, [], 1);
            
            dF{iStim} =  zscore(dF{iStim});
            
            dF{iStim} =  dF{iStim}(50:200+(round(300*30/info.nPlanes)), :); % 2 in of recording per stim
        
            [nFrames, nN] = size(dF{iStim});
            
            tt{iStim} = (1:nFrames)/(30/info.nPlanes) +tinit;
            
            tinit = max( tt{iStim} ) +50;
            
%             NdF = [];
%             if nFrames >0
%                 for iN = 1:nN
%                     bsl(:,iN) = calcium.baselineF(dF{iStim}(:, iN), 3 , 0);
%                 end
%                 NdF(iN, :) = zscore(dF{iStim}(:, iN)-bsl);
%                 NdF = npRet.my_conv(NdF, 2);
%                 figure;
%                 
%             end
%             figure;
        end
        
       for iStim = 1: numel(dat.Fcell)
            dF{iStim} = bsxfun(@rdivide, bsxfun(@times, dF{iStim}, stdi(:, iStim)'), mean(stdi, 2)');
            PlotDisplacedLFR(tt{iStim}, dF{iStim}, 5); hold on
            xlim([min(cat(2, tt{:})), max(cat(2, tt{:}))])
        end
    end
    
end
end

