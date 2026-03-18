function [allG, allR, allM] = plotROIcrops(imgG, imgR, imgM)

[nExp, nPlanes] = size(imgG);

allR = [];
allG = [];
allM = [];
nCells = 0;

for iExp = 1:nExp
    for iPlane = 1:nPlanes
        if ~isempty(imgG{iExp, iPlane})
        nCells = nCells + size(imgG,3);
        allR = cat(3, allR, imgR{iExp, iPlane});
        allG = cat(3, allG, imgG{iExp, iPlane});
        allM = cat(3, allM, imgM{iExp, iPlane});
        end
    end
end


nCells = size(allG,3);

nFigs  = ceil(nCells/16);

nCols = 8;
nRows = 4;
acc = 0;
for iF = 1: nFigs
    
    figure;
    
    for iCell = 1:16
        if acc < nCells
            acc = acc+1;
            panel = 1 + (iCell-1)*2;
            ax(panel) = subplot(nRows, nCols, panel);
            imagesc(imadjust(allG(:,:,acc))); axis image;
            colormap(ax(panel), gray); hold on
            contour(allM(:,:,acc), 'r','Linewidth', 0.5);
            formatAxes
            %            set(ax(panel), 'Xtick', [], 'Ytick', [])
            set(ax(panel), 'visible', 'off')
            
            ax(panel+1) = subplot(nRows, nCols, panel + 1);
            imagesc(imadjust(allR(:,:,acc))); axis image
            colormap(ax(panel+1), WhiteRed);
            formatAxes
            %            set(ax(panel+1), 'Xtick', [], 'Ytick', [], 'axes', 'off')
            set(ax(panel+1), 'visible', 'off')
            
            
        end
    end
end
figure;
ag = subplot(1,2,1);
imagesc(mean(allG, 3)); axis image
colormap(ag, gray);
formatAxes
set(ag, 'visible', 'off')
title('Mean presyn G')

ar = subplot(1,2,2);
imagesc(mean(allR, 3)); axis image
colormap(ar, WhiteRed);
formatAxes
set(ar , 'visible', 'off')
title('Mean presyn R')

end