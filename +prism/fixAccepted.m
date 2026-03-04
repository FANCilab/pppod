function [allMulti] = fixAccepted(G, R, rCC, zPx, mimg, starterYX)

neurons = [];
thisRoiPlane = [];
for iImg = 1:numel(rCC)
    thisZ = zPx(:,:,iImg);
    for icc = 1:numel(rCC{iImg})
        
        rCC{iImg}(icc).Centroid(3) =  thisZ(round(rCC{iImg}(icc).Centroid(2)), round(rCC{iImg}(icc).Centroid(1)));
        rCC{iImg}(icc).PixelList(:, 3) = NaN;
        for iPx = 1:size(rCC{iImg}(icc).PixelList(:,1), 1)
            rCC{iImg}(icc).PixelList(iPx, 3) = thisZ(rCC{iImg}(icc).PixelList(iPx,2), rCC{iImg}(icc).PixelList(iPx,1));
        end
    end
    neurons = [neurons; cat(1,rCC{iImg}(:).Centroid)];
    thisRoiPlane = [thisRoiPlane; ones(numel(rCC{iImg}),1)*iImg];
end

rCC = cat(1, rCC{:});

multiplets = findNeuronMultiplets_dev(neurons, 5, 20, [],[],[]);

[nx,ny, nz] = size(G);

for im = 1:numel(multiplets)
    
    for iD = 1: numel(multiplets{im})
        pix = sub2ind([nx ny],rCC(multiplets{im}(iD)).PixelList(:, 2), rCC(multiplets{im}(iD)).PixelList(:, 1));
        frameG = G(:,:, thisRoiPlane(multiplets{im}(iD)));
        frameR = R(:,:, thisRoiPlane(multiplets{im}(iD)));
        [multiG{im}(:,:,iD), multiR{im}(:,:,iD), multiM{im}(:,:,iD)] = s2pUtils.cutRoiImg({pix}, frameG, frameR, [], 35);
    end
end

allMulti.G = multiG;
allMulti.R = multiR;
allMulti.M = multiM;
allMulti.pos = neurons;

end