function Nframes = getNumberOfFrames(db)

nExp = numel(db.expts);

info = ppbox.infoPopulate(db.mouse_name, db.date, db.expts(1));


for iPlane = 1: info.nPlanes
    [root, refF, ~] = starter.getAnalysisRefs(db.mouse_name, db.date, db.expts, iPlane);
    
    if exist(fullfile(root, refF), 'file')
        load(fullfile(root, refF), 'dat');
        allFrames(iPlane, :) = [dat.ops.Nframes];
        clear dat
    end
end
Nframes(iPlane, :) = allFrames(iPlane, db.expID);

end


