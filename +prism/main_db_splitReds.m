%% create and save ImagingPlane file for each experiment;

clear;
% npRet.FR141_dbNpRF_1;
npRet.FR141_dbNpRF_2;
% npRet.FR140_dbNpRF_1;
% npRet.FR143_dbNpRF_1;


%% Split mCherry and dsRed fluorescence

%FR141_2
rectifySource = [1 1 0 0 0]; % rotate sources to make mCherry vertical
forceRatio = [0 0 0 0 0]; % force sources coeefficients to be the predicted ones
zLimit = Inf; 
useAllPx = [0 0 0 0 0]; %plane 20 1, plane 26 = 1 

%FR141_1
rectifySource = [0 1 0 0 0 0 0 0 0];
forceRatio = [0 0 0 0 0 0 1 1 1];
zLimit = 510;
useAllPx = [0 0 0 0 0 0 1 1 0];

%FR140
rectifySource = [1 0 0 0 0]; 
forceRatio = [0 0 0 0 0];
zLimit = Inf; 
useAllPx = [0 0 0 1 1];

%FR143
rectifySource = [0 0 0 0 0]; 
forceRatio = [0 0 0 0 0]; 
zLimit = Inf; 
useAllPx = [0 0 0 0 1];
%%
for iDb = 1:numel(db)
    regrOps(iDb).rectifySource = rectifySource(iDb);
    regrOps(iDb).forceRatio = forceRatio(iDb);
    regrOps(iDb).zLimit = zLimit;
    regrOps(iDb).useAllPx = useAllPx(iDb);
    
    prism.db_splitReds(db(iDb), regrOps(iDb))
end

%% process datasets with simultaneous pre and post neurons
clear;
corr.make_db_presyn_postsyn;
nDb = numel(db);

rectifySource = [1 1 1 0 1 0 1 0 0 0]; % rotate sources to make mCherry vertical
forceRatio = [0 0 0 0 0 0 0 0 0 0 ]; % force sources coeefficients to be the predicted ones
zLimit = Inf; 
useAllPx = [0 0 0 0 0 0 0 0 0 0 ]; %plane 20 1, plane 26 = 1 

%%

for iDb = 7 %numel(db)
    regrOps(iDb).rectifySource = rectifySource(iDb);
    regrOps(iDb).forceRatio = forceRatio(iDb);
    regrOps(iDb).zLimit = zLimit;
    regrOps(iDb).useAllPx = useAllPx(iDb);
    
    prism.db_splitReds(db(iDb), regrOps(iDb))
end


%%

iDb = 1; % 1 2 4 6 7 8
choices = s2pUtils.sort_measureROIRed(db(iDb));  


choices = s2pUtils.measureROIRed(db(iDb));   

%%
root = 'D:\OneDrive - University College London\Data\2P\';

if redFlag
    if exist(fullfile(saveDir, sprintf('ROI_presynRedOrCherry_%d.mat', db(1).starterID)))     
        load(fullfile(saveDir, sprintf('ROI_presynRedOrCherry_%d.mat', db(1).starterID)), 'choices');
    else
        choices = s2pUtils.sort_measureROIRed(db);   
    end
else
    if exist(fullfile(saveDir, sprintf('ROI_RedOrCherry_%d.mat', db(1).starterID)))     
        load(fullfile(saveDir, sprintf('ROI_RedOrCherry_%d.mat', db(1).starterID)), 'choices');
    else
        choices = s2pUtils.measureROIRed(db);   
    end
end


