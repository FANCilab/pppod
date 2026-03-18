function [expFolder, refF, refNeu, refSVD, refRaw] = getAnalysisRefs(mouse_name, date, expts, plane)

root = 'D:\OneDrive - University College London\Data\2P\';

if iscell(mouse_name)
    
    mouse_name = mouse_name{1};
    
    if nargin ==2
        date = date{1};
    elseif nargin >2 
       date = date{1};

        expts = cat(2, expts{:});
    end
     
end

if nargin <2
    expFolder= fullfile(root, mouse_name);
elseif nargin <3
    expFolder= fullfile(root, mouse_name, date);
else
    expStr =  sprintf('%d_', expts);
    expStr = expStr(1:end-1);
    expFolder = fullfile(root, mouse_name, date, expStr);
end

if nargin >3
    refRaw = ['F_', mouse_name,'_', date,'_', 'plane', num2str(plane), '.mat'];
    refF = ['F_', mouse_name,'_', date,'_', 'plane', num2str(plane), '_proc.mat'];
    refNeu = ['NEU_', mouse_name,'_', date,'_', 'plane', num2str(plane), '.mat'];
    refSVD = ['SVD_', mouse_name,'_', date,'_', 'plane', num2str(plane), '.mat'];
end

end