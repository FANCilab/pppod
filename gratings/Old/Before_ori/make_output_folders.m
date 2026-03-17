function out = make_output_folders(targetFolder)
% MAKE_OUTPUT_FOLDERS Create all folders used by the workflow.

    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end

    out = struct();
    out.root = targetFolder;
    out.responsiveness = fullfile(targetFolder, 'responsiveness');
    out.directionTimecourses = fullfile(targetFolder, 'direction_timecourses');
    out.sizeTimecourses = fullfile(targetFolder, 'size_timecourses');
    out.tfTimecourses = fullfile(targetFolder, 'tf_timecourses');
    out.sfTimecourses = fullfile(targetFolder, 'sf_timecourses');
    out.combinedTuning = fullfile(targetFolder, 'combined_tuning');
    out.pairwiseMatrices = fullfile(targetFolder, 'pairwise_response_matrices');

    folderNames = struct2cell(out);
    for iFolder = 1:numel(folderNames)
        folderPath = folderNames{iFolder};
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
    end
end
