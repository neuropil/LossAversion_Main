function [] = BatchFOOOF_TableGenerate()


% Navigate to Subject list
cd('Z:\LossAversion\Patient folders')
partCD = 'Z:\LossAversion\Patient folders\';
subDir = dir();
subDir2 = struct2table(subDir);
subDir3 = subDir2.name;
subjectLIST = subDir3(~ismember(subDir3,{'.','..','.DS_Store'}));


for ii = 1:length(subjectLIST)

    subjectID = subjectLIST{ii};
    % Check Processed folder
    partPath = strcat(partCD, subjectID , filesep, 'ProcessedEphys');
    if ~exist(partPath,'dir')
        continue
    else
        FOOOF_TABLE_GenerateJAT2(subjectID)
    end




end



end