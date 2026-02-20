%% combines all of the tables across participants into one big table for
% power and frequency WITH the block numbers 

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';

end % switch case

cd(partCD)

% Get directory of names you want
mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, 'CLASE');
dirCLASENames = string(dirNames(dirFilter));

% Empty Table Holder 
allPartBaseline = [];

for i = 1:length(dirCLASENames)

    tempPath = strcat(partCD, dirCLASENames(i),'\','Baseline');
    cd(tempPath)

    tempName = strcat(dirCLASENames(i), '_allFOOOF.mat');
    load(tempName) % loads as allBaselineTAB

    allPartBaseline = [allPartBaseline; allBaselineTAB];

    cd(partCD)

end % for / i