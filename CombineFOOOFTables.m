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
FOOOFBlockTab = [];

for i = 1:length(dirCLASENames)

    tempPath = strcat(partCD, dirCLASENames(i), '\Blocks');
    cd(tempPath)

    tempName = strcat(dirCLASENames(i), '_AllFOOOFBlocks.mat');
    load(tempName) % loads as allFOOOFBLocks

    FOOOFBlockTab = [FOOOFBlockTab; allFOOOFBlocks];

    cd(partCD)

end % for / i
