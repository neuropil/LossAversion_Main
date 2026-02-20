%% combines all of the tables across participants into one big table for
% power and frequency WITH the block numbers 

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\LH_Data\PSD\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\PSD\';

end % switch case

cd(partCD)

% Get directory of names you want
mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, 'CLASE');
dirCLASENames = string(dirNames(dirFilter));

% Empty Table Holder 
allPSDtab = [];

for i = 1:length(dirCLASENames)

    tempPath = strcat(partCD, dirCLASENames(i));
    cd(tempPath)

    tempName = strcat(dirCLASENames(i), '_PSDTab_v2.mat');
    load(tempName) % loads as allTABPSD

    allPSDtab = [allPSDtab; allTABPSD];

    cd(partCD)

end % for / i
