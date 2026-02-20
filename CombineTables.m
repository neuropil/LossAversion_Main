% combines all of the tables across participants into one big table for
% power and frequency without the block numbers 

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
allPartTables = [];

for i = 1:length(dirCLASENames)

tempPath = strcat(partCD, dirCLASENames(i),'\','Tables');
cd(tempPath)

tempName = strcat(dirCLASENames(i), '_allFreqTabs.mat');
load(tempName) % loads as allTables 

allPartTables = [allPartTables; allTables];

cd(partCD)

end % for / i 

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
allPartTablesBlocks = [];

for i = 1:length(dirCLASENames)

    tempPath = strcat(partCD, dirCLASENames(i),'\','Tables');
    cd(tempPath)

    tempName = strcat(dirCLASENames(i), '_allBlocksTab.mat');
    load(tempName) % loads as allTablesBlocks

    allPartTablesBlocks = [allPartTablesBlocks; allTablesBlocks];

    cd(partCD)

end % for / i

%% Combine tables for aperiodic 

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';

end % switch case

cd(partCD)

partIDs = {'CLASE009', 'CLASE018', 'CLASE019', 'CLASE022', 'CLASE029', ...
    'CLASE030', 'CLASE023', 'CLASE024', 'CLASE027', 'CLASE034', ...
    'CLASE001', 'CLASE007', 'CLASE008', 'CLASE026', 'CLASE031', 'CLASE035'};

allPartTAB = [];

for i = 1:length(partIDs)

    tempPartPath = string(strcat(partCD, partIDs(i), '\Aperiodic'));
    cd(tempPartPath)

    tempName = string(strcat(partIDs(i), '_aperiodic.mat'));
    load(tempName) % loads as allTAB

    allPartTAB = [allPartTAB; allTAB];

    cd(partCD)

end % for / i 
