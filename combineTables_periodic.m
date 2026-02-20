%%
partID = 'CLASE034';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';
        % saveCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';
        % saveCD = '';

end % switch case

% -- Save tables -- % 
% CD into the Tables folder 
savePath = strcat(partCD, partID, '\', 'Tables\');
cd(savePath)

% Concatenate all individual tables into one big one and save it as well as
% a csv file of the concatenated tables 
allTablesFreq = []; % Placeholder for the concatenated table

mdir = dir;
fdir = {mdir.name};
dirPTnames = fdir(contains(fdir,'periodFreq')); % List of part folder names



for i = 1:length(dirPTnames)

    load(dirPTnames{i}); % loads as partTabFreq 

    allTablesFreq = [allTablesFreq; partTabFreq];

end % for / i 

% Save 
saveName = strcat(partID, '_allPeriodicEpochTabs.mat');
% saveName = 'CLASE035_allPeriodicEpochTabs.mat';

% save table
save(saveName, "allTablesFreq");

% Save csv 
% tableSaveName = strcat(partID, '_allFreqTabs.csv');
% writetable(allTables, tableSaveName);
% writetable(allTables,'CLASE035_allFreqTabs.csv');


%%
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
allPartTablesFreq = [];

for i = 1:length(dirCLASENames)

tempPath = strcat(partCD, dirCLASENames(i),'\','Tables');
cd(tempPath)

tempName = strcat(dirCLASENames(i), '_allPeriodicEpochTabs.mat');
load(tempName) % loads as allTablesFreq 

allPartTablesFreq = [allPartTablesFreq; allTablesFreq];

cd(partCD)

end % for / i 