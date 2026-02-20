
partID = 'CLASE001';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';
        % saveCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';
        % saveCD = '';

end % switch case

partPath = strcat(partCD, partID);
cd(partPath)

mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, '.mat');
dirEpochNames = string(dirNames(dirFilter));

for i = 1:length(dirEpochNames)

    % Load file
    tempName = dirEpochNames(i);
    load(tempName) % loads as FOOOF_tab


    % Create repeating block/ trial number vectors
    trialNumVec = [repmat(1, 1, 81), repmat(2, 1, 81), repmat(3, 1, 81), ...
        repmat(4, 1, 81), repmat(5, 1, 81)]'; % trial number vector

    trialNumVec = table(trialNumVec, 'VariableNames', {'BlockNum'}); % make it a table

    FOOOF_tab = [FOOOF_tab, trialNumVec]; % Combine tables

    % Extract the rows from the table that are only Loss Aversion
    FOOOF_tab = FOOOF_tab(FOOOF_tab.LA,:);

    % Remove empty cells
    lastColumn = FOOOF_tab{:, 6};
    % Identify rows where the last column is empty
    emptyRows = cellfun(@isempty, lastColumn);
    % Remove those rows
    FOOOF_tab(emptyRows,:) = [];

    % Create subtables for each block
    subTab_Block1 = FOOOF_tab(FOOOF_tab.BlockNum == 1, :);
    subTab_Block2 = FOOOF_tab(FOOOF_tab.BlockNum == 2, :);
    subTab_Block3 = FOOOF_tab(FOOOF_tab.BlockNum == 3, :);
    subTab_Block4 = FOOOF_tab(FOOOF_tab.BlockNum == 4, :);
    subTab_Block5 = FOOOF_tab(FOOOF_tab.BlockNum == 5, :);


    % Run function on each subblock
    subBlock1 = processBlock(subTab_Block1, 1);
    subBlock2 = processBlock(subTab_Block2, 2);
    subBlock3 = processBlock(subTab_Block3, 3);
    subBlock4 = processBlock(subTab_Block4, 4);
    subBlock5 = processBlock(subTab_Block5, 5);

    %
    allBlocks = [subBlock1; subBlock2; subBlock3; subBlock4; subBlock5];

    % Temp name
    nameSplit = split(tempName, '_');
    partID = nameSplit(1);
    Hemi = nameSplit(2);
    BA = nameSplit(3);

    % Save name 
    saveName = strcat(partID, '_', Hemi, '_', BA, '_blocks.mat');

    % Repeat identifiers for length fo the table
    partID = repmat(partID, height(allBlocks),1);
    Hemi = repmat(Hemi, height(allBlocks),1);
    BA = repmat(BA, height(allBlocks),1);

    % Create ID table
    IDtab = table(partID, Hemi, BA, 'VariableNames',{'partID', 'Hemi', 'BrainArea'});

    % Combine freqTabs and ID table
    partTab = [IDtab, allBlocks];

    % Save loc
    savePath = strcat(partPath + "\Tables");
    cd(savePath)

    % save table
    save(saveName, "partTab");

    % CD back to fooof data
    cd(partPath)


end % for / i


%%

% partID = 'CLASE026';
% -- Save tables -- % 
% CD into the Tables folder 
% cd(savePath)

% Concatenate all individual tables into one big one and save it as well as
% a csv file of the concatenated tables 
allTablesBlocks = []; % Placeholder for the concatenated table

mdir = dir;
fdir = {mdir.name};
dirPTnames = fdir(contains(fdir,'blocks')); % List of part folder names

for i = 1:length(dirPTnames)

    load(dirPTnames{i}); % loads as partTab 

    allTablesBlocks = [allTablesBlocks; partTab];

end % for / i 

% Save 
% saveName = strcat(partID, '_allFreqTabs.mat');
saveName = 'CLASE035_allBlocksTab.mat';

% save table
save(saveName, "allTablesBlocks");

% Save csv 
% tableSaveName = strcat(partID, '_allFreqTabs.csv');
% writetable(allTables, tableSaveName);
% writetable(allTables,'CLASE035_allFreqTabs.csv');