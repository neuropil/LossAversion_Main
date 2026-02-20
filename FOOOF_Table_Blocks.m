function [] = FOOOF_Table_Blocks(partID)

% partID = 'CLASE001';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\Patient folders\';
        saveCD = 'Z:\LossAversion\LH_Data\PSD\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';
        % saveCD = 'Y:\LossAversion\LH_Data\PSD\';

end % switch case

% CD to participant folder
partPath = strcat(partCD, partID, '\Blocks');
cd(partPath)

% Get directory of names you want
mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, 'CLASE');
dirFOOOFNames = string(dirNames(dirFilter));

% LA Classifications
highLA = {'CLASE009', 'CLASE018', 'CLASE019', 'CLASE022', 'CLASE029', 'CLASE030'};
lowLA = {'CLASE023', 'CLASE024', 'CLASE027', 'CLASE034'};
neutralLA = {'CLASE001', 'CLASE007', 'CLASE008', 'CLASE026', 'CLASE031', 'CLASE035'};

if matches(partID, highLA)
    LAval = string('High');

elseif matches(partID, lowLA)
    LAval = string('Low');

else matches(partID, neutralLA)
    LAval = string('Neutral');
end % if else

allFOOOFBlocks = [];

for fi = 1:length(dirFOOOFNames)

    % Load temporary file
    tempName = dirFOOOFNames(fi);
    load(tempName) % Loads as FOOOFBlocks

    BlockNames = fieldnames(FOOOFBlocks);

    % Temp name
    nameSplit = split(tempName, '_');
    partID = nameSplit(1);
    Hemi = nameSplit(2);
    BA = nameSplit(3);


    for i = 1:length(BlockNames)

        tempBlockName = string(BlockNames{i});

        tempOffset = FOOOFBlocks.(BlockNames{i}).aperiodic_params(1);
        tempExponent = FOOOFBlocks.(BlockNames{i}).aperiodic_params(2);
        tempPSDTab = table(FOOOFBlocks.(BlockNames{i}).power_spectrum', FOOOFBlocks.(BlockNames{i}).freqs', ...
            'VariableNames', {'Power', 'Frequency'});

        tempBlockTab = table(partID, Hemi, BA, LAval, tempBlockName, tempOffset, tempExponent, {tempPSDTab}, ...
            'VariableNames', {'PartID', 'Hemi', 'BrainArea','LAVal', 'Block', 'Offset', 'Exponent', 'PSD'});


        allFOOOFBlocks = [allFOOOFBlocks; tempBlockTab];

    end % for / i

end % for / fi

% save name
tempSaveName = strcat(partID, '_AllFOOOFBlocks.mat');

% Save
save(tempSaveName, "allFOOOFBlocks");


end % function