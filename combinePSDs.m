function [] = combinePSDs(partID)

% partID = 'CLASE001';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\Patient folders\';
        saveCD = 'Z:\LossAversion\LH_Data\PSD\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\PSD\';
        % saveCD = 'Y:\LossAversion\LH_Data\PSD\';

end % switch case

% CD to participant folder
partPath = strcat(partCD, partID, '\');
cd(partPath)

% Get directory of names you want
mdir = dir;
dirNames = {mdir.name};
% dirFilter = contains(dirNames, 'CLASE');
dirFilter = contains(dirNames, '_v2');
dirPSDNames = string(dirNames(dirFilter));

allTABPSD = [];

highLA = {'CLASE009', 'CLASE018', 'CLASE019', 'CLASE022', 'CLASE029', 'CLASE030'};
lowLA = {'CLASE023', 'CLASE024', 'CLASE027', 'CLASE034'};
neutralLA = {'CLASE001', 'CLASE007', 'CLASE008', 'CLASE026', 'CLASE031', 'CLASE035'};

for i = 1:length(dirPSDNames)

     % Load temporary file
    tempName = dirPSDNames(i);
    load(tempName) % Loads as PSDtab 

    % Temp name
    nameSplit = split(tempName, '_');
    partID = nameSplit(1);
    Hemi = nameSplit(2);
    BA = nameSplit(3);

      if matches(partID, highLA)
        LAval = string('High');

    elseif matches(partID, lowLA)
        LAval = string('Low');

    else matches(partID, neutralLA)
        LAval = string('Neutral');

    end % if else

        tempTab = table(partID, Hemi, BA, {PSDtab}, LAval, ...
            'VariableNames', {'partID', 'Hemi', 'BrainArea', 'PSD', 'LAVal'});

        allTABPSD = [allTABPSD; tempTab];


end % for / i 

% Save name
saveName = strcat(partID, '_PSDTab_v2.mat');

% Save
save(saveName, "allTABPSD");

end % function 