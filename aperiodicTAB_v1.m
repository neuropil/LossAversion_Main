
function [] = aperiodicTAB_v1(partID)
% partID = 'CLASE009';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        % partCD = 'Z:\LossAversion\Patient folders\';
        saveCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';

    case 'LATERALHABENULA'
        % partCD = 'Y:\LossAversion\Patient folders\';
        saveCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';

end % switch case

savePath = strcat(saveCD, partID, '\Aperiodic');
cd(savePath)

% Get directory of names you want
mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, 'FOOOF_all');
dirEphysNames = string(dirNames(dirFilter));

allTAB = [];

highLA = {'CLASE009', 'CLASE018', 'CLASE019', 'CLASE022', 'CLASE029', 'CLASE030'};
lowLA = {'CLASE023', 'CLASE024', 'CLASE027', 'CLASE034'};
neutralLA = {'CLASE001', 'CLASE007', 'CLASE008', 'CLASE026', 'CLASE031', 'CLASE035'};

for i = 1:length(dirEphysNames)

    % Load temporary file
    tempName = dirEphysNames(i);
    load(tempName) % Loads as fooof_results

    tempOffset = fooof_results.aperiodic_params(1);
    tempExponent = fooof_results.aperiodic_params(2);

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

        tempTab = table(partID, Hemi, BA, tempOffset, tempExponent, LAval, ...
            'VariableNames', {'partID', 'Hemi', 'BrainArea', 'Offset', 'Exponent', 'LAVal'});

        allTAB = [allTAB; tempTab];

end % for / i
% Save name
saveName = strcat(partID, '_aperiodic.mat');

% Save
save(saveName, "allTAB");

end % function