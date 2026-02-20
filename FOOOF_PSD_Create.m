function [] = FOOOF_PSD_Create(partID)

% partID = 'CLASE001';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\Patient folders\';
        saveCD = 'Z:\LossAversion\LH_Data\PSD\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\Patient folders\';
        saveCD = 'Y:\LossAversion\LH_Data\PSD\';

end % switch case

% CD to participant folder
partPath = strcat(partCD, partID, '\', 'ProcessedEphys');
cd(partPath)

% Get directory of names you want
mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, 'allEphys');
dirEphysNames = string(dirNames(dirFilter));

for i = 1:length(dirEphysNames)
    % Load temporary file
    tempName = dirEphysNames(i);
    load(tempName) % Loads as allEphys

    % Average ephys across contacts
    if height(allEphys) > 1
        tempEphys = mean(allEphys);
    else
        tempEphys = allEphys;
    end % If else

    % Run PSD 
    % [psd, freqs] = pwelch(tempEphys,hamming(128), 64, 512, 500);
    [psd, freqs] = pwelch(tempEphys,hamming(256), 250, 512, 500);


    % Transpose, to make inputs row vectors
    freqs = freqs';
    psd = psd';

    % FOOOF settings
    settings = struct();  % Use defaults
    f_range = [1, 40];

    % Run FOOOF
    fooof_results = fooof(freqs, psd, f_range, settings, true);

    PSDtab = table(fooof_results.power_spectrum', fooof_results.freqs', ...
        'VariableNames', {'Power', 'Frequency'});
    
    % save name
    tempSaveName = erase(tempName, "allEphys.mat"); % remove EpochEphys from file name
    tempSaveName = tempSaveName + "PSD_v2"; % Add PSD to the end of the file name

    % CD to save location
    savePath = strcat(saveCD,partID);
    cd(savePath)

    % Save
    save(tempSaveName, "PSDtab");

    % CD back to participant folder
    cd(partPath)

end % for / i 

end % function 