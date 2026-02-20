
function [] = baseline_runFOOOF(partID)
% partID = 'CLASE001';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';

end % switch case

% CD to participant folder
partPath = strcat(partCD, partID, '\', 'Baseline');
cd(partPath)

% Get directory of names you want
mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, '_baseline');
dirEphysNames = string(dirNames(dirFilter));



for i = 1:length(dirEphysNames)

    % Load temporary file
    tempName = dirEphysNames(i);
    load(tempName) % Loads as baseline

    ephys = baseline.CombinedBaseline.CombinedEphys;

    % Average ephys across contacts
    if height(ephys) > 1
        tempEphys = mean(ephys);
    else
        tempEphys = ephys;
    end % If else

    % ---- Run PSD and FOOOF ---- %

    % PSD before running fooof - using pwelch as the PSD
    [psd, freqs] = pwelch(tempEphys,hamming(128), 64, 512, 500);

    % Transpose, to make inputs row vectors
    freqs = freqs';
    psd = psd';

    % FOOOF settings
    settings = struct();  % Use defaults
    f_range = [1, 40];

    % Run FOOOF
    fooof_results = fooof(freqs, psd, f_range, settings, true);

    % save name
    tempSaveName = erase(tempName, "_baseline.mat"); % remove EpochEphys from file name
    tempSaveName = tempSaveName + "_FOOOF_baseline.mat"; % Add PSD to the end of the file name

    % CD to save location
    cd(partPath)

    % Save
    save(tempSaveName, "fooof_results");


end % for / i

end % function 