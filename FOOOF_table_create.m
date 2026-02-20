function [] = FOOOF_table_create(partID)

% partID = 'CLASE001';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\Patient folders\';
        saveCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\Patient folders\';
        saveCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';

    case 'DESKTOP-I5CPDO7' % JAT WORK PC
        partCD = 'Y:\LossAversion\Patient folders\';
        saveCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';

end % switch case

% CD to participant folder
partPath = strcat(partCD, partID, '\', 'ProcessedEphys');
cd(partPath)

% Get directory of names you want
mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, 'EpochEphys');
dirEpochNames = string(dirNames(dirFilter));

for i = 1:length(dirEpochNames)

    % Load temporary file
    tempName = dirEpochNames(i);
    load(tempName,'EpochEphysTab')

    % Create new FOOOF table
    FOOOF_tab = EpochEphysTab(:, {'EpochID', 'LA', 'OutcomeGain', 'OutcomeLoss', 'OutcomeNeutral'}); % Use EpochEphysTab data

    FOOOF_values = cell(height(FOOOF_tab),1); % will replace with FOOOF output
    FOOOF_tab = addvars(FOOOF_tab, FOOOF_values, 'NewVariableNames','FOOOFoutput');

    for ei = 1:height(EpochEphysTab)

        % Temporary ephys - each epoch
        tempEphys = EpochEphysTab.Ephys{ei};

        % Mean Temp Ephys
        if height(tempEphys) > 1
            tempEphys = mean(tempEphys);
        else
            tempEphys = tempEphys;
        end % if else


        % ---- Run PSD and FOOOF ---- %

        if width(tempEphys) >= 128 % minimum number of samples for pwelch to work

            % PSD before running fooof - using pwelch as the PSD
            [psd, freqs] = pwelch(tempEphys,hamming(128), 64, 512, 500);

            % Transpose, to make inputs row vectors
            % freqs = freqs';
            % psd = psd';
            freqs = double(freqs(:));
            psd   = double(psd(:));

            % FOOOF settings
            settings = struct();  % Use defaults
            f_range = [1, 40];

            % Run FOOOF
            % fooof_results = fooof(freqs, psd, f_range, settings);
            fooof_results = fooof(freqs, psd, f_range, settings, true);

            % Add fooof_results to their row in the FOOOF_tab
            FOOOF_tab.FOOOFoutput{ei} = fooof_results;

        else
            continue
        end % if else


    end % for / ei

    % save name
    tempSaveName = erase(tempName, "EpochEphys.mat"); % remove EpochEphys from file name
    tempSaveName = tempSaveName + "FOOOF"; % Add PSD to the end of the file name

    % CD to save location
    savePath = strcat(saveCD,partID);
    cd(savePath)

    % Save
    save(tempSaveName, "FOOOF_tab");

    % CD back to participant folder
    cd(partPath)


end % for / i

end % function