function [] = FOOOF_table_create_JAT(partID)

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
        partCD = 'Z:\LossAversion\Patient folders\';
        saveCD = 'Z:\LossAversion\LH_Data\JAT_Data\';

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

        if width(tempEphys) >= 498 % minimum number of samples for pwelch to work % was 128

            % PSD before running fooof - using pwelch as the PSD
            % [psd, freqs] = pwelch(tempEphys,hamming(128), 64, 512, 500);

            % TEST
            % Fs = 500;
            % Nwin = 512;                  % ~1.0 Hz effective resolution (Fs/Nwin)
            % noverlap = round(0.75*Nwin); % 75% overlap
            % nfft = 1024;                 % (optional) denser frequency grid
            % [psd, freqs] = pwelch(tempEphys, hamming(Nwin), noverlap, nfft, Fs);

            % [psd, freqs] = pwelch(tempEphys, hamming(512), 384, 1024, 500);
            % 
            % logP_s = movmedian(log10(psd), 5);   % 5-bin median
            % pxx_s  = 10.^logP_s;

            Fs = 500;
            Nwin = 500;                 % use almost all data
            noverlap = 0;               % overlap doesn't help with one segment
            nfft = 1024;

            [psd,freqs] = pwelch(tempEphys, hamming(Nwin), noverlap, nfft, Fs);
            logP = log10(psd);
            logP_s = movmedian(logP, 3);     % 5-bin median (good for spiky junk)
            pxx_s  = 10.^logP_s;

            plot(freqs,pxx_s)
            xlim([0 100])

            % [psd, freqs] = pwelch(tempEphys,hamming(512), 128, 2048, 500);

            % Transpose, to make inputs row vectors
            % freqs = freqs';
            % psd = psd';
            freqs = double(freqs(:));
            psd   = double(psd(:));

            % FOOOF settings
            % settings = struct();  % Use defaults
            settings = struct(...
                'peak_width_limits', [2, 20], ...
                'max_n_peaks', 5, ...
                'min_peak_height', 0.0, ...
                'peak_threshold', 2.0, ...
                'aperiodic_mode', 'knee', ...
                'verbose', true);

            % settings_py = py.dict;
            % settings_py{'peak_width_limits'} = py.list({2.0, 20.0});   % <-- key fix
            % settings_py{'max_n_peaks'}       = int32(5);
            % settings_py{'min_peak_height'}   = 0.0;
            % settings_py{'peak_threshold'}    = 2.0;
            % settings_py{'aperiodic_mode'}    = 'knee';
            % settings_py{'verbose'}           = true;
            
            f_range = [2, 100]; % changed from 40
            % f_range_py = py.list({double(f_range(1)), double(f_range(2))});  % also do this

            % Run FOOOF
            % fooof_results = fooof(freqs, psd, f_range, settings);
            fooof_results = fooof(freqs, psd, f_range, settings, true);
            % fooof_results = fooof(freqs, psd, f_range_py, settings_py, true);


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