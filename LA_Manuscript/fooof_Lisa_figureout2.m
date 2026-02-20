function [allFoooftab] = fooof_Lisa_figureout2(subjectID)

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\Patient folders\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\Patient folders\';

    case 'DESKTOP-I5CPDO7' % JAT WORK PC
        partCD = 'Z:\LossAversion\Patient folders\';
        saveCD = 'Z:\LossAversion\LH_Data\JAT_FoooFtables\';


end % switch case

% CD to participant folder
partPath = strcat(partCD, subjectID , filesep, 'ProcessedEphys');
cd(partPath)

mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, 'EpochEphys');
dirEpochNames = string(dirNames(dirFilter));

% Table with Brain Area , Epoch type , fooof struct , subject , trial type

allFoooftab = table;

for i = 1:numel(dirEpochNames)

    % Load temporary file
    tempName = dirEpochNames(i);

    % Extract name
    tmpName = split(dirEpochNames{i},'_');

    useName = [tmpName{2},'_',tmpName{3}];

    load(tempName,'EpochEphysTab')

    % Create new FOOOF table
    FOOOF_tab = EpochEphysTab(:, {'EpochID', 'LA', 'OutcomeGain', 'OutcomeLoss', 'OutcomeNeutral'}); % Use EpochEphysTab data

    FOOOF_values = cell(height(FOOOF_tab),1); % will replace with FOOOF output
    FOOOF_tab = addvars(FOOOF_tab, FOOOF_values, 'NewVariableNames','FOOOFoutput');

    % epochLAtable = EpochEphysTab(FOOOF_tab.LA,:);

    % tmpROW = EpochEphysTab(find(FOOOF_tab.LA & matches(FOOOF_tab.EpochID,'outcome'),1,'first'),:);

    for eii = 1:height(EpochEphysTab)

            tempEphys = EpochEphysTab.Ephys{eii};

            if length(tempEphys(1,:)) < 499
                continue
            else
                fooofRes = testPSD_fooof(tempEphys);
                FOOOF_tab.FOOOFoutput{eii} = fooofRes;
            end
    end

    FOOOF_tab.brainA = repmat({useName},height(FOOOF_tab),1);

    allFoooftab = [allFoooftab ; FOOOF_tab];

end

% clean up empty columns
allFoooftab = allFoooftab(~cellfun(@(x) isempty(x), allFoooftab.FOOOFoutput, 'UniformOutput', true),:);

% subject ID
allFoooftab.SubID = repmat({subjectID},height(allFoooftab),1);

cd(saveCD)
save([subjectID,'_FooofTab.mat'],"allFoooftab")

end



function [fooof_results] = testPSD_fooof(inputLFP)

% inputLFP = detrend(inputLFP,'linear');   % remove slow trends
% inputLFP = inputLFP - mean(inputLFP);           % remove DC offset
% hpFilt = designfilt('highpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency',1,'SampleRate',500);
% inputLFP = filtfilt(hpFilt,inputLFP);
% [pxx,f] = pmtm(inputLFP,2.5,1024,500);
% plot(f,log10(pxx))
% xlim([1 100])

% Inputs you need:
%   x        : 1-second signal vector
%   Fs       : 500
%   f_range  : e.g. [3 45] for FOOOF
%   settings : MATLAB struct (your settings)

Fs = 500;
f_range = [3 45];

allPSD = zeros(1001,height(inputLFP));
for ii = 1:height(inputLFP)

    % --- preprocess (matches what improved your PSD) ---
    inputLFPi = detrend(inputLFP(ii,:),'linear');
    inputLFPi = inputLFPi - mean(inputLFPi);
    hp = designfilt('highpassiir','FilterOrder',4,...
        'HalfPowerFrequency',2,'SampleRate',Fs);
    inputLFPi = filtfilt(hp,inputLFPi);

    % --- PSD (multitaper recommended for 1-s) ---
    % nfft = 1024;
    % [psd, freqs] = pmtm(inputLFP, 3, nfft, Fs);

    [allPSD(:,ii), freqs] = pspectrum(inputLFPi,Fs,"power","FrequencyResolution",2,"Leakage",0.8,'FrequencyLimits',[2 55]);

    % --- restrict to f_range for plotting ---
    % idx = freqs >= f_range(1) & freqs <= f_range(2);
    % fplt = freqs(idx);
    % pplt = allPSD(idx);

end

meanPSD = mean(allPSD,2);

settings = struct( ...
    'peak_width_limits', [2, 10], ...   % narrower; wide peaks in 1 s = noise
    'max_n_peaks', 5, ...               % cap aggressively
    'min_peak_height', 0.15, ...        % raise floor
    'peak_threshold', 2, ...          % stricter than default 2
    'aperiodic_mode', 'fixed', ...      % knee not reliable here
    'verbose', true);

% --- run FOOOF (your wrapper expects full vectors + range) ---
fooof_results = fooof(freqs, meanPSD, f_range, settings, true);

% % --- overlay plot (log10 power is standard) ---
% figure; hold on;
% plot(fplt, log10(pplt), 'k', 'LineWidth', 1.5);                     % raw PSD
% plot(fooof_results.freqs, log10(fooof_results.ap_fit), 'b', 'LineWidth', 1.5);   % aperiodic fit
% plot(fooof_results.freqs, log10(fooof_results.fooofed_spectrum), 'r', 'LineWidth', 1.5); % full model
% plot(fooof_results.freqs, log10(fooof_results.power_spectrum - fooof_results.ap_fit), ...
%      'Color', [0.2 0.6 0.2], 'LineWidth', 1);                       % residual-ish (careful)
% 
% xlabel('Frequency (Hz)');
% ylabel('log_{10} Power');
% title('FOOOF fit overlay');
% legend({'PSD','Aperiodic fit','FOOOFed fit','(PSD - AP fit)'}, 'Location','northeast');
% xlim([f_range(1) f_range(2)]);
% grid on;

% --- print summary + peaks ---
disp('--- SETTINGS ---');
disp(settings);

disp('--- FOOOF RESULTS ---');
fprintf('Exponent: %.3f\n', fooof_results.aperiodic_params(2));
fprintf('Offset  : %.3f\n', fooof_results.aperiodic_params(1));
fprintf('R^2     : %.3f\n', fooof_results.r_squared);
fprintf('Error   : %.3f\n', fooof_results.error);

if isfield(fooof_results,'peak_params') && ~isempty(fooof_results.peak_params)
    peaks = fooof_results.peak_params;
    % columns: [CF, PW, BW] for fooof typically
    fprintf('\n--- PEAKS [CF  Amp  BW] ---\n');
    disp(peaks);
else
    disp('No peaks detected.');
end


end