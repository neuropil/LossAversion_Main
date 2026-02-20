

partID = 'CLASE009';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';
    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';
    case 'JAT'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';
end % switch case

% CD to participant folder for baseline
partPathBaseline = strcat(partCD, partID, '\', 'Baseline');
cd(partPathBaseline)

% Get directory of names you want for baseline
mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, '_baseline') & ~contains(dirNames, 'FOOOF_baseline');
dirEphysNames = string(dirNames(dirFilter));

for i = 1:length(dirEphysNames)

    % --- Load all files --- %

    % Load temp baseline file
    cd(partPathBaseline)
    tempBaselineName = dirEphysNames(i);
    load(tempBaselineName) % loads as baseline

    % CD to participant folder for trial FOOOF
    partPath = strcat(partCD, partID);
    cd(partPath)

    % Get directory names you want for FOOOF
    fdir = dir;
    fdirNames = {fdir.name};
    fdirFilter = contains(fdirNames, '.mat');
    fdirEphysNames = string(fdirNames(fdirFilter));

    % Extract subject, hemisphere, and area
    parts = split(tempBaselineName, '_');   % ["CLASE001","L","AH","baseline.mat"]
    subject = parts(1);
    hemi    = parts(2);
    area    = parts(3);

    % Build the prefix we want to match
    prefix = subject + "_" + hemi + "_" + area + "_";

    % Filter fdirEphysNames for the same hemisphere+area
    matchMask = startsWith(fdirEphysNames, prefix);
    matchedFiles = fdirEphysNames(matchMask);

    % Load FOOOF name that matches baseline name
    load(matchedFiles); % loads as FOOOF_tab

    % Create New Table
    ZscoreTab = FOOOF_tab(:, 1:5);
    ZscoreTab.Method1 = cell(height(ZscoreTab),1);
    ZscoreTab.Method2 = cell(height(ZscoreTab),1);
    ZscoreTab.Method3 = cell(height(ZscoreTab),1);

    % --- Combine all baseline spectrums --- %
    baselineSpec = [baseline.Break1.FOOOFOutput.fooofed_spectrum;...
        baseline.Break2.FOOOFOutput.fooofed_spectrum; ...
        baseline.Break3.FOOOFOutput.fooofed_spectrum; ...
        baseline.Break4.FOOOFOutput.fooofed_spectrum];

    freqs = baseline.Break1.FOOOFOutput.freqs; % Frequency
    freqsRound = round(freqs); % Rounded frequency

    % Baseline mean and SD for the whole spectrum
    baselineMean = mean(mean(baselineSpec));
    baselineSD = std(std(baselineSpec));

    % Baseline mean and SD for each frequency band
    thetaIDX = (freqs >= 1 & freqs <=8);
    alphaIDX = (freqs >= 8 & freqs <=12);
    betaIDX = (freqs >= 12 & freqs <=30);
    gammaIDX = (freqs >= 30);

    thetaBaselineMean = mean(mean(baselineSpec(:,thetaIDX)));
    alphaBaselineMean = mean(mean(baselineSpec(:,alphaIDX)));
    betaBaselineMean = mean(mean(baselineSpec(:,betaIDX)));
    gammaBaselineMean = mean(mean(baselineSpec(:,gammaIDX)));

    thetaBaselineSD = std(std(baselineSpec(:,thetaIDX)));
    alphaBaselineSD = std(std(baselineSpec(:,alphaIDX)));
    betaBaselineSD = std(std(baselineSpec(:,betaIDX)));
    gammaBaselineSD = std(std(baselineSpec(:,gammaIDX)));


    for j = 1:height(FOOOF_tab)

        if isempty(FOOOF_tab.FOOOFoutput{j,1})
            continue
        else

            tempFreqAll = FOOOF_tab.FOOOFoutput{j,1}.peak_params(:,1);
            tempPowerAll = FOOOF_tab.FOOOFoutput{j,1}.peak_params(:,2);

            % Empty cells
            method1Power = [];
            method2Power = [];
            method3Power = [];

            for jj = 1:height(tempFreqAll)

                tempFreq = tempFreqAll(jj);
                tempFreqRound = round(tempFreq);
                tempPower = tempPowerAll(jj);

                % --- Z-Score Method 1 --- %
                % Take the whole fooofed spectrum from the baseline and use that mean
                % and SD for the zscoring

                method_1 = (tempPower - baselineMean)/baselineSD;

                method1Power{jj,1} = method_1; % Add zscore to cell

                % --- Z-Score Method 2 --- %
                % Zscore power data based on what frequency band it was in

                if tempFreq <= 8 % theta
                    method_2 = (tempPower - thetaBaselineMean)/thetaBaselineSD;
                elseif  tempFreq >= 8 & tempFreq <=12 % alpha
                    method_2 = (tempPower - alphaBaselineMean)/alphaBaselineSD;
                elseif tempFreq >= 12 & tempFreq <=30 % beta
                    method_2 = (tempPower - betaBaselineMean)/betaBaselineSD;
                else tempFreq >= 30 % gamma
                    method_2 = (tempPower - gammaBaselineMean)/gammaBaselineSD;
                end

                method2Power{jj,1} = method_2; % Add zscore to cell

                % --- Z-Score Method 3 --- %
                % Zscore power data with the specific frequency it's releated to

                freqIDX = find(freqsRound == tempFreqRound); % Index into freqsRound to find what column has the same number as tempFreqRound

                tempBaselineSpec = baselineSpec(:,freqIDX);
                meanTempBaselineSpec = mean(tempBaselineSpec);
                SDTempBaselineSpec = std(tempBaselineSpec);

                method_3 = (tempPower - meanTempBaselineSpec)/SDTempBaselineSpec;

                method3Power{jj,1} = method_3; % Add zscore to cell

            end % for / jj

            % Add zscored data to table
            ZscoreTab.Method1{j} = method1Power;
            ZscoreTab.Method2{j} = method2Power;
            ZscoreTab.Method3{j} = method3Power;

        end % if else

    end % for / j

% Save 
cd('Y:\LossAversion\LH_Data\FOOOF_data\CLASE009\PeriodicZscore')

saveName = strcat(prefix, 'Zscore.mat');
save(saveName, "ZscoreTab");

end % for / i

