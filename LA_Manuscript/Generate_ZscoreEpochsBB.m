function [] = Generate_ZscoreEpochsBB()

% 1 ---------- LOOP THROUGH TABLE

% 2 ---------- FOR EACH ROW
%   i. ------- Get Hemisphere
%   ii. ------ Get Brain Area
%   iii. ----- Get Block ID
%   iv. ------ Find appropriate .mat file 

% 3 ---------- LOAD in BLOCK FILE

% 4 ---------- For each FOOOF
%   i -------- Collapse 5 second vectors into matrix for Block
%   ii ------- USE RESIDUAL

% 5 ---------- Z-score residual PSD for epoch

% 6 ---------- Re-extract z-scored power from Peaks

% 7 ---------- Add to table. 

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\Patient folders\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\Patient folders\';

    case 'DESKTOP-I5CPDO7' % JAT WORK PC
        preBlockLOC = 'Z:\LossAversion\LH_Data\JAT_BlockData';
        preZscoreLOC = 'Z:\LossAversion\LH_Data\JAT_FoooFtables\';
        postZscoreLOC = 'Z:\LossAversion\LH_Data\JAT_FoooFtablesZS';
end % switch case

cd(preBlockLOC)
bDir1 = dir('*.mat');
bDir2 = struct2table(bDir1);
bDir3 = bDir2.name;

cd(preZscoreLOC)

mDir1 = dir('*.mat');
mDir2 = struct2table(mDir1);
mDir3 = mDir2.name;

% Loop through each block file
for mmi = 1:height(mDir2)           % -------------------- 1
    mainSubFile = mDir3{mmi}; 
    cd(preZscoreLOC)
    load(mainSubFile,'allFoooftab');

    for tti = 1:height(allFoooftab) % -------------------- 2

        tmpROW = allFoooftab(tti,:);

        subID = tmpROW.SubID{1};               % Get sub ID
        brainArea = tmpROW.brainA{1};          % Get Brain Area
        blockID = tmpROW.BlockID(1);           % Get Block ID

        % CD to Block 
        cd(preBlockLOC)

        % Create Search name
        searchName = [subID,'_',brainArea,'_','blocks.mat'];

        if matches(searchName,bDir3)
            load(searchName,'preBlockFooof')

            % Get appropriate block matrix for all 3 fooofs 
            % 1. Fooof
            tmpBlockFOF = preBlockFooof.FoooF.(['B_',num2str(blockID)]);

            tmpBlMatFOOOF = zeros(5,982);
            for blmi = 1:5
                tmpBlMatFOOOF(blmi,:) = tmpBlockFOF{blmi}.power_spectrum - tmpBlockFOF{blmi}.ap_fit; 
            end

            % 2. JATspecP FIXED
            tmpBlockSPF_K = preBlockFooof.JATSpecP.(['B_',num2str(blockID)]);

            tmpBlMatSPF = zeros(5,1001);
            for blmi = 1:5
                tmpBlMatSPF(blmi,:) = tmpBlockSPF_K{blmi}.FIXED.residuals; 
            end

            % 3. JATspecP KNEE
            tmpBlMatSPK = zeros(5,1001);
            for blmi = 1:5
                tmpBlMatSPK(blmi,:) = tmpBlockSPF_K{blmi}.KNEE.residuals; 
            end

            % Z-score residuals for each block
            % 1. Fooof
            muIBI_Fooof = median(tmpBlMatFOOOF, 1, 'omitnan');
            sdIBI_Fooof = 1.4826 * mad(tmpBlMatFOOOF, 1, 1);
            sdIBI_Fooof(sdIBI_Fooof==0) = NaN;

            muIBI_Fooof = smoothdata(muIBI_Fooof, 'movmean', 7);
            sdIBI_Fooof = smoothdata(sdIBI_Fooof, 'movmean', 7);

            sdIBI_Fooof = max(sdIBI_Fooof, prctile(sdIBI_Fooof, 5));   % still floor it
            epoch_Fooof = tmpROW.FOOOFoutput{1}.power_spectrum - tmpROW.FOOOFoutput{1}.ap_fit; 
            zEpoch_Fooof = (epoch_Fooof - muIBI_Fooof) ./sdIBI_Fooof;

            freqs = tmpROW.FOOOFoutput{1}.freqs;
            peaks = tmpROW.FOOOFoutput{1}.peak_params;
            peaksN = [peaks , nan(height(peaks),1)];

            for ppi = 1:height(peaks)
                tmpHz = peaks(ppi,1);
                [~ , freqIND] = min(abs(freqs - tmpHz));
                zPkpwer = zEpoch_Fooof(freqIND);
                peaksN(ppi,4) = zPkpwer;
            end

            % Store the z-scored power from Peaks into the output table
            allFoooftab.FOOOFoutput{tti}.zPkpwer = peaksN;
            allFoooftab.FOOOFoutput{tti}.zScoreEpoch = zEpoch_Fooof;
            allFoooftab.FOOOFoutput{tti}.BlkMed = muIBI_Fooof;
            allFoooftab.FOOOFoutput{tti}.BlkMad = sdIBI_Fooof;

            % 2. JATspecP FIXED
            muIBI_spFIXED = median(tmpBlMatSPF, 1, 'omitnan');
            sdIBI_spFIXED = 1.4826 * mad(tmpBlMatSPF, 1, 1);
            sdIBI_spFIXED(sdIBI_spFIXED==0) = NaN;

            muIBI_spFIXED = smoothdata(muIBI_spFIXED, 'movmean', 9);
            sdIBI_spFIXED = smoothdata(sdIBI_spFIXED, 'movmean', 9);

            sdIBI_spFIXED = max(sdIBI_spFIXED, prctile(sdIBI_spFIXED, 5));   % still floor it
            epoch_spFIXED = transpose(tmpROW.SpecPARAMJ{1}.FIXED.residuals); 
            zEpoch_spFIXED = (epoch_spFIXED - muIBI_spFIXED) ./sdIBI_spFIXED;

            freqsF = tmpROW.SpecPARAMJ{1}.General.Freqs;
            peaksF = tmpROW.SpecPARAMJ{1}.FIXED.peaks.Summary;
            peaksNF = peaksF;
            peaksNF.AmpZS = nan(height(peaksNF),1);

            for ppi = 1:height(peaksNF)
                tmpHz = peaksNF.Hz(ppi);
                [~ , freqIND] = min(abs(freqsF - tmpHz));
                zPkpwer = zEpoch_spFIXED(freqIND);
                peaksNF.AmpZS(ppi) = zPkpwer;
            end

            % Store the z-scored power from Peaks into the output table
            allFoooftab.SpecPARAMJ{1}.FIXED.zPkpwer = peaksNF;
            allFoooftab.SpecPARAMJ{1}.FIXED.zScoreEpoch = zEpoch_spFIXED;
            allFoooftab.SpecPARAMJ{1}.FIXED.BlkMed = muIBI_spFIXED;
            allFoooftab.SpecPARAMJ{1}.FIXED.BlkMad = sdIBI_spFIXED;

            % 3. JATspecP KNEE
            muIBI_spKNEE = median(tmpBlMatSPK, 1, 'omitnan');
            sdIBI_spKNEE = 1.4826 * mad(tmpBlMatSPK, 1, 1);
            sdIBI_spKNEE(sdIBI_spKNEE==0) = NaN;

            muIBI_spKNEE = smoothdata(muIBI_spKNEE, 'movmean', 9);
            sdIBI_spKNEE = smoothdata(sdIBI_spKNEE, 'movmean', 9);

            sdIBI_spKNEE = max(sdIBI_spKNEE, prctile(sdIBI_spKNEE, 5));   % still floor it
            epoch_spKNEE = transpose(tmpROW.SpecPARAMJ{1}.KNEE.residuals);
            zEpoch_spKNEE = (epoch_spKNEE - muIBI_spKNEE) ./sdIBI_spKNEE;

            freqsK = tmpROW.SpecPARAMJ{1}.General.Freqs;
            peaksK = tmpROW.SpecPARAMJ{1}.KNEE.peaks.Summary;
            peaksNK = peaksK;
            peaksNK.AmpZS = nan(height(peaksNK),1);

            for ppi = 1:height(peaksNK)
                tmpHz = peaksNK.Hz(ppi);
                [~ , freqIND] = min(abs(freqsK - tmpHz));
                zPkpwer = zEpoch_spKNEE(freqIND);
                peaksNK.AmpZS(ppi) = zPkpwer;
            end

            % Store the z-scored power from Peaks into the output table
            allFoooftab.SpecPARAMJ{1}.KNEE.zPkpwer = peaksNK;
            allFoooftab.SpecPARAMJ{1}.KNEE.zScoreEpoch = zEpoch_spKNEE;
            allFoooftab.SpecPARAMJ{1}.KNEE.BlkMed = muIBI_spKNEE;
            allFoooftab.SpecPARAMJ{1}.KNEE.BlkMad = sdIBI_spKNEE;

            % Z_trial = (R_trial - muIBI) ./ sdIBI;
            % Z_trial = (R_trial - muIBI) ./ sdIBI;     % [nF x nT]
        else
            test = 1;

        end
    end % LOOP THROUGH ROWS
    cd(postZscoreLOC)
    namePARTS = split(mainSubFile,'.');
    saveNAME = [namePARTS{1},'_ZS.',namePARTS{2}];
    save(saveNAME,'allFoooftab');
    disp(saveNAME)
end % LOOP THROUGH SUBJECTS


% Inputs:
%   freqs   [nF x 1] frequency vector (3–55 Hz, same for baseline + trials)
%   R_ibi   [nF x nB] residual spectra for baseline (IBI 1s windows)
%   R_trial [nF x nT] residual spectra for trials (1s windows)
%
% Residual definition assumed:
%   R = log10(PSD) - log10(aperiodic_fit_knee)

% ===== #3: baseline mean & SD per frequency, then z-score trials =====
% muIBI  = mean(R_ibi,  2, 'omitnan');      % [nF x 1]
% sdIBI  = std (R_ibi,  0, 'omitnan');      % [nF x 1]
% sdIBI(sdIBI==0) = NaN;
% 
% muIBI = median(R_ibi, 2, 'omitnan');
% sdIBI = 1.4826 * mad(R_ibi, 1, 2);
% sdIBI(sdIBI==0) = NaN;

end