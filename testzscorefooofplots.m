close 
subplot(3,1,1)
plot(epoch_Fooof)
title('Trial Residual')
subplot(3,1,2)
plot(muIBI_Fooof)
title('MAD IBI Residual')
subplot(3,1,3)
plot(zEpoch_Fooof)
title('Z-scored Trial Residual')

% hold on
% zEpoch_Fooof_s = smoothdata(zEpoch_Fooof,1,'movmean',5);
% z_max = max(zEpoch_Fooof_s,[],1);
% plot(z_max)

freqs = tmpROW.FOOOFoutput{1}.freqs;
peaks = tmpROW.FOOOFoutput{1}.peak_params;
peaksN = [peaks , nan(height(peaks),1)];

for ppi = 1:height(peaks)

    tmpHz = peaks(ppi,1);
    [~ , freqIND] = min(abs(freqs - tmpHz));
    zPkpwer = zEpoch_Fooof(freqIND); 
    peaksN(ppi,4) = zPkpwer;

end

%%

close 
subplot(3,1,1)
plot(epoch_spKNEE)
title('Trial Residual')
subplot(3,1,2)
plot(muIBI_spKNEE)
title('MAD IBI Residual')
subplot(3,1,3)
plot(zEpoch_spKNEE)
title('Z-scored Trial Residual')


%% 

preZscoreLOC = 'Z:\LossAversion\LH_Data\JAT_FoooFtables\';
cd(preZscoreLOC)
mDir1 = dir('*.mat');
mDir2 = struct2table(mDir1);
mDir3 = mDir2.name;

for mmmii = 2:length(mDir3)
    subjectID = mDir3{mmmii}(1:8); % Extract subject ID from filename
    FOOOF_TABLE_GenerateJAT2(subjectID)
end