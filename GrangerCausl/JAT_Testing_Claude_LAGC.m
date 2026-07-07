%% LOCATION


% X:\LossAversion\LH_Data\JAT_TrialData
cd('X:\LossAversion\LH_Data\JAT_TrialData')

%%


test_load_region_epochs('CLASE001_L_AC_TrialDATA.mat')   % expect PASS now

%%
results = Run_GC_Pipeline('/path/to/Subject_Data', '/path/to/GC_Results', ...
    'subjects', {'CLASE001'}, 'events', {'CHOICE'}, 'nPerm', 20);