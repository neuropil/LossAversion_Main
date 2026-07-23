%% LOCATION


% X:\LossAversion\LH_Data\JAT_TrialData
% cd('X:\LossAversion\LH_Data\JAT_TrialData')

% THIS IS THE AO NAS Drive
cd('Y:\LossAversion\LH_Data\JAT_TrialData')

%%


test_load_region_epochs('CLASE001_L_AH_TrialDATA.mat')   % expect PASS now

%%

path2sub = 'Y:\LossAversion\LH_Data\JAT_TrialData';
path2GC = 'Y:\LossAversion\LH_Data\GC_RESULTS';



%%

subjectLIST = {'CLASE001','CLASE007','CLASE008','CLASE009','CLASE018','CLASE019','CLASE022',...
    'CLASE023','CLASE024','CLASE026','CLASE027','CLASE029','CLASE030','CLASE031','CLASE034',...
    'CLASE035'};

for ii = 1:length(subjectLIST)

    results = Run_GC_Pipeline(path2sub, path2GC, ...
        'subjects', subjectLIST(ii), 'events', {'CHOICE','OUTCOME'}, 'nPerm', 20);

end

%%

summary = Plot_GC_Summary(subjResult)

%%

oooo = Plot_GC_Contrast2(subjResult, 'CHOICE','OUTCOME')
oooo.omnibus.p        % did the network change at all? (FWER, max-statistic)
oooo.node.dOut, oooo.node.pOut, oooo.node.hOut   % per-region outgoing influence change
oooo.node.dIn,  oooo.node.pIn,  oooo.node.hIn    % per-region incoming
% o.type.intraL / .intraR / .inter        % by connection class
oooo.edge.pFWER                            % per-edge, FWER-corrected