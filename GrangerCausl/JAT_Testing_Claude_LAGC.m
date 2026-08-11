%% LOCATION


% X:\LossAversion\LH_Data\JAT_TrialData
% cd('X:\LossAversion\LH_Data\JAT_TrialData')

% THIS IS THE AO NAS Drive
cd('Y:\LossAversion\LH_Data\JAT_TrialData')

%%


test_GC_Extensions()   % expect PASS now

%%

path2sub = 'Y:\LossAversion\LH_Data\JAT_TrialData';
path2GC = 'Y:\LossAversion\LH_Data\GC_RESULTS';


%%

results = Run_GC_Pipeline(path2sub, path2GC, 'mode','spectral', 'maxOrder',60, 'nPerm',0);
% inspect the orders/windows it reports, then:
results = Run_GC_Pipeline(path2sub, path2GC, 'mode','spectral', ...
    'order', 5, 'winSamples', 1500, 'nPerm', 100);




%%

out = Diagnose_GC_Floor('CLASE001_GC.mat', 'dataDir', path2sub, 'nNull', 200);

%%

out = Verify_GC_Cube('CLASE001_GC.mat', path2sub, 'trials', 1:10)

%%

results = Run_GC_Pipeline(path2sub, path2GC, 'mode','spectral', 'maxOrder', 60);
Plot_Spectral_GC(fullfile(outDir,'CLASE001_GC.mat'));

%%

% subjectLIST = {'CLASE001','CLASE007','CLASE008','CLASE009','CLASE018','CLASE019','CLASE022',...
%     'CLASE023','CLASE024','CLASE026','CLASE027','CLASE029','CLASE030','CLASE031','CLASE034',...
%     'CLASE035'};
% 
% for ii = 1:length(subjectLIST)
% 
%     results = Run_GC_Pipeline(path2sub, path2GC, ...
%         'subjects', subjectLIST(ii), 'events', {'CHOICE','OUTCOME'}, 'nPerm', 20);
% 
% end

%%

% summary = Plot_GC_Summary(subjResult)

%%

% oooo = Plot_GC_Contrast2(subjResult, 'CHOICE','OUTCOME')
% oooo.omnibus.p        % did the network change at all? (FWER, max-statistic)
% oooo.node.dOut, oooo.node.pOut, oooo.node.hOut   % per-region outgoing influence change
% oooo.node.dIn,  oooo.node.pIn,  oooo.node.hIn    % per-region incoming
% % o.type.intraL / .intraR / .inter        % by connection class
% oooo.edge.pFWER                            % per-edge, FWER-corrected

%%

TR = Run_GC_TimeResolved(path2sub, 'CLASE001', ...
        'winLen',250, 'step',25, 'tZero',501, 'baseline',[1 375], ...
        'spectral',true, 'commonTrials',true, 'nPerm',0, ...
        'epochOpts',{'winSamples',1500});
%%
Plot_GC_TimeResolved(TR, 'view','band', 'events',{'CHOICE'}, ...
    'baseline',true, 'topK',6, 'bands',struct('theta',[4 8],'alpha',[8 12]),...
    'ylim', [-0.01 0.06]);

Plot_GC_TimeResolved(TR, 'view','band', 'events',{'OUTCOME'}, ...
    'baseline',true, 'topK',6, 'bands',struct('theta',[4 8],'alpha',[8 12]),...
    'ylim', [-0.01 0.06]);

%%
close all
TR = Run_GC_TimeResolved(path2sub, 'CLASE019', ...
    'winLen',250, 'step',25, 'tZero',501, 'baseline',[1 375], ...
    'spectral',true, 'commonTrials',true, 'nPerm',0, ...
    'epochOpts',{'winSamples',1500});
%%
Plot_GC_TimeResolved(TR, 'view','band', 'events',{'CHOICE'}, ...
    'baseline',true, 'topK',6, 'bands',struct('theta',[4 8],'alpha',[8 12]),...
    'ylim', [-0.01 0.06]);

Plot_GC_TimeResolved(TR, 'view','band', 'events',{'OUTCOME'}, ...
    'baseline',true, 'topK',6, 'bands',struct('theta',[4 8],'alpha',[8 12]),...
    'ylim', [-0.01 0.06]);