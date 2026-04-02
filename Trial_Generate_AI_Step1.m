function [trialTableOUT] = Trial_Generate_AI_Step1(tempPtID , Hemi , BrainArea , contactNUMs)

PCname = getenv('COMPUTERNAME');

switch PCname
    case 'DESKTOP-I5CPDO7'

        matNWB_25 = {'CLASE001','CLASE007','CLASE008','CLASE009','CLASE018','CLASE019','CLASE022',...
            'CLASE023','CLASE024','CLASE026'};

        matNWB_27 = {'CLASE027','CLASE029','CLASE030','CLASE031','CLASE034','CLASE035'};

        if matches(tempPtID,matNWB_25)
            nwbMatCD = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.5.0.0';
            remNWB = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.7.0\matnwb-2.7.0';
        elseif matches(tempPtID,matNWB_27)
            nwbMatCD = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.7.0\matnwb-2.7.0';
            remNWB = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.5.0.0';
        end
        rmpath(genpath(remNWB));
        addpath(genpath(nwbMatCD));

        synologyCD = 'Z:\LossAversion\Patient folders'; % Synology path
        NLXEventCD = 'C:\Users\Admin\Documents\Github\NLX-Event-Viewer\NLX_IO_Code'; % NLX event reader path
        addpath(NLXEventCD); % add NLX files to path

        epoch_dir = 'Z:\LossAversion\LH_Data\JAT_ADAtest\epochs';
end

%% Inputs 
% tempPtID = 'CLASE035'; % Patient ID
% Hemi = 'L';
% BrainArea = 'OF';

%%
% [cleanVolts] = CleanEphys_v2('CLASE035', 'LOF', 6, 1:7);

% ---- Load In Files ---- %

% NWB % 
% Load in NWB File 
tmpPTpath = strcat(synologyCD,'\', tempPtID,'\'); % Patient path on synology

paths = [];
paths.NWBdata = [strcat(tmpPTpath, 'NWB-processing','\', 'NWB_Data', '\')]; % Path to PtID NWB data
cd(paths.NWBdata) % CD to pt NWB data

nwbdir = dir; % NWB Dir
nwbdirNames = {nwbdir.name}; % Names of files in NWB folder
nwbdirFilter = contains(nwbdirNames, 'filter'); % Find the NWB files that have 'Filter' in the name
tempLAname = string(nwbdirNames(nwbdirFilter)); % String the NWB file that has 'filter' in it

tmp_LA = nwbRead(tempLAname); % Load the filter NWB file

% Voltage data for all macrowires and their channels 
ma_data = tmp_LA.processing.get('ecephys').nwbdatainterface.get...
    ('LFP').electricalseries.get('MacroWireSeries').data.load;

[cleanVolts] = CleanEphys_V3(tmp_LA, ma_data, BrainArea, 6, contactNUMs);

% Load in timestamps
ma_timestamps = tmp_LA.processing.get('ecephys').nwbdatainterface.get...
    ('LFP').electricalseries.get('MacroWireSeries').timestamps.load;

% Downsample timedata
ma_timestampsDS = downsample(ma_timestamps, 8); % this is downstampled by a factor of 8

% Load in Behavorial Time Table % 
paths.BehavTab = [strcat(tmpPTpath, 'Behavioral-data\EventBehavior\')]; % path to behavioral table 
cd(paths.BehavTab)

behavDir = dir;
behavDirNames = {behavDir.name};
behavDirMat = contains(behavDirNames, 'v2');
tmpBehavName = string(behavDirNames(behavDirMat));
load(tmpBehavName,'eventTABLE')           % Loads as eventTABLE 

% LA_Behavior data file
paths.LABehav = [strcat(tmpPTpath, '\LA_M\')];
cd(paths.LABehav)
matFILE = dir('*.mat');
matFILEu = struct2table(matFILE);
labEHAVf = matFILEu.name;

load(labEHAVf,'fullbehTAB');


%%
% Indicies of TRIAL starts (i.e., CHOICESHOW)
sanityCheck = unique(eventTABLE.TrialEvID(eventTABLE.TrialEvNum == 1));

try
    assert(isscalar(sanityCheck) && matches(sanityCheck, 'choiceShow'), ...
        'SanityCheck:UnexpectedValue', ...
        'Expected a single value of ''choiceShow'' but got: %s', ...
        strjoin(sanityCheck', ', '))
catch ME
    error(ME.identifier, 'Error in %s (line %d): %s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message)
end

% Find where each block starts (TrialEvNum == 1)
trialStarts = find(eventTABLE.TrialEvNum == 1);
nBlocks     = numel(trialStarts);

% Preallocate
TrialID = zeros(height(eventTABLE), 1);

for i = 1:nBlocks
    if i < nBlocks
        idx = trialStarts(i):trialStarts(i+1)-1;
    else
        idx = trialStarts(i):numel(eventTABLE.TrialEvNum);
    end
    TrialID(idx) = i;
end

eventTABLE.TrialID = TrialID;

trialSTARTS_inds = find(eventTABLE.TrialEvNum == 1);

fs = 500;
% trialTABLE = table;
epoch_types  = {'gamble',       'cue',        'outcome'};
event_prefixes = {'choiceShow', 'respWindowS',  'outDispS'};   % prefix-matched against TrialEvID
% ---- Epoch windows [start, end] in seconds relative to event onset ----
win.gamble  = [-0.2,  1.5];   % choiceShow
win.cue     = [-0.2,  1.5];   % respWindS (response window start)
win.outcome = [-0.2,  1.5];   % outDispS
windows      = {win.gamble, win.cue, win.outcome};

onset_us = struct();
for e = 1:numel(epoch_types)
    pfx  = event_prefixes{e};
    % Match rows whose TrialEvID starts with the prefix
    mask = strncmp(eventTABLE.TrialEvID, pfx, numel(pfx));
    sub_tbl = eventTABLE(mask, :);
    % Sort by TrialID so order matches behavioral table
    sub_tbl = sortrows(sub_tbl, 'TrialID');
    onset_us.(epoch_types{e}) = sub_tbl.TrialEvTm;   % [n_trials x 1] µs Unix
    fprintf('Epoch %-8s : %d trials found\n', epoch_types{e}, height(sub_tbl));
end
n_bipolar = height(cleanVolts);
neural2use = cleanVolts;
% ---- Epoch loop ----
for e = 1:numel(epoch_types)
    etype    = epoch_types{e};
    win_sec  = windows{e};
    onsets   = onset_us.(etype);   % [n_trials x 1] µs

    n_trials  = numel(onsets);
    win_samp  = round(win_sec * fs);            % samples rel. to onset
    n_samp    = diff(win_samp) + 1;
    t_axis    = linspace(win_sec(1), win_sec(2), n_samp);

    % Preallocate  [n_bipolar x n_time x n_trials]
    epoch_data = nan(n_bipolar, n_samp, n_trials);


    % CHECK THAT IT IS EPOCH and NOT TRIAL START
    for tr = 1:n_trials
        % Convert µs onset → sample index in neural recording
        % onset_sec   = (onsets(tr) - rec_start_us) / 1e6;
        % onset_samp  = round(onset_sec * fs) + 1;  % 1-based

        [~ , tmpTrialONSET] = min(abs(ma_timestampsDS - onsets(tr)));

        idx_start = tmpTrialONSET + win_samp(1);
        idx_end   = tmpTrialONSET + win_samp(2);

        % Boundary check
        % if idx_start < 1 || idx_end > size(neural, 2)
        %     fprintf('    [%s %s %s] Trial %d out of bounds — skipping\n', ...
        %         subCL, rgn, etype, tr);
        %     continue
        % end

        epoch_data(:, :, tr) = neural2use(1:n_bipolar, idx_start:idx_end);
    end
    % MAKE BASELINE 250ms BEFORE CHOICESHOW for EACH TRIAL EPOCH
    % ---- Baseline correct (z-score to pre-stimulus window) ----
    baseline_idx = t_axis >= -0.5 & t_axis < 0;
    if sum(baseline_idx) > 1
        epoch_data = baseline_correct(epoch_data, baseline_idx);
    else
        warning('No baseline samples found for %s %s %s — skipping BC', ...
            subCL, rgn, etype);
    end

    % ---- Assemble output struct ----
    out.data           = epoch_data;     % [n_bipolar x n_time x n_trials]
    out.t_axis         = t_axis;
    out.fs             = fs;
    out.subject_CL     = tempPtID;
    out.region         = BrainArea;
    out.hemisphere     = BrainArea(1);      % {n_bipolar x 1} cell
    out.contact_labels = contactNUMs; % {n_bipolar x 1} cell  e.g. 'LAMY'
    out.epoch          = etype;
    out.n_trials       = n_trials;
    out.n_bipolar      = n_bipolar;
    out.beh            = fullbehTAB;

    out_name = fullfile(epoch_dir, ...
        sprintf('%s_%s_%s.mat', tempPtID, BrainArea, etype));
    save(out_name, 'out', '-v7.3');
    fprintf('  Saved: %s\n', out_name);
end  % epoch loop


















































allTrialTables = cell(numel(trialSTARTS_inds),1);
allTrialChanns = cell(numel(trialSTARTS_inds),1);
allTrialTimes = cell(numel(trialSTARTS_inds),1);

for ii = 1:length(trialSTARTS_inds) % FOR EACH TRIAL

    % Index of start time for choice show
    tmpTStartInd = trialSTARTS_inds(ii);

    % Raw time of choice show
    tmpStartTS = eventTABLE.TrialEvTm(tmpTStartInd);
    % FIND this TS in ma_timestamps as INDEX
    [~ , choiceSHOWind] = min(abs(ma_timestampsDS - tmpStartTS));

    % Trial start will encompass ITI - 1000 ms
    trialStartI_NLX = choiceSHOWind - 500;

    tmpTrialTab = eventTABLE(eventTABLE.TrialID == ii,:);
    % tmpTrialTabRS = tmpTrialTab(~matches(tmpTrialTab.TrialEvID,'choiceShow'),:);

    tmpTrialTabRSC = removevars(tmpTrialTab,["Trials","OffsetSecs"]);

    for ttid = 1:height(tmpTrialTabRSC)

        [~, tmpEpochTS_NLX] = min(abs(ma_timestampsDS - tmpTrialTabRSC.TrialEvTm(ttid)));
        tmpTrialTabRSC.AbsolINDEX(ttid) = tmpEpochTS_NLX;
        tmpTrialTabRSC.RelINDEX(ttid) = (tmpEpochTS_NLX - trialStartI_NLX) + 1;

    end


    if ii ~= length(trialSTARTS_inds)

        % End of trial - one sample prior to start of next trial
        endOfTrialrow = height(tmpTrialTabRSC) + 1;

        tmpTrialTabRSC.Blocks(endOfTrialrow) = tmpTrialTabRSC.Blocks(1);
        tmpTrialTabRSC.TrialEvNum(endOfTrialrow) = 6;
        tmpTrialTabRSC.TrialEvID{endOfTrialrow} = 'endOfTrial';
        tmpTrialTabRSC.TrialID(endOfTrialrow) = tmpTrialTabRSC.TrialID(1);

        % Get next trial start with 1 sample offset
        upcomingTrialStartTS = trialSTARTS_inds(ii + 1);

        upcomingTmpStartTS = eventTABLE.TrialEvTm(upcomingTrialStartTS);

        [~, upcomingTrialStartI_NLX] = min(abs(ma_timestampsDS - upcomingTmpStartTS));

        endOfcurrentTRIAL_TS = upcomingTrialStartI_NLX - 1;

        tmpTrialTabRSC.TrialEvTm(endOfTrialrow) = ma_timestampsDS(endOfcurrentTRIAL_TS);
        tmpTrialTabRSC.AbsolINDEX(endOfTrialrow) = endOfcurrentTRIAL_TS;
        tmpTrialTabRSC.RelINDEX(endOfTrialrow) = (endOfcurrentTRIAL_TS - trialStartI_NLX) + 1;

        % IF LAST trial of BLOCK
        if tmpTrialTabRSC.RelINDEX(endOfTrialrow) > 5000

            tmpLASTindex = tmpTrialTabRSC.AbsolINDEX(endOfTrialrow - 1);
            tmpLASTindexShift = tmpLASTindex + 125;
            tmpTrialTabRSC.TrialEvTm(endOfTrialrow) = ma_timestampsDS(tmpLASTindexShift);
            tmpTrialTabRSC.AbsolINDEX(endOfTrialrow) = tmpLASTindexShift;
            tmpTrialTabRSC.RelINDEX(endOfTrialrow) = (tmpLASTindexShift - trialStartI_NLX) + 1;

        end

    else % LAST TRIAL
        % End of trial - LAST event + 250 ms
        endOfTrialrow = height(tmpTrialTabRSC) + 1;

        tmpTrialTabRSC.Blocks(endOfTrialrow) = tmpTrialTabRSC.Blocks(1);
        tmpTrialTabRSC.TrialEvNum(endOfTrialrow) = 6;
        tmpTrialTabRSC.TrialEvID{endOfTrialrow} = 'endOfTrial';
        tmpTrialTabRSC.TrialID(endOfTrialrow) = tmpTrialTabRSC.TrialID(1);

        tmpLASTindex = tmpTrialTabRSC.AbsolINDEX(endOfTrialrow - 1);

        tmpLASTindexShift = tmpLASTindex + 125; % 125 samples = 250 ms

        tmpTrialTabRSC.TrialEvTm(endOfTrialrow) = ma_timestampsDS(tmpLASTindexShift);
        tmpTrialTabRSC.AbsolINDEX(endOfTrialrow) = tmpLASTindexShift;
        tmpTrialTabRSC.RelINDEX(endOfTrialrow) = (tmpLASTindexShift - trialStartI_NLX) + 1;

    end

    allTrialTables{ii} = tmpTrialTabRSC;

    % Get trial voltage from CHANNELS
    allTrialChanns{ii} = cleanVolts(:, trialStartI_NLX:...
        tmpTrialTabRSC.AbsolINDEX(end));
    allTrialTimes{ii} = ma_timestampsDS(trialStartI_NLX:...
        tmpTrialTabRSC.AbsolINDEX(end));

end

trialTableOUT = table(allTrialTables,allTrialChanns,allTrialTimes,'VariableNames',...
    {'TrialTables','TrialChanVolts','TrialTimes'});

% Save

% saveLoc = [strcat('Y:\LossAversion\LH_Data\FOOOF_data\', tempPtID, '\', 'Baseline')];
% cd(saveLoc)
% 
% saveName = [strcat(tempPtID, '_',  Hemi, '_', BrainArea, '_', 'baseline.mat')];
% 
% save(saveName, "baseline");

end



function data_bc = baseline_correct(data, bl_idx)
    % data   : [n_chan x n_time x n_trials]
    % bl_idx : logical index into time dimension (pre-stimulus samples)
    bl_mean = mean(data(:, bl_idx, :), 2);   % [n_chan x 1 x n_trials]
    bl_std  = std( data(:, bl_idx, :), 0, 2);
    bl_std(bl_std == 0) = 1;                 % guard against flat baseline
    data_bc = (data - bl_mean) ./ bl_std;
end