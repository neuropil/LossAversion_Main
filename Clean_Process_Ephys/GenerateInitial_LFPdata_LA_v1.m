function [] = GenerateInitial_LFPdata_LA_v1(subjectID , hemi2use, ...
    thresh2use , contacts2use , nlxBrainArea , saveBrainArea)

% Inputs
% subjectID = 'CLASE035'; % Patient ID
% Hemi = 'L';
% BrainArea = 'LOF';

% Load in ephys that is artifact rejected and bipolar refrenced

% [cleanVolts] = CleanEphys_v2('CLASE035', 'LOF', 6, 2:7);

% [cleanVolts] = CleanEphys_v2('PartID', 'ShortBrainAreaName', Std, num of contacts);
% Num of contacts  = 1:3 or 1

PCname = getenv('COMPUTERNAME');

switch PCname
    case 'DLPFC' % laptop
        nwbMatCD = 'C:\Users\Lisa\Documents\MATLAB'; % NWB_read path
        synologyCD = 'Z:\LossAversion\Patient folders'; % Synology path
        NLXEventCD = 'E:\GitKraken\NLX-Event-Viewer\NLX_IO_Code'; % NLX event reader path
        addpath('E:\GitKraken\NLX-Event-Viewer\NLX_IO_Code'); % add NLX files to path
        addpath(genpath([strcat(nwbMatCD,'\matnwb-2.5.0.0')])); % add nwb_read and subfolders to path
    case 'LATERALHABENULA' % lab computer
        nwbMatCD = 'C:\Users\Lisa\Documents\MATLAB';
        synologyCD = 'Y:\LossAversion\Patient folders'; % Synology path
        NLXEventCD = 'Z:\GitKraken\NLX-Event-Viewer\NLX_IO_Code'; % NLX event reader path
        addpath('Z:\GitKraken\NLX-Event-Viewer\NLX_IO_Code'); % add NLX files to path
        addpath(genpath([strcat(nwbMatCD,'\matnwb-2.5.0.0')])); % add nwb_read and subfolders to path
    case 'DESKTOP-I5CPDO7'

        matNWB_25 = {'CLASE001','CLASE007','CLASE008','CLASE009','CLASE018','CLASE019','CLASE022',...
            'CLASE023','CLASE024','CLASE026'};

        matNWB_27 = {'CLASE027','CLASE029','CLASE030','CLASE031','CLASE034','CLASE035'};

        if matches(subjectID,matNWB_25)
            nwbMatCD = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.5.0.0';
            remNWB = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.7.0\matnwb-2.7.0';
        elseif matches(subjectID,matNWB_27)
            nwbMatCD = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.7.0\matnwb-2.7.0';
            remNWB = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.5.0.0';
        end
        rmpath(genpath(remNWB));
        addpath(genpath(nwbMatCD));


        synologyCD = 'Z:\LossAversion\Patient folders'; % Synology path
        NLXEventCD = 'C:\Users\Admin\Documents\Github\NLX-Event-Viewer\NLX_IO_Code'; % NLX event reader path
        addpath(NLXEventCD); % add NLX files to path

end

% ---- Load In Files ---- %
% NWB %
% Load in NWB File
tmpPTpath = strcat(synologyCD,'\', subjectID,'\'); % Patient path on synology

paths = [];
paths.NWBdata = [strcat(tmpPTpath, 'NWB-processing','\', 'NWB_Data', '\')]; % Path to PtID NWB data
cd(paths.NWBdata) % CD to pt NWB data

nwbdir = dir; % NWB Dir
nwbdirNames = {nwbdir.name}; % Names of files in NWB folder
nwbdirFilter = contains(nwbdirNames, 'filter'); % Find the NWB files that have 'Filter' in the name
tempLAname = string(nwbdirNames(nwbdirFilter)); % String the NWB file that has 'filter' in it

tmp_LA = nwbRead(tempLAname); % Load the filter NWB file

% Load in timestamps
ma_timestamps = tmp_LA.processing.get('ecephys').nwbdatainterface.get...
    ('LFP').electricalseries.get('MacroWireSeries').timestamps.load;

% Downsample timedata
ma_timestampsDS = downsample(ma_timestamps, 8); % this is downstampled by a factor of 8

% LFP data
ma_data = tmp_LA.processing.get('ecephys').nwbdatainterface.get...
    ('LFP').electricalseries.get('MacroWireSeries').data.load;

% CLEAN
[cleanVolts] = CleanEphys_V3(tmp_LA, ma_data, nlxBrainArea, thresh2use, contacts2use);


% Load in Behavorial Time Table %
paths.BehavTab = [strcat(tmpPTpath, 'Behavioral-data\EventBehavior\')]; % path to behavioral table
cd(paths.BehavTab)

behavDir = dir;
behavDirNames = {behavDir.name};
behavDirMat = contains(behavDirNames, 'v2');
tmpBehavName = string(behavDirNames(behavDirMat));
load(tmpBehavName , 'eventTABLE')           % Loads as eventTABLE

% Load behavioral participant data %
paths.BehFolder = [strcat(tmpPTpath,'Behavioral-data', '\')];
cd(paths.BehFolder)

behDIR = dir; % behavioral folder DIR
behDIRNames = {behDIR.name}; % Names of files in behavioral folder
behDIRmat = contains(behDIRNames, 'mat'); % Find the files that have ".mat" in the name
tempBehavName = string(behDIRNames(behDIRmat));

eventTab = load(tempBehavName);         % load behavioral time table file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% EPOCH EPHYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- EPOCHS ---- %
% Seperate ephys per epochs
% --- 1. Locate event rows in eventTABLE --------------------------------- %
% idxChoice  = find(eventTABLE.TrialEvNum == 1);   % gamble onset
% idxRespShow   = find(eventTABLE.TrialEvNum == 2);   % response window open
% idxRespExit = find(eventTABLE.TrialEvNum == 3);   % response window close
% idxOutShow    = find(eventTABLE.TrialEvNum == 4);   % outcome onset
% idxOutExit  = find(eventTABLE.TrialEvNum == 5);   % outcome offset
% 
% nTrials    = length(idxChoice);
% epochEphys = cell(nTrials * 3, 1);
% epochRepeat = repmat({'choice'; 'response'; 'outcome'}, nTrials, 1);
% 
% for ci = 1:nTrials
%     choiceTime   = eventTABLE.TrialEvTm(idxChoice(ci));
%     respShowTime = eventTABLE.TrialEvTm(idxRespShow(ci));
%     respExitTime = eventTABLE.TrialEvTm(idxRespExit(ci));
%     outShowTime  = eventTABLE.TrialEvTm(idxOutShow(ci));
%     outExitTime  = eventTABLE.TrialEvTm(idxOutExit(ci));
% 
%     [~, choiceLoc]   = min(abs(ma_timestampsDS - choiceTime));
%     [~, respShowLoc] = min(abs(ma_timestampsDS - respShowTime));
%     [~, respExitLoc] = min(abs(ma_timestampsDS - respExitTime));
%     [~, outShowLoc]  = min(abs(ma_timestampsDS - outShowTime));
%     [~, outExitLoc]  = min(abs(ma_timestampsDS - outExitTime));
% 
%     row = (ci - 1) * 3 + 1;
%     epochEphys{row,   1} = cleanVolts(:, choiceLoc   : respShowLoc - 1);
%     epochEphys{row+1, 1} = cleanVolts(:, respShowLoc : respExitLoc - 1);
%     epochEphys{row+2, 1} = cleanVolts(:, outShowLoc  : outExitLoc  - 1);
% end
% 
% Epoch_EphysTAB = cell2table([epochRepeat, epochEphys], ...
%     "VariableNames", ["EpochID", "Ephys"]);

% --- 1. Locate event rows in eventTABLE --------------------------------- %
idxChoiceShow  = find(eventTABLE.TrialEvNum == 1);   % gamble onset
idxRespWindS   = find(eventTABLE.TrialEvNum == 2);   % response window open
idxRespWindEnd = find(eventTABLE.TrialEvNum == 3);   % response window close
idxOutDispS    = find(eventTABLE.TrialEvNum == 4);   % outcome onset
idxOutDispEnd  = find(eventTABLE.TrialEvNum == 5);   % outcome offset

% Sanity check: all event types must have the same trial count
nTrials = length(idxChoiceShow);
assert(length(idxRespWindS)   == nTrials, 'Mismatch: respWindS count');
assert(length(idxRespWindEnd) == nTrials, 'Mismatch: respWindEnd count');
assert(length(idxOutDispS)    == nTrials, 'Mismatch: outDispS count');
assert(length(idxOutDispEnd)  == nTrials, 'Mismatch: outDispEnd count');

fprintf('[Epoch] Found %d trials across 5 event types.\n', nTrials);

% --- 2. Pre-allocate outputs -------------------------------------------- %
epochChoice   = cell(nTrials, 1);   % {nCh × nSamples} per trial
epochResponse = cell(nTrials, 1);
epochOutcome  = cell(nTrials, 1);

epochLabels   = repmat({'choice'; 'response'; 'outcome'}, nTrials, 1);
trialNums     = reshape(repmat(1:nTrials, 3, 1), [], 1);   % [1 1 1 2 2 2 ...]
nSamplesLog   = zeros(nTrials * 3, 1);   % sample count per row (diagnostic)

misalignLog = false(nTrials, 3);   % cols: [choice, response, outcome]

% --- 3. Extract epochs -------------------------------------------------- %
for ci = 1:nTrials

    % --- 3a. Retrieve event timestamps (from eventTABLE) ---------------- %
    t_choiceShow  = eventTABLE.TrialEvTm(idxChoiceShow(ci));
    t_respWindS   = eventTABLE.TrialEvTm(idxRespWindS(ci));
    t_respWindEnd = eventTABLE.TrialEvTm(idxRespWindEnd(ci));
    t_outDispS    = eventTABLE.TrialEvTm(idxOutDispS(ci));
    t_outDispEnd  = eventTABLE.TrialEvTm(idxOutDispEnd(ci));

    % --- 3b. Map timestamps → nearest sample indices in ma_timestampsDS - %
    [~, loc_choiceShow]  = min(abs(ma_timestampsDS - t_choiceShow));
    [~, loc_respWindS]   = min(abs(ma_timestampsDS - t_respWindS));
    [~, loc_respWindEnd] = min(abs(ma_timestampsDS - t_respWindEnd));
    [~, loc_outDispS]    = min(abs(ma_timestampsDS - t_outDispS));
    [~, loc_outDispEnd]  = min(abs(ma_timestampsDS - t_outDispEnd));

    % --- 3c. Validate ordering (catch corrupted/misaligned trials) ------- %
    % --- 3c & 3d. Validate ordering, slice, or NaN-fill if misaligned --- %
    row = (ci - 1) * 3 + 1;

    % Choice epoch: choiceShow → respWindS
    if loc_choiceShow < loc_respWindS
        epochChoice{ci} = cleanVolts(:, loc_choiceShow : loc_respWindS - 1);
    else
        warning('[Epoch] Trial %d MISALIGNED: choiceShow >= respWindS. Filling choice epoch with NaN.', ci);
        epochChoice{ci} = NaN(size(cleanVolts, 1), 1);   % [nCh × 1] NaN placeholder
        misalignLog(ci, 1) = true;
    end

    % Response epoch: respWindS → respWindEnd
    if loc_respWindS < loc_respWindEnd
        epochResponse{ci} = cleanVolts(:, loc_respWindS : loc_respWindEnd - 1);
    else
        warning('[Epoch] Trial %d MISALIGNED: respWindS >= respWindEnd. Filling response epoch with NaN.', ci);
        epochResponse{ci} = NaN(size(cleanVolts, 1), 1);
        misalignLog(ci, 2) = true;
    end

    % Outcome epoch: outDispS → outDispEnd
    if loc_outDispS < loc_outDispEnd
        epochOutcome{ci} = cleanVolts(:, loc_outDispS : loc_outDispEnd - 1);
    else
        warning('[Epoch] Trial %d MISALIGNED: outDispS >= outDispEnd. Filling outcome epoch with NaN.', ci);
        epochOutcome{ci} = NaN(size(cleanVolts, 1), 1);
        misalignLog(ci, 3) = true;
    end

    % Log sample counts (NaN epochs log as 1 — the placeholder width)
    nSamplesLog(row)   = size(epochChoice{ci},   2);
    nSamplesLog(row+1) = size(epochResponse{ci}, 2);
    nSamplesLog(row+2) = size(epochOutcome{ci},  2);

    % --- 3d. Slice cleanVolts (non-overlapping: end = next onset - 1) --- %
    %   cleanVolts : [nChannels × nTotalSamples]
    epochChoice{ci}   = cleanVolts(:, loc_choiceShow  : loc_respWindS   - 1);
    epochResponse{ci} = cleanVolts(:, loc_respWindS   : loc_respWindEnd - 1);
    epochOutcome{ci}  = cleanVolts(:, loc_outDispS    : loc_outDispEnd  - 1);

    % Log sample counts for diagnostics
    row = (ci - 1) * 3 + 1;
    nSamplesLog(row)   = size(epochChoice{ci},   2);
    nSamplesLog(row+1) = size(epochResponse{ci}, 2);
    nSamplesLog(row+2) = size(epochOutcome{ci},  2);

end % for ci

% --- 4. Assemble output table ------------------------------------------- %
epochEphysAll = [epochChoice; epochResponse; epochOutcome];

% Re-sort so rows interleave as: trial1/choice, trial1/response,
% trial1/outcome, trial2/choice ... (matches epochLabels / trialNums order)
sortIdx       = reshape(reshape(1:nTrials*3, nTrials, 3)', [], 1);
epochEphysAll = epochEphysAll(sortIdx);

Epoch_EphysTAB = table( ...
    epochLabels, ...
    trialNums, ...
    nSamplesLog, ...
    epochEphysAll, ...
    'VariableNames', {'EpochID', 'TrialNum', 'nSamples', 'Ephys'});

% Epoch_Ephys = [epochRepeat, epochEphys];        % Ephys by epochs
% Epoch_EphysTAB = cell2table(Epoch_Ephys, "VariableNames",["EpochID" "Ephys"]);

% ---- Behavior ---- %

% Trial type info
% riskyloss < 0 = gain/loss trial :: either gain X or lose Y
% riskyloss == 0 = gain only :: either gain X or lose 0
% choice 1 = gamble, 0 = alternative

% Check trial
checkIndex = eventTab.subjdata.cs.ischecktrial;

% Loss aversion (LA) or risk aversion (RA)
% Gain/loss trials - this measures loss aversion
LA_trials = eventTab.subjdata.cs.riskyLoss < 0 & ~checkIndex;
% Gain only trials - this measures risk aversion
% RA_trials = eventTab.subjdata.cs.riskyLoss == 0 & ~checkIndex;

% Outcome results
outcomeLoss = eventTab.subjdata.cs.outcome < 0 & ~checkIndex;
outcomeGain = eventTab.subjdata.cs.outcome > 0 & ~checkIndex;
outcomeNeutral = eventTab.subjdata.cs.outcome == 0 & ~checkIndex;

% Repeat results
LA_rep = num2cell(reshape(repmat(LA_trials', 3, 1), [], 1));
outcomeLoss_rep = num2cell(reshape(repmat(outcomeLoss', 3, 1), [], 1));
outcomeGain_rep = num2cell(reshape(repmat(outcomeGain', 3, 1), [], 1));
outcomeNeutral_rep = num2cell(reshape(repmat(outcomeNeutral', 3, 1), [], 1));

% Create behavior table
tmpBehavTab = cell2table([LA_rep, outcomeLoss_rep, outcomeGain_rep, outcomeNeutral_rep],...
    "VariableNames", ["LA" "OutcomeLoss" "OutcomeGain" "OutcomeNeutral"]);

EpochEphysTab = [Epoch_EphysTAB, tmpBehavTab]; % save this variable

% Save
saveLoc = [strcat(synologyCD,'\', subjectID, '\', 'ProcessedEphys')];
cd(saveLoc)

saveName = [strcat(subjectID, '_',  hemi2use, '_', saveBrainArea, '_', 'EpochEphys.mat')];

save(saveName, "EpochEphysTab");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ALL EPHYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Index where blocks are 
idxBlock1 = find(eventTABLE.Blocks(:) == 1);
idxBlock1_First = idxBlock1(1);
idxBlock1_Last = idxBlock1(end);

idxBlock2 = find(eventTABLE.Blocks(:) == 2);
idxBlock2_First = idxBlock2(1);
idxBlock2_Last = idxBlock2(end);

idxBlock3 = find(eventTABLE.Blocks(:) == 3);
idxBlock3_First = idxBlock3(1);
idxBlock3_Last = idxBlock3(end);

idxBlock4 = find(eventTABLE.Blocks(:) == 4);
idxBlock4_First = idxBlock4(1);
idxBlock4_Last = idxBlock4(end);

idxBlock5 = find(eventTABLE.Blocks(:) == 5);
idxBlock5_First = idxBlock5(1);
idxBlock5_Last = idxBlock5(end);

% Get out time stamps for each behavioral marker
block1_FirstTime = eventTABLE.TrialEvTm(idxBlock1_First);
block1_LastTime = eventTABLE.TrialEvTm(idxBlock1_Last);

block2_FirstTime = eventTABLE.TrialEvTm(idxBlock2_First);
block2_LastTime = eventTABLE.TrialEvTm(idxBlock2_Last);

block3_FirstTime = eventTABLE.TrialEvTm(idxBlock3_First);
block3_LastTime = eventTABLE.TrialEvTm(idxBlock3_Last);

block4_FirstTime = eventTABLE.TrialEvTm(idxBlock4_First);
block4_LastTime = eventTABLE.TrialEvTm(idxBlock4_Last);

block5_FirstTime = eventTABLE.TrialEvTm(idxBlock5_First);
block5_LastTime = eventTABLE.TrialEvTm(idxBlock5_Last);

% Get out the location of those time stamps for the behavior marker
[~, block1_FirstLoc] = min(abs(ma_timestampsDS - block1_FirstTime));
[~, block1_LastLoc] = min(abs(ma_timestampsDS - block1_LastTime));

[~, block2_FirstLoc] = min(abs(ma_timestampsDS - block2_FirstTime));
[~, block2_LastLoc] = min(abs(ma_timestampsDS - block2_LastTime));

[~, block3_FirstLoc] = min(abs(ma_timestampsDS - block3_FirstTime));
[~, block3_LastLoc] = min(abs(ma_timestampsDS - block3_LastTime));

[~, block4_FirstLoc] = min(abs(ma_timestampsDS - block4_FirstTime));
[~, block4_LastLoc] = min(abs(ma_timestampsDS - block4_LastTime));

[~, block5_FirstLoc] = min(abs(ma_timestampsDS - block5_FirstTime));
[~, block5_LastLoc] = min(abs(ma_timestampsDS - block5_LastTime));

% Get out ephys for each behavioral event / screen seen
block1_ephys = cleanVolts(:, block1_FirstLoc:block1_LastLoc);
block2_ephys = cleanVolts(:, block2_FirstLoc:block2_LastLoc);
block3_ephys = cleanVolts(:, block3_FirstLoc:block3_LastLoc);
block4_ephys = cleanVolts(:, block4_FirstLoc:block4_LastLoc);
block5_ephys = cleanVolts(:, block5_FirstLoc:block5_LastLoc);

allEphys = [block1_ephys block2_ephys block3_ephys block4_ephys block5_ephys];

% Save 
saveLoc = [strcat(synologyCD,'\', subjectID, '\', 'ProcessedEphys')];
cd(saveLoc)

saveName = [strcat(subjectID, '_',  hemi2use, '_', saveBrainArea, '_', 'allEphys.mat')];

save(saveName, "allEphys");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- EPHYS ---- % 

% idxChoice = find(eventTABLE.TrialEvNum(:) == 1);        % Find rows where choice show are
% idxOutExit = find(eventTABLE.TrialEvNum(:) == 5);
% 
% trialNum = num2cell(transpose(1:135));
% ephysTEMP = double.empty;
% 
% for i = 1:length(trialNum)
% 
%     % Get out time stamps for each behavioral marker
%     choiceTime = eventTABLE.TrialEvTm(idxChoice(i));
%     outExitTime = eventTABLE.TrialEvTm(idxOutExit(i));
% 
%     % Get out the location of those time stamps for the behavior marker
%     [~, choiceEphysLoc] = min(abs(ma_timestampsDS - choiceTime));
%     [~, outExitEphysLoc] = min(abs(ma_timestampsDS - outExitTime));
% 
%     % Pre Choice Buffer
% 
%     % Equation:
%     % sampling rate * number of seconds (x) = number of samples needed
%     % EG: 500 * 0.25 = 125
% 
%     tempPreTimeBeh = choiceEphysLoc - 500;          % One second buffer
% 
%     % Get out ephys for each trial
% 
%     tmpTrialEphys = cleanVolts(:, tempPreTimeBeh:outExitEphysLoc);
% 
%     ephysTEMP{i,1} = tmpTrialEphys;
% 
% end % for / i
% 
% ephysTAB = cell2table([trialNum, ephysTEMP], "VariableNames",["trialNum" "ephys"]);

% --- 1. Parameters ------------------------------------------------------- %
fs = 500;
preBufferSec     = 1.0;                          % pre-choice buffer (s)
preBufferSamples = round(fs * preBufferSec);     % in samples

% --- 2. Locate event rows ------------------------------------------------ %
idxChoiceShow = find(eventTABLE.TrialEvNum == 1);   % gamble onset
idxOutDispEnd = find(eventTABLE.TrialEvNum == 5);   % outcome offset

nTrials = length(idxChoiceShow);
assert(length(idxOutDispEnd) == nTrials, ...
    'Mismatch: choiceShow and outDispEnd counts differ');

fprintf('[WholeTrial] Extracting %d trials with %.1fs pre-choice buffer.\n', ...
    nTrials, preBufferSec);

% --- 3. Pre-allocate ----------------------------------------------------- %
trialEphys  = cell(nTrials, 1);
nSamplesLog = zeros(nTrials, 1);
clippedLog  = false(nTrials, 1);   % flag trials where buffer was clipped

% --- 4. Extract ---------------------------------------------------------- %
for ci = 1:nTrials

    % Retrieve timestamps
    t_choiceShow = eventTABLE.TrialEvTm(idxChoiceShow(ci));
    t_outDispEnd = eventTABLE.TrialEvTm(idxOutDispEnd(ci));

    % Map to nearest sample indices
    [~, loc_choiceShow] = min(abs(ma_timestampsDS - t_choiceShow));
    [~, loc_outDispEnd] = min(abs(ma_timestampsDS - t_outDispEnd));

    % Apply pre-choice buffer with boundary guard
    startLoc = loc_choiceShow - preBufferSamples;
    if startLoc < 1
        warning('[WholeTrial] Trial %d: pre-buffer clipped to recording start.', ci);
        startLoc    = 1;
        clippedLog(ci) = true;
    end

    % Validate ordering
    assert(startLoc < loc_outDispEnd, ...
        'Trial %d: startLoc not before outDispEnd', ci);

    % Slice
    trialEphys{ci}  = cleanVolts(:, startLoc : loc_outDispEnd);
    nSamplesLog(ci) = size(trialEphys{ci}, 2);

end % for ci

% --- 5. Assemble table --------------------------------------------------- %
ephysTAB = table( ...
    (1:nTrials)',   ...
    nSamplesLog,    ...
    trialEphys,     ...
    'VariableNames', {'TrialNum', 'nSamples', 'Ephys'});

% ---- BEHAVIOR ---- %

% Trial type info
% riskyloss < 0 = gain/loss trial :: either gain X or lose Y
% riskyloss == 0 = gain only :: either gain X or lose 0
% choice 1 = gamble, 0 = alternative

% Check trial
checkIndex = eventTab.subjdata.cs.ischecktrial;

% Loss aversion (LA) or risk aversion (RA)
% Gain/loss trials - this measures loss aversion
LA_trials = eventTab.subjdata.cs.riskyLoss < 0 & ~checkIndex;
% Gain only trials - this measures risk aversion
% RA_trials = eventTab.subjdata.cs.riskyLoss == 0 & ~checkIndex;

% Outcome results
outcomeLoss = eventTab.subjdata.cs.outcome < 0 & ~checkIndex;
outcomeGain = eventTab.subjdata.cs.outcome > 0 & ~checkIndex;
outcomeNeutral = eventTab.subjdata.cs.outcome == 0 & ~checkIndex;

BehaviorTab = table( LA_trials, outcomeLoss, outcomeGain, outcomeNeutral, ...
    'VariableNames',["LA" "OutcomeLoss" "OutcomeGain" "OutcomeNeutral"]);

% combine tables 

TrialEphysTab = [ephysTAB, BehaviorTab]; % save this variable 

% Save 
saveLoc = [strcat(synologyCD,'\', subjectID, '\', 'ProcessedEphys')];
cd(saveLoc)

saveName = [strcat(subjectID, '_',  hemi2use, '_', saveBrainArea, '_', 'Trial.mat')];

save(saveName, "TrialEphysTab");


end