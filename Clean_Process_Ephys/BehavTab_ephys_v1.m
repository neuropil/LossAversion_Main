% This will create a table called EpochEphysTab which has each epoch in the experiment
% bipolar refrenced and artifact rejected. It also has all of the behavior
% info added to the table. Save this table in the ProcessedEphys folder in
% each participant's file. 

% Inputs 
tempPtID = 'CLASE035'; % Patient ID
Hemi = 'L';
BrainArea = 'LOF';

%% Load in ephys that is artifact rejected and bipolar refrenced 

[cleanVolts] = CleanEphys_v2('CLASE035', 'LOF', 6, 2:7);
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
end

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
load(tmpBehavName)           % Loads as eventTABLE 

% Load behavioral participant data %
paths.BehFolder = [strcat(tmpPTpath,'Behavioral-data', '\')];
cd(paths.BehFolder)

behDIR = dir; % behavioral folder DIR
behDIRNames = {behDIR.name}; % Names of files in behavioral folder
behDIRmat = contains(behDIRNames, 'mat'); % Find the files that have ".mat" in the name
tempBehavName = string(behDIRNames(behDIRmat));

eventTab = load(tempBehavName);         % load behavioral time table file


% ---- EPOCHS ---- %
% Seperate ephys per epochs

idxChoice = find(eventTABLE.TrialEvNum(:) == 1);        % Find rows where choice show are
idxRespShow = find(eventTABLE.TrialEvNum(:) == 2);      % Rows where response window shows 
idxRespExit = find(eventTABLE.TrialEvNum(:) == 3);
idxOutShow = find(eventTABLE.TrialEvNum(:) == 4);
idxOutExit = find(eventTABLE.TrialEvNum(:) == 5);

epochRepeat = repmat({'choice', 'response', 'outcome'}, 1,135)';

epochEphys = double.empty; 

for ci = 1:length(idxChoice)

    % Get out time stamps for each behavioral marker
    choiceTime = eventTABLE.TrialEvTm(idxChoice(ci));
    respShowTime = eventTABLE.TrialEvTm(idxRespShow(ci));
    respExitTime = eventTABLE.TrialEvTm(idxRespExit(ci));
    outShowTime = eventTABLE.TrialEvTm(idxOutShow(ci));
    outExitTime = eventTABLE.TrialEvTm(idxOutExit(ci));

    % Get out the location of those time stamps for the behavior marker
    [~, choiceEphysLoc] = min(abs(ma_timestampsDS - choiceTime));
    [~, respShowEphysLoc] = min(abs(ma_timestampsDS - respShowTime));
    [~, respExitEphysLoc] = min(abs(ma_timestampsDS - respExitTime));
    [~, outShowEphysLoc] = min(abs(ma_timestampsDS - outShowTime));
    [~, outExitEphysLoc] = min(abs(ma_timestampsDS - outExitTime));
    
    % Get out ephys for each behavioral event / screen seen
    choiceEphys = cleanVolts(:, choiceEphysLoc:respShowEphysLoc);
    respEphys = cleanVolts(:, respShowEphysLoc:respExitEphysLoc);
    outEphys = cleanVolts(:, outShowEphysLoc:outExitEphysLoc);


    for ii = 1 

    % Add Ephys to epoch location 
    epochEphys2{ii,1} = choiceEphys;        % Choice
    epochEphys2{(ii+1),1} = respEphys;      % Response
    epochEphys2{(ii+2),1} = outEphys;       % Outcome

    epochEphys = [epochEphys; epochEphys2];

    epochEphys2 = double.empty;

    end % for / ii 

end % for / ci 

Epoch_Ephys = [epochRepeat, epochEphys];        % Ephys by epochs 
Epoch_EphysTAB = cell2table(Epoch_Ephys, "VariableNames",["EpochID" "Ephys"]);

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
saveLoc = [strcat(synologyCD,'\', tempPtID, '\', 'ProcessedEphys')];
cd(saveLoc)

saveName = [strcat(tempPtID, '_',  Hemi, '_', BrainArea, '_', 'EpochEphys.mat')];

save(saveName, "EpochEphysTab");