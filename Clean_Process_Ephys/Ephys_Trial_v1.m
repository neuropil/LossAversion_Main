% This will create a table called TrialEphysTab which has each trial in the experiment
% bipolar refrenced and artifact rejected. It also has all of the behavior
% info added to the table. Save this table in the ProcessedEphys folder in
% each participant's file. 

%% Inputs 

tempPtID = 'CLASE035'; % Patient ID
Hemi = 'L';
BrainArea = 'LOF';

%%

[cleanVolts] = CleanEphys_v2('CLASE035', 'LOF', 6, 2:7);

tempPtID = 'CLASE022'; % Patient ID
Hemi = 'L';
BrainArea = 'PH';

%%

[cleanVolts] = CleanEphys_v2('CLASE022', 'LPH', 6, 1:4);

%
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

% ---- EPHYS ---- % 

idxChoice = find(eventTABLE.TrialEvNum(:) == 1);        % Find rows where choice show are
idxOutExit = find(eventTABLE.TrialEvNum(:) == 5);

trialNum = num2cell(transpose(1:135));
ephysTEMP = double.empty;

for i = 1:length(trialNum)

    % Get out time stamps for each behavioral marker
    choiceTime = eventTABLE.TrialEvTm(idxChoice(i));
    outExitTime = eventTABLE.TrialEvTm(idxOutExit(i));

    % Get out the location of those time stamps for the behavior marker
    [~, choiceEphysLoc] = min(abs(ma_timestampsDS - choiceTime));
    [~, outExitEphysLoc] = min(abs(ma_timestampsDS - outExitTime));

    % Pre Choice Buffer

    % Equation:
    % sampling rate * number of seconds (x) = number of samples needed
    % EG: 500 * 0.25 = 125

    tempPreTimeBeh = choiceEphysLoc - 500;          % One second buffer

    % Get out ephys for each trial

    tmpTrialEphys = cleanVolts(:, tempPreTimeBeh:outExitEphysLoc);

    ephysTEMP{i,1} = tmpTrialEphys;

end % for / i

ephysTAB = cell2table([trialNum, ephysTEMP], "VariableNames",["trialNum" "ephys"]);

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
saveLoc = [strcat(synologyCD,'\', tempPtID, '\', 'ProcessedEphys')];
cd(saveLoc)

saveName = [strcat(tempPtID, '_',  Hemi, '_', BrainArea, '_', 'Trial.mat')];

save(saveName, "TrialEphysTab");