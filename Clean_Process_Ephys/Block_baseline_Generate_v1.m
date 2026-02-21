function [baseline] = Block_baseline_Generate_v1(tempPtID , Hemi , BrainArea , contactNUMs)

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

end

%% Inputs 
% tempPtID = 'CLASE035'; % Patient ID
% Hemi = 'L';
% BrainArea = 'OF';

%%
% [cleanVolts] = CleanEphys_v2('CLASE035', 'LOF', 6, 1:7);

[cleanVolts] = CleanEphys_v2(tempPtID, BrainArea, 6, contactNUMs);

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

%%
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
% block1_FirstTime = eventTABLE.TrialEvTm(idxBlock1_First);
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

% Time 
break1Time = (block2_FirstLoc - block1_LastLoc)/500;
break2Time = (block3_FirstLoc - block2_LastLoc)/500;
break3Time = (block4_FirstLoc - block3_LastLoc)/500;
break4Time = (block5_FirstLoc - block4_LastLoc)/500;

totalTime = break1Time + break2Time + break3Time + break4Time;

% Get out ephys for each behavioral event / screen seen
break1_ephys = cleanVolts(:, block1_LastLoc:block2_FirstLoc);
break2_ephys = cleanVolts(:, block2_LastLoc:block3_FirstLoc);
break3_ephys = cleanVolts(:, block3_LastLoc:block4_FirstLoc);
break4_ephys = cleanVolts(:, block4_LastLoc:block5_FirstLoc);

% 5 second Pre-block baselines
preBlock1_INDs = block1_FirstLoc - 2500:(block1_FirstLoc-1);
midBlk1a2 = round(median(block1_LastLoc:block2_FirstLoc));
preBlock2_INDs = midBlk1a2 - 1250:(midBlk1a2 + 1250)-1;
midBlk2a3 = round(median(block2_LastLoc:block3_FirstLoc));
preBlock3_INDs = midBlk2a3 - 1250:(midBlk2a3 + 1250)-1;
midBlk3a4 = round(median(block3_LastLoc:block4_FirstLoc));
preBlock4_INDs = midBlk3a4 - 1250:(midBlk3a4 + 1250)-1;
midBlk4a5 = round(median(block4_LastLoc:block5_FirstLoc));
preBlock5_INDs = midBlk4a5 - 1250:(midBlk4a5 + 1250)-1;

preBlock1 = cleanVolts(:,preBlock1_INDs);
preBlock2 = cleanVolts(:,preBlock2_INDs);
preBlock3 = cleanVolts(:,preBlock3_INDs);
preBlock4 = cleanVolts(:,preBlock4_INDs);
preBlock5 = cleanVolts(:,preBlock5_INDs);

allEphys = [break1_ephys break2_ephys break3_ephys break4_ephys];

% set up data to save 
baseline = struct;
% Break 1 
baseline.Break1.Break1Time = break1Time;
baseline.Break1.Break1Ephys = break1_ephys;
% Break 2 
baseline.Break2.Break2Time = break2Time;
baseline.Break2.Break2Ephys = break2_ephys;
% Break 3 
baseline.Break3.Break3Time = break3Time;
baseline.Break3.Break3Ephys = break3_ephys;
% Break 4 
baseline.Break4.Break4Time = break4Time;
baseline.Break4.Break4Ephys = break4_ephys;
% Combined
baseline.CombinedBaseline.CombinedEphys = allEphys;
baseline.CombinedBaseline.CombinedTime = totalTime;

baseline.PreBlock1.Indices = preBlock1_INDs;
baseline.PreBlock1.LFP     = preBlock1;
baseline.PreBlock2.Indices = preBlock2_INDs;
baseline.PreBlock2.LFP     = preBlock2;
baseline.PreBlock3.Indices = preBlock3_INDs;
baseline.PreBlock3.LFP     = preBlock3;
baseline.PreBlock4.Indices = preBlock4_INDs;
baseline.PreBlock4.LFP     = preBlock4;
baseline.PreBlock5.Indices = preBlock5_INDs;
baseline.PreBlock5.LFP     = preBlock5;
baseline.WholeLFP = cleanVolts;

%% Save

% saveLoc = [strcat('Y:\LossAversion\LH_Data\FOOOF_data\', tempPtID, '\', 'Baseline')];
% cd(saveLoc)
% 
% saveName = [strcat(tempPtID, '_',  Hemi, '_', BrainArea, '_', 'baseline.mat')];
% 
% save(saveName, "baseline");

end