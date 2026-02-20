
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
        saveCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';
end

%% Inputs 
tempPtID = 'CLASE035'; % Patient ID
Hemi = 'L';
BrainArea = 'LOF';
%%
[cleanVolts] = CleanEphys_v2('CLASE035', 'LOF', 6, 2:7);

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

% Get out ephys for each behavioral event / screen seen. And Average over
% contacts 

if height(cleanVolts) > 1
    block1_ephys = mean(cleanVolts(:, block1_FirstLoc:block1_LastLoc));
    block2_ephys = mean(cleanVolts(:, block2_FirstLoc:block2_LastLoc));
    block3_ephys = mean(cleanVolts(:, block3_FirstLoc:block3_LastLoc));
    block4_ephys = mean(cleanVolts(:, block4_FirstLoc:block4_LastLoc));
    block5_ephys = mean(cleanVolts(:, block5_FirstLoc:block5_LastLoc));

else
    block1_ephys = cleanVolts(:, block1_FirstLoc:block1_LastLoc);
    block2_ephys = cleanVolts(:, block2_FirstLoc:block2_LastLoc);
    block3_ephys = cleanVolts(:, block3_FirstLoc:block3_LastLoc);
    block4_ephys = cleanVolts(:, block4_FirstLoc:block4_LastLoc);
    block5_ephys = cleanVolts(:, block5_FirstLoc:block5_LastLoc);

end % if else

%% FOOOF 

% Set up struct 
FOOOFBlocks = struct;

% Run PSD before running FOOOF 
[psd_Block1, freqs_Block1] = pwelch(block1_ephys,hamming(128), 64, 512, 500);
[psd_Block2, freqs_Block2] = pwelch(block2_ephys,hamming(128), 64, 512, 500);
[psd_Block3, freqs_Block3] = pwelch(block3_ephys,hamming(128), 64, 512, 500);
[psd_Block4, freqs_Block4] = pwelch(block4_ephys,hamming(128), 64, 512, 500);
[psd_Block5, freqs_Block5] = pwelch(block5_ephys,hamming(128), 64, 512, 500);

% Transpose, to make inputs row vectors
psd_Block1 = psd_Block1';
freqs_Block1 = freqs_Block1';

psd_Block2 = psd_Block2';
freqs_Block2 = freqs_Block2';

psd_Block3 = psd_Block3';
freqs_Block3 = freqs_Block3';

psd_Block4 = psd_Block4';
freqs_Block4 = freqs_Block4';

psd_Block5 = psd_Block5';
freqs_Block5 = freqs_Block5';

% FOOOF %
% FOOOF settings
settings = struct();  % Use defaults
f_range = [1, 40];

fooof_results_Block1 = fooof(freqs_Block1, psd_Block1, f_range, settings, true);
fooof_results_Block2 = fooof(freqs_Block2, psd_Block2, f_range, settings, true);
fooof_results_Block3 = fooof(freqs_Block3, psd_Block3, f_range, settings, true);
fooof_results_Block4 = fooof(freqs_Block4, psd_Block4, f_range, settings, true);
fooof_results_Block5 = fooof(freqs_Block5, psd_Block5, f_range, settings, true);

% Add results to struct 
FOOOFBlocks.Block1 = fooof_results_Block1;
FOOOFBlocks.Block2 = fooof_results_Block2;
FOOOFBlocks.Block3 = fooof_results_Block3;
FOOOFBlocks.Block4 = fooof_results_Block4;
FOOOFBlocks.Block5 = fooof_results_Block5;

%% Save 
 % save name
saveName = [strcat(tempPtID, '_',  Hemi, '_', BrainArea, '_', 'FOOOFByBlocks.mat')];

 % CD to save location
 savePath = strcat(saveCD,tempPtID, '\Blocks');
 cd(savePath)

 save(saveName, "FOOOFBlocks");
