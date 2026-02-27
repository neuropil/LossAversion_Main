function [trialTableOUT] = Trial_Generate_v1(tempPtID , Hemi , BrainArea , contactNUMs)

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
load(tmpBehavName,'eventTABLE')           % Loads as eventTABLE 

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

% trialTABLE = table;

allTrialTables = cell(numel(trialSTARTS_inds),1);
allTrialChanns = cell(numel(trialSTARTS_inds),1);
allTrialTimes = cell(numel(trialSTARTS_inds),1);

for ii = 1:length(trialSTARTS_inds) % FOR EACH TRIAL

    tmpTStartInd = trialSTARTS_inds(ii);

    tmpStartTS = eventTABLE.TrialEvTm(tmpTStartInd);

    [~, trialStartI_NLX] = min(abs(ma_timestampsDS - tmpStartTS));

    tmpTrialTab = eventTABLE(eventTABLE.TrialID == ii,:);
    tmpTrialTabRS = tmpTrialTab(~matches(tmpTrialTab.TrialEvID,'choiceShow'),:);

    tmpTrialTabRSC = removevars(tmpTrialTabRS,["Trials","OffsetSecs"]);

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
    allTrialChanns{ii} = cleanVolts(:, tmpTrialTabRSC.AbsolINDEX(1):...
        tmpTrialTabRSC.AbsolINDEX(end));
    allTrialTimes{ii} = ma_timestampsDS(tmpTrialTabRSC.AbsolINDEX(1):...
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