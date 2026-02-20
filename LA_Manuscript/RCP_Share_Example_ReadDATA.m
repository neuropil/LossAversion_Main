%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% USE AS EXAMPLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% ******************************** IMPORTANT **************************
%%%%%%%%%% IF YOU USE MATLAB AND MATNWB
%%% USE MATLAB 2024a
%%% USE MATNWB 2.5.0
%%%%% ******************************** IMPORTANT **************************

%%%% - NWB --- Load
nwbLOC = ['DIRECTORY FOR NWB FILES'];
cd(nwbLOC)

%%%%%%%%%%%%%%%%%%%%%%%% Load filtered LFP data %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 500 Hz, high-pass filtered 0.5 Hz, with Notch filter
tmpNwbdir = dir('*.nwb');
tmpNwballname = {tmpNwbdir.name};
tmpNwbFilter = tmpNwballname{contains(tmpNwballname,'filter')};

%%%%%%%%%%%%%%%%%%%%%%%% Load RAW LFP data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4000 Hz, no filter
tmpNwbRaw = tmpNwballname{contains(tmpNwballname,'raw')};

%%%%%%%%%%%%%%%%%%%%%%%% PICK data to LOAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toLOAD = 'filt'; % 'raw'

switch toLOAD
    case 'filt'
        nwbDATASel = nwbRead(tmpNwbFilter);

        % Loads in Macrowire timestamps
        LFP_timestamps = nwbDATASel.processing.get('ecephys').nwbdatainterface.get...
            ('LFP').electricalseries.get('MacroWireSeries').timestamps.load;
        % Voltage data for all macrowires and their channels
        LFP_data = nwbDATASel.processing.get('ecephys').nwbdatainterface.get...
            ('LFP').electricalseries.get('MacroWireSeries').data.load;
        LFP_dataD = double(LFP_data);
        %%%%% IMPORTANT --- time stamps are sampled at 4kHz need to
        %%%%% be downsampled for FILTERED data.
        LFP_time = downsample(LFP_timestamps,8);

    case 'raw'
        nwbDATASel = nwbRead(tmpNwbRaw);

        % Loads in Macrowire timestamps
        LFP_timestamps = nwbDATASel.acquisition.get('MacroWireSeries').timestamps.load;
        % Voltage data for all macrowires and their channels
        LFP_data = nwbDATASel.acquisition.get('MacroWireSeries').data.load;
        LFP_dataD = double(LFP_data);

        LFP_time = LFP_timestamps;
end

%%%%% Channel ID --- USE 'MA' to get Macro contacts
chanLabels = cellstr(nwbDATASel.general_extracellular_ephys_electrodes.vectordata.get('label').data.load()); %use MA only!
MAchan = find(contains(chanLabels,'MA_'));
chanID = cellstr(nwbDATASel.general_extracellular_ephys_electrodes.vectordata.get('location').data.load());
hemisphere = cellstr(nwbDATASel.general_extracellular_ephys_electrodes.vectordata.get('hemisph').data.load());
shortBnames = cellstr(nwbDATASel.general_extracellular_ephys_electrodes.vectordata.get('shortBAn').data.load());
%%%% Use MA index to obtain appropriate channels
%%%%%%%% LONG names from Epilepsy
chanID2 = chanID(MAchan);
%%%%%%%% Hemisphere
chanHemi = hemisphere(MAchan);
%%%%%%%% ShortName from Epilepsy
chanSname = shortBnames(MAchan);


%%%%%%%%%%%%%%%%%%%%% BEHAVIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(['MAIN DIRECTORY with Behavior and MNI files'])
tmpMATBehdir = dir('*.mat');
tmpMATBehFl = {tmpMATBehdir.name};
tmpBEHavior = tmpMATBehFl{contains(tmpMATBehFl,'Behavior')};
load(tmpBEHavior,'eventTABLE','behavInfo')

% Use TrialEvTm to locate sample of interest in the LFP_time vector
% SAMPLE index for the START of Trial 1
[~,sampleINDEX] = min(abs(LFP_time - eventTABLE.TrialEvTm(1)));

% BEHAVIOR of interest
% Trial TYPE: 1) Loss , 2) Gain, 3) Neutral
% Gain/loss trials - this measures loss aversion
LA_trials = behavInfo.RiskyLossVal < 0 & ~checkIndex;
% Gain only trials - this measures risk aversion
RA_trials = behavInfo.RiskyLossVal == 0 & ~checkIndex;
gamble_trials = behavInfo.ChoiceSelected == 1 & ~checkIndex;
alternative_trials = behavInfo.ChoiceSelected == 0 & ~checkIndex;


outcomeLoss = eventTab.subjdata.cs.outcome < 0 & ~checkIndex;
outcomeNeutral = eventTab.subjdata.cs.outcome == 0 & ~checkIndex;
outcomeGain = eventTab.subjdata.cs.outcome > 0 & ~checkIndex;

%%%%%%%%%%%%%%%%%%%% MNI LABELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpBEHavior = tmpMATBehFl{contains(tmpMATBehFl,'MNI')};
load(tmpBEHavior,'tmpCSV')
% TABLE with channel IDS and MNI/Freesurfer labels








































