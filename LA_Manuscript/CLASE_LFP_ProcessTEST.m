function [] = CLASE_LFP_ProcessTEST(subjectID)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% USE AS EXAMPLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% - Behavior
behaviorLOC = [mainDIR_clase , filesep , 'Patient folders' , filesep ,...
    subjectID,'\Behavioral-data\EventBehavior'];
cd(behaviorLOC)
tmpMATdir = dir('*.mat');
tmpMATFl = {tmpMATdir.name};
load(tmpMATFl{1},'eventTABLE')



%%%% - NWB
nwbLOC = [mainDIR_clase , filesep , 'Patient folders' , filesep ,...
    subjectID,'\NWB-processing\NWB_Data'];
cd(nwbLOC)

tmpNwbdir = dir('*.nwb');
tmpNwbFl = {tmpNwbdir.name};
tmpNwbFl2 = tmpNwbFl{contains(tmpNwbFl,'filter')};
nwbDATA = nwbRead(tmpNwbFl2);

% Loads in Macrowire timestamps
LFP_timestamps = nwbDATA.processing.get('ecephys').nwbdatainterface.get...
    ('LFP').electricalseries.get('MacroWireSeries').timestamps.load;
% Voltage data for all macrowires and their channels
LFP_data = nwbDATA.processing.get('ecephys').nwbdatainterface.get...
    ('LFP').electricalseries.get('MacroWireSeries').data.load;
% To get sampling frequency info you get it from the description
% LFP_sessionInfo = tmp_LA.processing.get('ecephys').nwbdatainterface.get...
%     ('LFP').electricalseries.get('MacroWireSeries');
LFP_dataD = double(LFP_data);
LFP_time = downsample(LFP_timestamps,8);

% channel ID
chanLabels = cellstr(nwbDATA.general_extracellular_ephys_electrodes.vectordata.get('label').data.load()); %use MA only!
MAchan = find(contains(chanLabels,'MA_'));
chanID = cellstr(nwbDATA.general_extracellular_ephys_electrodes.vectordata.get('location').data.load());
hemisphere = cellstr(nwbDATA.general_extracellular_ephys_electrodes.vectordata.get('hemisph').data.load());
shortBnames = cellstr(nwbDATA.general_extracellular_ephys_electrodes.vectordata.get('shortBAn').data.load());
%%% INDEX and add to VAR
chanID2 = chanID(MAchan);
chanHemi = hemisphere(MAchan);
chanSname = shortBnames(MAchan);

%%%%% PUT WIRE and CONTACT SELECTION INTO FUNCTION
[accBPconP , bipLabsACC] = selectWireContacts('anterior cingulate' , chanSname , 3 ,...
    chanID2, LFP_dataD);

% --- Setup ---
% Fs = 500;         % Your sampling rate (Hz) - change this to match your data
% theta_band = [4 8];  % Frequency range for theta (Hz)
% 
% timeUU = 0:1/Fs:(length(accBPconP)-1)/500;
% 
% % Compute power spectral density (PSD) using Welch's method (or other method)
% window_size = Fs;  % Use 1 second window for power analysis (adjust as needed)
% noverlap = floor(window_size/2); % 50% overlap
% nfft = window_size; % use the whole window for the FFT
% 
% [psd_wire2, f] = pwelch(accBPconP, window_size, noverlap, nfft, Fs);

%%% QUESTIONS:







end








































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SUB FUNCTIONS

function [wireBPconP , bipLabsWireU] = selectWireContacts(brainName , shortNameS , contactPair ,...
    allchannels, alldata)

% WIRE 1
% brainName = 'amygdala';
% contactPair = 3;

outWire = matches(allchannels,brainName);
outData = alldata(outWire,:);
numchanWire = height(outData);
shortNameF = shortNameS(outWire);
shortName = shortNameF{1};
% Create channel labels for bipolar data
bipLabsWire = cell(numchanWire - 1, 1);
for i = 1:(numchanWire - 1)
    labString = [shortName,'%d - ',shortName,'%d'];
    % bipLabsWire{i} = sprintf('AMY%d - AMY%d', i + 1, i); % e.g., 'Ch2 - Ch1'
    bipLabsWire{i} = sprintf(labString, i + 1, i);
end
% Bipolar referecing
wireDataBP = diff(outData,1,1);
wireBPconP = wireDataBP(contactPair,:);

bipLabsWireU = bipLabsWire{contactPair};

end