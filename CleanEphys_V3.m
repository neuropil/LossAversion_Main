function [cleanVolts] = CleanEphys_V3(tmp_LA, ma_data, shortBA, std_thresh, numCon)
% CleanEphys_V3  Artifact rejection and bipolar re-referencing for iEEG data.
%
% Inputs:
%   tmp_LA      - NWB data structure containing electrode metadata
%   ma_data     - Raw voltage matrix [electrodes x samples] (macro contacts only)
%   shortBA     - Short brain area name string (e.g., 'LAMY')
%   std_thresh  - Standard deviation threshold for artifact rejection (e.g., 6)
%   numCon      - Contact indices within the brain area (e.g., 1:3)
%
% Output:
%   cleanVolts  - Artifact-rejected, bipolar re-referenced voltage matrix
%                 [numCon-1 x samples] if bipolar; [numCon x samples] if single contact
%
% Notes:
%   Bipolar re-referencing produces N-1 channels from N contacts (adjacent pairs).

%% Extract electrode metadata

% Short brain area labels
shortBAname = tmp_LA.general_extracellular_ephys_electrodes.vectordata.get('shortBAn').data.load;

% Electrode type labels — filter out microwire (MI) contacts, keep macro (MA)
electrodeLabels    = cellstr(tmp_LA.general_extracellular_ephys_electrodes.vectordata.get('label').data.load);
isMicrowire        = contains(electrodeLabels, 'MI');
shortBAnameMA      = cellstr(shortBAname);
shortBAnameMA(isMicrowire) = [];

%% Select channels for target brain area

baFlag    = contains(shortBAnameMA, shortBA);
baVoltRaw = ma_data(baFlag, :);

%% Artifact rejection

% Extract requested contacts
tempVolt = baVoltRaw(numCon, :);

for i = 1:length(numCon)
    contact = double(tempVolt(i, :));

    % Compute threshold based on mean ± (std_thresh * SD)
    threshold     = mean(contact) + std(contact) * std_thresh;
    overThreshold = abs(contact) > threshold;
    fracOver      = sum(overThreshold) / length(contact);

    % Reject contact (set to NaN) if >4% of samples exceed threshold
    if fracOver > 0.04
        tempVolt(i, :) = NaN;
    end
end

artVolt = tempVolt;

%% Bipolar re-referencing

if length(numCon) > 1
    % Subtract adjacent channel pairs: channel(i+1) - channel(i)
    nPairs    = height(artVolt) - 1;
    baVoltBI  = zeros(nPairs, width(artVolt));

    for bi = 1:nPairs
        baVoltBI(bi, :) = artVolt(bi+1, :) - artVolt(bi, :);
    end

    cleanVolts = baVoltBI;
else
    cleanVolts = double(artVolt);
end

end
