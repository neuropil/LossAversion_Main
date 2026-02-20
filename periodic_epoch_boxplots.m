
cd('Y:\LossAversion\LH_Data\FOOOF_data\PeriodicComponents')
load('PerodicByEpoch.mat'); % loads as allPartTablesFreq

% Only loss aversion trials 
LAtab = allPartTablesFreq(allPartTablesFreq.LA == 1,:);

% LA group tables 
lowLA = LAtab(LAtab.LAval == 'Low',:);
highLA = LAtab(LAtab.LAval == 'High',:);

%% Boxplots
% One brain area 
BA = 'OF';
% subTab = LAtab(matches(LAtab.BrainArea, BA), :); 

% Only for amygdala 
% subTab = LAtab(matches(LAtab.BrainArea, BA) & ...
%     matches(LAtab.Hemi, 'R'), :);

% Orbital frontal 
subTab = LAtab(ismember(LAtab.BrainArea, {'MOF', 'LOF'}), :);

% Define the desired frequency order
desiredFreqOrder = {'delta', 'theta', 'alpha', 'beta', 'gamma'};  

% Convert to categorical with specified order, then sort
UNfreq = unique(categorical(subTab.FreqLabel, desiredFreqOrder));
UNfreq = sort(UNfreq);

tiledlayout(1, length(UNfreq))

for i = 1:length(UNfreq)
    tempFreq = char(UNfreq(i));  % convert categorical to char
    tempFreqTab = subTab(matches(subTab.FreqLabel, tempFreq), :);

    % LAcats = {'Low', 'Neutral', 'High'};
        LAcats = {'Low', 'High'};

    % yMax = max(tempFreqTab.Power);

    nexttile
    boxchart(categorical(tempFreqTab.LAval, LAcats), tempFreqTab.Power, ...
        "BoxFaceAlpha", 0, 'MarkerStyle','none')
    ylim([0 1.3])
    title(tempFreq)
end

%% Swarmchart 

tiledlayout(1, length(UNfreq))

for i = 1:length(UNfreq)
    tempFreq = char(UNfreq(i));  % convert categorical to char
    tempFreqTab = subTab(matches(subTab.FreqLabel, tempFreq), :);

    LAcats = {'Low', 'High'};

    yMax = max(tempFreqTab.Power);

    nexttile
    swarmchart(categorical(tempFreqTab.LAval, LAcats), tempFreqTab.Power)
    % ylim([0 1.3])
    title(tempFreq)
end
