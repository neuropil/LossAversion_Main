cd('Y:\LossAversion\LH_Data\FOOOF_data\PeriodicComponents')
load('PeriodicByBlocks.mat');
%%

% One brain area 
BA = 'AH';
subTab = allTablesBlocks_Periodic(matches(allTablesBlocks_Periodic.BrainArea, BA), :); 

% Only for amygdala 
% subTab = allTablesBlocks_Periodic(matches(allTablesBlocks_Periodic.BrainArea, BA) & ...
%     matches(allTablesBlocks_Periodic.Hemi, 'L'), :);

% Orbital frontal 
% subTab = allTablesBlocks_Periodic(ismember(allTablesBlocks_Periodic.BrainArea, {'MOF', 'LOF'}), :);


% Define the desired frequency order
desiredFreqOrder = {'delta', 'theta', 'alpha', 'beta', 'gamma'};  

% Convert to categorical with specified order, then sort
UNfreq = unique(categorical(subTab.Frequency, desiredFreqOrder));
UNfreq = sort(UNfreq);

tiledlayout(1, length(UNfreq))

for i = 1:length(UNfreq)
    tempFreq = char(UNfreq(i));  % convert categorical to char
    tempFreqTab = subTab(matches(subTab.Frequency, tempFreq), :);

    LAcats = {'Low', 'Neutral', 'High'};

    yMax = max(tempFreqTab.Power);

    nexttile
    boxchart(categorical(tempFreqTab.LAVal, LAcats), tempFreqTab.Power, ...
        "BoxFaceAlpha", 0, 'MarkerStyle','none')
    % ylim([0 1.3])
    title(tempFreq)
end


%% Box charts with epoch types too 

BA = 'AMY';

subTab = allTablesBlocks_Periodic(matches(allTablesBlocks_Periodic.BrainArea, BA), :);

UNfreq = unique(subTab.Frequency);

tiledlayout(1, length(UNfreq))

for i = 1:length(UNfreq)

    tempFreq = UNfreq{i}; % temporary frequency name 

    tempFreqTab = subTab(matches(subTab.Frequency, tempFreq), :); % temp freq tab 

    LAcats = {'Low', 'Neutral', 'High'}; 
    epochCats = {'Choice', 'Response', 'OutGG', 'OutGL', 'OutAN'};

    yMax = max(tempFreqTab.Power);

    nexttile
    boxchart( ...
        categorical(tempFreqTab.LAVal, LAcats), ...
        tempFreqTab.Power, ...
        'GroupByColor', categorical(tempFreqTab.Epoch, epochCats))

    hold on 
    ylim([0 yMax])
    hold on
    title(tempFreq)
    % hold on 
    % legend(epochCats, 'Location', 'bestoutside')
    hold on 

end % for / i

hold on 
    legend(epochCats, 'Location', 'bestoutside')


%% plots each epoch on an individual row 

BA = 'AMY';
subTab = allTablesBlocks_Periodic(matches(allTablesBlocks_Periodic.BrainArea, BA), :);
UNfreq = unique(subTab.Frequency);

% Define Epochs to plot in separate rows
epochRows = {'Choice', 'Response', 'OutGG', 'OutGL', 'OutAN'};
numRows = length(epochRows);
numCols = length(UNfreq);

% Create tiled layout with 5 rows
tiledlayout(numRows, numCols)

LAcats = {'Low', 'Neutral', 'High'};

for i = 1:numCols
    tempFreq = UNfreq{i};
    tempFreqTab = subTab(matches(subTab.Frequency, tempFreq), :);
    % yMax = max(tempFreqTab.Power);

    for j = 1:numRows
        epochName = epochRows{j};
        epochTab = tempFreqTab(matches(tempFreqTab.Epoch, epochName), :);

        nexttile((j - 1) * numCols + i) % Calculate tile index
        boxchart( ...
            categorical(epochTab.LAVal, LAcats), ...
            epochTab.Power, ...
            'GroupByColor', categorical(epochTab.Epoch, epochRows))

        title([tempFreq ' - ' epochName])
        ylim([0 1.25])
        hold on
    end
end

% Add legend outside
legend(epochRows, 'Location', 'bestoutside')

%%
subTab_AMY_R = tempTab(matches(tempTab.BrainArea,'AMY') & ...
    matches(tempTab.Hemi, 'R'), :);


