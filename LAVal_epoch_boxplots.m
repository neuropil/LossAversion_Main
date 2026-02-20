cd('Y:\LossAversion\LH_Data\FOOOF_data\PeriodicComponents')
load('PerodicByEpoch.mat')

%%
LAtab = allPartTablesFreq(allPartTablesFreq.LA == 1, :); % LA trials only 

% Recode 'LOF' and 'MOF' to 'OF'
ofIDX = ismember(LAtab{:,"BrainArea"}, {'LOF', 'MOF'}); % index 
LAtab{ofIDX, "BrainArea"} = {'OF'};

% Frequency Tab 
betaTab = LAtab(strcmp(LAtab{:, "FreqLabel"}, 'beta'), :);
thetaTab = LAtab(strcmp(LAtab{:, "FreqLabel"}, 'theta'), :);

%% Brain area tab 
BA = 'OF';
% subTab = betaTab(strcmp(betaTab{:, "BrainArea"}, BA), :);
subTab = thetaTab(strcmp(thetaTab{:, "BrainArea"}, BA), :);

% Remove outliers % 
% Identify outliers
isout = isoutlier(subTab.Power, 'quartiles');

% Replace outliers with NaN
data_with_nan = subTab.Power;
data_with_nan(isout) = NaN;

% Combine Condition and LAval into one label
groupLabels = strcat(subTab.EpochID, "-", subTab.LAval);

% Convert to categorical with specified order
groupCat = categorical(groupLabels, ...
    {'response-Low', 'response-High', 'outcome-Low', 'outcome-High'}, ...
    'Ordinal', true);

%
% Unique group labels
uniqueGroups = categories(groupCat);

% Create the boxchart
% figure;
% % boxchart(groupCat, subTab.Power);
% boxchart(groupCat, data_with_nan);
% xlabel('LA Value and Epoch');
% ylabel('Power');
% % title('Power grouped by LA Value and Condition');



%% Create figure and hold for multiple boxcharts
figure;
hold on;

% Loop through each group and plot separately
for i = 1:numel(uniqueGroups)
    thisGroup = uniqueGroups{i};
    idx = groupCat == thisGroup;

    % Choose color: red for 'High', blue otherwise
    if contains(thisGroup, 'High')
        color = 'r';
    else
        color = 'b';
    end

    % Plot boxchart for this group
    % boxchart(repmat(i, sum(idx), 1), data_with_nan(idx), 'BoxFaceColor', color);
       boxchart(repmat(i, sum(idx), 1), data_with_nan(idx), 'BoxFaceAlpha', 0, 'BoxEdgeColor',color, 'WhiskerLineColor',color, 'MarkerColor',color);
end
% Set x-axis labels
xticks(1:numel(uniqueGroups));
xticklabels({'Response', 'Response', 'Outcome', 'Outcome'});
% xlabel('LA Value and Epoch');
ylabel('Power');
ylim([0 1.2])
title('Orbitofrontal');
hold off;



%% Plot every brain area on the same plot by what epoch and LA type 
% Filter for 'Low' LAval and 'response' EpochID
lowLAtab = betaTab(betaTab.LAval == 'High', :);
tempTab = lowLAtab(matches(lowLAtab.EpochID, 'response'), :);

% Identify outliers
isout2 = isoutlier(tempTab.Power, 'quartiles');

% Replace outliers with NaN
data_with_nan2 = tempTab.Power;
data_with_nan2(isout2) = NaN;

color = 'r';

% Create boxchart
figure;
boxchart(categorical(tempTab.BrainArea), data_with_nan2, ...
    'BoxFaceAlpha', 0, ...
    'BoxEdgeColor', color, ...
    'WhiskerLineColor', color, ...
    'MarkerColor', color);

xlabel('Brain Area');
ylabel('Power');
ylim([0 1.8]);
title('Response Epoch');















