
cd('Y:\LossAversion\LH_Data\FOOOF_data') % Lab computer
cd('Z:\LossAversion\LH_Data\FOOOF_data') % Laptop 
% load('aperiodic_HighVsLow_Neutral.mat') % loads as allPartTAB 
% AllFOOOFBlocks might be what i used to get more data points 

% Use this so you get more data points 
allPartTAB = load('FOOOF_blocks_allPart.mat'); 
allPartTAB = allPartTAB.FOOOFBlockTab;

%% ---- Amydala ---- % 
maxYOffset = max(allPartTAB.Offset);
minYOffset = min(allPartTAB.Offset);

maxYExp = max(allPartTAB.Exponent);
minYExp = min(allPartTAB.Exponent);

% Subtable %
subTab_AMY_R = allPartTAB(matches(allPartTAB.BrainArea,'AMY') & ...
    matches(allPartTAB.Hemi, 'R'), :);
subTab_AMY_L = allPartTAB(matches(allPartTAB.BrainArea, 'AMY') & ...
    matches(allPartTAB.Hemi, 'L'), :);

% Right 
tiledlayout(1,2)
nexttile
boxchart(categorical(subTab_AMY_R.LAVal), subTab_AMY_R.Offset, "BoxFaceAlpha", 0)
hold on 
scatter(categorical(subTab_AMY_R.LAVal), subTab_AMY_R.Offset, 30, 'b', "filled")
hold on 
minYOffset = round(minYOffset, 2);
maxYOffset = round(maxYOffset, 2);
ylim([minYOffset maxYOffset])% use round function 
ylabel('Right Offset')
% axis square

nexttile
boxchart(categorical(subTab_AMY_R.LAVal), subTab_AMY_R.Exponent, 'BoxFaceAlpha', 0)
hold on 
scatter(categorical(subTab_AMY_R.LAVal), subTab_AMY_R.Exponent, 30, 'blue', 'filled')
hold on 
ylim([minYExp maxYExp])
ylabel('Right Exponent')

figure
% Left 
tiledlayout(1,2)
nexttile
boxchart(categorical(subTab_AMY_L.LAVal), subTab_AMY_L.Offset, 'BoxFaceAlpha', 0)
hold on 
scatter(categorical(subTab_AMY_L.LAVal), subTab_AMY_L.Offset, 50, 'b', "filled")
hold on 
ylim([minYOffset maxYOffset])
ylabel('Left Offset')

nexttile
boxchart(categorical(subTab_AMY_L.LAVal), subTab_AMY_L.Exponent, 'BoxFaceAlpha', 0)
hold on 
scatter(categorical(subTab_AMY_L.LAVal), subTab_AMY_L.Exponent, 50, 'blue', 'filled')
hold on 
ylim([minYExp maxYExp])
ylabel('Left Exponent')






% ---- Left and Right combined ----- %
% maxYOffset = max(allPartTAB.Offset);
% minYOffset = min(allPartTAB.Offset);
% 
% maxYExp = max(allPartTAB.Exponent);
% minYExp = min(allPartTAB.Exponent);
% 
% % AMY 
% subTab_AMY = allPartTAB(matches(allPartTAB.BrainArea,'AMY'),:);
% 
% tiledlayout(1,2)
% nexttile
% boxchart(categorical(subTab_AMY.LAVal), subTab_AMY.Offset)
% hold on 
% NeuColor = repmat([0 0 1], sum(matches(subTab_AMY.LAVal, 'Neutral')), 1);
% HighColor  = repmat([1 0 0], sum(matches(subTab_AMY.LAVal, 'High')), 1);
% LowColor = repmat([0 1 0], sum(matches(subTab_AMY.LAVal, 'Low')), 1);
% ScatColor = [HighColor; LowColor; NeuColor];
% scatter(categorical(subTab_AMY.LAVal), subTab_AMY.Offset, 50, ScatColor, "filled")
% 
% hold on 
% ylim([minYOffset maxYOffset])
% ylabel('Offset')
% 
% nexttile
% boxchart(categorical(subTab_AMY.LAVal), subTab_AMY.Exponent)
% hold on 
% scatter(categorical(subTab_AMY.LAVal), subTab_AMY.Exponent, 50, 'blue', 'filled')
% hold on 
% ylim([minYExp maxYExp])
% ylabel('Exponent')

% AH
subTab_AH = allPartTAB(matches(allPartTAB.BrainArea,'AH'),:);

tiledlayout(1,2)
nexttile
boxchart(categorical(subTab_AH.LAVal), subTab_AH.Offset, 'BoxFaceAlpha', 0)
hold on 
scatter(categorical(subTab_AH.LAVal), subTab_AH.Offset, 50, 'blue', 'filled')
hold on
ylim([minYOffset maxYOffset])
ylabel('Offset')

nexttile
boxchart(categorical(subTab_AH.LAVal), subTab_AH.Exponent,  'BoxFaceAlpha', 0)
hold on 
scatter(categorical(subTab_AH.LAVal), subTab_AH.Exponent, 50, 'blue', 'filled')
hold on 
ylim([minYExp maxYExp])
ylabel('Exponent')

% PH 
subTab_PH = allPartTAB(matches(allPartTAB.BrainArea,'PH'),:);

tiledlayout(1,2)
nexttile
boxchart(categorical(subTab_PH.LAVal), subTab_PH.Offset,  'BoxFaceAlpha', 0)
hold on 
scatter(categorical(subTab_PH.LAVal), subTab_PH.Offset, 50, 'blue','filled')
hold on
ylim([minYOffset maxYOffset])
ylabel('Offset')

nexttile
boxchart(categorical(subTab_PH.LAVal), subTab_PH.Exponent,'BoxFaceAlpha', 0)
hold on 
scatter(categorical(subTab_PH.LAVal), subTab_PH.Exponent, 50, 'blue','filled')
hold on 
ylim([minYExp maxYExp])
ylabel('Exponent')

% MOF & LOF 
subTab_OF = allPartTAB(ismember(allPartTAB.BrainArea, {'MOF', 'LOF'}), :);

tiledlayout(1,2)
nexttile
boxchart(categorical(subTab_OF.LAVal), subTab_OF.Offset, 'BoxFaceAlpha', 0)
hold on 
scatter(categorical(subTab_OF.LAVal), subTab_OF.Offset, 50, 'blue', 'filled')
hold on 
ylim([minYOffset maxYOffset])
ylabel('Offset')

nexttile
boxchart(categorical(subTab_OF.LAVal), subTab_OF.Exponent, 'BoxFaceAlpha', 0)
hold on 
scatter(categorical(subTab_OF.LAVal), subTab_OF.Exponent, 50, 'blue', 'filled')
hold on
ylim([minYExp maxYExp])
ylabel('Exponent')

