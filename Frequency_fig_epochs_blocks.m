
cd('Y:\LossAversion\LH_Data\FOOOF_data')
load('FOOOF_allBlockTables.mat'); % Loads as allPartTablesBlocks
%%
Hemi = 'R';
BrainArea = 'AMY';
BlockNum = '5';

subTab_choice = allPartTablesBlocks(matches(allPartTablesBlocks.Hemi,Hemi) & ...
    matches(allPartTablesBlocks.BrainArea, BrainArea) & ...
    matches(allPartTablesBlocks.Epoch, 'Choice') & ...
    matches(string(allPartTablesBlocks.BlockNum), BlockNum),:);

subTab_Resp = allPartTablesBlocks(matches(allPartTablesBlocks.Hemi,Hemi) & ...
    matches(allPartTablesBlocks.BrainArea, BrainArea) & ...
    matches(allPartTablesBlocks.Epoch, 'Response') & ...
    matches(string(allPartTablesBlocks.BlockNum), BlockNum),:);

subTab_OutGG = allPartTablesBlocks(matches(allPartTablesBlocks.Hemi,Hemi) & ...
    matches(allPartTablesBlocks.BrainArea, BrainArea) & ...
    matches(allPartTablesBlocks.Epoch, 'OutGG') & ...
    matches(string(allPartTablesBlocks.BlockNum), BlockNum),:);

subTab_OutGL = allPartTablesBlocks(matches(allPartTablesBlocks.Hemi,Hemi) & ...
    matches(allPartTablesBlocks.BrainArea, BrainArea) & ...
    matches(allPartTablesBlocks.Epoch, 'OutGL') & ...
    matches(string(allPartTablesBlocks.BlockNum), BlockNum),:);

subTab_OutAN = allPartTablesBlocks(matches(allPartTablesBlocks.Hemi,Hemi) & ...
    matches(allPartTablesBlocks.BrainArea, BrainArea) & ...
    matches(allPartTablesBlocks.Epoch, 'OutAN') & ...
    matches(string(allPartTablesBlocks.BlockNum), BlockNum),:);

unBands = unique(allPartTablesBlocks.Frequency);

for i = 1:length(unBands)

    xTick = i;
    yTickChoice = mean(subTab_choice.Power(matches(subTab_choice.Frequency,unBands{i}),:));
    scatter(xTick,yTickChoice, 100, "blue", 'filled')

    hold on % Error bars
    choiceSTD = std(subTab_choice.Power(matches(subTab_choice.Frequency,unBands{i}),:));
    choiceError = errorbar(xTick, yTickChoice, choiceSTD, 'vertical');
    choiceError.Color = 'blue'; % This sets the color of the error bars and the line

    hold on
    yTickResp = mean(subTab_Resp.Power(matches(subTab_Resp.Frequency,unBands{i}),:));
    scatter(xTick,yTickResp,100, "magenta", 'filled')

    hold on % Error bars 
    respSTD = std(subTab_Resp.Power(matches(subTab_Resp.Frequency,unBands{i}),:));
    respError = errorbar(xTick, yTickResp, respSTD, 'vertical');
    respError.Color = 'magenta';

    hold on
    yTickOutGG = mean(subTab_OutGG.Power(matches(subTab_OutGG.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGG,100, "green", 'filled')

    hold on % Error bars 
    outGGSTD = std(subTab_OutGG.Power(matches(subTab_OutGG.Frequency,unBands{i}),:));
    outGGError = errorbar(xTick, yTickOutGG, outGGSTD, 'vertical');
    outGGError.Color = 'green';

    hold on 
    yTickOutGL = mean(subTab_OutGL.Power(matches(subTab_OutGL.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGL,100, "red", 'filled')

    hold on % Error bars 
    outGLSTD = std(subTab_OutGL.Power(matches(subTab_OutGL.Frequency,unBands{i}),:));
    outGLError = errorbar(xTick, yTickOutGL, outGLSTD, 'vertical');
    outGLError.Color = 'red';

    hold on 
    yTickOutAN = mean(subTab_OutAN.Power(matches(subTab_OutAN.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutAN,100, "cyan", 'filled')

    hold on % Error bars 
    outANSTD = std(subTab_OutAN.Power(matches(subTab_OutAN.Frequency,unBands{i}),:));
    outANError = errorbar(xTick, yTickOutAN, outANSTD, 'vertical');
    outANError.Color = 'cyan';

end % for / i 

xticks(1:4)
xticklabels(unBands)
xlim([0 5])
title('RAMY')