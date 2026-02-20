
cd('Y:\LossAversion\LH_Data\FOOOF_data')
allTablesPart = load("FOOOF_alltables.mat");        % Load table
allTablesPart = allTablesPart.allPartTables;        % Rewrite variable

%%
Hemi = 'R';
BrainArea = 'AMY';

subTab_choice = allTablesPart(matches(allTablesPart.Hemi,Hemi) & ...
    matches(allTablesPart.BrainArea, BrainArea) & ...
    matches(allTablesPart.Epoch, 'Choice'),:);

    % choiceSTD = std(subTab_choice.Power);

subTab_Resp = allTablesPart(matches(allTablesPart.Hemi,Hemi) & ...
    matches(allTablesPart.BrainArea, BrainArea) & ...
    matches(allTablesPart.Epoch, 'Response'),:);

    % respSTD = std(subTab_Resp.Power);

subTab_OutGG = allTablesPart(matches(allTablesPart.Hemi,Hemi) & ...
    matches(allTablesPart.BrainArea, BrainArea) & ...
    matches(allTablesPart.Epoch, 'OutGG'),:);

    outGGSTD = std(subTab_OutGG.Power);

subTab_OutGL = allTablesPart(matches(allTablesPart.Hemi,Hemi) & ...
    matches(allTablesPart.BrainArea, BrainArea) & ...
    matches(allTablesPart.Epoch, 'OutGL'),:);

    outGLSTD = std(subTab_OutGL.Power);

subTab_OutAN = allTablesPart(matches(allTablesPart.Hemi,Hemi) & ...
    matches(allTablesPart.BrainArea, BrainArea) & ...
    matches(allTablesPart.Epoch, 'OutAN'),:);

    outANSTD = std(subTab_OutAN.Power);

unBands = unique(allTablesPart.Frequency);

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





%% Don't need to use below. Rewrote it so it's dynamic and more efficient 
%% allTablesPart = readtable("FOOOF_alltables_combined.csv");
allTablesPart = load("FOOOF_alltables.mat");

allTablesPart = allTablesPart.allPartTables;

%% LAH 
subTabLAH_choice = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'Choice'),:);

    LAH_choiceSTD = std(subTabLAH_choice.Power);

subTabLAH_Resp = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'Response'),:);

subTabLAH_OutGG = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'OutGG'),:);

subTabLAH_OutGL = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'OutGL'),:);

subTabLAH_OutAN = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'OutAN'),:);

unBands = unique(subTabLAH_choice.Frequency);

for i = 1:length(unBands)

    xTick = i;
    yTickChoice = mean(subTabLAH_choice.Power(matches(subTabLAH_choice.Frequency,unBands{i}),:));
    scatter(xTick,yTickChoice, 100, "blue", 'filled')

    % add error bar with the same color as the scatter 
    hold on
    h = errorbar(xTick, yTickChoice, LAH_choiceSTD, 'vertical');
    h.Color = 'blue';  % This sets the color of the error bars and the line


    hold on
    yTickResp = mean(subTabLAH_Resp.Power(matches(subTabLAH_Resp.Frequency,unBands{i}),:));
    scatter(xTick,yTickResp,100, "magenta", 'filled')

    hold on
    yTickOutGG = mean(subTabLAH_OutGG.Power(matches(subTabLAH_OutGG.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGG,100, "green", 'filled')

    hold on 
    yTickOutGL = mean(subTabLAH_OutGL.Power(matches(subTabLAH_OutGL.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGL,100, "red", 'filled')

    hold on 
    yTickOutAN = mean(subTabLAH_OutAN.Power(matches(subTabLAH_OutAN.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutAN,100, "cyan", 'filled')

end % for / i 

xticks(1:4)
xticklabels(unBands)
xlim([0 5])
title('LAH')


%% LPH 
subTabLPH_choice = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'Choice'),:);

subTabLPH_Resp = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'Response'),:);

subTabLPH_OutGG = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'OutGG'),:);

subTabLPH_OutGL = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'OutGL'),:);

subTabLPH_OutAN = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'OutAN'),:);

unBands = unique(subTabLPH_choice.Frequency);

for i = 1:length(unBands)

    xTick = i;
    yTickChoice = mean(subTabLPH_choice.Power(matches(subTabLPH_choice.Frequency,unBands{i}),:));
    scatter(xTick,yTickChoice, 100, "blue", 'filled')
    
    hold on
    yTickResp = mean(subTabLPH_Resp.Power(matches(subTabLPH_Resp.Frequency,unBands{i}),:));
    scatter(xTick,yTickResp,100, "magenta", 'filled')

    hold on
    yTickOutGG = mean(subTabLPH_OutGG.Power(matches(subTabLPH_OutGG.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGG,100, "green", 'filled')

    hold on 
    yTickOutGL = mean(subTabLPH_OutGL.Power(matches(subTabLPH_OutGL.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGL,100, "red", 'filled')

    hold on 
    yTickOutAN = mean(subTabLPH_OutAN.Power(matches(subTabLPH_OutAN.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutAN,100, "cyan", 'filled')

end % for / i 

xticks(1:4)
xticklabels(unBands)
xlim([0 5])
title('LPH')


%%
% RAH 
subTabRAH_choice = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'Choice'),:);

subTabRAH_Resp = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'Response'),:);

subTabRAH_OutGG = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'OutGG'),:);

subTabRAH_OutGL = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'OutGL'),:);

subTabRAH_OutAN = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'OutAN'),:);

unBands = unique(subTabRAH_choice.Frequency);

for i = 1:length(unBands)

    xTick = i;
    yTickChoice = mean(subTabRAH_choice.Power(matches(subTabRAH_choice.Frequency,unBands{i}),:));
    scatter(xTick,yTickChoice, 100, "blue", 'filled')
    
    hold on
    yTickResp = mean(subTabRAH_Resp.Power(matches(subTabRAH_Resp.Frequency,unBands{i}),:));
    scatter(xTick,yTickResp,100, "magenta", 'filled')

    hold on
    yTickOutGG = mean(subTabRAH_OutGG.Power(matches(subTabRAH_OutGG.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGG,100, "green", 'filled')

    hold on 
    yTickOutGL = mean(subTabRAH_OutGL.Power(matches(subTabRAH_OutGL.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGL,100, "red", 'filled')

    hold on 
    yTickOutAN = mean(subTabRAH_OutAN.Power(matches(subTabRAH_OutAN.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutAN,100, "cyan", 'filled')

end % for / i 

xticks(1:4)
xticklabels(unBands)
xlim([0 5])
title('RAH')

%%
% RPH 
subTabRPH_choice = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'Choice'),:);

subTabRPH_Resp = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'Response'),:);

subTabRPH_OutGG = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'OutGG'),:);

subTabRPH_OutGL = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'OutGL'),:);

subTabRPH_OutAN = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'PH') & ...
    matches(allTablesPart.Epoch, 'OutAN'),:);

unBands = unique(subTabRPH_choice.Frequency);

for i = 1:length(unBands)

    xTick = i;
    yTickChoice = mean(subTabRPH_choice.Power(matches(subTabRPH_choice.Frequency,unBands{i}),:));
    scatter(xTick,yTickChoice, 100, "blue", 'filled')
    
    hold on
    yTickResp = mean(subTabRPH_Resp.Power(matches(subTabRPH_Resp.Frequency,unBands{i}),:));
    scatter(xTick,yTickResp,100, "magenta", 'filled')

    hold on
    yTickOutGG = mean(subTabRPH_OutGG.Power(matches(subTabRPH_OutGG.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGG,100, "green", 'filled')

    hold on 
    yTickOutGL = mean(subTabRPH_OutGL.Power(matches(subTabRPH_OutGL.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGL,100, "red", 'filled')

    hold on 
    yTickOutAN = mean(subTabRPH_OutAN.Power(matches(subTabRPH_OutAN.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutAN,100, "cyan", 'filled')

end % for / i 

xticks(1:4)
xticklabels(unBands)
xlim([0 5])
title('RPH')
>>>>>>> origin/main


%% LAMY 
subTabLAmy_choice = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'Choice'),:);

subTabLAmy_Resp = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'Response'),:);

subTabLAmy_OutGG = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'OutGG'),:);

subTabLAmy_OutGL = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'OutGL'),:);

subTabLAmy_OutAN = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'OutAN'),:);

unBands = unique(subTabLAmy_choice.Frequency);

for i = 1:length(unBands)

    xTick = i;
    yTickChoice = mean(subTabLAmy_choice.Power(matches(subTabLAmy_choice.Frequency,unBands{i}),:));
    scatter(xTick,yTickChoice, 100, "blue", 'filled')
    
    hold on
    yTickResp = mean(subTabLAmy_Resp.Power(matches(subTabLAmy_Resp.Frequency,unBands{i}),:));
    scatter(xTick,yTickResp,100, "magenta", 'filled')

    hold on
    yTickOutGG = mean(subTabLAmy_OutGG.Power(matches(subTabLAmy_OutGG.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGG,100, "green", 'filled')

    hold on 
    yTickOutGL = mean(subTabLAmy_OutGL.Power(matches(subTabLAmy_OutGL.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGL,100, "red", 'filled')

    hold on 
    yTickOutAN = mean(subTabLAmy_OutAN.Power(matches(subTabLAmy_OutAN.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutAN,100, "cyan", 'filled')

end % for / i 

xticks(1:4)
xticklabels(unBands)
xlim([0 5])
title('LAMY')

%% RAMY 
subTabRAmy_choice = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'Choice'),:);

subTabRAmy_Resp = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'Response'),:);

subTabRAmy_OutGG = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'OutGG'),:);

subTabRAmy_OutGL = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'OutGL'),:);

subTabRAmy_OutAN = allTablesPart(matches(allTablesPart.Hemi,'R') & ...
    matches(allTablesPart.BrainArea, 'AMY') & ...
    matches(allTablesPart.Epoch, 'OutAN'),:);

unBands = unique(subTabRAmy_choice.Frequency);

for i = 1:length(unBands)

    xTick = i;
    yTickChoice = mean(subTabRAmy_choice.Power(matches(subTabRAmy_choice.Frequency,unBands{i}),:));
    scatter(xTick,yTickChoice, 100, "blue", 'filled')
    
    hold on
    yTickResp = mean(subTabRAmy_Resp.Power(matches(subTabRAmy_Resp.Frequency,unBands{i}),:));
    scatter(xTick,yTickResp,100, "magenta", 'filled')

    hold on
    yTickOutGG = mean(subTabRAmy_OutGG.Power(matches(subTabRAmy_OutGG.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGG,100, "green", 'filled')

    hold on 
    yTickOutGL = mean(subTabRAmy_OutGL.Power(matches(subTabRAmy_OutGL.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutGL,100, "red", 'filled')

    hold on 
    yTickOutAN = mean(subTabRAmy_OutAN.Power(matches(subTabRAmy_OutAN.Frequency,unBands{i}),:));
    scatter(xTick,yTickOutAN,100, "cyan", 'filled')

end % for / i 

xticks(1:4)
xticklabels(unBands)
xlim([0 5])
title('RAMY')

%% 










%% With JAT 
subTabRAmy_choice = allTablesPart(matches(allTablesPart.Hemi,'L') & ...
    matches(allTablesPart.BrainArea, 'AH') & ...
    matches(allTablesPart.Epoch, 'Choice'),:);

unBands = unique(subTabRAmy_choice.Frequency);

for i = 1:length(unBands)

    xTick = i;
    yTickChoice = mean(subTabRAmy_choice.Power(matches(subTabRAmy_choice.Frequency,unBands{i}),:));
    scatter(xTick,yTickChoice, 100, "black", 'filled')
    hold on 
    

end % for / i 

xticks(1:4)
xticklabels(unBands)
xlim([0 5])









