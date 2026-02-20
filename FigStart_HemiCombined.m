cd('Y:\LossAversion\LH_Data\PSD')
% use the '_v2' version for the updated pwelch inputs
load("AllPartPSDTab_v2.mat")
%% Subtable for brain area 

BA = 'OF';

subTab_BA_Neutral = allPSDtab(matches(allPSDtab.BrainArea,BA) ...
    & matches(allPSDtab.LAVal, 'Neutral'),:);
subTab_BA_High = allPSDtab(matches(allPSDtab.BrainArea,BA) ...
    & matches(allPSDtab.LAVal, 'High'),:);
subTab_BA_Low = allPSDtab(matches(allPSDtab.BrainArea,BA) ...
    & matches(allPSDtab.LAVal, 'Low'),:);

%% Orbital frontal brain area 

subTab_BA_Neutral = allPSDtab(ismember(allPSDtab.BrainArea, {'MOF', 'LOF'}) ...
    & matches(allPSDtab.LAVal, 'Neutral'), :);

subTab_BA_High = allPSDtab(ismember(allPSDtab.BrainArea, {'MOF', 'LOF'}) ...
    & matches(allPSDtab.LAVal, 'High'), :);

subTab_BA_Low = allPSDtab(ismember(allPSDtab.BrainArea, {'MOF', 'LOF'}) ...
    & matches(allPSDtab.LAVal, 'Low'), :);


    %% --- Left and Right combined ---- %%
% Neutral 
Neutral = [];

for i = 1:height(subTab_BA_Neutral)

    tempPSD = subTab_BA_Neutral.PSD{i}.Power;

    Neutral = [Neutral tempPSD];

end % for / i 

% High 
High = [];

for i = 1:height(subTab_BA_High)

    tempPSD = subTab_BA_High.PSD{i}.Power;

    High = [High tempPSD];

end % for / i 

% Low 
Low = [];

for i = 1:height(subTab_BA_Low)

    tempPSD = subTab_BA_Low.PSD{i}.Power;

    Low = [Low tempPSD];

end % for / i 

% Means %
NeutralMean = mean(Neutral, 2); 
LowMean = mean(Low, 2);
HighMean = mean(High, 2);

% STD %
NeutralSTD = std(Neutral, 0, 2);
HighSTD = std(High, 0, 2); 
LowSTD = std(Low, 0, 2);

% Upper / Lower 
upperNeuPSD = (NeutralMean + NeutralSTD)'; % Make row vector 
lowerNeuPSD = (NeutralMean - NeutralSTD)';

upperHighPSD = (HighMean + HighSTD)';
lowerHighPSD = (HighMean - HighSTD)';

upperLowPSD = (LowMean + LowSTD)';
lowerLowPSD = (LowMean - LowSTD)';

% Fxx
fxx = subTab_BA_Neutral.PSD{1}.Frequency'; % Make sure it's a row vector and not a column vector 
xMin = min(fxx);
xMax = max(fxx);

xDATA = [fxx fliplr(fxx)];
yDATANeutral = [lowerNeuPSD fliplr(upperNeuPSD)];
yDATAHigh = [lowerHighPSD fliplr(upperHighPSD)];
yDATALow = [lowerLowPSD fliplr(upperLowPSD)];

%% Plot 
% plot(fxx, NeutralMean, 'Color','green', 'LineWidth', 2)
% hold on 
% plot(fxx, HighMean, 'Color', 'red', 'LineWidth', 2)
% hold on 
% plot(fxx, LowMean, 'Color', 'blue', 'LineWidth', 2)
% hold on 
% xlim([xMin xMax])
% % ylim([0 5.5])
% yticks([2 3.75 5.5])
% 
% Neu1 = patch(xDATA, yDATANeutral,'g');
% Neu1.FaceAlpha = 0.2;
% Neu1.EdgeAlpha = 0;
% 
% High1 = patch(xDATA, yDATAHigh,'r');
% High1.FaceAlpha = 0.2;
% High1.EdgeAlpha = 0;
% 
% Low1 =  patch(xDATA, yDATALow,'b');
% Low1.FaceAlpha = 0.2;
% Low1.EdgeAlpha = 0; 
% 
% axis square

%% Plot 


% Plot lines
n1 = plot(fxx, NeutralMean, 'Color','green', 'LineWidth', 2)
hold on 
h1 = plot(fxx, HighMean, 'Color', 'red', 'LineWidth', 2)
hold on 
l1 = plot(fxx, LowMean, 'Color', 'blue', 'LineWidth', 2)
hold on 
ylabel('LFP Power')
xlabel('Frequency (Hz)')
%xlim([xMin xMax])
% ylim([0 5.5])
%yticks([2 3.5 5.5])

ylim([0 6])
yticks(0:1.5:5.5)
xlim([1 40])
xticks([1,5:5:40])



Neu1 = patch(xDATA, yDATANeutral,'g');
Neu1.FaceAlpha = 0.2;
Neu1.EdgeAlpha = 0;

High1 = patch(xDATA, yDATAHigh,'r');
High1.FaceAlpha = 0.2;
High1.EdgeAlpha = 0;

Low1 =  patch(xDATA, yDATALow,'b');
Low1.FaceAlpha = 0.2;
Low1.EdgeAlpha = 0; 

axis square

hold on 
title('Orbital Frontal')
 
% Add legend
legend([n1 h1 l1], {'Neutral', 'High', 'Low'}, 'Location', 'northeast')
