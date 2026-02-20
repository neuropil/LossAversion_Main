cd('Y:\LossAversion\LH_Data\PSD')
% use the '_v2' version for the updated pwelch inputs
%%
% AMYGDALA %
subTab_AMY = allPSDtab(matches(allPSDtab.BrainArea,'AMY'),:);

subTab_AMY_Neutral_R = subTab_AMY(matches(subTab_AMY.LAVal, 'Neutral') ...
    & matches(subTab_AMY.Hemi, 'R'),:);
subTab_AMY_Neutral_L = subTab_AMY(matches(subTab_AMY.LAVal, 'Neutral') ...
    & matches(subTab_AMY.Hemi, 'L'),:);
subTab_AMY_High_R = subTab_AMY(matches(subTab_AMY.LAVal, 'High') ...
    & matches(subTab_AMY.Hemi, 'R'), :);
subTab_AMY_High_L = subTab_AMY(matches(subTab_AMY.LAVal, 'High') ...
    & matches(subTab_AMY.Hemi, 'L'), :);
subTab_AMY_Low_R = subTab_AMY(matches(subTab_AMY.LAVal, 'Low') ...
    & matches(subTab_AMY.Hemi, 'R'), :);
subTab_AMY_Low_L = subTab_AMY(matches(subTab_AMY.LAVal, 'Low') ...
    & matches(subTab_AMY.Hemi, 'L'), :);

% Neutral 
AMYNeutral_R = [];

for i = 1:height(subTab_AMY_Neutral_R)

    tempPSD = subTab_AMY_Neutral_R.PSD{i}.Power;

    AMYNeutral_R = [AMYNeutral_R tempPSD];

end % for / i 

AMYNeutral_L = [];

for i = 1:height(subTab_AMY_Neutral_L)

    tempPSD = subTab_AMY_Neutral_L.PSD{i}.Power;

    AMYNeutral_L = [AMYNeutral_L tempPSD];

end % for / i 


% High 
AMYHigh_R = [];

for i = 1:height(subTab_AMY_High_R)

    tempPSD = subTab_AMY_High_R.PSD{i}.Power;

    AMYHigh_R = [AMYHigh_R tempPSD];

end % for / i 


AMYHigh_L = [];

for i = 1:height(subTab_AMY_High_L)

    tempPSD = subTab_AMY_High_L.PSD{i}.Power;

    AMYHigh_L = [AMYHigh_L tempPSD];

end % for / i 


% Low 
AMYLow_R = [];

for i = 1:height(subTab_AMY_Low_R)

    tempPSD = subTab_AMY_Low_R.PSD{i}.Power;

    AMYLow_R = [AMYLow_R tempPSD];

end % for / i 

AMYLow_L = [];

for i = 1:height(subTab_AMY_Low_L)

    tempPSD = subTab_AMY_Low_L.PSD{i}.Power;

    AMYLow_L = [AMYLow_L tempPSD];

end % for / i 


% Means %
% Right 
AMYNeutralMean_R = mean(AMYNeutral_R, 2); 
AMYHighMean_R = mean(AMYHigh_R, 2);
AMYLowMean_R = mean(AMYLow_R, 2);

% Left
AMYNeutralMean_L = mean(AMYNeutral_L, 2); 
AMYHighMean_L = mean(AMYHigh_L, 2);
AMYLowMean_L = mean(AMYLow_L, 2);

% STD %
% Right 
AMYNeutralSTD_R = std(AMYNeutral_R, 0, 2);
AMYHighSTD_R = std(AMYHigh_R, 0, 2);
AMYLowSTD_R = std(AMYLow_R, 0, 2);
 
% Left
AMYNeutralSTD_L = std(AMYNeutral_L, 0, 2);
AMYHighSTD_L = std(AMYHigh_L, 0, 2); 
AMYLowSTD_L = std(AMYLow_L, 0, 2);

% Upper / Lower %
% Right
upperNeuPSD_R = (AMYNeutralMean_R + AMYNeutralSTD_R)'; % Make row vector 
lowerNeuPSD_R = (AMYNeutralMean_R - AMYNeutralSTD_R)';

upperHighPSD_R = (AMYHighMean_R + AMYHighSTD_R)';
lowerHighPSD_R = (AMYHighMean_R - AMYHighSTD_R)';

upperLowPSD_R = (AMYLowMean_R + AMYLowSTD_R)';
lowerLowPSD_R = (AMYLowMean_R - AMYLowSTD_R)';

% Left
upperNeuPSD_L = (AMYNeutralMean_L + AMYNeutralSTD_L)'; % Make row vector 
lowerNeuPSD_L = (AMYNeutralMean_L - AMYNeutralSTD_L)';

upperHighPSD_L = (AMYHighMean_L + AMYHighSTD_L)';
lowerHighPSD_L = (AMYHighMean_L - AMYHighSTD_L)';

upperLowPSD_L = (AMYLowMean_L + AMYLowSTD_L)';
lowerLowPSD_L = (AMYLowMean_L - AMYLowSTD_L)';

%%
% Plot % 

% Fxx
fxx = subTab_AMY_Neutral_R.PSD{1}.Frequency'; % Make sure it's a row vector and not a column vector 
% xMin = min(fxx);
% xMax = max(fxx);

% Right 
xDATA_R = [fxx fliplr(fxx)];
yDATANeutral_R = [lowerNeuPSD_R fliplr(upperNeuPSD_R)];
yDATAHigh_R = [lowerHighPSD_R fliplr(upperHighPSD_R)];
yDATALow_R = [lowerLowPSD_R fliplr(upperLowPSD_R)];

% Plot 
nR = plot(fxx, AMYNeutralMean_R, 'Color','green', 'LineWidth', 2)
hold on 
hR= plot(fxx, AMYHighMean_R, 'Color', 'red', 'LineWidth', 2)
hold on 
lR = plot(fxx, AMYLowMean_R, 'Color', 'blue', 'LineWidth', 2)
hold on 

ylabel('LFP Power')
xlabel('Frequency (Hz)')
% xlim([xMin xMax])
% yticks([2 3.5 5.5])
ylim([0 6])
yticks(0:1.5:5.5)
xlim([1 40])
xticks([1,5:5:40])


Neu_R = patch(xDATA_R, yDATANeutral_R,'g');
Neu_R.FaceAlpha = 0.2;
Neu_R.EdgeAlpha = 0;

High_R = patch(xDATA_R, yDATAHigh_R,'r');
High_R.FaceAlpha = 0.2;
High_R.EdgeAlpha = 0;

Low_R =  patch(xDATA_R, yDATALow_R,'b');
Low_R.FaceAlpha = 0.2;
Low_R.EdgeAlpha = 0;

axis square

hold on
title('Right Amygdala')

% Add legend 
legend([nR hR lR], {'Neutral', 'High', 'Low'}, 'Location', 'northeast')


hold off 





%% Left 
figure;

% Left 
xDATA_L = [fxx fliplr(fxx)];
yDATANeutral_L = [lowerNeuPSD_L fliplr(upperNeuPSD_L)];
yDATAHigh_L = [lowerHighPSD_L fliplr(upperHighPSD_L)];
yDATALow_L = [lowerLowPSD_L fliplr(upperLowPSD_L)];

% Plot 
nL = plot(fxx, AMYNeutralMean_L, 'Color','green', 'LineWidth', 2)
hold on 
hL = plot(fxx, AMYHighMean_L, 'Color', 'red', 'LineWidth', 2)
hold on 
lL = plot(fxx, AMYLowMean_L, 'Color', 'blue', 'LineWidth', 2)
hold on 
ylabel('LFP Power')
xlabel('Frequency (Hz)')

ylim([0 6])
yticks(0:1.5:5.5)
xlim([1 40])
xticks([1,5:5:40])

% xlim([xMin xMax])
% yticks([2 3.5 5.5])

Neu_L = patch(xDATA_L, yDATANeutral_L,'g');
Neu_L.FaceAlpha = 0.2;
Neu_L.EdgeAlpha = 0;

High_L = patch(xDATA_L, yDATAHigh_L,'r');
High_L.FaceAlpha = 0.2;
High_L.EdgeAlpha = 0;

Low_L =  patch(xDATA_L, yDATALow_L,'b');
Low_L.FaceAlpha = 0.2;
Low_L.EdgeAlpha = 0;

axis square

hold on 
title('Left Amygdala')

% Add legend 
legend([nL hL lL], {'Neutral', 'High', 'Low'}, 'Location', 'northeast')





            %% --- Left and Right combined ---- %%
% Neutral 
AMYNeutral = [];

for i = 1:height(subTab_AMY_Neutral_R)

    tempPSD = subTab_AMY_Neutral_R.PSD{i}.Power;

    AMYNeutral = [AMYNeutral tempPSD];

end % for / i 

% High 
AMYHigh = [];

for i = 1:height(subTab_AMY_High_R)

    tempPSD = subTab_AMY_High_R.PSD{i}.Power;

    AMYHigh = [AMYHigh tempPSD];

end % for / i 

% Low 
AMYLow = [];

for i = 1:height(subTab_AMY_Low)

    tempPSD = subTab_AMY_Low.PSD{i}.Power;

    AMYLow = [AMYLow tempPSD];

end % for / i 

% Means %
AMYNeutralMean = mean(AMYNeutral, 2); 
AMYLowMean = mean(AMYLow, 2);
AMYHighMean = mean(AMYHigh, 2);

% STD %
AMYNeutralSTD = std(AMYNeutral, 0, 2);
AMYHighSTD = std(AMYHigh, 0, 2); 
AMYLowSTD = std(AMYLow, 0, 2);

% Upper / Lower 
upperNeuPSD = (AMYNeutralMean + AMYNeutralSTD)'; % Make row vector 
lowerNeuPSD = (AMYNeutralMean - AMYNeutralSTD)';

upperHighPSD = (AMYHighMean + AMYHighSTD)';
lowerHighPSD = (AMYHighMean - AMYHighSTD)';

upperLowPSD = (AMYLowMean + AMYLowSTD)';
lowerLowPSD = (AMYLowMean - AMYLowSTD)';

% Fxx
fxx = subTab_AMY_Neutral_R.PSD{1}.Frequency'; % Make sure it's a row vector and not a column vector 

xDATA_L = [fxx fliplr(fxx)];
yDATANeutral_L = [lowerNeuPSD fliplr(upperNeuPSD)];
yDATAHigh = [lowerHighPSD fliplr(upperHighPSD)];
yDATALow = [lowerLowPSD fliplr(upperLowPSD)];

%% Plot 
plot(fxx, AMYNeutralMean, 'Color','green')
hold on 
plot(fxx, AMYHighMean, 'Color', 'red')
hold on 
plot(fxx, AMYLowMean, 'Color', 'blue')

Neu1 = patch(xDATA_L, yDATANeutral_L,'g');
Neu1.FaceAlpha = 0.3;

High1 = patch(xDATA_L, yDATAHigh,'r');
High1.FaceAlpha = 0.3;

Low1 =  patch(xDATA_L, yDATALow,'b');
Low1.FaceAlpha = 0.3





%% Code from JAT %% 
% you've created and plotted a meanPSD variable
% create a stdPSD variable

% make sure these are row vectors
upperPSD = meanPSD + stdPSD
lowerPSD = meanPSD - stdPSD

% assuming you are fxx as the x-axis
% assuming that fxx is a row vector (not a column vector)
xDATA_L = [fxx fliplr(fxx)];
yDATA = [lowerPSD fliplr(upperPSD)];

Neu1 = patch(xDATA_L, yDATA,'r')
Neu1.FaceAlpha = 0.3


% Old Plot 
plot(subTab_AMY_Neutral_R.PSD{1}.Frequency, mean(AMYNeutral, 2), 'Color','green')
hold on 
plot(subTab_AMY_Neutral_R.PSD{1}.Frequency, mean(AMYHigh, 2), 'Color','red')
hold on 
plot(subTab_AMY_Neutral_R.PSD{1}.Frequency, mean(AMYLow, 2), 'Color','blue')