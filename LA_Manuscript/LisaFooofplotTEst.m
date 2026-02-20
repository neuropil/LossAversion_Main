%%

for i = 1:3
    hold on
    switch i
        case 1
            plot(fooof_results.freqs,fooof_results.ap_fit,'k-')
        case 2
            plot(fooof_results.freqs,fooof_results.power_spectrum,'g-')
        case 3
            plot(fooof_results.freqs,fooof_results.fooofed_spectrum,'r-')
    end
end
legend('ap_fit','psd','fooof')

%%

% load in file
cd('Y:\LossAversion\LH_Data\JAT_Data\CLASE001')

load('CLASE001_L_AH_FOOOF.mat')


% color2use = [1 0 0 0.6;   % choice   (red)
%              0 1 0 0.6;   % response (green)
%              0 0 0 0.6];  % outcome  (black)

color2useNA = brewermap(3,'Dark2');
color2useWA = [color2useNA , repmat(0.3,3,1)];
color2useWAmu = [color2useNA , repmat(0.8,3,1)];

for ii = 1:height(FOOOF_tab)

    hold on
    % Use stock Frequency x axis
    x_axisP = FOOOF_tab.FOOOFoutput{1,:}.freqs;

    if isempty(FOOOF_tab.FOOOFoutput{ii,:})
        continue
    else
        y_axisP = FOOOF_tab.FOOOFoutput{ii,:}.fooofed_spectrum;

        switch FOOOF_tab.EpochID{ii}
            case 'choice'
                uCC = 1;
            case 'response'
                uCC = 2;
            case 'outcome'
                uCC = 3;
        end
        plot(x_axisP , y_axisP,'Color',color2useWA(uCC,:))
    end

end

legend({'choice','response','outcome'})

%%

laONLY_tab = FOOOF_tab(FOOOF_tab.LA,:);

x_axisP = laONLY_tab.FOOOFoutput{1,:}.freqs;

choiceMAT = nan(sum(matches(laONLY_tab.EpochID,'choice')),39);
responseMAT = nan(sum(matches(laONLY_tab.EpochID,'response')),39);
outcomeMAT = nan(sum(matches(laONLY_tab.EpochID,'outcome')),39);

choiceI = 1;
responseI = 1;
outcomeI = 1;

for ii = 1:height(laONLY_tab)
    
    if isempty(laONLY_tab.FOOOFoutput{ii,:})
        continue
    else
        tmpPow = laONLY_tab.FOOOFoutput{ii,:}.fooofed_spectrum;

        switch laONLY_tab.EpochID{ii}
            case 'choice'
                choiceMAT(choiceI,:) = tmpPow;
                choiceI = choiceI + 1;
            case 'response'
                responseMAT(responseI,:) = tmpPow;
                responseI = responseI + 1;                
            case 'outcome'
                outcomeMAT(outcomeI,:) = tmpPow;
                outcomeI = outcomeI + 1;
        end
    end
end

% Calculate the mean power for each condition across trials
choiceCL = choiceMAT(any(~isnan(choiceMAT),2),:);
responseCL = responseMAT(any(~isnan(responseMAT),2),:);
outcomeCL = outcomeMAT(any(~isnan(outcomeMAT),2),:);

meanChoice = mean(choiceCL, 1);
meanResponse = mean(responseCL, 1);
meanOutcome = mean(outcomeCL, 1);

% Plot the mean power for each condition
figure;
hold on;
% plot(x_axisP, choiceCL, 'Color', color2useWA(1,:), 'LineWidth', 1.5);
% plot(x_axisP, responseCL, 'Color', color2useWA(2,:), 'LineWidth', 1.5);
% plot(x_axisP, outcomeCL, 'Color', color2useWA(3,:), 'LineWidth', 1.5);

plot(x_axisP, meanChoice, 'Color', color2useWAmu(1,:), 'LineWidth', 4);
plot(x_axisP, meanResponse, 'Color', color2useWAmu(2,:), 'LineWidth', 4);
plot(x_axisP, meanOutcome, 'Color', color2useWAmu(3,:), 'LineWidth', 4);

xlabel('Frequency (Hz)');
ylabel('Mean Power');
title('Mean Power Spectrum by Condition');
legend({'Choice', 'Response', 'Outcome'});
hold off;



%%

color2useNA = brewermap(3,'Dark2');
color2useWA = [color2useNA , repmat(0.3,3,1)];
color2useWAmu = [color2useNA , repmat(0.8,3,1)];

for ii = 1:numel(epMap)

    hold on
    % Use stock Frequency x axis
    x_axisP = newFooof{1}.freqs;

    if isempty(newFooof{ii})
        continue
    else
        y_axisP = newFooof{ii}.fooofed_spectrum;

        switch epMap{ii}
            case 'choice'
                uCC = 1;
            case 'response'
                uCC = 2;
            case 'outcome'
                uCC = 3;
        end
        plot(x_axisP , y_axisP,'Color',color2useWA(uCC,:))
    end

end

legend({'choice','response','outcome'})




%%

x_axisP = newFooof{1}.freqs;

choiceMAT = nan(sum(matches(epMap,'choice')),793);
responseMAT = nan(sum(matches(epMap,'response')),793);
outcomeMAT = nan(sum(matches(epMap,'outcome')),793);

choiceI = 1;
responseI = 1;
outcomeI = 1;

for ii = 1:numel(epMap)
    
    if isempty(newFooof{ii})
        continue
    else
        tmpPow = newFooof{ii}.fooofed_spectrum;

        switch epMap{ii}
            case 'choice'
                choiceMAT(choiceI,:) = tmpPow;
                choiceI = choiceI + 1;
            case 'response'
                responseMAT(responseI,:) = tmpPow;
                responseI = responseI + 1;                
            case 'outcome'
                outcomeMAT(outcomeI,:) = tmpPow;
                outcomeI = outcomeI + 1;
        end
    end
end

% Calculate the mean power for each condition across trials
choiceCL = choiceMAT(any(~isnan(choiceMAT),2),:);
responseCL = responseMAT(any(~isnan(responseMAT),2),:);
outcomeCL = outcomeMAT(any(~isnan(outcomeMAT),2),:);

meanChoice = mean(choiceCL, 1);
meanResponse = mean(responseCL, 1);
meanOutcome = mean(outcomeCL, 1);

% Plot the mean power for each condition
figure;
hold on;
% plot(x_axisP, choiceCL, 'Color', color2useWA(1,:), 'LineWidth', 1.5);
% plot(x_axisP, responseCL, 'Color', color2useWA(2,:), 'LineWidth', 1.5);
% plot(x_axisP, outcomeCL, 'Color', color2useWA(3,:), 'LineWidth', 1.5);

plot(x_axisP, meanChoice, 'Color', color2useWAmu(1,:), 'LineWidth', 4);
plot(x_axisP, meanResponse, 'Color', color2useWAmu(2,:), 'LineWidth', 4);
plot(x_axisP, meanOutcome, 'Color', color2useWAmu(3,:), 'LineWidth', 4);

xlabel('Frequency (Hz)');
ylabel('Mean Power');
title('Mean Power Spectrum by Condition');
legend({'Choice', 'Response', 'Outcome'});
hold off;