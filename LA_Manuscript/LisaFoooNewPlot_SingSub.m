function [] = LisaFoooNewPlot_SingSub(Subject,TrialType)

close all

% load in file
cd('Z:\LossAversion\LH_Data\JAT_FoooFtables')

fileName = [Subject,'_FooofTab.mat'];
load(fileName,'allFoooftab')

%%%%% CREATE LIST OF BRAIN REGIONS

brainList = allFoooftab.brainA;
% Convert to string for consistency
brainList = string(brainList);

% Get unique regions (stable preserves order of appearance)
uniqRegions = unique(brainList, 'stable');

% Display menu
fprintf('\nSelect brain region:\n');
for ii = 1:numel(uniqRegions)
    fprintf('  %d) %s\n', ii, uniqRegions(ii));
end
fprintf('  %d) All regions\n', numel(uniqRegions) + 1);

% Get user selection
sel = input('Enter selection number: ');

% Parse selection
if sel == numel(uniqRegions) + 1
    selectedRegions = uniqRegions;     % ALL
else
    selectedRegions = uniqRegions(sel); % single region
end

color2useNA = brewermap(3,'Dark2');
% color2useWA = [color2useNA , repmat(0.3,3,1)];
color2useWAmu = [color2useNA , repmat(0.8,3,1)];

% Convert strings to categoricals
allFoooftab.brainA = categorical(allFoooftab.brainA);
allFoooftab.EpochID = categorical(allFoooftab.EpochID);

% Temporay table with Only LA | brain = L_AH

nGroups = numel(selectedRegions);          % or height(G), etc.

nCols = ceil(sqrt(nGroups));
nRows = ceil(nGroups / nCols);

t = tiledlayout(nRows, nCols, ...
    "TileSpacing","compact", ...
    "Padding","compact");

for brainII = 1:length(selectedRegions)

    tmpSubBR_SP_OLA = allFoooftab(allFoooftab.LA &...
        allFoooftab.brainA == selectedRegions(brainII),:);

    x_axisP = tmpSubBR_SP_OLA.FOOOFoutput{1}.freqs;

    choiceMAT = nan(sum(tmpSubBR_SP_OLA.EpochID == "choice"),793);
    responseMAT = nan(sum(tmpSubBR_SP_OLA.EpochID == "response"),793);
    outcomeMAT = nan(sum(tmpSubBR_SP_OLA.EpochID == "outcome"),793);

    choiceI = 1;
    responseI = 1;
    outcomeI = 1;

    for ii = 1:height(tmpSubBR_SP_OLA)

        epochFOOOF = tmpSubBR_SP_OLA.FOOOFoutput{ii};

        if isempty(epochFOOOF)
            continue
        else
            tmpPow = epochFOOOF.fooofed_spectrum;

            switch tmpSubBR_SP_OLA.EpochID(ii)
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
    nexttile(t, brainII);
    hold on;
    % plot(x_axisP, choiceCL, 'Color', color2useWA(1,:), 'LineWidth', 1.5);
    % plot(x_axisP, responseCL, 'Color', color2useWA(2,:), 'LineWidth', 1.5);
    % plot(x_axisP, outcomeCL, 'Color', color2useWA(3,:), 'LineWidth', 1.5);

    plot(x_axisP, meanChoice, 'Color', color2useWAmu(1,:), 'LineWidth', 4);
    plot(x_axisP, meanResponse, 'Color', color2useWAmu(2,:), 'LineWidth', 4);
    plot(x_axisP, meanOutcome, 'Color', color2useWAmu(3,:), 'LineWidth', 4);

    xlim([3 45])

    xlabel('Frequency (Hz)');
    ylabel('Mean Power');
    title(selectedRegions(brainII),'Interpreter','none');
    legend({'Choice', 'Response', 'Outcome'});
    hold off;


end



end
