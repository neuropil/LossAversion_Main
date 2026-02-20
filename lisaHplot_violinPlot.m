% Load file
load('PerodicByEpoch.mat')

%%

% allPartTablesFreq.Properties.VariableNames;

lowLAtab = allPartTablesFreq(matches(allPartTablesFreq.LAval,'Low'),:); % Low table
hihLAtab = allPartTablesFreq(matches(allPartTablesFreq.LAval,'High'),:); % High table 
%%
% BrainRegion
brainREG = 'OF';
% Frequency band
freq2use = 'theta';
% Hemi 
Hemi = 'R';


allX = [];         % Store x-tick positions
allLabels = {};    % Store corresponding labels


subCOUNT = 1;
for iila = 1:2

    switch iila 
        case 1
            useLAtab = lowLAtab;
            % offset = 0; % First group starts at 0
            baseX = 0; % Start lowLAtab at x = 1, 2, 3, ...
            plotColor = "blue";
        case 2
            useLAtab = hihLAtab;
            % offset = 0.5; % Second group offset by 0.5 to separate plots
            baseX = 4.5; % Start hihLAtab far to the right (e.g., x = 101, 102, ...)
            plotColor = "red";
    end


    uniSUBS = unique(useLAtab.partID); % Subject ID's 

    for iisu = 1:length(uniSUBS)

        tmpSUBtab = useLAtab(matches(useLAtab.partID,uniSUBS{iisu}),:); % Temporary subject table

        if ismember(brainREG, {'AMY'})

            brainBANDtab = tmpSUBtab(matches(tmpSUBtab.BrainArea, brainREG) &...
                matches(tmpSUBtab.FreqLabel, freq2use) & ...
                matches(tmpSUBtab.Hemi, Hemi),:);

        elseif ismember(brainREG, {'OF'})

            brainBANDtab = tmpSUBtab(ismember(tmpSUBtab.BrainArea, {'MOF', 'LOF'}), :);

        else
            brainBANDtab = tmpSUBtab(matches(tmpSUBtab.BrainArea, brainREG) &...
                matches(tmpSUBtab.FreqLabel, freq2use),:); % Brain area and freq band for one participant
        end

        % Skip plotting if brainBANDtab is empty
        % if isempty(brainBANDtab)
        %     continue
        % end


        hold on
        % xDATA = ones(height(brainBANDtab),1)*subCOUNT;
        % xDATA = ones(height(brainBANDtab),1)*iisu*iila;
        % xDATA = ones(height(brainBANDtab),1)*(iisu + offset); % Offset x-position
        % xDATA = ones(height(brainBANDtab),1)*(baseX + iisu); % Separate x-position for each group
        xPOS = baseX + iisu; % Unique x-position
        xDATA = ones(height(brainBANDtab),1) * xPOS;
        yDATA = brainBANDtab.Power;

        % figure
        % swP = swarmchart(xDATA,yDATA,50,"red","filled");
        swP = swarmchart(xDATA, yDATA, 50, plotColor, "filled");
        swP.MarkerFaceAlpha = 0.4;
        swP.MarkerEdgeColor = "none";
        swP.XJitter = "rand";
        swP.XJitterWidth = 0.2;

        hold on 

        vnP = violinplot(xDATA, yDATA, Orientation= "vertical", DensityDirection= "positive");
        vnP.EdgeColor = 'none';
        vnP.FaceColor = plotColor; % Set violin plot color

        hold on


        % Store x-tick position and label
        allX(end+1) = xPOS;
        allLabels{end+1} = uniSUBS{iisu};


    end

  
end



% Set x-axis ticks and labels after plotting
xticks(allX)
xticklabels(allLabels)
xtickangle(45) % Optional: rotate labels for readability

ylim([-0.05 2.5])
ylabel('Power')

title(strcat(brainREG, ' | ', freq2use))