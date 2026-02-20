
cd('Y:\LossAversion\LH_Data\FOOOF_data')
load("FOOOF_allBlockTables.mat"); % Loads as allPartTablesBlocks
unBA = unique(allPartTablesBlocks.BrainArea);
aov = struct;

for i = 1:length(unBA)
    figure;
    tempBA = unBA{i};
    BATable = allPartTablesBlocks(matches(allPartTablesBlocks.BrainArea, tempBA),:);

    unFreq = unique(BATable.Frequency);

    tiledlayout(4,1,"TileSpacing","tight","Padding","compact")
    colors = 'bgrc';
    for fi = 1:length(unFreq)
        
        tempFreq = unFreq{fi};
        FreqTab = BATable(matches(BATable.Frequency, tempFreq),:);
        
        nexttile

        for bi = 1:5
          blockMean = mean(FreqTab.Power(ismember(FreqTab.BlockNum, bi)));
          blockSTD = std(FreqTab.Power(ismember(FreqTab.BlockNum, bi)));
          
          scatter(bi,blockMean, 60, colors(fi), "filled")
          hold on 
          line([bi bi], [blockMean-blockSTD blockMean+blockSTD], 'color', 'k' )
        end % for / bi 
        xlim([0.5 5.5])
        xticks(1:5)
        xlabel('Blocks')
        title(tempFreq)

        ylim([0 1.0])
        ylabel('Power')

        % aov = anova(BATable, "Power ~ BlockNum + Frequency + BlockNum:Frequency");

    end % for / fi 
sgtitle(tempBA)

aov.(tempBA) = anova(BATable, "Power ~ BlockNum + Frequency + BlockNum:Frequency"); % create a struct for the ANOVA 
end % for / i 


%%
% Posthoc from anova 
m = multcompare(aov, ["BlockNum", "Frequency"]);
mAH = multcompare(aov.AH, ["BlockNum", "Frequency"]);
mAMY = multcompare(aov.AMY, ["BlockNum", "Frequency"]);
mLOF = multcompare(aov.LOF, ["BlockNum", "Frequency"]);
mMOF = multcompare(aov.MOF, ["BlockNum", "Frequency"]);
mPH = multcompare(aov.PH, ["BlockNum", "Frequency"]);

% Get only the significant post hoc values 
np = m(m.pValue <= 0.05, :);
sigmAH = mAH(mAH.pValue <= 0.05, :);
sigmAMY = mAMY(mAMY.pValue <= 0.05, :);
sigmLOF = mLOF(mLOF.pValue <= 0.05, :);
sigmMOF = mMOF(mMOF.pValue <= 0.05, :);
sigmPH = mPH(mPH.pValue <= 0.05, :);