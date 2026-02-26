function [] = LAM_AperiodicParmCompare(inputU, foooftype)

% Compare Offset and Exponent 
% Brain Region | LA type | Epoch

% Define brain regions, LA types, and epochs for comparison

% JAT knee

cd('Z:\LossAversion\LH_Data\JAT_FoooFtablesZS')
load('CLASE007_FooofTab_ZS.mat','allFoooftab')

allFoooftab.EpochID = categorical(allFoooftab.EpochID);
% allFoooftab.EpochID == 'choice
allFoooftab.brainA = categorical(allFoooftab.brainA);
% allFoooftab.brainA == 'L_AH';
lahtab = allFoooftab(allFoooftab.brainA == 'L_AMY' &...
    allFoooftab.LA & allFoooftab.EpochID ~= 'response',:);

lah_data_exp = cell(1,3);
lah_data_off = cell(1,3);
chC = 1;
reC = 1;
ouC = 1;

switch foooftype
    case 'ff'

    case 'sf'

        useC = 'FIXED';

    case 'sk'

        useC = 'KNEE';

end

for cii = 1:height(lahtab)

    tmpRow = lahtab.EpochID(cii);
    switch tmpRow
        case 'choice'
            lah_data_exp{1,1}(chC) = lahtab.SpecPARAMJ{cii}.(useC).AfitParms.Z_Exp;
            lah_data_off{1,1}(chC) = lahtab.SpecPARAMJ{cii}.(useC).AfitParms.Z_Off;
            chC = chC + 1;
        case 'response'
            lah_data_exp{1,2}(reC) = lahtab.SpecPARAMJ{cii}.(useC).AfitParms.Z_Exp;
            lah_data_off{1,2}(reC) = lahtab.SpecPARAMJ{cii}.(useC).AfitParms.Z_Off;
            reC = reC + 1;
        case 'outcome'
            lah_data_exp{1,3}(ouC) = lahtab.SpecPARAMJ{cii}.(useC).AfitParms.Z_Exp;
            lah_data_off{1,3}(ouC) = lahtab.SpecPARAMJ{cii}.(useC).AfitParms.Z_Off;
            ouC = ouC + 1;
    end

    lahtab.EXP(cii) = lahtab.SpecPARAMJ{cii}.(useC).AfitParms.Z_Exp;
    lahtab.OFF(cii) = lahtab.SpecPARAMJ{cii}.(useC).AfitParms.Z_Off;
end

lahtab = lahtab(~(lahtab.EXP > 10 | lahtab.EXP < -10),:);
lahtab = lahtab(~(lahtab.OFF > 10 | lahtab.OFF < -10),:);


% if matches(inputU,'OFF')
%     makePlot(lah_data_off , 'OFF')
% else
%     makePlot(lah_data_exp , 'EXP')
% end
% 
% lah_data_exp2 = lah_data_exp(1,3);

lahtab = removevars(lahtab , ["LA","OutcomeGain","OutcomeLoss",...
    "OutcomeNeutral","FOOOFoutput","RAWlfp","SpecPARAMJ",...
    "TrialID","brainA","SubID"]);

useIdx = ismember(lahtab.EpochID,{'choice','outcome'});
lahtab2 = lahtab(useIdx,:);

if matches(inputU,'OFF')
    makePlot(lahtab , 'OFF')
else
    makePlot(lahtab , 'EXP')
end

% blockMeans = groupsummary(lahtab2,{'BlockID','EpochID'},'mean','EXP');
% lahtab3 = removevars(lahtab2,"BlockID")
% blockMeans = groupsummary(lahtab2, {'EpochID'}, 'mean', 'EXP');
% blockMeans2 = removevars(blockMeans,"GroupCount");
% wide = unstack(blockMeans2,'mean_EXP','EpochID');
% 
% pairedDiff = wide.outcome - wide.choice;
% 
% obsDiff = mean(pairedDiff);
% 
% bootDiff = bootstrp(5000,@mean,pairedDiff);
% 
% ci = prctile(bootDiff,[2.5 97.5]);
% p = mean(abs(bootDiff) >= abs(obsDiff));

choiceVals  = lahtab2.EXP(lahtab2.EpochID =='choice');
outcomeVals = lahtab2.EXP(lahtab2.EpochID == 'outcome');

nBoot = 5000;
bootDiff = zeros(nBoot,1);

for i = 1:nBoot
    bootChoice  = choiceVals(randi(numel(choiceVals),numel(choiceVals),1));
    bootOutcome = outcomeVals(randi(numel(outcomeVals),numel(outcomeVals),1));
    
    bootDiff(i) = mean(bootOutcome) - mean(bootChoice);
end

obsDiff = mean(outcomeVals) - mean(choiceVals);

ci = prctile(bootDiff,[2.5 97.5])
p  = mean(abs(bootDiff) >= abs(obsDiff))


end







function [] = makePlot(lahtab2 , titleIN)



virCMAP = colormap(viridis(3));
close all

for piil = 1:2

    if piil == 1
        yData = lahtab2.(titleIN)(lahtab2.EpochID == 'choice');
    else
        yData = lahtab2.(titleIN)(lahtab2.EpochID == 'outcome');
    end

    xData = ones(numel(yData),1)*piil;
    % yData = choice;

    hold on
    s1 = scatter(xData,yData,40,virCMAP(piil,:),'filled');
    s1.XJitter = "rand";
    s1.XJitterWidth = 0.25;

    tmpMedian = median(yData);

    % Plot the median line for each category
    line([piil - 0.2, piil + 0.2], [tmpMedian, tmpMedian], 'Color','k', 'LineWidth', 2);
    text(piil + 0.25,tmpMedian,num2str(round(tmpMedian,2)))

end

xlim([0.7 2.3])
xticks(1:2)
xticklabels({'Choice','Outcome'})

if matches(titleIN , 'OFF')
    ylabel('Z-scored Aperiodic offset [all trials] L-AH Clase001','Interpreter','none')
else
    ylabel('Z-scored Aperiodic exponent [all trials] L-AH Clase001','Interpreter','none')
end

% aovDATA = transpose([expORoff{1}, expORoff{2}, expORoff{3}]);
% grpDATA = transpose([repmat({'Choice'},1,numel(expORoff{1})), repmat({'Response'},1,numel(expORoff{2})),...
%     repmat({'Outcome'},1,numel(expORoff{3}))]);
% grpDATA = categorical(grpDATA);

% Perform Kruskal-Wallis test for the three conditions
aov = anova(lahtab2.EpochID,lahtab2.(titleIN));
multcompare(aov);
aovSTATS = stats(aov);
pVALUE = round(aovSTATS.pValue(1),3);
title(['p-value ',num2str(pVALUE)])



end