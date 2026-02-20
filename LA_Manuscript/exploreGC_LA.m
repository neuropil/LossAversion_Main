function [datapvalTAB2] = exploreGC_LA(subjectID, compareC)

getPCname = getenv('COMPUTERNAME');

switch getPCname
    case 'DESKTOP-I5CPDO7'  % Work PC

        cd('D:\Dropbox\LisaLossGC')

end

close all
% load("CLASE007_L_AMY_L_AH_Choice.mat","gc_Choice_dataT")
% load("CLASE007_L_AMY_L_AH_Outcome.mat","gc_Outcome_dataT")
dirLIST = dir('*.mat');
dirTAB = struct2table(dirLIST);
dirNames = string(dirTAB.name);
nameParts = split(dirNames,'_');
subjectSall = nameParts(:,1);
subjectLIST = dirNames(subjectSall == string(subjectID));

fileNameTab = getFileNameTab(subjectLIST);

switch compareC
    case 'hemi' % Loop through condition and compare L v R

        condU = unique(fileNameTab.Condition);

        histoComp = cell(3,2);
        datapvalTAB = zeros(3,4);
        for cii = 1:length(condU)

            leftCond = fileNameTab.FileName(fileNameTab.Hemi == "L" & fileNameTab.Condition == condU(cii));
            load(leftCond,['gc_' , convertStringsToChars(condU(cii)) , '_dataT'])
            gc_LeftData = eval(['gc_' , convertStringsToChars(condU(cii)), '_dataT']);

            rightCond = fileNameTab.FileName(fileNameTab.Hemi == "R" & fileNameTab.Condition == condU(cii));
            load(rightCond,['gc_' , convertStringsToChars(condU(cii)) , '_dataT'])
            gc_RightData = eval(['gc_' , convertStringsToChars(condU(cii)), '_dataT']);

            [fromAMY_STATS , toAMY_STATS] = histoplotStat_side(gc_LeftData , gc_RightData);

            histoComp{cii,1} = fromAMY_STATS;
            histoComp{cii,2} = toAMY_STATS;

            datapvalTAB(cii,1) = fromAMY_STATS.pval;
            datapvalTAB(cii,2) = fromAMY_STATS.LRmed_diff;
            datapvalTAB(cii,3) = toAMY_STATS.pval;
            datapvalTAB(cii,4) = toAMY_STATS.LRmed_diff;

        end

        datapvalTAB2 = array2table(datapvalTAB,'VariableNames',{'From_AMY_pval',...
            'From_AMY_LRdiff','TO_AMY_pval','To_AMY_LRdiff'});
        datapvalTAB2.Condition = condU;

    case 'condition' % Loop through hemi and compare conditions

        hemiLOOP = ["L","R"];
        for hii = 1:2

            if any(fileNameTab.Hemi == hemiLOOP(hii))

                choiceCond = fileNameTab.FileName(fileNameTab.Hemi == hemiLOOP(hii) & fileNameTab.Condition == "Choice");
                load(choiceCond,'gc_Choice_dataT')
                gc_CHOICE = eval('gc_Choice_dataT');

                responCond = fileNameTab.FileName(fileNameTab.Hemi == hemiLOOP(hii) & fileNameTab.Condition == "Response");
                load(responCond,'gc_Response_dataT')
                gc_RESPONSE = eval('gc_Response_dataT');

                outcomCond = fileNameTab.FileName(fileNameTab.Hemi == hemiLOOP(hii) & fileNameTab.Condition == "Outcome");
                load(outcomCond,'gc_Outcome_dataT')
                gc_OUTCOME = eval('gc_Outcome_dataT');

                [fromAMY_STATS , toAMY_STATS] = histoplotStat_cond(gc_CHOICE , gc_RESPONSE , gc_OUTCOME);


                test = 1;

            else
                continue
            end

        end
    case 'epoch' % Loop through hemi and condition and compare Region TO vs FROM



end

test = 1;
% close all




end






function [fileNameTab] = getFileNameTab(subjectLIST)

hemiA = strings(height(subjectLIST),1);
brain1A = strings(height(subjectLIST),1);
brain2A = strings(height(subjectLIST),1);
conditA = strings(height(subjectLIST),1);

for ii = 1:height(subjectLIST)

    tmpFile = subjectLIST{ii};
    fileParts = split(tmpFile,{'_','.'});

    hemiA(ii) = fileParts(2);
    brain1A(ii) = fileParts(3);
    brain2A(ii) = fileParts(5);
    conditA(ii) = fileParts(6);

end

fileNameTab = table(subjectLIST,hemiA,brain1A,brain2A,conditA,'VariableNames',...
    {'FileName','Hemi','Brain1','Brain2','Condition'});

end



function [FROM_AMY_stats , TO_AMY_stats] = histoplotStat_side(leftCond , rightCond)


tiledlayout(2,1)

nexttile
% leftCond = cond1.Properties.VariableNames;
% rightCond = cond2.Properties.VariableNames;

h1 = histogram(leftCond.("AVG FROM LAMY TO LAH"));
hold on
h2 = histogram(rightCond.("AVG FROM RAMY TO RAH"));

h1.Normalization = "probability";
h1.BinWidth = 0.002;
h2.Normalization = "probability";
h2.BinWidth = 0.002;
h1.EdgeColor = "white";
h2.EdgeColor = "white";
redColor = hex2rgb("#C53B29");
greyColor = [0.4 0.4 0.4];
colororder([greyColor;redColor])
title("Compare FROM AMY TO LAH between L and R : Choice")
subtitle("Between L and R : Choice")
xlabel("GC")
ylabel("Probability")
leg = legend([h1,h2]);
leg.Box = "off";
leg.IconColumnWidth = 9;
leg.String = {'Left', 'Right'};

xlim([0 0.25])

ax = gca;
ax.Box = "off";
ax.LineWidth = 2;
ax.TitleHorizontalAlignment = "left";

% STATS
[~,fromApval,fromAstats] = kstest2(leftCond.("AVG FROM LAMY TO LAH"),...
    rightCond.("AVG FROM RAMY TO RAH"));

FROM_LR_Med_diff = median(leftCond.("AVG FROM LAMY TO LAH"),'omitnan') -...
    median(rightCond.("AVG FROM RAMY TO RAH"),'omitnan');

FROM_AMY_stats.pval = fromApval;
FROM_AMY_stats.stats = fromAstats;
FROM_AMY_stats.LRmed_diff = FROM_LR_Med_diff;

nexttile
h1 = histogram(leftCond.("AVG FROM LAH TO LAMY"));
hold on
h2 = histogram(rightCond.("AVG FROM RAH TO RAMY"));

h1.Normalization = "probability";
h1.BinWidth = 0.002;
h2.Normalization = "probability";
h2.BinWidth = 0.002;
h1.EdgeColor = "white";
h2.EdgeColor = "white";
redColor = hex2rgb("#C53B29");
greyColor = [0.4 0.4 0.4];
colororder([greyColor;redColor])
title("Compare FROM AH TO AMY between L and R : Choice")
subtitle("Between L and R : Choice")
xlabel("GC")
ylabel("Probability")
leg = legend([h1,h2]);
leg.Box = "off";
leg.IconColumnWidth = 9;
leg.String = {'Left', 'Right'};

xlim([0 0.25])

ax = gca;
ax.Box = "off";
ax.LineWidth = 2;
ax.TitleHorizontalAlignment = "left";

TO_LR_Med_diff = median(leftCond.("AVG FROM LAH TO LAMY"),'omitnan') -...
    median(rightCond.("AVG FROM RAH TO RAMY"),'omitnan');

[~,toApval,toAstats] = kstest2(leftCond.("AVG FROM LAH TO LAMY"),...
    rightCond.("AVG FROM RAH TO RAMY"));
TO_AMY_stats.pval = toApval;
TO_AMY_stats.stats = toAstats;
TO_AMY_stats.LRmed_diff = TO_LR_Med_diff;





end









function [FROM_AMY_stats , TO_AMY_stats] = histoplotStat_cond(CHOICE , RESPONSE , OUTCOME)


tiledlayout(2,1)

nexttile

if any(contains(CHOICE.Properties.VariableNames,'AVG FROM LAH'))
    fromAMY = 'AVG FROM LAMY TO LAH';
    toAMY = 'AVG FROM LAH TO LAMY';
    hemiID = 'L';
else
    fromAMY = 'AVG FROM RAMY TO RAH';
    toAMY = 'AVG FROM RAH TO RAMY';
    hemiID = 'R';
end
 

choice_FROM_AMY = histogram(CHOICE.(fromAMY));
hold on
response_FROM_AMY = histogram(RESPONSE.(fromAMY));
outcome_FROM_AMY = histogram(OUTCOME.(fromAMY));

choice_FROM_AMY.Normalization = "probability";
choice_FROM_AMY.BinWidth = 0.005;
response_FROM_AMY.Normalization = "probability";
response_FROM_AMY.BinWidth = 0.005;
outcome_FROM_AMY.Normalization = "probability";
outcome_FROM_AMY.BinWidth = 0.005;

choice_FROM_AMY.EdgeColor = "white";
response_FROM_AMY.EdgeColor = "white";
outcome_FROM_AMY.EdgeColor = "white";

redColor = hex2rgb("#C53B29");
greyColor = [0.4 0.4 0.4];
blueColor = [0, 180, 216]/255;

colororder([greyColor;redColor;blueColor])
title(['Compare FROM AMY between Choice v Response v Outcome | ', hemiID])
subtitle(['Between Epochs : ',hemiID])
xlabel("GC")
ylabel("Probability")
leg = legend([choice_FROM_AMY,response_FROM_AMY,outcome_FROM_AMY]);
leg.Box = "off";
leg.IconColumnWidth = 9;
leg.String = {'Choice', 'Response','Outcome'};

xlim([0 0.15])

ax = gca;
ax.Box = "off";
ax.LineWidth = 2;
ax.TitleHorizontalAlignment = "left";

groupID = [repmat({'C'},135,1); repmat({'R'},135,1);repmat({'O'},135,1)];
dataIN = [CHOICE.(fromAMY) ; RESPONSE.(fromAMY) ; OUTCOME.(fromAMY)];

% STATS
[FROM_AMY_pval,FROM_AMY_tabstst,FROM_AMY_stats] = kruskalwallis(dataIN,groupID,"off");

FROM_AMY_stats.pval = FROM_AMY_pval;
FROM_AMY_stats.stats = FROM_AMY_stats;
FROM_AMY_stats.metatable = FROM_AMY_tabstst;

nexttile
choice_TO_AMY = histogram(CHOICE.(toAMY));
hold on
response_TO_AMY = histogram(RESPONSE.(toAMY));
outcome_TO_AMY = histogram(OUTCOME.(toAMY));

choice_TO_AMY.Normalization = "probability";
choice_TO_AMY.BinWidth = 0.005;
response_TO_AMY.Normalization = "probability";
response_TO_AMY.BinWidth = 0.005;
outcome_TO_AMY.Normalization = "probability";
outcome_TO_AMY.BinWidth = 0.005;

choice_TO_AMY.EdgeColor = "white";
response_TO_AMY.EdgeColor = "white";
outcome_TO_AMY.EdgeColor = "white";

redColor = hex2rgb("#C53B29");
greyColor = [0.4 0.4 0.4];
blueColor = [0, 180, 216]/255;

colororder([greyColor;redColor;blueColor])
title(['Compare TO AMY between Choice v Response v Outcome | ', hemiID])
subtitle(['Between Epochs : ',hemiID])
xlabel("GC")
ylabel("Probability")
leg = legend([choice_TO_AMY,response_TO_AMY,outcome_TO_AMY]);
leg.Box = "off";
leg.IconColumnWidth = 9;
leg.String = {'Choice', 'Response','Outcome'};

ax = gca;
ax.Box = "off";
ax.LineWidth = 2;
ax.TitleHorizontalAlignment = "left";

xlim([0 0.15])

dataIN2 = [CHOICE.(toAMY) ; RESPONSE.(toAMY) ; OUTCOME.(toAMY)];

% STATS
[TO_AMY_pval,TO_AMY_tabstst,TO_AMY_stats] = kruskalwallis(dataIN2,groupID,"off");

TO_AMY_stats.pval = TO_AMY_pval;
TO_AMY_stats.stats = TO_AMY_stats;
TO_AMY_stats.metatable = TO_AMY_tabstst;




end