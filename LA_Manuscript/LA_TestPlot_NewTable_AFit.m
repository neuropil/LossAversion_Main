    close all
for ii = 1:3

    switch ii
        case 1
            table2use = FOallTABLE;
        case 2
            table2use = FNallTABLE;
        case 3
            table2use = KNallTABLE;
    end


    conds = {'Z_Exp','Z_Off'};
    figure
    tiledlayout(2,1)
    for i = 1:2
        nexttile

        table2use.EpochID = categorical(table2use.EpochID);
        table2use.BArea = categorical(table2use.BArea);

        BA_COND1 = table2use.EpochID == 'choice' & table2use.BArea == 'PH';
        s1 = scatter(ones(sum(BA_COND1),1), table2use.(conds{i})(BA_COND1),30,'k','filled');
        s1.XJitter = "rand";
        s1.XJitterWidth = 0.2;

        hold on
        BA_COND2 = table2use.EpochID == 'response' & table2use.BArea == 'PH';
        s2 = scatter(ones(sum(BA_COND2),1)*2, table2use.(conds{i})(BA_COND2),30,'r','filled');

        s2.XJitter = "rand";
        s2.XJitterWidth = 0.2;

        BA_COND3 = table2use.EpochID == 'outcome' & table2use.BArea == 'PH';
        s3 = scatter(ones(sum(BA_COND3),1)*3, table2use.(conds{i})(BA_COND3),30,'g','filled');

        s3.XJitter = "rand";
        s3.XJitterWidth = 0.2;

        scatter(1,mean(table2use.(conds{i})(BA_COND1)),70,'red','filled')
        scatter(2,mean(table2use.(conds{i})(BA_COND2)),70,'black','filled')
        scatter(3,mean(table2use.(conds{i})(BA_COND3)),70,'black','filled')

    end
end