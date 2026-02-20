% load data
cd('C:\Users\Admin\Documents\Github\LossAversionManuscript')

dataFile = readtable('estimation_results_20250820.csv');

cmap2use = colormap(plasma);
close
cmapLOCS = round(linspace(1,256,3));
trimCmap = cmap2use(cmapLOCS,:);
trimCmap(3,:) = cmap2use(235,:);

%%


close all
data2use = {'lambda','rho','mu'};
data2useSE = {'lambdaSE','rhoSE','muSE'};

tiledlayout(1,3,"TileSpacing","tight","Padding","compact")

dataTableSort = sortrows(dataFile,"lambda","descend");
lossAverseLINE = find(dataTableSort.lossaverse,1,'last') + 0.5;
gainlossNeutralLINE = find(dataTableSort.gainlossneutral,1,'last') + 0.5;

barALPHA = 0.7;
lineALPHA = 0.6;

for iii = 1:3

    switch iii
        case 1 % lambda

            nexttile
            % Data
            barLL = barh(dataTableSort.(data2use{iii}));
            barLL.FaceColor = trimCmap(iii,:);
            barLL.FaceAlpha = barALPHA;
            barLL.EdgeColor = "none";

            yticks(1:numel(dataTableSort.(data2use{iii})))
            yticklabels(dataTableSort.subjectIDs)

            box off

            xticks(0:0.5:5)
            xlim([0 5])
            xlabel('lambda | Loss aversion')
            ylabel('Subject Number')

            % Data SE
            yLINEs = [1:17;1:17];
            xLINEs = [transpose(dataTableSort.(data2use{iii}) - dataTableSort.(data2useSE{iii})) ;...
                transpose(dataTableSort.(data2use{iii}))];

            hold on
            line(xLINEs,yLINEs,'Color','k','LineWidth',1)

            % Lambda line
            xline(1,'--','')

            % Lambda mean
            meanLam = xline(mean(dataTableSort.(data2use{iii})),'-','mean');
            meanLam.Color = trimCmap(iii,:);
            meanLam.Alpha = lineALPHA;
            meanLam.LineWidth = 4;
            meanLam.LabelVerticalAlignment = "bottom";
            meanLam.LabelOrientation = "horizontal";

            % Loss lines
            yline(lossAverseLINE,'k-')
            yline(gainlossNeutralLINE,'k-')

            title('')
            st = subtitle('Loss aversion');
            st.Units = 'normalized';
            st.Position(1) = 0;     
            st.HorizontalAlignment = 'left';

            text(3.5,15.5,'Loss averse')
            text(3.5,11,'Gain loss neutral')
            text(3.5,4.5,'Gain seeking')


        case 2 % mu

            nexttile
            % Data
            barLL = barh(dataTableSort.(data2use{iii}));
            barLL.FaceColor = trimCmap(iii,:);
            barLL.FaceAlpha = barALPHA;
            barLL.EdgeColor = "none";

            yticks(1:numel(dataTableSort.(data2use{iii})))
            yticklabels(dataTableSort.subjectIDs)

            box off

            xticks(0:0.5:2.5)
            xlim([0 2.5])
            xlabel('rho | Risk attitude')
            ylabel('Subject Number')

            % Data SE
            yLINEs = [1:17;1:17];
            xLINEs = [transpose(dataTableSort.(data2use{iii}) - dataTableSort.(data2useSE{iii})) ;...
                transpose(dataTableSort.(data2use{iii}))];

            hold on
            line(xLINEs,yLINEs,'Color','k','LineWidth',1)

            % Lambda line
            xline(1,'--','')

            % Lambda mean
            meanLam = xline(mean(dataTableSort.(data2use{iii})),'-','mean');
            meanLam.Color = trimCmap(iii,:);
            meanLam.Alpha = lineALPHA;
            meanLam.LineWidth = 4;
            meanLam.LabelVerticalAlignment = "bottom";
            meanLam.LabelOrientation = "horizontal";

            % Loss lines
            yline(lossAverseLINE,'k-')
            yline(gainlossNeutralLINE,'k-')

            title('')
            st = subtitle('Risk attitude');
            st.Units = 'normalized';
            st.Position(1) = 0;     
            st.HorizontalAlignment = 'left';

        case 3 % cc

            nexttile
            % Data
            barLL = barh(dataTableSort.(data2use{iii}));
            barLL.FaceColor = trimCmap(iii,:);
            barLL.FaceAlpha = barALPHA;
            barLL.EdgeColor = "none";

            yticks(1:numel(dataTableSort.(data2use{iii})))
            yticklabels(dataTableSort.subjectIDs)

            box off

            xticks(0:20:100)
            xlim([0 100])
            xlabel('mu | Choice consistency')
            ylabel('Subject Number')

            % Data SE
            yLINEs = [1:17;1:17];
            xLINEs = [transpose(dataTableSort.(data2use{iii}) - dataTableSort.(data2useSE{iii})) ;...
                transpose(dataTableSort.(data2use{iii}))];

            hold on
            line(xLINEs,yLINEs,'Color','k','LineWidth',1)

            %  mean
            meanLam = xline(mean(dataTableSort.(data2use{iii})),'-','mean');
            meanLam.Color = trimCmap(iii,:);
            meanLam.Alpha = lineALPHA;
            meanLam.LineWidth = 4;
            meanLam.LabelVerticalAlignment = "bottom";
            meanLam.LabelOrientation = "horizontal";

            % Loss lines
            yline(lossAverseLINE,'k-')
            yline(gainlossNeutralLINE,'k-')

            title('')
            st = subtitle('Choice consistency');
            st.Units = 'normalized';
            st.Position(1) = 0;     
            st.HorizontalAlignment = 'left';

    end


end



