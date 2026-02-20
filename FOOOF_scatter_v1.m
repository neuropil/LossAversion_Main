
partID = 'CLASE007';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';
        % saveCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';
        % saveCD = '';

end % switch case


% CD to participant folder
partPath = strcat(partCD, partID);
cd(partPath)

% do stuff for the length of files in the folder

% Get directory of names you want
mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, '.mat');
dirEpochNames = string(dirNames(dirFilter));

for i = 1:length(dirEpochNames)

    % Load temporary file
    tempName = dirEpochNames(i);
    load(tempName)


    % ---- Pull out epochs for all trial types ---- %

    % Find the rows where the EpochID column is "Outcome"
    rowsWithChoice = strcmp(FOOOF_tab.EpochID, "choice");
    rowsWithResponse = strcmp(FOOOF_tab.EpochID, "response");
    rowsWithOutcome = strcmp(FOOOF_tab.EpochID, "outcome");

    % Extract those rows from the table and create a new table
    tempChoice = FOOOF_tab(rowsWithChoice, :); % choice epoch
    tempResp = FOOOF_tab(rowsWithResponse, :); % Response epoch
    tempOut = FOOOF_tab(rowsWithOutcome, :); % tmpOut is now the table that has all the outcome epochs in it

    % Extract the rows from the table that are only Loss Aversion
    tempChoice = tempChoice(tempChoice.LA,:);
    tempResp = tempResp(tempResp.LA,:);
    tempOut = tempOut(tempOut.LA,:);

    % Remove the empty cells from tempResp
    lastColumn = tempResp{:, end}; % Extract the last column
    % Identify rows where the last column is empty
    emptyRows = cellfun(@isempty, lastColumn);
    % Remove those rows
    tempResp(emptyRows, :) = [];

    % Create GG/GL/AN outcome tables
    tempOutGG = tempOut(tempOut.OutcomeGain,:);
    tempOutGL = tempOut(tempOut.OutcomeLoss,:);
    tempOutAN = tempOut(tempOut.OutcomeNeutral, :);

    % Cells to hold the power and frequency to plot later
    choicePower = [];
    choiceFreq = [];
    choiceOffset = [];
    choiceExp = [];

    respPower = [];
    respFreq = [];
    respOffset = [];
    respExp = [];

    outGGPower = [];
    outGGFreq = [];
    outGGOffset = [];
    outGGExp = [];

    outGLPower = [];
    outGLFreq = [];
    outGLOffset = [];
    outGLExp = [];

    outANPower = [];
    outANFreq = [];
    outANOffset = [];
    outANExp = [];

    % ---- pull out data ---- %

    % Choice
    for ci = 1:height(tempChoice)
        tempChoicePower = tempChoice.FOOOFoutput{ci,1}.peak_params(:,2);
        tempChoiceFreq = tempChoice.FOOOFoutput{ci,1}.peak_params(:,1);
        tempChoiceOffset = tempChoice.FOOOFoutput{ci,1}.aperiodic_params(:,1);
        tempChoiceExp = tempChoice.FOOOFoutput{ci,1}.aperiodic_params(:,2);

        choicePower = [choicePower; tempChoicePower];
        choiceFreq = [choiceFreq; tempChoiceFreq];
        choiceOffset = [choiceOffset; tempChoiceOffset];
        choiceExp = [choiceExp; tempChoiceExp];

    end % for / ci

    % Response
    for ri = 1:height(tempResp)
        tempRespPower = tempResp.FOOOFoutput{ri,1}.peak_params(:,2);
        tempRespFreq = tempResp.FOOOFoutput{ri,1}.peak_params(:,1);
        tempRespOffset = tempResp.FOOOFoutput{ri,1}.aperiodic_params(:,1);
        tempRespExp = tempResp.FOOOFoutput{ri,1}.aperiodic_params(:,2);

        respPower = [respPower; tempRespPower];
        respFreq = [respFreq; tempRespFreq];
        respOffset = [respOffset; tempRespOffset];
        respExp = [tempRespExp; tempRespExp];

    end % for / ri

    % Outcome GG
    for gg = 1:height(tempOutGG)
        tempOutGGPower = tempOutGG.FOOOFoutput{gg,1}.peak_params(:,2);
        tempOutGGFreq = tempOutGG.FOOOFoutput{gg,1}.peak_params(:,1);
        tempOutGGOffset = tempOutGG.FOOOFoutput{gg,1}.aperiodic_params(:,1);
        tempOutGGExp = tempOutGG.FOOOFoutput{gg,1}.aperiodic_params(:,2);

        outGGPower = [outGGPower; tempOutGGPower];
        outGGFreq = [outGGFreq; tempOutGGFreq];
        outGGOffset = [outGGOffset; tempOutGGOffset];
        outGGExp = [outGGExp; tempOutGGExp];

    end % for / gg

    % Outcome GL
    for gl = 1:height(tempOutGL)
        tempOutGLPower = tempOutGL.FOOOFoutput{gl,1}.peak_params(:,2);
        tempOutGLFreq = tempOutGL.FOOOFoutput{gl,1}.peak_params(:,1);
        tempOutGLOffset = tempOutGL.FOOOFoutput{gl,1}.aperiodic_params(:,1);
        tempOutGLExp = tempOutGL.FOOOFoutput{gl,1}.aperiodic_params(:,2);

        outGLPower = [outGLPower; tempOutGLPower];
        outGLFreq = [outGLFreq; tempOutGLFreq];
        outGLOffset = [outGLOffset; tempOutGLOffset];
        outGLExp = [outGLExp; tempOutGLExp];

    end % for / gl

    % Outcome AN
    for an = 1:height(tempOutAN)
        tempOutANPower = tempOutAN.FOOOFoutput{an,1}.peak_params(:,2);
        tempOutANFreq = tempOutAN.FOOOFoutput{an,1}.peak_params(:,1);
        tempOutANOffset = tempOutAN.FOOOFoutput{an,1}.aperiodic_params(:,1);
        tempOutANExp = tempOutAN.FOOOFoutput{an,1}.aperiodic_params(:,2);

        outANPower = [outANPower; tempOutANPower];
        outANFreq = [outANFreq; tempOutANFreq];
        outANOffset = [outANOffset; tempOutANOffset];
        outANExp = [outANExp; tempOutANExp];

    end % for / an


    % Plotting % 

    % Temporary plot name 
    tempNamePlot = erase(tempName, "_FOOOF.mat");

    tiledlayout(1,2)
    nexttile
    scatter(outGGFreq,outGGPower, 'green', 'filled')
    hold on
    scatter(outGLFreq,outGLPower,'red','filled')
    hold on
    scatter(outANFreq,outANPower, 'cyan', 'filled')
    hold on 
    title([tempNamePlot + "-outcome"])
    

    % outcome all same color
    nexttile
    scatter(choiceFreq,choicePower, 'blue','filled')
    hold on
    scatter(respFreq,respPower, 'magenta','filled')
    hold on
    scatter(outGGFreq,outGGPower, 'green', 'filled')
    hold on
    scatter(outGLFreq,outGLPower,'green','filled')
    hold on
    scatter(outANFreq,outANPower, 'green', 'filled')
    hold on 
    title(tempNamePlot)


    % ----- Aperiodic ----- %

    % Create Exponent labels
    ChoiceLabel = repmat({'Choice'}, size(choiceExp));
    ResponseLabel = repmat({'Response'}, size(respExp));
    OutGGLabel = repmat({'OutGG'}, size(outGGExp));
    OutGLLabel = repmat({'OutGL'}, size(outGLExp));
    OutANlabel = repmat({'OutAN'}, size(outANExp));

    % Exponent
    figure;
    tiledlayout(1,2)
    nexttile
    boxchart(categorical([ChoiceLabel; ResponseLabel; OutGGLabel; OutGLLabel; OutANlabel]), ...
        [choiceExp; respExp; outGGExp; outGLExp; outANExp])
    hold on 
    title(tempNamePlot + "-Exponent")

    % Create Offset labels
    ChoiceLabelOffset = repmat({'Choice'}, size(choiceOffset));
    RespLabelOffset = repmat({'Response'}, size(respOffset));
    OutGGLabelOffset = repmat({'OutGG'}, size(outGGOffset));
    OutGLLabelOffset = repmat({'OutGL'}, size(outGLOffset));
    OutANlabelOffset = repmat({'OutAN'}, size(outANOffset));

    % Offset
    nexttile
    boxchart(categorical([ChoiceLabelOffset; RespLabelOffset; OutGGLabelOffset; OutGLLabelOffset; OutANlabelOffset]), ...
        [choiceOffset; respOffset; outGGOffset; outGLOffset; outANOffset])
    hold on 
    title(tempNamePlot + "-Offset")


    % Save 
    savePath = strcat(partPath + "\Figures");
    cd(savePath)

    saveas(figure(1), [tempNamePlot + ".jpg"]);
    saveas(figure(2), [tempNamePlot + "-aperiodic.jpg"]);
        
    close all 
    
    cd(partPath)


end % for / i