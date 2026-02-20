
% Create 'Tables' folder in FOOOF participant folder before running script
partID = 'CLASE035';

PCname = getenv('COMPUTERNAME');
switch PCname
    case 'DLPFC'
        partCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';
        % saveCD = 'Z:\LossAversion\LH_Data\FOOOF_data\';

    case 'LATERALHABENULA'
        partCD = 'Y:\LossAversion\LH_Data\FOOOF_data\';
        % saveCD = '';

end % switch case

partPath = strcat(partCD, partID);
cd(partPath)


mdir = dir;
dirNames = {mdir.name};
dirFilter = contains(dirNames, '.mat');
dirEpochNames = string(dirNames(dirFilter));

for i = 1:length(dirEpochNames)

    tempName = dirEpochNames(i);
    load(tempName)

    % Extract the rows from the table that are only Loss Aversion
    FOOOF_tab = FOOOF_tab(FOOOF_tab.LA,:);

    % Remove empty cells
    lastColumn = FOOOF_tab{:, end};
    % Identify rows where the last column is empty
    emptyRows = cellfun(@isempty, lastColumn);
    % Remove those rows
    FOOOF_tab(emptyRows,:) = [];

    % Find the rows where the EpochID column is "Outcome"
    rowsWithChoice = strcmp(FOOOF_tab.EpochID, "choice");
    rowsWithResponse = strcmp(FOOOF_tab.EpochID, "response");
    rowsWithOutcome = strcmp(FOOOF_tab.EpochID, "outcome");

    % Extract those rows from the table and create a new table
    tempChoice = FOOOF_tab(rowsWithChoice, :); % choice epoch
    tempResp = FOOOF_tab(rowsWithResponse, :); % Response epoch
    tempOut = FOOOF_tab(rowsWithOutcome, :); % tmpOut is now the table that has all the outcome epochs in it

    % Create GG/GL/AN outcome tables
    tempOutGG = tempOut(tempOut.OutcomeGain,:);
    tempOutGL = tempOut(tempOut.OutcomeLoss,:);
    tempOutAN = tempOut(tempOut.OutcomeNeutral, :);

    % Cells to hold the power and frequency
    choicePower = [];
    choiceFreq = [];

    respPower = [];
    respFreq = [];

    outGGPower = [];
    outGGFreq = [];

    outGLPower = [];
    outGLFreq = [];

    outANPower = [];
    outANFreq = [];

    % ---- Pull out power and frequency ---- %

    % Choice
    for ci = 1:height(tempChoice)
        tempChoicePower = tempChoice.FOOOFoutput{ci,1}.peak_params(:,2);
        tempChoiceFreq = tempChoice.FOOOFoutput{ci,1}.peak_params(:,1);

        choicePower = [choicePower; tempChoicePower];
        choiceFreq = [choiceFreq; tempChoiceFreq];

    end % for / ci

    % Response
    for ri = 1:height(tempResp)
        tempRespPower = tempResp.FOOOFoutput{ri,1}.peak_params(:,2);
        tempRespFreq = tempResp.FOOOFoutput{ri,1}.peak_params(:,1);

        respPower = [respPower; tempRespPower];
        respFreq = [respFreq; tempRespFreq];

    end % for / ri

    % Outcome GG
    for gg = 1:height(tempOutGG)
        tempOutGGPower = tempOutGG.FOOOFoutput{gg,1}.peak_params(:,2);
        tempOutGGFreq = tempOutGG.FOOOFoutput{gg,1}.peak_params(:,1);

        outGGPower = [outGGPower; tempOutGGPower];
        outGGFreq = [outGGFreq; tempOutGGFreq];

    end % for / gg

    % Outcome GL
    for gl = 1:height(tempOutGL)
        tempOutGLPower = tempOutGL.FOOOFoutput{gl,1}.peak_params(:,2);
        tempOutGLFreq = tempOutGL.FOOOFoutput{gl,1}.peak_params(:,1);

        outGLPower = [outGLPower; tempOutGLPower];
        outGLFreq = [outGLFreq; tempOutGLFreq];

    end % for / gl

    % Outcome AN
    for an = 1:height(tempOutAN)
        tempOutANPower = tempOutAN.FOOOFoutput{an,1}.peak_params(:,2);
        tempOutANFreq = tempOutAN.FOOOFoutput{an,1}.peak_params(:,1);

        outANPower = [outANPower; tempOutANPower];
        outANFreq = [outANFreq; tempOutANFreq];

    end % for / an



    % -- Choice -- %
    choiceLab = string(repmat({'Choice'}, size(choicePower)));
    choiceTemp = [choicePower choiceFreq];

    thetaIDX = (choiceTemp(:,2) >= 1 & choiceTemp(:,2) <= 8);
    alphaIDX = (choiceTemp(:,2) >= 8 & choiceTemp(:,2) <= 12);
    betaIDX = (choiceTemp(:,2) >= 12 & choiceTemp(:,2) <= 30);
    gammaIDX = (choiceTemp(:,2) >= 30);

    ChoiceArray = cell(size(thetaIDX));
    ChoiceArray(thetaIDX) = {'theta'};
    ChoiceArray(alphaIDX) = {'alpha'};
    ChoiceArray(betaIDX) = {'beta'};
    ChoiceArray(gammaIDX) = {'gamma'};

    choiceTab = table(choiceLab, choiceTemp(:,1), ChoiceArray, 'VariableNames',{'Epoch', 'Power', 'Frequency'});

    % -- Response -- %
    responseLab = repmat({'Response'}, size(respFreq));

    if length(respPower) > 1
        responseTemp = [respPower respFreq];

        thetaIDX = (responseTemp(:,2) >= 1 & responseTemp(:,2) <= 8);
        alphaIDX = (responseTemp(:,2) >= 8 & responseTemp(:,2) <= 12);
        betaIDX = (responseTemp(:,2) >= 12 & responseTemp(:,2) <= 30);
        gammaIDX = (responseTemp(:,2) >= 30);

        RespArray = cell(size(thetaIDX));
        RespArray(thetaIDX) = {'theta'};
        RespArray(alphaIDX) = {'alpha'};
        RespArray(betaIDX) = {'beta'};
        RespArray(gammaIDX) = {'gamma'};

        RespTab = table(responseLab, responseTemp(:,1), RespArray, 'VariableNames',{'Epoch', 'Power', 'Frequency'});

    else
        responseTemp = [respPower respFreq];
        RespArray = cell(size(responseTemp));
        RespTab = table(responseLab, responseTemp, RespArray, 'VariableNames',{'Epoch', 'Power', 'Frequency'});
    end % if else

    % -- Outcome GG -- %
    outGGLab = repmat({'OutGG'}, size(outGGPower));
    outGGTemp = [outGGPower outGGFreq];

    thetaIDX = (outGGTemp(:,2) >= 1 & outGGTemp(:,2) <= 8);
    alphaIDX = (outGGTemp(:,2) >= 8 & outGGTemp(:,2) <= 12);
    betaIDX = (outGGTemp(:,2) >= 12 & outGGTemp(:,2) <= 30);
    gammaIDX = (outGGTemp(:,2) >= 30);

    OutGGArray = cell(size(thetaIDX));
    OutGGArray(thetaIDX) = {'theta'};
    OutGGArray(alphaIDX) = {'alpha'};
    OutGGArray(betaIDX) = {'beta'};
    OutGGArray(gammaIDX) = {'gamma'};

    OutGGTab = table(outGGLab, outGGTemp(:,1), OutGGArray, 'VariableNames',{'Epoch', 'Power', 'Frequency'});

    % -- Outcome GL -- %
    outGLLab = repmat({'OutGL'}, size(outGLPower));
    outGLTemp = [outGLPower outGLFreq];

    thetaIDX = (outGLTemp(:,2) >= 1 & outGLTemp(:,2) <= 8);
    alphaIDX = (outGLTemp(:,2) >= 8 & outGLTemp(:,2) <= 12);
    betaIDX = (outGLTemp(:,2) >= 12 & outGLTemp(:,2) <= 30);
    gammaIDX = (outGLTemp(:,2) >= 30);

    OutGLArray = cell(size(thetaIDX));
    OutGLArray(thetaIDX) = {'theta'};
    OutGLArray(alphaIDX) = {'alpha'};
    OutGLArray(betaIDX) = {'beta'};
    OutGLArray(gammaIDX) = {'gamma'};

    OutGLTab = table(outGLLab, outGLTemp(:,1), OutGLArray, 'VariableNames',{'Epoch', 'Power', 'Frequency'});

    % -- Outcome AN -- %
    outANLab = repmat({'OutAN'}, size(outANPower));
    outANTemp = [outANPower outANFreq];

    thetaIDX = (outANTemp(:,2) >= 1 & outANTemp(:,2) <= 8);
    alphaIDX = (outANTemp(:,2) >= 8 & outANTemp(:,2) <= 12);
    betaIDX = (outANTemp(:,2) >= 12 & outANTemp(:,2) <= 30);
    gammaIDX = (outANTemp(:,2) >= 30);

    OutANArray = cell(size(thetaIDX));
    OutANArray(thetaIDX) = {'theta'};
    OutANArray(alphaIDX) = {'alpha'};
    OutANArray(betaIDX) = {'beta'};
    OutANArray(gammaIDX) = {'gamma'};

    OutANTab = table(outANLab, outANTemp(:,1), OutANArray, 'VariableNames',{'Epoch', 'Power', 'Frequency'});


    % -- Combine frequency tables -- %
    freqTabs = vertcat(choiceTab, RespTab, OutGGTab, OutGLTab, OutANTab);

    % Temp name
    nameSplit = split(tempName, '_');
    partID = nameSplit(1);
    Hemi = nameSplit(2);
    BA = nameSplit(3);

    % Repeat identifiers for length fo the table
    partID = repmat(partID, height(freqTabs),1);
    Hemi = repmat(Hemi, height(freqTabs),1);
    BA = repmat(BA, height(freqTabs),1);

    % Create ID table
    IDtab = table(partID, Hemi, BA, 'VariableNames',{'partID', 'Hemi', 'BrainArea'});

    % Combine freqTabs and ID table
    partTab = [IDtab, freqTabs];

    % Save loc
    savePath = strcat(partPath + "\Tables");
    cd(savePath) 

    % save table 
    save(tempName, "partTab");

    % CD back to fooof data 
    cd(partPath)

end % for / i

%%
% partID = 'CLASE026';
% -- Save tables -- % 
% CD into the Tables folder 
% cd(savePath)

% Concatenate all individual tables into one big one and save it as well as
% a csv file of the concatenated tables 
allTables = []; % Placeholder for the concatenated table

mdir = dir;
fdir = {mdir.name};
dirPTnames = fdir(contains(fdir,'CLASE')); % List of part folder names

for i = 1:length(dirPTnames)

    load(dirPTnames{i}); % loads as partTab 

    allTables = [allTables; partTab];

end % for / i 

% Save 
% saveName = strcat(partID, '_allFreqTabs.mat');
saveName = 'CLASE035_allFreqTabs.mat';

% save table
save(saveName, "allTables");

% Save csv 
% tableSaveName = strcat(partID, '_allFreqTabs.csv');
% writetable(allTables, tableSaveName);
writetable(allTables,'CLASE035_allFreqTabs.csv');




