% partID = 'CLASE001';

function [] = FOOOF_freqTAB_Create(partID)

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
    load(tempName) % loads as FOOOF_tab

    % Find the rows where the EpochID column is "Outcome"
    rowsWithChoice = strcmp(FOOOF_tab.EpochID, "choice");
    rowsWithResponse = strcmp(FOOOF_tab.EpochID, "response");
    rowsWithOutcome = strcmp(FOOOF_tab.EpochID, "outcome");

    % Extract those rows from the table and create a new table
    tempChoice = FOOOF_tab(rowsWithChoice, :); % choice epoch
    tempResp = FOOOF_tab(rowsWithResponse, :); % Response epoch
    tempOut = FOOOF_tab(rowsWithOutcome, :); % tmpOut is now the table that has all the outcome epochs in it

    % Add new columns to tables
    tempChoice.TrialID = (1:135)';
    tempChoice.Frequency = zeros(height(tempChoice), 1);
    tempChoice.Power = zeros(height(tempChoice), 1);

    tempResp.TrialID = (1:135)';
    tempResp.Frequency = zeros(height(tempResp), 1);
    tempResp.Power = zeros(height(tempResp), 1);

    tempOut.TrialID = (1:135)';
    tempOut.Frequency = zeros(height(tempOut), 1);
    tempOut.Power =  zeros(height(tempOut), 1);

     % Identify rows where the last column is empty ----- STOPPED HERE
     % CHECK THIS CODE 
    emptyRowsOut = cellfun(@isempty, tempOut.FOOOFoutput);
    % Remove those rows
    tempOut(emptyRowsOut,:) = [];

    % Create GG/GL/AN outcome tables
    tempOutGG = tempOut(tempOut.OutcomeGain,:);
    tempOutGL = tempOut(tempOut.OutcomeLoss,:);
    tempOutAN = tempOut(tempOut.OutcomeNeutral, :);


    % Create empty tables to concatenate
    tempChoiceTAB = table([], [], [], [], [], [], [], [], [], ...
        'VariableNames', {'EpochID', 'LA', 'OutcomeGain', 'OutcomeLoss', ...
        'OutcomeNeutral', 'FOOOFoutput', 'TrialID', 'Frequency', 'Power'});

    tempRespTAB = table([], [], [], [], [], [], [], [], [], ...
        'VariableNames', {'EpochID', 'LA', 'OutcomeGain', 'OutcomeLoss', ...
        'OutcomeNeutral', 'FOOOFoutput', 'TrialID', 'Frequency', 'Power'});

    tempOutGGTAB = table([], [], [], [], [], [], [], [], [], ...
        'VariableNames', {'EpochID', 'LA', 'OutcomeGain', 'OutcomeLoss', ...
        'OutcomeNeutral', 'FOOOFoutput', 'TrialID', 'Frequency', 'Power'});

    tempOutGLTAB = table([], [], [], [], [], [], [], [], [], ...
        'VariableNames', {'EpochID', 'LA', 'OutcomeGain', 'OutcomeLoss', ...
        'OutcomeNeutral', 'FOOOFoutput', 'TrialID', 'Frequency', 'Power'});

    tempOutANTAB = table([], [], [], [], [], [], [], [], [], ...
        'VariableNames', {'EpochID', 'LA', 'OutcomeGain', 'OutcomeLoss', ...
        'OutcomeNeutral', 'FOOOFoutput', 'TrialID', 'Frequency', 'Power'});

    %%% Choice %%%
    for ci = 1:height(tempChoice)
        tempChoicePower = tempChoice.FOOOFoutput{ci,1}.peak_params(:,2);
        tempChoiceFreq = tempChoice.FOOOFoutput{ci,1}.peak_params(:,1);

        tempChoiceInner = repmat(tempChoice(ci,:), length(tempChoiceFreq), 1);

        tempChoiceInner.Frequency = tempChoiceFreq;
        tempChoiceInner.Power = tempChoicePower;

        tempChoiceTAB = [tempChoiceTAB ; tempChoiceInner];

    end % for / ci

    % Frequecy labeling
    thetaIDX = (tempChoiceTAB.Frequency >= 1 & tempChoiceTAB.Frequency <= 8);
    alphaIDX = (tempChoiceTAB.Frequency >= 8 & tempChoiceTAB.Frequency <= 12);
    betaIDX = (tempChoiceTAB.Frequency >= 12 & tempChoiceTAB.Frequency <= 30);
    gammaIDX = (tempChoiceTAB.Frequency >= 30);

    ChoiceArray = cell(size(thetaIDX));
    ChoiceArray(thetaIDX) = {'theta'};
    ChoiceArray(alphaIDX) = {'alpha'};
    ChoiceArray(betaIDX) = {'beta'};
    ChoiceArray(gammaIDX) = {'gamma'};

    tempChoiceTAB.FreqLabel = ChoiceArray;


    %%% Response %%%
    % Identify rows where the last column is empty
    emptyRows = cellfun(@isempty, tempResp.FOOOFoutput);
    % Remove those rows
    tempResp(emptyRows,:) = [];

    for ri = 1:height(tempResp)

        tempRespPower = tempResp.FOOOFoutput{ri,1}.peak_params(:,2);
        tempRespFreq = tempResp.FOOOFoutput{ri,1}.peak_params(:,1);

        tempRespInner = repmat(tempResp(ri, :), length(tempRespFreq), 1);

        tempRespInner.Frequency = tempRespFreq;
        tempRespInner.Power = tempRespPower;

        tempRespTAB = [tempRespTAB;tempRespInner];

    end % for / ri

    % Frequecy labeling
    thetaIDX = (tempRespTAB.Frequency >= 1 & tempRespTAB.Frequency <= 8);
    alphaIDX = (tempRespTAB.Frequency >= 8 & tempRespTAB.Frequency <= 12);
    betaIDX = (tempRespTAB.Frequency >= 12 & tempRespTAB.Frequency <= 30);
    gammaIDX = (tempRespTAB.Frequency >= 30);

    RespArray = cell(size(thetaIDX));
    RespArray(thetaIDX) = {'theta'};
    RespArray(alphaIDX) = {'alpha'};
    RespArray(betaIDX) = {'beta'};
    RespArray(gammaIDX) = {'gamma'};

    tempRespTAB.FreqLabel = RespArray;




    %%% OutGG %%%

    for gg = 1:height(tempOutGG)
        tempOutGGPower = tempOutGG.FOOOFoutput{gg,1}.peak_params(:,2);
        tempOutGGFreq = tempOutGG.FOOOFoutput{gg,1}.peak_params(:,1);

        tempOutGGInner = repmat(tempOutGG(gg,:), length(tempOutGGFreq), 1);

        tempOutGGInner.Frequency = tempOutGGFreq;
        tempOutGGInner.Power = tempOutGGPower;

        tempOutGGTAB = [tempOutGGTAB ; tempOutGGInner];

    end % for / gg

    % Frequecy labeling
    thetaIDX = (tempOutGGTAB.Frequency >= 1 & tempOutGGTAB.Frequency <= 8);
    alphaIDX = (tempOutGGTAB.Frequency >= 8 & tempOutGGTAB.Frequency <= 12);
    betaIDX = (tempOutGGTAB.Frequency >= 12 & tempOutGGTAB.Frequency <= 30);
    gammaIDX = (tempOutGGTAB.Frequency >= 30);

    OutGGArray = cell(size(thetaIDX));
    OutGGArray(thetaIDX) = {'theta'};
    OutGGArray(alphaIDX) = {'alpha'};
    OutGGArray(betaIDX) = {'beta'};
    OutGGArray(gammaIDX) = {'gamma'};

    tempOutGGTAB.FreqLabel = OutGGArray;


    %%% OutGL %%%

    for gl = 1:height(tempOutGL)
        tempOutGLPower = tempOutGL.FOOOFoutput{gl,1}.peak_params(:,2);
        tempOutGLFreq = tempOutGL.FOOOFoutput{gl,1}.peak_params(:,1);

        tempOutGLInner = repmat(tempOutGL(gl,:), length(tempOutGLFreq), 1);

        tempOutGLInner.Frequency = tempOutGLFreq;
        tempOutGLInner.Power = tempOutGLPower;

        tempOutGLTAB = [tempOutGLTAB ; tempOutGLInner];

    end % for / gg

    % Frequecy labeling
    thetaIDX = (tempOutGLTAB.Frequency >= 1 & tempOutGLTAB.Frequency <= 8);
    alphaIDX = (tempOutGLTAB.Frequency >= 8 & tempOutGLTAB.Frequency <= 12);
    betaIDX = (tempOutGLTAB.Frequency >= 12 & tempOutGLTAB.Frequency <= 30);
    gammaIDX = (tempOutGLTAB.Frequency >= 30);

    OutGLArray = cell(size(thetaIDX));
    OutGLArray(thetaIDX) = {'theta'};
    OutGLArray(alphaIDX) = {'alpha'};
    OutGLArray(betaIDX) = {'beta'};
    OutGLArray(gammaIDX) = {'gamma'};

    tempOutGLTAB.FreqLabel = OutGLArray;



    %%% OutAN %%% 

    for an = 1:height(tempOutAN)
        tempOutANPower = tempOutAN.FOOOFoutput{an,1}.peak_params(:,2);
        tempOutANFreq = tempOutAN.FOOOFoutput{an,1}.peak_params(:,1);

        tempOutANInner = repmat(tempOutAN(an,:), length(tempOutANFreq), 1);

        tempOutANInner.Frequency = tempOutANFreq;
        tempOutANInner.Power = tempOutANPower;

        tempOutANTAB = [tempOutANTAB ; tempOutANInner];

    end % for / gg

     % Frequecy labeling
    thetaIDX = (tempOutANTAB.Frequency >= 1 & tempOutANTAB.Frequency <= 8);
    alphaIDX = (tempOutANTAB.Frequency >= 8 & tempOutANTAB.Frequency <= 12);
    betaIDX = (tempOutANTAB.Frequency >= 12 & tempOutANTAB.Frequency <= 30);
    gammaIDX = (tempOutANTAB.Frequency >= 30);

    OutANArray = cell(size(thetaIDX));
    OutANArray(thetaIDX) = {'theta'};
    OutANArray(alphaIDX) = {'alpha'};
    OutANArray(betaIDX) = {'beta'};
    OutANArray(gammaIDX) = {'gamma'};

    tempOutANTAB.FreqLabel = OutANArray;


      % -- Combine frequency tables -- %
    freqTabs = vertcat(tempChoiceTAB, tempRespTAB, tempOutGGTAB, tempOutGLTAB, tempOutANTAB);

   

    % Add LA Val % 
    % LA Classifications
highLA = {'CLASE009', 'CLASE018', 'CLASE019', 'CLASE022', 'CLASE029', 'CLASE030'};
lowLA = {'CLASE023', 'CLASE024', 'CLASE027', 'CLASE034'};
neutralLA = {'CLASE001', 'CLASE007', 'CLASE008', 'CLASE026', 'CLASE031', 'CLASE035'};

if matches(partID, highLA)
    LAval = string('High');

elseif matches(partID, lowLA)
    LAval = string('Low');

else matches(partID, neutralLA)
    LAval = string('Neutral');
end % if else

freqTabs.LAval = repmat(LAval, height(freqTabs), 1); % add laval to table

    % Save %

    % Temp name
    nameSplit = split(tempName, '_');
    partID = nameSplit(1);
    Hemi = nameSplit(2);
    BA = nameSplit(3);

    % Save name 
    saveName = strcat(partID, '_', Hemi, '_', BA, '_', 'periodFreq');
    
    % Repeat identifiers for length fo the table
    partID = repmat(partID, height(freqTabs),1);
    Hemi = repmat(Hemi, height(freqTabs),1);
    BA = repmat(BA, height(freqTabs),1);

    % Create ID table
    IDtab = table(partID, Hemi, BA, 'VariableNames',{'partID', 'Hemi', 'BrainArea'});

    % Combine freqTabs and ID table
    partTabFreq = [IDtab, freqTabs];

    % Save loc
    savePath = strcat(partPath + "\Tables");
    cd(savePath) 

    % save table 
    save(saveName, "partTabFreq");

    % CD back to fooof data 
    cd(partPath)
    
end % for / i

end % function 