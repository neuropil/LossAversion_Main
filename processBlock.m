function [epochTab] = processBlock(subBlock, blockNum)

% Create subtables for each epoch
subTab_Choice = subBlock(matches(subBlock.EpochID, 'choice'), :);
subTab_Response = subBlock(matches(subBlock.EpochID, 'response'), :);
subTab_Outcome = subBlock(matches(subBlock.EpochID, 'outcome'), :);

% Create GG/GL/AN outcome tables
subTab_OutGG = subTab_Outcome(subTab_Outcome.OutcomeGain, :);
subTab_OutGL = subTab_Outcome(subTab_Outcome.OutcomeLoss, :);
subTab_OutAN = subTab_Outcome(subTab_Outcome.OutcomeNeutral, :);

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
for ci = 1:height(subTab_Choice)
    tempChoicePower = subTab_Choice.FOOOFoutput{ci,1}.peak_params(:,2);
    tempChoiceFreq = subTab_Choice.FOOOFoutput{ci,1}.peak_params(:,1);

    choicePower = [choicePower; tempChoicePower];
    choiceFreq = [choiceFreq; tempChoiceFreq];

end % for / ci

% Response
for ri = 1:height(subTab_Response)
    tempRespPower = subTab_Response.FOOOFoutput{ri,1}.peak_params(:,2);
    tempRespFreq = subTab_Response.FOOOFoutput{ri,1}.peak_params(:,1);

    respPower = [respPower; tempRespPower];
    respFreq = [respFreq; tempRespFreq];

end % for / ri

% Outcome GG
for gg = 1:height(subTab_OutGG)
    tempOutGGPower = subTab_OutGG.FOOOFoutput{gg,1}.peak_params(:,2);
    tempOutGGFreq = subTab_OutGG.FOOOFoutput{gg,1}.peak_params(:,1);

    outGGPower = [outGGPower; tempOutGGPower];
    outGGFreq = [outGGFreq; tempOutGGFreq];

end % for / gg

% Outcome GL
for gl = 1:height(subTab_OutGL)
    tempOutGLPower = subTab_OutGL.FOOOFoutput{gl,1}.peak_params(:,2);
    tempOutGLFreq = subTab_OutGL.FOOOFoutput{gl,1}.peak_params(:,1);

    outGLPower = [outGLPower; tempOutGLPower];
    outGLFreq = [outGLFreq; tempOutGLFreq];

end % for / gl

% Outcome AN
for an = 1:height(subTab_OutAN)
    tempOutANPower = subTab_OutAN.FOOOFoutput{an,1}.peak_params(:,2);
    tempOutANFreq = subTab_OutAN.FOOOFoutput{an,1}.peak_params(:,1);

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

if length(respPower) >= 1
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

if length(outGGPower) >= 1
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

else
    outGGTemp = [outGGPower outGGFreq];

    OutGGArray = cell(size(outGGTemp));
    OutGGTab = table(outGGLab, outGGTemp, OutGGArray, 'VariableNames',{'Epoch', 'Power', 'Frequency'});
end


% -- Outcome GL -- %
outGLLab = repmat({'OutGL'}, size(outGLPower));

if length(outGLPower) >= 1
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

else
    outGLTemp = [outGLPower outGLFreq];

    OutGLArray = cell(size(outGLTemp));
    OutGLTab = table(outGLLab, outGLTemp, OutGLArray, 'VariableNames',{'Epoch', 'Power', 'Frequency'});

end % if else


% -- Outcome AN -- %
outANLab = repmat({'OutAN'}, size(outANPower));

if length(outANPower) >= 1
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

else

    outANTemp = [outANPower outANFreq];

    OutANArray = cell(size(outANTemp));
    OutANTab = table(outANLab, outANTemp, OutANArray, 'VariableNames',{'Epoch', 'Power', 'Frequency'});
end % if else
% -- Combine frequency tables -- %
freqTab = vertcat(choiceTab, RespTab, OutGGTab, OutGLTab, OutANTab);

blockLabel = repmat(blockNum, height(freqTab), 1);
blockTab = table(blockLabel, 'VariableNames',{'BlockNum'});

% Combine all tables %
epochTab = horzcat(freqTab, blockTab);


end % function