tmpSub = 'Z:\LossAversion\Patient folders\CLASE007';
cd(tmpSub)

behavioralloc = [tmpSub , filesep , 'Behavioral-data\EventBehavior'];

cd(behavioralloc)

load('clase007_BehEvTable_v2.mat')

% Use TrialEvTm to find location in TimeVec


%% 

mainDir = 'Z:\LossAversion\Patient folders';
cd(mainDir)
matdlist = dir();
matdlist2 = struct2table(matdlist);
matdlist3 = matdlist2.name;
matdlist4 = matdlist3(~ismember(matdlist3,{'.','..','.DS_Store'}));


for matI = 1:length(matdlist4)

    % CD to behavior dir
    behavdir = [mainDir , filesep , matdlist4{matI}, filesep , 'Behavioral-data'];

    cd(behavdir)

    % Load raw
    matt1 = dir('*.mat');
    matt1b = struct2table(matt1);
    matt1c = matt1b.name;

    load(matt1c,'subjdata');

    % CD to TTl table

    behavdir2 = [mainDir , filesep , matdlist4{matI}, filesep ,...
        'Behavioral-data' , filesep , 'EventBehavior'];

    cd(behavdir2)

    matt2 = dir('*.mat');
    matt2b = struct2table(matt2);
    matt2c = matt2b.name;

    if length(matt2) > 1
       tmpMatb = matt2c{contains(matt2c,'v2')};
       load(tmpMatb,'eventTABLE')
    else
        load(matt2c,'eventTABLE')
    end

    % Fix Event table

    if matches(matdlist4{matI},'CLASE006')
        tmpREP = repelem((1:135)', 5);
        tmpREP2 = tmpREP(2:end);
        eventTABLE.Trials = tmpREP2;
    else
        eventTABLE.Trials = repelem((1:135)', 5);
    end
    eventTABLE = removevars(eventTABLE,"OffsetSecs");

    % Combine tables

    if matches(matdlist4{matI},'CLASE006')
        trialIDn = transpose(2:135);
        trialStartTm = eventTABLE.TrialEvTm(matches(eventTABLE.TrialEvID,'choiceShow'));

        behavInfo = table(trialIDn,trialStartTm , subjdata.cs.riskyGain(2:135) , subjdata.cs.riskyLoss(2:135) ,...
            subjdata.cs.alternative(2:135) , subjdata.cs.ischecktrial(2:135) , subjdata.cs.choice(2:135),...
            subjdata.cs.outcome(2:135) , 'VariableNames',{'TrialID','TrialStartEvent','RiskyGainVal',...
            'RiskyLossVal','AlternativeVal','IsCheckTrial', 'ChoiceSelected','Outcome'});
    else
        trialIDn = transpose(1:135);
        trialStartTm = eventTABLE.TrialEvTm(matches(eventTABLE.TrialEvID,'choiceShow'));

        behavInfo = table(trialIDn, trialStartTm , subjdata.cs.riskyGain , subjdata.cs.riskyLoss ,...
            subjdata.cs.alternative , subjdata.cs.ischecktrial , subjdata.cs.choice,...
            subjdata.cs.outcome , 'VariableNames',{'TrialID','TrialStartEvent','RiskyGainVal',...
            'RiskyLossVal','AlternativeVal','IsCheckTrial', 'ChoiceSelected','Outcome'});
    end

    subDiri = [mainDir , filesep , matdlist4{matI} , filesep , 'RCPshare'];

    if ~exist(subDiri,'dir')
        mkdir(subDiri)
    end

    cd(subDiri)

    saveFileName = ['BehaviorTTL_',matdlist4{matI},'.mat'];

    save(saveFileName,'eventTABLE','behavInfo');

    disp([matdlist4{matI} , ' Done'])



end