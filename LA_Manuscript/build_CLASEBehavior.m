function [] = build_CLASEBehavior(subjectID)

pcNAME = getenv("COMPUTERNAME");
switch pcNAME
    case 'DESKTOP-I5CPDO7'
        mainDIR_clase = 'Z:\LossAversion';
    case 'otherwise'
end

%%%%% - Behavior
behaviorLOC1 = [mainDIR_clase , filesep , 'Patient folders' , filesep ,...
    subjectID,'\Behavioral-data\EventBehavior'];
cd(behaviorLOC1)
tmpMATdir = dir('*.mat');
tmpMATfile = {tmpMATdir.name};
load(tmpMATfile{1},'eventTABLE')

procTABLE = eventTABLE;

behaviorLOC2 = [mainDIR_clase , filesep , 'Patient folders' , filesep ,...
    subjectID,'\Behavioral-data'];
cd(behaviorLOC2)
tmpMATdir2 = dir('*.mat');
tmpMATfile2 = {tmpMATdir2.name};
load(tmpMATfile2{1},'subjdata');

matBEH = subjdata;

matBEH.cs.triBlock = transpose(matBEH.cs.triBlock);
matBEH.cs.subjectIndex = cellstr(matBEH.cs.subjectIndex);

matTAB1 = struct2table(matBEH.cs);

if ~isfield(matBEH.ts,'studystop')
    matTMPst = rmfield(matBEH.ts,["blockStart","studystart","blockBreakStart"]);
else
    matTMPst = rmfield(matBEH.ts,["blockStart","studystart","blockBreakStart","studystop"]);
    matMetaST.studystart = matBEH.ts.studystart;
end
matTMPstF = fieldnames(matTMPst);
for mmfi = 1:length(matTMPstF)
    matTMPst.(matTMPstF{mmfi}) = transpose(matTMPst.(matTMPstF{mmfi}));
end

tmpFIELDS = fieldnames(matTMPst);
rowNUMS = structfun(@(x) height(x), matTMPst);

if ~all(ismember(rowNUMS,135))
    rmINDEX = ~ismember(rowNUMS,135);
    rm2FIELD = tmpFIELDS{rmINDEX};
    matTMPst = rmfield(matTMPst,rm2FIELD);
end

matTAB2 = struct2table(matTMPst);

matMetaST.blockStart = matBEH.ts.blockStart;
matMetaST.blockBreakStart = matBEH.ts.blockBreakStart;
matMetaST.studystop = matBEH.ts.studystart;
matMetaST.params = matBEH.params;
if isfield(matBEH.ts,'studystop')
    matMetaST.payoutTrial = matBEH.payoutTrial;
end

fullbehTAB = [matTAB1 , matTAB2];

newBehDir = [mainDIR_clase , '\Patient folders\' , subjectID , '\LA_M'];

if ~exist(newBehDir,"dir")
    mkdir(newBehDir)
end

cd(newBehDir)

behavFILEname = [subjectID , '_BehTimePerformTAB.mat'];

save(behavFILEname,"fullbehTAB","matMetaST","procTABLE")

end