function [] = processLA_behav_v3(subjID , ttlStyle, ttlID, NWBdir, NWBname , behDIRsave)

arguments

    % Required
    subjID (1,1) string
    ttlStyle (1,1) double
    % Optional
    ttlID (1,:) char = 'NA'
    NWBdir (1,1) string = "NA"
    NWBname (1,1) string = "NA"
    behDIRsave (1,1) string = "NA"

end

% NWB Directory
if matches(ttlID,'NA')
    ttlIDu = '(0x0001)';
else
    ttlIDu = ttlID;
end

% NWB Directory
if matches(NWBdir,"NA")
    NWBdirU = uigetdir;
else
    NWBdirU = NWBdir;
end

% Behavior save Directory
if matches(behDIRsave,"NA")
    BHsavedirU = uigetdir;
else
    BHsavedirU = behDIRsave;
end

% NWB file name
if matches(NWBname,"NA")
    NWBfname = uigetfile(NWBname);
else
    NWBfname = NWBname;
end

cd(NWBdirU)

nwbCHECK = which('nwbtest.m');
if isempty(nwbCHECK)
    % find the way to bring up documents folder or search for folder
    matNWB = uigetdir;
    addpath(genpath(matNWB));
end

% Change directory to NWB location
cd(NWBdirU)
% Load NWB file
tmpNWB_LA = nwbRead(NWBfname);

% behavioral timestamps data are in microseconds
eventStamps = tmpNWB_LA.acquisition.get('events').timestamps.load;
eventSids = tmpNWB_LA.acquisition.get('events').data.load;
eventIDcs = cellstr(eventSids);
% BLOCK START TTL
% Total = 5
% Each Trial has 5 TTLS
% initial display of choice options % choice
% start of response window % respWindowS
% end of response window/response % respWindowE
% start of outcome display % outDispS
% end of outcome display % outDispE
% Total of 135 Trials x 5 (5 blocks of 27 trials)
% TRIALS = 675
% TOTAL = 675;


% Clean TTL id codes
switch ttlStyle
    case 1 % NO HEX flags
        % CHECK 1: Check for all TTL events
        hexFlagsTTL = eventIDcs(contains(eventIDcs,'TTL Input'));
        hexOnly = extractBetween(hexFlagsTTL,'(',')');
        decFhex = hex2dec(hexOnly);
        decFhex2 = decFhex(decFhex ~= 0);

        tempTab = tabulate(categorical(decFhex2));
        if height(tempTab) < 5
            ckTT1 = false;
        else
            ckTT1 = true;
        end

        % CHECK 2: Check for 5 blocks
        [~,peakLOCS,~] = findpeaks(diff(eventStamps),'MinPeakDistance',100,'MinPeakHeight',10000000);
        peakLOCS = peakLOCS(peakLOCS > 10);

        blockINDst = zeros(4,2);

        for bi = 1:5
            if bi == 1
                if length(peakLOCS) > 4
                    blockINDst(bi,1) = peakLOCS(bi) + 1;
                    blockINDst(bi,2) = peakLOCS(bi + 1);
                else
                    blockINDst(bi,1) = 1;
                    blockINDst(bi,2) = peakLOCS(bi);
                end
            elseif bi == 5
                if length(peakLOCS) > 4
                    blockINDst(bi,1) = peakLOCS(bi) + 1;
                    blockINDst(bi,2) = blockINDst(bi,1) + 269;
                else
                    blockINDst(bi,1) = peakLOCS(bi - 1) + 1;
                    blockINDst(bi,2) = blockINDst(bi,1) + 269;
                end
            else
                if length(peakLOCS) > 4
                    blockINDst(bi,1) = peakLOCS(bi) + 1;
                    blockINDst(bi,2) = peakLOCS(bi + 1);
                else
                    blockINDst(bi,1) = peakLOCS(bi - 1) + 1;
                    blockINDst(bi,2) = peakLOCS(bi);
                end
            end
        end

        if any(reshape(blockINDst > numel(eventStamps),numel(blockINDst),1))
            ckTT2 = false;
        else
            ckTT2 = true;
        end

        % CHECK 3: Check the each block has 135 elements
        ttlCOUNT = zeros(5,1);
        eventBlocksAll = cell(5,1);
        eventBlockIDs = cell(5,1);
        for tTC = 1:5

            testBlockID = eventIDcs(blockINDst(tTC,1):blockINDst(tTC,2));
            blockNUMS = blockINDst(tTC,1):blockINDst(tTC,2);

            if matches(ttlIDu,'XXXXX')
                ttlCOUNT(tTC) = sum(~contains(testBlockID,'(0x0000)'));
                tmpTTLcount = testBlockID(~contains(testBlockID,'(0x0000)'));
                eventBlockIDs{tTC,1} = transpose(blockNUMS(~contains(testBlockID,'(0x0000)')));
            else
                ttlCOUNT(tTC) = sum(contains(testBlockID,ttlIDu));
                tmpTTLcount = testBlockID(contains(testBlockID,ttlIDu));
                eventBlockIDs{tTC,1} = transpose(blockNUMS(contains(testBlockID,ttlIDu)));
            end
            eventBlocksAll{tTC,1} = tmpTTLcount;

            

        end

        if any(ttlCOUNT ~= 135)
            ckTT3 = false;
        else
            ckTT3 = true;
        end


        [eventTABLE] = getTTLevTab(ckTT1, ckTT2, ckTT3, ttlCOUNT, eventBlockIDs, eventStamps);


    case 2 % MIXED/MISSING HEX flags

        % Determine which hex flags are present and which are missing


    case 3 % HEX flags complete

        [eventTABLE] = getNewTTLs(eventIDcs,eventStamps);

end


cd(BHsavedirU)
saveNAME = [char(subjID) , '_BehEvTable_v2.mat'];
save(saveNAME , "eventTABLE");


end



function [newEVENTs] = getNewEvents(newEvts)

if length(newEvts) > 675

    newEvtsd = [0 ; diff(newEvts)];
    [~,pIND] = findpeaks(newEvtsd,"MinPeakHeight", 10000000);
    if length(newEvtsd) - pIND(end) == 135
        startINDn = length(newEvtsd) - 675 + 1;
        newEVENTs = newEvts(startINDn:length(newEvtsd));
    end


else
    disp('MISCOUNT TTL!!!!!!!!!')
    return
end


end




function [outTable] = getNewTTLs(newEvts , newEvtTS)

hexFlagsTTL = newEvts(contains(newEvts,'TTL Input'));
hexOnly = extractBetween(hexFlagsTTL,'(',')');
decFhex = hex2dec(hexOnly);
decFhex2 = decFhex(decFhex ~= 0);
tempTab = tabulate(categorical(decFhex2));
finTab = tempTab(ismember(cell2mat(tempTab(:,2)),135),:);
finInds = str2double(finTab(:,1));

% uniqueIDS = cellfun(@(x) str2double(x), finTab(:,1));
tmpblockIND = repmat(transpose(1:5),135,1);

if height(finTab) == 5 && all(cell2mat(finTab(:,2)) == 135)

    alltrials = [];
    trialepID = {};
    trialepNum = [];
    allblocks = [];
    newEvts2use = zeros(675,1);
    for bbi = 1:5
        tmpBlck = 135;

        decFhexIND = decFhex2 == finInds(bbi);
        tmpTS2use = newEvtTS(decFhexIND);

        %         tmpTS2use = eventTS(eventS{bbi});
        newEvts2use(tmpblockIND == bbi) = tmpTS2use;
        trialepNumi = repmat(transpose(1:5),27,1);
        alltrialsi = transpose(1:tmpBlck);
        trialepIDi = repmat(transpose({'choiceShow','respWindowS','respWindowE','outDispS',...
            'outDispE'}),27,1);

        blockTi = repmat(bbi,tmpBlck,1);
        alltrials = [alltrials ; alltrialsi];
        trialepID = [trialepID ; trialepIDi];
        trialepNum = [trialepNum ; trialepNumi];
        allblocks = [allblocks ; blockTi];


    end

    %     newEvts2use = [newEvts2use ; tmpTS2use];

end

% compute offset
offsetCk = [diff(newEvts2use/1000000) ; nan];


outTable = table(allblocks, alltrials, trialepNum, trialepID, newEvts2use,...
    offsetCk,'VariableNames',{'Blocks','Trials','TrialEvNum','TrialEvID',...
    'TrialEvTm','OffsetSecs'});



end





function [outTable] = getTTLevTab(check1, check2, check3, blockTot, eventS, eventTS)
% Check 1 = TTL hexes
% Check 2 = 5 blocks
% Check 3 = 135 per block

% if ~check1 && check2 && ~check3

alltrials = [];
trialepID = {};
trialepNum = [];
allblocks = [];
newEvts2use = [];
for bbi = 1:height(blockTot)
    tmpBlck = blockTot(bbi);

    tmpTS2use = eventTS(eventS{bbi});

    if tmpBlck ~= 135

        if tmpBlck < 135

            trialepNum0 = repmat(transpose(1:5),27,1);
            trialepNumi = trialepNum0(2:end);
            alltrialsi = transpose(2:135);
            trialepID0 = repmat(transpose({'choiceShow','respWindowS','respWindowE','outDispS',...
                'outDispE'}),27,1);
            trialepIDi = trialepID0(2:end);
        elseif tmpBlck > 135
            tsDiff = diff(tmpTS2use)/1000000;
            if tsDiff(1) < 1.5 || tsDiff(1) > 2.3
                trialepNumi = repmat(transpose(1:5),27,1);
                alltrialsi = transpose(1:135);
                trialepIDi = repmat(transpose({'choiceShow','respWindowS','respWindowE','outDispS',...
                    'outDispE'}),27,1);
                tmpTS2use = tmpTS2use(2:end);
                tmpBlck = 135;
            end



        end
    else
        trialepNumi = repmat(transpose(1:5),27,1);
        alltrialsi = transpose(1:tmpBlck);
        trialepIDi = repmat(transpose({'choiceShow','respWindowS','respWindowE','outDispS',...
            'outDispE'}),27,1);
    end
    blockTi = repmat(bbi,tmpBlck,1);
    alltrials = [alltrials ; alltrialsi];
    trialepID = [trialepID ; trialepIDi];
    trialepNum = [trialepNum ; trialepNumi];
    allblocks = [allblocks ; blockTi];
    newEvts2use = [newEvts2use ; tmpTS2use];

end
% end
alltrials2 = zeros(size(alltrials));
curVal = 0;
for tti = 1:length(alltrials)
    
    if tti == 1
        alltrials2(tti) = alltrials(tti);
        curVal = curVal + alltrials2(tti);
    elseif tti ~= length(alltrials)

        tmpCheck = alltrials(tti) - alltrials(tti - 1);
        if tmpCheck == 1
            alltrials2(tti) = curVal + 1;
            curVal = alltrials2(tti);
        else
            alltrials2(tti) = (alltrials(tti - 1) + tmpCheck) + curVal;
            curVal = alltrials2(tti);
        end
    else
         alltrials2(tti) = curVal + 1;
    end

end

% compute offset
offsetCk = [diff(newEvts2use/1000000) ; nan];
trialNumSet = reshape(repmat(1:135,5,1),675,1);
trialNumSetu = trialNumSet(alltrials2);

outTable = table(allblocks, alltrials, trialNumSetu, trialepNum, trialepID, newEvts2use,...
    offsetCk,'VariableNames',{'Blocks','Trials','TrialiNum','TrialEvNum','TrialEvID',...
    'TrialEvTm','OffsetSecs'});



end




