function [] = LA_Build_allTable_BAND_V1()


cd('Z:\LossAversion\LH_Data\JAT_FoooFtablesZS')

matDIR = dir('*.mat');
matDIR2 = struct2table(matDIR);
matDIR3 = matDIR2.name;

FNallTABLE = table;
FOallTABLE = table;
KNallTABLE = table;

for ii = 1:length(matDIR3)

    load(matDIR3{ii},'allFoooftab')

    BAonly = {};
    HMonly = {};
    FN_Z_Band = table;
    FO_Z_Band = table;
    KN_Z_Band = table;

    for ii2 = 1:height(allFoooftab)

        % Band parameters
        tmpROW = allFoooftab(ii2,:);
        remCols = ["FOOOFoutput","RAWlfp","SpecPARAMJ"];
        keepTab = removevars(tmpROW,remCols);

        BAonlytmp = extractAfter(char(allFoooftab.brainA(ii2)),'_');
        HMonlytmp = extractBefore(char(allFoooftab.brainA(ii2)),'_');

        % FN
        tmpFN = tmpROW.SpecPARAMJ{1}.FIXED.zPkpwer;
        bandTabFN = loopThruBands(tmpFN);

        if istable(bandTabFN)
            FN_keepTab = keepTab(repelem(1:height(keepTab), height(bandTabFN)), :);
            FN_Z_Bcomb = [FN_keepTab , bandTabFN];

            FN_Z_Bcomb.HMonly = repmat({HMonlytmp},height(FN_Z_Bcomb),1);
            FN_Z_Bcomb.BAonly = repmat({BAonlytmp},height(FN_Z_Bcomb),1);

            FN_Z_Band = [FN_Z_Band ; FN_Z_Bcomb];
        end

        % FO
        tmpFO = tmpROW.FOOOFoutput{1}.zPkpwer;
        tmpFO2 = array2table(tmpFO,'VariableNames',{'Hz','X','Y','AmpZS'});
        bandTabFO = loopThruBands(tmpFO2);

        if istable(bandTabFO)
            FO_keepTab = keepTab(repelem(1:height(keepTab), height(bandTabFO)), :);
            FO_Z_Bcomb = [FO_keepTab , bandTabFO];

            FO_Z_Bcomb.HMonly = repmat({HMonlytmp},height(FO_Z_Bcomb),1);
            FO_Z_Bcomb.BAonly = repmat({BAonlytmp},height(FO_Z_Bcomb),1);

            FO_Z_Band = [FO_Z_Band ; FO_Z_Bcomb];
        end

        % KN
        tmpKN = tmpROW.SpecPARAMJ{1}.KNEE.zPkpwer;
        bandTabKN = loopThruBands(tmpKN);

        if istable(bandTabKN)
            KN_keepTab = keepTab(repelem(1:height(keepTab), height(bandTabKN)), :);
            KN_Z_Bcomb = [KN_keepTab , bandTabKN];

            KN_Z_Bcomb.HMonly = repmat({HMonlytmp},height(KN_Z_Bcomb),1);
            KN_Z_Bcomb.BAonly = repmat({BAonlytmp},height(KN_Z_Bcomb),1);

            KN_Z_Band = [KN_Z_Band ; KN_Z_Bcomb];
        end

        disp(ii2)
    end

    FNallTABLE = [FNallTABLE ; FN_Z_Band];
    % tmpFOtab = [allFoooftab2 , table(FO_Z_Exp,FO_Z_Off,BAonly,HMonly,...
    %     'VariableNames',{'Z_Exp','Z_Off','BArea','HEMI'})];
    FOallTABLE = [FOallTABLE ; FO_Z_Band];
    KNallTABLE = [KNallTABLE ; KN_Z_Band];
end

FNallTABLE = renamevars(FNallTABLE,["HMonly","BAonly"],["HEMI","BArea"]);
FOallTABLE = renamevars(FOallTABLE,["HMonly","BAonly"],["HEMI","BArea"]);
KNallTABLE = renamevars(KNallTABLE,["HMonly","BAonly"],["HEMI","BArea"]);
test = 1;

















end



function [outTable] = loopThruBands(FN_Bands)

FN_Bands = FN_Bands(~(FN_Bands.Hz < 4),:);

outBand = cell(height(FN_Bands),1);
zPower = zeros(height(FN_Bands),1);

for i = 1:height(FN_Bands)
    if FN_Bands.Hz(i) > 4 && FN_Bands.Hz(i) < 8
        outBand{i} = 'theta';
        zPower(i) = FN_Bands.AmpZS(i);
    elseif FN_Bands.Hz(i) > 8 && FN_Bands.Hz(i) < 13
        outBand{i} = 'beta';
        zPower(i) = FN_Bands.AmpZS(i);
    elseif FN_Bands.Hz(i) > 13 && FN_Bands.Hz(i) < 30
        outBand{i} = 'beta';
        zPower(i) = FN_Bands.AmpZS(i);
    elseif FN_Bands.Hz(i) > 30 && FN_Bands.Hz(i) < 75
        outBand{i} = 'gamma';
        zPower(i) = FN_Bands.AmpZS(i);
    end
end

if isempty(outBand)
    outTable = nan;
else
    outBand = categorical(outBand);

    outTable = table(outBand,zPower,'VariableNames',{'Band','Zpower'});
end
end