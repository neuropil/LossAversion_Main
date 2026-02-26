function [] = LA_Build_allTable_AFP_V1()


cd('Z:\LossAversion\LH_Data\JAT_FoooFtablesZS')

matDIR = dir('*.mat');
matDIR2 = struct2table(matDIR);
matDIR3 = matDIR2.name;

FNallTABLE = table;
FOallTABLE = table;
KNallTABLE = table;

for ii = 1:length(matDIR3)

    load(matDIR3{ii},'allFoooftab')

    FN_Z_Exp = zeros(height(allFoooftab),1);
    FN_Z_Off = zeros(height(allFoooftab),1);
    FO_Z_Exp = zeros(height(allFoooftab),1);
    FO_Z_Off = zeros(height(allFoooftab),1);
    KN_Z_Exp = zeros(height(allFoooftab),1);
    KN_Z_Off = zeros(height(allFoooftab),1);

    BAonly = cell(height(allFoooftab),1);
    HMonly = cell(height(allFoooftab),1);

    for ii2 = 1:height(allFoooftab)

        % Afit parameters
        FN_Z_Exp(ii2) = allFoooftab.SpecPARAMJ{ii2}.FIXED.AfitParms.Z_Exp;
        FN_Z_Off(ii2) = allFoooftab.SpecPARAMJ{ii2}.FIXED.AfitParms.Z_Off;
        FO_Z_Exp(ii2) = allFoooftab.FOOOFoutput{ii2}.AfitParms.Z_Exp;
        FO_Z_Off(ii2) = allFoooftab.FOOOFoutput{ii2}.AfitParms.Z_Off;
        KN_Z_Exp(ii2) = allFoooftab.SpecPARAMJ{ii2}.KNEE.AfitParms.Z_Exp;
        KN_Z_Off(ii2) = allFoooftab.SpecPARAMJ{ii2}.KNEE.AfitParms.Z_Off;

        BAonly{ii2} = extractAfter(char(allFoooftab.brainA(1)),'_');
        HMonly{ii2} = extractBefore(char(allFoooftab.brainA(1)),'_');
        % % Bands
        % % Loop through FOOOF type
        % for ffi = 1:3
        %     switch ffi
        %         case 1
        %             FN_Bands = allFoooftab.SpecPARAMJ{ii2}.FIXED.zPkpwer;
        % 
        %         case 2
        % 
        % 
        % 
        % 
        %     end
        % end % End of band for loop
    end

    allFoooftab2 = removevars(allFoooftab,["FOOOFoutput","RAWlfp","SpecPARAMJ",...
        "brainA"]);

    tmpFNtab = [allFoooftab2 , table(FN_Z_Exp,FN_Z_Off,BAonly,HMonly,...
        'VariableNames',{'Z_Exp','Z_Off','BArea','HEMI'})];

    FNallTABLE = [FNallTABLE ; tmpFNtab];

    tmpFOtab = [allFoooftab2 , table(FO_Z_Exp,FO_Z_Off,BAonly,HMonly,...
        'VariableNames',{'Z_Exp','Z_Off','BArea','HEMI'})];

    FOallTABLE = [FOallTABLE ; tmpFOtab];

    tmpKNtab = [allFoooftab2 , table(KN_Z_Exp,KN_Z_Off,BAonly,HMonly,...
        'VariableNames',{'Z_Exp','Z_Off','BArea','HEMI'})];

    KNallTABLE = [KNallTABLE ; tmpKNtab];

end


test = 1;

















end



% function [theta, alpha, beta, gamma] = loopThruBands(FN_Bands)
% 
% 
% for i = 1:height(FN_Bands)
% 
% 
% 
% 
% 
% 
% end
% 
% 
% end