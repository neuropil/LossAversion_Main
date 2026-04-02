function [] = BatchGenerateInitialFiles_LA_v1()


% saveLoc = 'Z:\LossAversion\LH_Data\JAT_BlockData';
cd('C:\Users\Admin\Documents\Github\LossAversion_Main\LA_Manuscript')
conAllsubs = readtable('ContactNumbersLA_UPDATE.xlsx');

for cii = 21:height(conAllsubs)

    tmpRow = conAllsubs(cii,:);

    tempPtID = tmpRow.SubjectCL{1};
    Hemi = tmpRow.HemiS{1};
    BrainArea = tmpRow.HSBA{1};
    saveBAname = tmpRow.nSBA{1};
    conNumsTi = tmpRow.ContactNums{1};

    conNums = parseCONs(conNumsTi);

    GenerateInitial_LFPdata_LA_v1(tempPtID , Hemi, ...
        6 , conNums , BrainArea , saveBAname);


end



end


function [numSS] = parseCONs(conNumsTi)

if contains(conNumsTi,',')
    % outNum = num2cell(conNumsTi);
    outNumF = extractBefore(conNumsTi,',');
    outNumL = extractAfter(conNumsTi,',');
    numSS = [str2double(outNumF) , str2double(outNumL)];
elseif contains(conNumsTi,';')
    % outNum = num2cell(conNumsTi);
    % numSS = str2double(outNum{1}) : str2double(outNum{3});
    outNumF = extractBefore(conNumsTi,';');
    outNumL = extractAfter(conNumsTi,';');
    numSS = str2double(outNumF) : str2double(outNumL);
else
    numSS = str2double(conNumsTi);
end
disp([newline num2str(numSS) newline]); disp([conNumsTi newline])

end