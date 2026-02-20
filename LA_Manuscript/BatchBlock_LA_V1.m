function [] = BatchBlock_LA_V1()


saveLoc = 'Z:\LossAversion\LH_Data\JAT_BlockData';
cd('C:\Users\Admin\Documents\Github\LossAversion_Main\LossAversionManuscript')
conAllsubs = readtable('ContactNumbersLA.xlsx');

for cii = 56:height(conAllsubs)

    tmpRow = conAllsubs(cii,:);
 
    tempPtID = tmpRow.SubjectCL{1};
    Hemi = tmpRow.HemiS{1};
    BrainArea = tmpRow.HSBA{1};
    conNumsTi = tmpRow.ContactNums{1};

    conNums = parseCONs(conNumsTi);

    [baselineBlock] = Block_baseline_Generate_v1(tempPtID , Hemi , BrainArea , conNums);

    saveName = [tmpRow.SubjectCL{1},'_',tmpRow.HemiS{1},'_',tmpRow.nSBA{1},'_blocks.mat'];
    cd(saveLoc)
    save(saveName,'baselineBlock');


end



end


function [numSS] = parseCONs(conNumsTi)

if contains(conNumsTi,',')
    outNum = num2cell(conNumsTi);
    numSS = [str2double(outNum{1}) , str2double(outNum{3})];
elseif contains(conNumsTi,';')
    outNum = num2cell(conNumsTi);
    numSS = str2double(outNum{1}) : str2double(outNum{3});
else
    numSS = str2double(conNumsTi);
end


end