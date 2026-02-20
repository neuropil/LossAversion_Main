mainDir = 'Z:\LossAversion\Patient folders';
cd(mainDir)
matdlist = dir();
matdlist2 = struct2table(matdlist);
matdlist3 = matdlist2.name;
matdlist4 = matdlist3(~ismember(matdlist3,{'.','..','.DS_Store'}));


for matI = 1:length(matdlist4)

    % CD to behavior dir
    nwbdir = [mainDir , filesep , matdlist4{matI}, filesep , 'NWB-Processing\NWB_Data'];

    if ~exist(nwbdir,'dir')
        disp('missing')
        keyboard
    end

    cd(nwbdir)

    tmpNwbdir = dir('*.nwb');
    tmpNwbFl = {tmpNwbdir.name};
    tmpNwbFl2 = tmpNwbFl{contains(tmpNwbFl,'filter')};

    tmpNW = nwbRead(tmpNwbFl2)

    disp(['success read' , num2str(matI)])

    % 
    % % Load raw
    % matt1 = dir('*.csv');
    % matt1b = struct2table(matt1);
    % matt1c = matt1b.name;
    % 
    % tmpCSV = readtable(matt1c);
    % 
    % % Fix Event table
    % 
    % tmpCSV = removevars(tmpCSV,["Coord_x","Coord_y","Coord_z","Interpolation",...
    %     "LocationType","Radius","OrigCoord_x","OrigCoord_y","OrigCoord_z",...
    %     "SurfaceElectrode","DistanceShifted","SurfaceType","VertexNumber",...
    %     "Sphere_x","Sphere_y","Sphere_z","T1A","T1R","T1S","MRVoxel_I",...
    %     "MRVoxel_J","MRVoxel_K"]);
    % 
    % subDiri = [mainDir , filesep , matdlist4{matI} , filesep , 'RCPshare'];
    % 
    % if ~exist(subDiri,'dir')
    %     mkdir(subDiri)
    % end
    % 
    % cd(subDiri)
    % 
    % saveFileName = ['MNI_Labels_',matdlist4{matI},'.mat'];
    % 
    % save(saveFileName,'tmpCSV');
    % 
    % disp([matdlist4{matI} , ' Done'])



end