%% Clear
clear;
clc;

%% Add path
addpath('/home/kevin/Desktop/Lisa_Data/Functions/');

%% Set parameters
fs = 500;
epochEvents = {'choice';'response';'outcome'};
nSamp = floor((100*fs)./1000);
nIter = 100;

%% Load the subject directory
subDir = '/home/kevin/Desktop/Lisa_Data/Subject_Data/';
subjects = dir(subDir);
subjects = string({subjects.name}');
subjects = subjects(3:end);

%% Loop over subjects
for ii = 10:10%1:length(subjects) - 8 was the first one with channel names

    %% Create the path to load the data
    subData = fullfile(subDir,subjects(ii,1),'ProcessedEphys');

    %% Get only the epoch ephys
    channel_data = dir(fullfile(subData,'*EpochEphys.mat'));
    channel_data = string({channel_data.name}');

    %% Loop over files
    for jj = 1:length(channel_data)

        %% Create the file path
        file = fullfile(subData,channel_data(jj,1));

        %% Get the electrode IDs
        electrode = extractBefore(channel_data(jj,1),'_EpochEphys.mat'); % subject 10 is messed up - has the wrong data
        electrode = extractAfter(electrode,append(subjects(ii,1),'_'));

        %% Load the data
        data = load(file);
        data = data.EpochEphysTab;
        clear file

        %% Grab the EpochIDs
        epochID = string(data.EpochID);

        %% Get number of contacts
        nContact = size(data.Ephys{1},1);

        %% Save wire and channel names
        for kk = 1:nContact
            wire{jj,1}{kk,1} = electrode;
            wireFull{jj,1}{kk,1} = append(electrode,'-',num2str(kk));
        end
        clear nContact

        %% Loop over epochEvents to allocate data
        for kk = 1:length(epochEvents)

            %% String match
            I = strcmp(epochID,epochEvents{kk,1});

            %% Get the corresponding data
            fullData.(epochEvents{kk,1}).LFP(jj,:) = data.Ephys(I);

            %% Clear
            clear I;

        end

    end

    %% Loop to do GC
    for jj = 1:length(epochEvents)

        %% Get the number of epochs
        % nEpochs = size(fullData.(epochEvents{jj,1}).LFP,2);

        %% Specify number of epochs
        nEpochs = 135;

        %% Channel names
        chNames = vertcat(wireFull{:});

        %% Save the individual channel names
        fullData.(epochEvents{jj,1}).chanNames = chNames;

        %% Note channel interactions
        for kk = 1:length(chNames)
            for ll = 1:length(chNames)
                if kk ~= ll
                    fullData.(epochEvents{jj,1}).chanInfo{kk,ll} = string(append('FROM ',chNames{ll},' | TO ',chNames{kk}));
                end
            end
        end

        %% Loop over epochs to compute GC
        for kk = 1:nEpochs

            %% Print
            fprintf('Analyzing event %s epoch %d..\n',epochEvents{jj,1},kk)

            %% Ignore epoch if any part of the epoch is less than 100 ms
            if size(fullData.choice.LFP{1,kk},2) < nSamp || size(fullData.response.LFP{1,kk},2) < nSamp || size(fullData.outcome.LFP{1,kk},2) < nSamp
                continue
            end

            %% Get the data
            GC_data = fullData.(epochEvents{jj,1}).LFP(:,kk);
            GC_data = cell2mat(GC_data);

            %% Normalize the epoch
            GC_data = normalize(GC_data,2,"zscore");

            %% Perform GC
            for ll = 1:size(GC_data,1)
                for mm = 1:size(GC_data,1)
                    if ll ~= mm

                        %% Create the data for GC
                        reduced = GC_data(ll,:);
                        full = [GC_data(ll,:);GC_data(mm,:)];

                        %% Perform autoregression
                        [R_reduced,~] = Autoregressive_Process_V1(reduced,nSamp);
                        [R_full,~] = Autoregressive_Process_V1(full,nSamp);

                        %% Calculatye the variances
                        var_reduced = var(R_reduced,0,2);
                        var_full = var(R_full,0,2);

                        %% Calculate the GC value
                        fullData.(epochEvents{jj,1}).GC(ll,mm,kk) = log(var_reduced./var_full(1,1));

                        %% Calculate random permutations from sender
                        for iter = 1:nIter
                            ix = randperm(size(GC_data,2));
                            full_stat = [GC_data(ll,:);GC_data(mm,ix)];
                            [R_full_stat,~] = Autoregressive_Process_V1(full_stat,nSamp);
                            var_full_stat = var(R_full_stat,0,2);
                            gc_stat(iter,1) = log(var_reduced./var_full_stat(1,1));
                            clear ix full_stat R_full_stat var_full_stat
                        end

                        %% Calculate if significant
                        sig_val = prctile(gc_stat,95);
                        if fullData.(epochEvents{jj,1}).GC(ll,mm,kk) >= sig_val
                            fullData.(epochEvents{jj,1}).GC_Sig(ll,mm,kk) = 1;
                        else
                            fullData.(epochEvents{jj,1}).GC_Sig(ll,mm,kk) = 0;
                        end

                        %% Clear
                        clear sig_val gc_stat reduced full var_reduced var_full R_full R_reduced

                    end
                end
            end

            %% Clear
            clear GC_data

        end

        %%
        %x = 0;

    end

    %% Save the data
    save_path = fullfile(subDir,subjects(ii,1),'GC_Results');
    if ~isfolder(save_path)
        mkdir(save_path)
    end
    save_file = fullfile(save_path,'gc_analysis.mat');
    save(save_file,'fullData','-v7.3');

    %% Clear
    clear channel_data chNames data electrode epochID fullData save_file save_path subData wire wireFull
    
end