%% step1_epoch_data.m
% Epochs continuous data around each trigger type per subject per region
% Output: saves one struct per subject per region into /epochs/

% clear; clc;
saveLoc = 'Z:\LossAversion\LH_Data\JAT_TrialData';
cd('C:\Users\Admin\Documents\Github\LossAversion_Main\LA_Manuscript')
conAllsubs = readtable('ContactNumbersLA.xlsx');

for cii = 1:height(conAllsubs)

    tmpRow = conAllsubs(cii,:);
 
    tempPtID = tmpRow.SubjectCL{1};
    Hemi = tmpRow.HemiS{1};
    BrainArea = tmpRow.HSBA{1};
    conNumsTi = tmpRow.ContactNums{1};

    conNums = parseCONs(conNumsTi);

    [subjectTrialInfo] = Trial_Generate_AI_Step1(tempPtID , Hemi , BrainArea , conNums);

    saveName = [tmpRow.SubjectCL{1},'_',tmpRow.HemiS{1},'_',tmpRow.nSBA{1},'_TrialDATA.mat'];
    cd(saveLoc)
    save(saveName,'subjectTrialInfo');


end

%% ---- USER SETTINGS ----
% project_dir   = '/path/to/project_root';
% data_dir      = fullfile(project_dir, 'data');
% epoch_dir     = fullfile(project_dir, 'epochs');
% if ~exist(epoch_dir, 'dir'), mkdir(epoch_dir); end
% 
% fs            = 500;                    % sampling frequency (Hz)
% subjects      = {'sub-01', 'sub-02'};   % add all subjects
% regions       = {'AMY', 'OFC', 'ACC'};  % add all regions
% 
% % Trigger codes — edit to match your actual codes
% TRIG_GAMBLE   = 10;   % gamble options appear on screen
% TRIG_CUE      = 20;   % cue to respond
% TRIG_OUTCOME  = 30;   % outcome screen
% 
% % Epoch windows [start, end] in SECONDS relative to trigger onset
% % Adjust based on your task timing
% win.gamble    = [-0.2,  1.5];   % pre-stim baseline + gamble viewing
% win.cue       = [-0.2,  1.5];   % pre-cue baseline + response window
% win.outcome   = [-0.2,  2.0];   % pre-outcome baseline + feedback
% 
% epoch_types   = {'gamble', 'cue', 'outcome'};
% trig_codes    = [TRIG_GAMBLE, TRIG_CUE, TRIG_OUTCOME];
% windows       = {win.gamble, win.cue, win.outcome};

%% ---- MAIN LOOP ----
for s = 1:numel(subjects)
    sub = subjects{s};
    sub_dir = fullfile(data_dir, sub);

    % --- Load triggers ---
    % Expected: struct with fields trig.timestamps (samples) and trig.codes
    trig_file = fullfile(sub_dir, sprintf('%s_triggers.mat', sub));
    trig      = load(trig_file);           % loads 'trig' struct

    % --- Load behavior ---
    % Expected: struct with nTrials x 1 fields
    beh_file  = fullfile(sub_dir, sprintf('%s_behavior.mat', sub));
    beh       = load(beh_file);            % loads 'beh' struct

    for r = 1:numel(regions)
        rgn = regions{r};

        % --- Load neural data ---
        % Expected: matrix [n_channels x n_samples]
        neural_file = fullfile(sub_dir, sprintf('%s_region-%s.mat', sub, rgn));
        if ~exist(neural_file, 'file')
            fprintf('Skipping %s %s — file not found\n', sub, rgn);
            continue
        end
        nd      = load(neural_file);
        neural  = nd.data;          % [n_channels x n_samples]
        n_chan  = size(neural, 1);

        % --- Epoch loop ---
        for e = 1:numel(epoch_types)
            etype    = epoch_types{e};
            code     = trig_codes(e);
            win_sec  = windows{e};

            % Find trigger onset samples for this code
            trial_onsets = trig.timestamps(trig.codes == code);
            n_trials     = numel(trial_onsets);

            win_samp     = round(win_sec * fs);   % [start_samp, end_samp]
            n_samp       = diff(win_samp) + 1;
            t_axis       = linspace(win_sec(1), win_sec(2), n_samp); % time in sec

            % Preallocate: [channels x time x trials]
            epoch_data   = nan(n_chan, n_samp, n_trials);

            for tr = 1:n_trials
                idx_start = trial_onsets(tr) + win_samp(1);
                idx_end   = trial_onsets(tr) + win_samp(2);

                % Boundary check
                if idx_start < 1 || idx_end > size(neural, 2)
                    fprintf('  Trial %d out of bounds — skipping\n', tr);
                    continue
                end

                epoch_data(:, :, tr) = neural(:, idx_start:idx_end);
            end

            % --- Baseline correct (z-score to pre-stim window) ---
            baseline_idx = t_axis >= win_sec(1) & t_axis < 0;
            epoch_data   = baseline_correct(epoch_data, baseline_idx);

            % --- Save ---
            out.data      = epoch_data;   % [chan x time x trials]
            out.t_axis    = t_axis;
            out.fs        = fs;
            out.subject   = sub;
            out.region    = rgn;
            out.epoch     = etype;
            out.n_trials  = n_trials;
            out.beh       = beh;          % attach behavior to each file

            out_name = fullfile(epoch_dir, ...
                sprintf('%s_region-%s_epoch-%s.mat', sub, rgn, etype));
            save(out_name, 'out', '-v7.3');
            fprintf('Saved: %s\n', out_name);
        end
    end
end

%% ---- HELPER: baseline correction ----
function data_bc = baseline_correct(data, bl_idx)
    % data: [chan x time x trials]
    % bl_idx: logical index into time dimension
    bl_mean  = mean(data(:, bl_idx, :), 2);   % [chan x 1 x trials]
    bl_std   = std(data(:, bl_idx, :), 0, 2); % [chan x 1 x trials]
    bl_std(bl_std == 0) = 1;                  % avoid divide by zero
    data_bc  = (data - bl_mean) ./ bl_std;    % z-score to baseline
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