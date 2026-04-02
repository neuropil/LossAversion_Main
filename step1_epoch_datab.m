%% step1_epoch_data.m
% Epochs continuous iEEG data using the channel roster table and eventTABLE.
%
% DATA STRUCTURES ASSUMED
% -----------------------
%  chanTable  : table with columns
%                 SubjectCL   – subject ID matching neural file naming (e.g. 'CLASE001')
%                 SubjectMW   – alternate ID (e.g. 'MW_9')
%                 SBA         – macro-region label (e.g. 'AMY', 'OFC', 'AH')
%                 HemiS       – hemisphere ('L' or 'R')
%                 HnSBA       – hemisphere+region label (e.g. 'LAMY')
%                 ContactNums – string of contact indices, semicolon-separated (e.g. '1;3')
%
%  eventTABLE : table (675 x 7) with columns
%                 Blocks, Trials, TrialEvNum, TrialEvID (string), TrialEvTm (µs Unix),
%                 OffsetSecs, TrialID
%                 Epoch-defining events:
%                   'choiceShow...' → gamble onset
%                   'respWindS...'  → response window / cue onset
%                   'outDispS...'   → outcome display onset
%
%  Neural file : .mat per subject per region, variable 'data' [n_chan x n_samples]
%                plus scalar 'rec_start_us' (recording start time in µs Unix)
%                so that sample k corresponds to time: rec_start_us + (k-1)/fs * 1e6
%
% OUTPUT
% ------
%  One .mat per subject × region × epoch type, saved to epoch_dir.
%  Each file contains struct 'out' with fields:
%    data      [n_bipolar x n_time x n_trials]  (baseline z-scored)
%    t_axis    [1 x n_time]  seconds rel. to epoch onset
%    fs, subject, region, hemisphere, epoch, n_trials
%    contacts  cell array of contact-pair labels
%    beh       struct with per-trial behavioral variables

clear; clc;

%% ════════════════════════════════════════════════════════════
%  USER SETTINGS  — edit these paths and parameters
%% ════════════════════════════════════════════════════════════
project_dir = '/path/to/project_root';
data_dir    = fullfile(project_dir, 'data');
epoch_dir   = fullfile(project_dir, 'epochs');
if ~exist(epoch_dir, 'dir'), mkdir(epoch_dir); end

% Channel roster file  (the table shown in Image 1)
chan_table_file = fullfile(project_dir, 'chanTable.mat');  % variable name: chanTable

% Behavioral / event file  (contains eventTABLE and beh struct)
beh_file = fullfile(project_dir, 'behaviorData.mat');      % variables: eventTABLE, beh

fs = 500;   % sampling frequency (Hz)

% ---- Epoch windows [start, end] in seconds relative to event onset ----
win.gamble  = [-0.2,  1.5];   % choiceShow
win.cue     = [-0.2,  1.5];   % respWindS (response window start)
win.outcome = [-0.2,  1.5];   % outDispS

epoch_types  = {'gamble',       'cue',        'outcome'};
event_prefixes = {'choiceShow', 'respWindowE',  'outDispS'};   % prefix-matched against TrialEvID
% NOTE: cue epoch locked to respWindEnd (end of response window = decision made)
windows      = {win.gamble,     win.cue,      win.outcome};

%% ════════════════════════════════════════════════════════════
%  LOAD SHARED TABLES
%% ════════════════════════════════════════════════════════════
ct  = load(chan_table_file);   chanTable  = ct.chanTable;
bd  = load(beh_file);
eventTABLE = bd.eventTABLE;
beh        = bd.beh;           % behavioral struct — attached to every epoch file

% Unique subjects in the channel roster
subjects = unique(chanTable.SubjectCL, 'stable');

%% ════════════════════════════════════════════════════════════
%  EXTRACT EPOCH ONSET TIMES (µs) PER EPOCH TYPE
%  One onset per trial, sorted by TrialID
%% ════════════════════════════════════════════════════════════
onset_us = struct();
for e = 1:numel(epoch_types)
    pfx  = event_prefixes{e};
    % Match rows whose TrialEvID starts with the prefix
    mask = strncmp(eventTABLE.TrialEvID, pfx, numel(pfx));
    sub_tbl = eventTABLE(mask, :);
    % Sort by TrialID so order matches behavioral table
    sub_tbl = sortrows(sub_tbl, 'TrialID');
    onset_us.(epoch_types{e}) = sub_tbl.TrialEvTm;   % [n_trials x 1] µs Unix
    fprintf('Epoch %-8s : %d trials found\n', epoch_types{e}, height(sub_tbl));
end

%% ════════════════════════════════════════════════════════════
%  MAIN LOOP: subjects → rows in chanTable → epoch types
%% ════════════════════════════════════════════════════════════
for s = 1:numel(subjects)
    subCL  = subjects{s};
    % All channel-table rows for this subject
    sub_rows = chanTable(strcmp(chanTable.SubjectCL, subCL), :);
    subMW    = sub_rows.SubjectMW{1};   % e.g. 'MW_9'

    fprintf('\n▶ Subject %s (%s)  —  %d channel rows\n', subCL, subMW, height(sub_rows));

    % ---- Neural file (one file per subject per region) ----
    % Unique regions for this subject
    regions_sub = unique(sub_rows.SBA, 'stable');

    for r = 1:numel(regions_sub)
        rgn = regions_sub{r};

        % All rows for this subject+region
        rgn_rows  = sub_rows(strcmp(sub_rows.SBA, rgn), :);
        hemi_list = rgn_rows.HemiS;        % 'L' or 'R' per row

        % ---- Load neural data ----
        % File convention: <SubjectMW>_region-<SBA>.mat
        % Variables inside: data [n_chan x n_samples], rec_start_us (scalar µs)
        neural_file = fullfile(data_dir, subMW, ...
            sprintf('%s_region-%s.mat', subMW, rgn));
        if ~exist(neural_file, 'file')
            fprintf('  Skipping %s %s — neural file not found\n', subCL, rgn);
            continue
        end
        nd            = load(neural_file);
        neural        = nd.data;          % [n_chan x n_samples]
        rec_start_us  = nd.rec_start_us;  % scalar: recording start time in µs

        % ---- Parse contact pairs from ContactNums ----
        % ContactNums is a semicolon-separated string like '1;3' meaning
        % bipolar pair between contacts 1 and 3.
        % Each row in rgn_rows is one bipolar channel already stored in 'neural'.
        % We use the row index (within this region) as the channel index.
        n_bipolar = height(rgn_rows);
        contact_labels = rgn_rows.HnSBA;   % e.g. 'LAMY', used as channel label

        % ---- Epoch loop ----
        for e = 1:numel(epoch_types)
            etype    = epoch_types{e};
            win_sec  = windows{e};
            onsets   = onset_us.(etype);   % [n_trials x 1] µs

            n_trials  = numel(onsets);
            win_samp  = round(win_sec * fs);            % samples rel. to onset
            n_samp    = diff(win_samp) + 1;
            t_axis    = linspace(win_sec(1), win_sec(2), n_samp);

            % Preallocate  [n_bipolar x n_time x n_trials]
            epoch_data = nan(n_bipolar, n_samp, n_trials);

            for tr = 1:n_trials
                % Convert µs onset → sample index in neural recording
                onset_sec   = (onsets(tr) - rec_start_us) / 1e6;
                onset_samp  = round(onset_sec * fs) + 1;  % 1-based

                idx_start = onset_samp + win_samp(1);
                idx_end   = onset_samp + win_samp(2);

                % Boundary check
                if idx_start < 1 || idx_end > size(neural, 2)
                    fprintf('    [%s %s %s] Trial %d out of bounds — skipping\n', ...
                        subCL, rgn, etype, tr);
                    continue
                end

                epoch_data(:, :, tr) = neural(1:n_bipolar, idx_start:idx_end);
            end

            % ---- Baseline correct (z-score to pre-stimulus window) ----
            baseline_idx = t_axis >= win_sec(1) & t_axis < 0;
            if sum(baseline_idx) > 1
                epoch_data = baseline_correct(epoch_data, baseline_idx);
            else
                warning('No baseline samples found for %s %s %s — skipping BC', ...
                    subCL, rgn, etype);
            end

            % ---- Assemble output struct ----
            out.data           = epoch_data;     % [n_bipolar x n_time x n_trials]
            out.t_axis         = t_axis;
            out.fs             = fs;
            out.subject_CL     = subCL;
            out.subject_MW     = subMW;
            out.region         = rgn;
            out.hemisphere     = hemi_list;      % {n_bipolar x 1} cell
            out.contact_labels = contact_labels; % {n_bipolar x 1} cell  e.g. 'LAMY'
            out.epoch          = etype;
            out.n_trials       = n_trials;
            out.n_bipolar      = n_bipolar;
            out.beh            = beh;

            out_name = fullfile(epoch_dir, ...
                sprintf('%s_region-%s_epoch-%s.mat', subCL, rgn, etype));
            save(out_name, 'out', '-v7.3');
            fprintf('  Saved: %s\n', out_name);
        end  % epoch loop
    end  % region loop
end  % subject loop

fprintf('\nDone.\n');

%% ════════════════════════════════════════════════════════════
%  HELPER: baseline z-score correction
%% ════════════════════════════════════════════════════════════
function data_bc = baseline_correct(data, bl_idx)
    % data   : [n_chan x n_time x n_trials]
    % bl_idx : logical index into time dimension (pre-stimulus samples)
    bl_mean = mean(data(:, bl_idx, :), 2);   % [n_chan x 1 x n_trials]
    bl_std  = std( data(:, bl_idx, :), 0, 2);
    bl_std(bl_std == 0) = 1;                 % guard against flat baseline
    data_bc = (data - bl_mean) ./ bl_std;
end