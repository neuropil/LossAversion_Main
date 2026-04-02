%% step2_extract_features.m
% For each epoch file, extract band power per trial per channel
% Output: feature struct saved to /features/

clear; clc;

project_dir  = 'Z:\LossAversion\LH_Data\JAT_ADAtest';
epoch_dir    = fullfile(project_dir, 'epochs');
feat_dir     = fullfile(project_dir, 'features');
if ~exist(feat_dir, 'dir'), mkdir(feat_dir); end

fs = 500;

% Frequency bands of interest
bands.theta  = [4   8];
bands.alpha  = [8  12];
bands.beta   = [13 30];
bands.hfa    = [70 150];   % high frequency activity
band_names   = fieldnames(bands);

% Time windows for feature extraction (seconds), per epoch type
% These are the POST-STIMULUS windows you care about
feat_win.gamble  = [0.0  1.0];   % gamble viewing period
feat_win.cue     = [0.0  0.8];   % pre-response window
feat_win.outcome = [0.0  1.0];   % outcome processing

epoch_files = dir(fullfile(epoch_dir, '*.mat'));

for f = 1:numel(epoch_files)
    load(fullfile(epoch_dir, epoch_files(f).name), 'out');

    etype     = out.epoch;
    t         = out.t_axis;
    data      = out.data;         % [chan x time x trials]
    n_chan    = size(data, 1);
    n_trials  = size(data, 3);

    % Time index for feature window
    fw        = feat_win.(etype);
    t_idx     = t >= fw(1) & t <= fw(2);

    feat      = struct();
    feat.subject  = out.subject_CL;
    feat.region   = out.region;
    feat.epoch    = etype;
    feat.n_trials = n_trials;
    feat.n_chan   = n_chan;
    feat.beh      = out.beh;

    % --- Extract power per band ---
    % Result: [n_chan x n_trials] for each band
    for b = 1:numel(band_names)
        bname    = band_names{b};
        brange   = bands.(bname);

        band_power = nan(n_chan, n_trials);
        for tr = 1:n_trials
            seg = data(:, t_idx, tr);          % [chan x t_win]
            for ch = 1:n_chan
                band_power(ch, tr) = mean_band_power(seg(ch,:), fs, brange);
            end
        end
        feat.(bname) = band_power;             % [n_chan x n_trials]
    end

    % --- Also store channel-mean for quick modeling ---
    for b = 1:numel(band_names)
        bname = band_names{b};
        feat.([bname '_mean']) = mean(feat.(bname), 1)'; % [n_trials x 1]
    end

    out_name = fullfile(feat_dir, ...
        sprintf('%s_%s_%s_features.mat', ...
        out.subject_CL, out.region, etype));
    save(out_name, 'feat', '-v7.3');
    fprintf('Features saved: %s\n', out_name);
end

%% ---- HELPER: mean band power via pwelch ----
function mp = mean_band_power(sig, fs, band)
    % sig: [1 x time]
    [pxx, f] = pwelch(sig, hamming(round(fs*0.25)), [], [], fs);
    f_idx    = f >= band(1) & f <= band(2);
    mp       = mean(pxx(f_idx));
end