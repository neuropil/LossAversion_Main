%% ========================================================================
%  FULL LFP PIPELINE FOR LOSS AVERSION DATA
%  Adapted from Kirkby et al. methods for block-structured task data
%
%  Pipeline stages:
%   1.  Load data & define inputs
%   2.  Preprocessing (filter, downsample, re-reference)
%   3.  ICA-based artifact removal
%   4.  Z-score across full session
%   5.  Segment into 10s epochs & assign block/trial labels
%   6.  Compute coherence matrices per segment
%   7.  PCA + ICA to extract ICNs (full session)
%   8.  Summarize ICN activity per block
%   9.  PRIMARY ANALYSIS: Ridge regression → predict lambda per block (n=5)
%   10. SECONDARY ANALYSIS: Logistic regression → predict trial choice (n=135)
%   11. Visualization
%   12. Save results
%
%  Behavioral inputs required:
%   PRIMARY:   lambda_blocks [n_blocks × 1]  loss aversion coeff per block
%   SECONDARY: trial_onsets  [n_trials × 1]  trial onset in samples (at fs_orig)
%              trial_choices [n_trials × 1]  binary accept(1)/reject(0) per trial
%
%  Assumptions:
%   - Continuous LFP recording at original sampling rate
%   - Event markers for block start/end times
%   - L and R hemisphere electrodes kept as separate nodes
%   - Brain regions: STC, AMY, HPC, OFC, CIN, INS (bilateral)
%   - FastICA toolbox installed (https://research.ics.aalto.fi/ica/fastica/)
% =========================================================================

clear; clc;

%% ========================================================================
%  PARAMETERS
% =========================================================================

% Sampling rates
fs_orig = 2000;               % original acquisition rate (Hz) — adjust to your system
fs      = 512;                % target rate after downsampling (Hz)

% Filtering
bp_low      = 0.5;            % bandpass low cutoff (Hz)
bp_high     = 256;            % bandpass high cutoff (Hz)
bp_ord      = 8;              % Chebyshev Type I filter order
notch_freqs = [60,120,180,240]; % line noise frequencies (Hz)
notch_bw    = 4;              % notch bandwidth (Hz)
notch_ord   = 5;              % Butterworth notch filter order

% Segmentation
seg_dur = 10;                 % segment duration (seconds)
seg_len = seg_dur * fs;       % segment length (samples at fs)

% Frequency bands [low, high] Hz
freq_bands = [4  8;           % theta
              8  13;          % alpha
              13 30;          % beta
              30 70];         % gamma
band_names = {'theta','alpha','beta','gamma'};
n_bands    = size(freq_bands, 1);

% Task structure
n_blocks       = 5;           % number of task blocks
n_trials_block = 27;          % trials per block
n_trials_total = n_blocks * n_trials_block; % 135 total trials

% Coherence
n_surrogates = 20;            % phase-randomized surrogates for noise floor

% Regression — shared
n_lambda_vals = 100;          % regularization path length
alpha_en      = 0.1;          % elastic net alpha (0=ridge, 1=lasso)
n_shuffles    = 100;          % null distribution permutations

% Secondary analysis: neural feature window around trial onset
% ICN state averaged over a window preceding each trial
pre_trial_win_sec  = 5;       % seconds before trial onset
post_trial_win_sec = 0;       % seconds after trial onset (onset only)
% Together these define which 10s segments contribute to each trial's
% neural feature vector

%% ========================================================================
%  ELECTRODE METADATA
%  Populate with your actual implantation information
% =========================================================================

% chans.label  : {n_chans × 1} channel names  e.g. {'LAM1','LAM2','RAM1',...}
% chans.region : {n_chans × 1} region labels  e.g. {'AMY','AMY','AMY',...}
% chans.hemi   : {n_chans × 1} hemisphere     e.g. {'L','L','R',...}
% chans.lead   : [n_chans × 1] lead index     e.g. [1,1,1,2,2,2,...]
%
% L and R kept as distinct nodes — preserves ipsilateral, contralateral,
% and interhemispheric coherence pairs (all potentially relevant for
% loss aversion)

chans.label  = {};
chans.region = {};
chans.hemi   = {};
chans.lead   = [];

n_chans = length(chans.label);

%% ========================================================================
%  STAGE 1: LOAD DATA
% =========================================================================

fprintf('=== STAGE 1: Loading data ===\n');

% ---- Neural data --------------------------------------------------------
% raw_lfp: [n_chans × n_timepoints] continuous recording at fs_orig
% raw_lfp = your_load_function('lfp_data.mat');

% ---- Block markers ------------------------------------------------------
% block_times: [n_blocks × 2] [start, end] in samples at fs_orig
%   OR in seconds — convert below if needed
% block_times = your_load_function('block_markers.mat');
% If in seconds: block_times = round(block_times_sec .* fs_orig);

% ---- PRIMARY behavioral input -------------------------------------------
% lambda_blocks: [n_blocks × 1]
%   Loss aversion coefficient (lambda) fit across all 27 trials per block
%   using prospect theory utility function:
%     U(gamble) = p * V(gain) + (1-p) * lambda * V(loss)
%   Fit using logistic regression or maximum likelihood per block
% lambda_blocks = your_load_function('lambda_blocks.mat');

% ---- SECONDARY behavioral inputs ----------------------------------------
% trial_onsets:  [n_trials_total × 1] trial onset times in samples at fs_orig
% trial_choices: [n_trials_total × 1] binary: 1=accept gamble, 0=reject
% trial_block:   [n_trials_total × 1] block membership (1–5) per trial
%
% trial_onsets  = your_load_function('trial_onsets.mat');
% trial_choices = your_load_function('trial_choices.mat');
% trial_block   = your_load_function('trial_block.mat');

fprintf('  Data loaded. Channels: %d, Timepoints: %d\n', ...
        size(raw_lfp,1), size(raw_lfp,2));

%% ========================================================================
%  STAGE 2: PREPROCESSING
% =========================================================================

fprintf('=== STAGE 2: Preprocessing ===\n');

%% 2a. Bandpass filter — Chebyshev Type I, 8th order, 0.5–256 Hz
fprintf('  Bandpass filtering [%.1f–%.1f Hz]...\n', bp_low, bp_high);
[b_bp, a_bp] = cheby1(bp_ord, 0.5, ...
                       [bp_low bp_high]/(fs_orig/2), 'bandpass');
lfp_filt = filtfilt(b_bp, a_bp, raw_lfp')';

%% 2b. Downsample to 512 Hz
fprintf('  Downsampling %d → %d Hz...\n', fs_orig, fs);
ds_factor = fs_orig / fs;
if mod(ds_factor, 1) ~= 0
    error(['Downsampling factor (%.2f) must be integer. ' ...
           'Adjust fs_orig or fs.'], ds_factor);
end
lfp_ds       = downsample(lfp_filt', ds_factor)';
n_timepoints = size(lfp_ds, 2);

%% 2c. Notch filters — Butterworth, 5th order, ±2 Hz around line frequencies
fprintf('  Notch filtering at [%s] Hz...\n', num2str(notch_freqs));
lfp_notch = lfp_ds;
for nf = 1:length(notch_freqs)
    f0    = notch_freqs(nf);
    f_lo  = (f0 - notch_bw/2) / (fs/2);
    f_hi  = (f0 + notch_bw/2) / (fs/2);
    [b_n, a_n] = butter(notch_ord, [f_lo f_hi], 'stop');
    lfp_notch = filtfilt(b_n, a_n, lfp_notch')';
end

%% 2d. Re-reference to common average within each lead
fprintf('  Re-referencing within leads...\n');
lfp_reref    = lfp_notch;
unique_leads = unique(chans.lead);
for l = 1:length(unique_leads)
    idx           = chans.lead == unique_leads(l);
    lfp_reref(idx,:) = lfp_reref(idx,:) - mean(lfp_reref(idx,:), 1);
end

%% ========================================================================
%  STAGE 3: ICA-BASED ARTIFACT REMOVAL
%  Requires FastICA toolbox: https://research.ics.aalto.fi/ica/fastica/
%  For 15-min sessions treat as single segment (paper used ~1hr segments)
% =========================================================================

fprintf('=== STAGE 3: ICA artifact removal ===\n');
lfp_clean = ica_artifact_removal(lfp_reref);

%% ========================================================================
%  STAGE 4: Z-SCORE ACROSS FULL SESSION
%  Must be done BEFORE segmentation to preserve between-block variance
%  Using within-block z-scoring would normalize out the signal of interest
% =========================================================================

fprintf('=== STAGE 4: Z-scoring across full session ===\n');
lfp_z = zscore(lfp_clean, 0, 2);   % [n_chans × n_timepoints]
fprintf('  Done. Shape: [%d × %d]\n', size(lfp_z,1), size(lfp_z,2));

%% ========================================================================
%  STAGE 5: SEGMENTATION & LABELING
%  Segment continuously — trial boundaries ignored (10s >> 6s trial)
%  Label each segment by block membership and nearest trial
% =========================================================================

fprintf('=== STAGE 5: Segmentation and labeling ===\n');

% Convert block marker times to downsampled sample indices
block_samples_ds = round(block_times ./ ds_factor);   % [n_blocks × 2]

% Convert trial onset times to downsampled indices
trial_onsets_ds = round(trial_onsets ./ ds_factor);   % [n_trials × 1]

% Build contiguous 10s segments
n_segs     = floor(n_timepoints / seg_len);
seg_onset  = ((0:n_segs-1)' .* seg_len) + 1;          % [n_segs × 1]
seg_offset = seg_onset + seg_len - 1;
seg_mid    = seg_onset + floor(seg_len/2);

% Assign block label to each segment (based on midpoint)
seg_block  = zeros(n_segs, 1);   % 0 = inter-block
for b = 1:n_blocks
    in_block = seg_mid >= block_samples_ds(b,1) & ...
               seg_mid <= block_samples_ds(b,2);
    seg_block(in_block) = b;
end

% For each trial, find which 10s segment contains its onset
% Used in secondary analysis to align neural features to choices
trial_seg_idx = zeros(n_trials_total, 1);
for t = 1:n_trials_total
    % Find segment whose window covers the trial onset
    % Use segment starting up to pre_trial_win_sec before onset
    earliest = trial_onsets_ds(t) - pre_trial_win_sec * fs;
    seg_candidates = find(seg_onset >= earliest & ...
                          seg_onset <= trial_onsets_ds(t));
    if ~isempty(seg_candidates)
        % Use the segment immediately preceding trial onset
        trial_seg_idx(t) = seg_candidates(end);
    else
        % Fallback: segment containing trial onset
        trial_seg_idx(t) = find(seg_onset <= trial_onsets_ds(t) & ...
                                seg_offset >= trial_onsets_ds(t), 1, 'last');
    end
end

% Report
fprintf('  Total segments: %d (%.1f min)\n', n_segs, n_segs*seg_dur/60);
for b = 1:n_blocks
    ns = sum(seg_block==b);
    fprintf('  Block %d: %d segments (%.1f min)\n', b, ns, ns*seg_dur/60);
end
fprintf('  Inter-block: %d segments\n', sum(seg_block==0));
fprintf('  Trials with valid segment alignment: %d / %d\n', ...
        sum(trial_seg_idx>0), n_trials_total);

%% ========================================================================
%  STAGE 6: COHERENCE MATRICES
%  Computed per 10s segment, all electrode pairs, 4 frequency bands
%  Phase-randomized surrogate subtraction removes amplitude co-variation
%  noise floor (Srinath & Ray 2014)
% =========================================================================

fprintf('=== STAGE 6: Computing coherence matrices ===\n');

pairs   = nchoosek(1:n_chans, 2);
n_pairs = size(pairs, 1);
fprintf('  %d pairs × %d segments × %d bands\n', n_pairs, n_segs, n_bands);

% Store coherence: [n_pairs × n_segs × n_bands]
coh_matrices = zeros(n_pairs, n_segs, n_bands, 'single');

win = hann(seg_len);   % non-overlapping Hanning window = full segment

for s = 1:n_segs
    if mod(s, 20) == 0
        fprintf('  Segment %d / %d\n', s, n_segs);
    end
    seg_data = lfp_z(:, seg_onset(s):seg_offset(s));   % [n_chans × seg_len]
    
    for p = 1:n_pairs
        x = seg_data(pairs(p,1), :)';
        y = seg_data(pairs(p,2), :)';
        
        % Original coherence
        [Cxy, f] = mscohere(x, y, win, 0, [], fs);
        
        % Noise floor via phase-randomized surrogates
        Cxy_surr = zeros(length(f), n_surrogates);
        for sur = 1:n_surrogates
            Cxy_surr(:,sur) = mscohere(phase_randomize(x), ...
                                        phase_randomize(y), win, 0, [], fs);
        end
        Cxy_corr = max(Cxy - mean(Cxy_surr, 2), 0);   % floor at zero
        
        % Band-average
        for band = 1:n_bands
            f_idx = f >= freq_bands(band,1) & f < freq_bands(band,2);
            coh_matrices(p, s, band) = mean(Cxy_corr(f_idx));
        end
    end
end
fprintf('  Coherence computation complete.\n');

%% ========================================================================
%  STAGE 7: ICN EXTRACTION — PCA + ICA ON FULL SESSION
%  Use ALL segments (including inter-block) for stable decomposition
%  Marchenko-Pastur threshold determines number of significant PCs
%  FastICA then separates into maximally independent networks
% =========================================================================

fprintf('=== STAGE 7: ICN extraction (PCA + ICA) ===\n');

icn_weights = cell(n_bands, 1);   % spatial maps: {band}[n_pairs × n_ics]
icn_proj    = cell(n_bands, 1);   % projections: {band}[n_segs × n_ics]
n_ics       = zeros(n_bands, 1);

for band = 1:n_bands
    fprintf('  [%s band]\n', band_names{band});
    
    M = double(coh_matrices(:,:,band));   % [n_pairs × n_segs]
    
    %% PCA — dimensionality reduction and orthogonalization
    [coeff, score, latent] = pca(M');
    % coeff:  [n_pairs × n_pairs]  loadings
    % score:  [n_segs  × n_pairs]  scores
    % latent: [n_pairs × 1]        eigenvalues
    
    % Marchenko-Pastur significance threshold
    q        = size(M,2) / size(M,1);     % n_segs / n_pairs
    sigma2   = mean(latent);
    lambda_mp = sigma2 * (1 + sqrt(1/q))^2;
    n_sig    = max(sum(latent > lambda_mp), 2);   % at least 2
    fprintf('    Significant PCs: %d (MP threshold: %.4f)\n', ...
            n_sig, lambda_mp);
    
    sig_scores = score(:, 1:n_sig);    % [n_segs × n_sig]
    sig_coeff  = coeff(:, 1:n_sig);    % [n_pairs × n_sig]
    
    %% ICA — maximally independent sources from significant PCs
    % fastica() input: [n_features × n_observations]
    try
        [icaweights, ~, ~] = fastica(sig_scores', ...
                                      'approach', 'symm', ...
                                      'g',        'tanh', ...
                                      'verbose',  'off');
    catch ME
        error(['FastICA failed: %s\n' ...
               'Install toolbox from: ' ...
               'https://research.ics.aalto.fi/ica/fastica/'], ME.message);
    end
    
    n_ic = size(icaweights, 1);
    n_ics(band) = n_ic;
    fprintf('    ICNs extracted: %d\n', n_ic);
    
    % Map ICA weights back to electrode-pair space
    icn_map             = sig_coeff * icaweights';   % [n_pairs × n_ics]
    icn_weights{band}   = icn_map;
    
    % Project full session coherence onto ICNs: [n_segs × n_ics]
    icn_proj{band}      = (icn_map' * M)';
end

%% ========================================================================
%  STAGE 8: SUMMARIZE ICN ACTIVITY PER BLOCK
%  Variance of ICN projection within each block captures transient
%  peaks in ICN activity (following paper's 1-min sliding window approach,
%  adapted to block structure)
% =========================================================================

fprintf('=== STAGE 8: ICN activity per block ===\n');

max_ics       = max(n_ics);
icn_block_var = zeros(n_blocks, max_ics, n_bands);

for band = 1:n_bands
    proj = icn_proj{band};   % [n_segs × n_ics]
    for b = 1:n_blocks
        idx = seg_block == b;
        if sum(idx) < 2
            warning('Block %d, band %s: fewer than 2 segments', ...
                    b, band_names{band});
            continue;
        end
        n_ic = n_ics(band);
        icn_block_var(b, 1:n_ic, band) = var(proj(idx,:), 0, 1);
    end
end

%% ========================================================================
%  STAGE 9: PRIMARY ANALYSIS
%  Ridge regression: predict lambda per block (n=5)
%  LOO cross-validation appropriate for n=5
%  Null distribution from 100 row-permuted shuffles
% =========================================================================

fprintf('=== STAGE 9: PRIMARY — Block-level ridge regression ===\n');

%% Build feature matrix [n_blocks × n_features]
X_block = reshape(icn_block_var, n_blocks, []);

% Remove zero-variance features (absent regions / unused ICNs)
active_cols   = var(X_block, 0, 1) > 0;
X_block       = X_block(:, active_cols);
fprintf('  Feature matrix: [%d blocks × %d features]\n', ...
        size(X_block,1), size(X_block,2));

% Standardize features (zero mean, unit variance)
[X_block_norm, X_block_mu, X_block_sig] = zscore(X_block);

%% Leave-one-out cross-validation partition
cv_loo = cvpartition(n_blocks, 'LeaveOut');

%% Elastic net (alpha=0.1 ≈ ridge) with LOO CV
[B_block, FI_block] = lasso(X_block_norm, lambda_blocks, ...
                             'Alpha',     alpha_en, ...
                             'CV',        cv_loo, ...
                             'NumLambda', n_lambda_vals);

idx_min       = FI_block.IndexMinMSE;
B_best_block  = B_block(:, idx_min);
int_best_block = FI_block.Intercept(idx_min);

% In-sample R²
y_pred_block = X_block_norm * B_best_block + int_best_block;
SS_res = sum((lambda_blocks - y_pred_block).^2);
SS_tot = sum((lambda_blocks - mean(lambda_blocks)).^2);
R2_block = 1 - SS_res / SS_tot;
fprintf('  R² = %.4f\n', R2_block);
fprintf('  Non-zero coefficients: %d / %d\n', ...
        sum(B_best_block ~= 0), length(B_best_block));

%% Null distribution: permute block labels 100x
fprintf('  Computing null distribution (%d shuffles)...\n', n_shuffles);
R2_null_block = zeros(n_shuffles, 1);
for shuf = 1:n_shuffles
    perm_idx  = randperm(n_blocks);
    X_shuf    = X_block_norm(perm_idx, :);
    [B_s, FI_s] = lasso(X_shuf, lambda_blocks, ...
                         'Alpha',     alpha_en, ...
                         'CV',        cv_loo, ...
                         'NumLambda', n_lambda_vals);
    i_s  = FI_s.IndexMinMSE;
    yhat = X_shuf * B_s(:,i_s) + FI_s.Intercept(i_s);
    R2_null_block(shuf) = 1 - sum((lambda_blocks - yhat).^2) / SS_tot;
end
p_block = mean(R2_null_block >= R2_block);
fprintf('  p-value vs null: %.4f\n', p_block);

%% ========================================================================
%  STAGE 10: SECONDARY ANALYSIS
%  Logistic regression: predict trial-level accept/reject (n=135)
%  Uses ICN state in window preceding each trial onset
%  5-fold stratified CV (preserves class balance across folds)
%  Null distribution from 100 trial-label permutations
%
%  Neural feature per trial: mean ICN projection variance in the
%  pre_trial_win_sec window of segments preceding trial onset
% =========================================================================

fprintf('=== STAGE 10: SECONDARY — Trial-level logistic regression ===\n');

%% Build trial-level feature matrix [n_trials × n_features]
% For each trial, average ICN projection over the preceding window of segments

max_ics_total = sum(n_ics);
X_trial = zeros(n_trials_total, max_ics_total * n_bands);

col = 1;
for band = 1:n_bands
    proj = icn_proj{band};   % [n_segs × n_ics]
    n_ic = n_ics(band);
    
    for t = 1:n_trials_total
        s_idx = trial_seg_idx(t);
        if s_idx == 0
            % No valid segment found — leave as zeros
            col_end = col + n_ic - 1;
            X_trial(t, col:col_end) = 0;
            continue;
        end
        
        % Collect segments within pre-trial window
        win_start = trial_onsets_ds(t) - pre_trial_win_sec * fs;
        win_segs  = find(seg_onset >= win_start & ...
                         seg_onset <= trial_onsets_ds(t) & ...
                         seg_block == trial_block(t));
        
        if isempty(win_segs)
            win_segs = s_idx;   % fallback to single segment
        end
        
        % Mean ICN projection variance over window segments
        col_end = col + n_ic - 1;
        X_trial(t, col:col_end) = mean(var(proj(win_segs,:), 0, 1), 1);
    end
    col = col + n_ic;
end

% Remove zero-variance features
active_trial = var(X_trial, 0, 1) > 0;
X_trial      = X_trial(:, active_trial);
fprintf('  Feature matrix: [%d trials × %d features]\n', ...
        size(X_trial,1), size(X_trial,2));

% Standardize
[X_trial_norm, ~, ~] = zscore(X_trial);

%% 5-fold stratified CV — preserves accept/reject balance across folds
cv_trial = cvpartition(trial_choices, 'KFold', 5, 'Stratify', true);

%% Logistic elastic net
[B_trial, FI_trial] = lassoglm(X_trial_norm, trial_choices, ...
                                 'binomial', ...
                                 'Alpha',     alpha_en, ...
                                 'CV',        cv_trial, ...
                                 'NumLambda', n_lambda_vals);

idx_trial     = FI_trial.IndexMinDeviance;
B_best_trial  = B_trial(:, idx_trial);
int_best_trial = FI_trial.Intercept(idx_trial);

% Predicted probabilities and accuracy
log_odds     = X_trial_norm * B_best_trial + int_best_trial;
p_accept     = 1 ./ (1 + exp(-log_odds));
y_pred_trial = double(p_accept >= 0.5);
acc_trial    = mean(y_pred_trial == trial_choices);
fprintf('  Accuracy = %.3f (chance = %.3f)\n', acc_trial, ...
        max(mean(trial_choices), 1-mean(trial_choices)));
fprintf('  Non-zero coefficients: %d / %d\n', ...
        sum(B_best_trial ~= 0), length(B_best_trial));

% AUC
[~, ~, ~, AUC] = perfcurve(trial_choices, p_accept, 1);
fprintf('  AUC = %.3f\n', AUC);

%% Null distribution: permute trial choice labels 100x
fprintf('  Computing null distribution (%d shuffles)...\n', n_shuffles);
acc_null  = zeros(n_shuffles, 1);
AUC_null  = zeros(n_shuffles, 1);
for shuf = 1:n_shuffles
    y_shuf = trial_choices(randperm(n_trials_total));
    cv_s   = cvpartition(y_shuf, 'KFold', 5, 'Stratify', true);
    [B_s, FI_s] = lassoglm(X_trial_norm, y_shuf, 'binomial', ...
                             'Alpha',     alpha_en, ...
                             'CV',        cv_s, ...
                             'NumLambda', n_lambda_vals);
    i_s    = FI_s.IndexMinDeviance;
    lo_s   = X_trial_norm * B_s(:,i_s) + FI_s.Intercept(i_s);
    p_s    = 1 ./ (1 + exp(-lo_s));
    acc_null(shuf) = mean(double(p_s >= 0.5) == y_shuf);
    [~,~,~, AUC_null(shuf)] = perfcurve(y_shuf, p_s, 1);
end
p_trial_acc = mean(acc_null  >= acc_trial);
p_trial_auc = mean(AUC_null  >= AUC);
fprintf('  p-value (accuracy) vs null: %.4f\n', p_trial_acc);
fprintf('  p-value (AUC)      vs null: %.4f\n', p_trial_auc);

%% ========================================================================
%  STAGE 11: VISUALIZATION
% =========================================================================

fprintf('=== STAGE 11: Visualization ===\n');

%% Figure 1: Mean coherence matrices per band
figure('Name','Coherence Matrices','Position',[50 50 1400 320]);
for band = 1:n_bands
    subplot(1, n_bands, band);
    C_sq = squareform(mean(coh_matrices(:,:,band), 2, 'omitnan'));
    imagesc(C_sq); axis square; colorbar;
    title(band_names{band}); xlabel('Channel'); ylabel('Channel');
end
sgtitle('Mean Band Coherence (Full Session)');

%% Figure 2: ICN block variance heatmaps
figure('Name','ICN Block Variance','Position',[50 430 1400 320]);
for band = 1:n_bands
    subplot(1, n_bands, band);
    imagesc(squeeze(icn_block_var(:, 1:n_ics(band), band))');
    xlabel('Block'); ylabel('ICN index');
    title(band_names{band}); colorbar;
    xticks(1:n_blocks); xticklabels(arrayfun(@(x) sprintf('B%d',x), ...
           1:n_blocks, 'UniformOutput', false));
end
sgtitle('ICN Projection Variance per Block');

%% Figure 3: PRIMARY — block regression results
figure('Name','Primary: Block Regression','Position',[50 810 800 380]);
subplot(1,2,1);
scatter(lambda_blocks, y_pred_block, 100, 'filled', 'MarkerFaceColor',[0.2 0.5 0.9]);
hold on;
xlims = xlim; plot(xlims, xlims, 'k--');
xlabel('Observed \lambda'); ylabel('Predicted \lambda');
title(sprintf('Block regression  R² = %.3f', R2_block)); axis square;

subplot(1,2,2);
histogram(R2_null_block, 20, 'FaceColor',[0.7 0.7 0.7]); hold on;
xline(R2_block, 'r-', 'LineWidth', 2);
xlabel('R² (null distribution)'); ylabel('Count');
title(sprintf('p = %.3f (n=%d shuffles)', p_block, n_shuffles));
sgtitle('PRIMARY: Lambda ~ ICN variance (block level)');

%% Figure 4: SECONDARY — trial logistic regression results
figure('Name','Secondary: Trial Logistic Regression','Position',[900 810 900 380]);

subplot(1,3,1);
% ROC curve
[tpr, fpr, ~, ~] = perfcurve(trial_choices, p_accept, 1);
plot(fpr, tpr, 'b-', 'LineWidth', 2); hold on;
plot([0 1],[0 1],'k--');
xlabel('False Positive Rate'); ylabel('True Positive Rate');
title(sprintf('ROC  AUC = %.3f', AUC)); axis square;

subplot(1,3,2);
histogram(acc_null, 20, 'FaceColor',[0.7 0.7 0.7]); hold on;
xline(acc_trial, 'r-', 'LineWidth', 2);
xlabel('Accuracy (null)'); ylabel('Count');
title(sprintf('Accuracy: %.3f  p = %.3f', acc_trial, p_trial_acc));

subplot(1,3,3);
histogram(AUC_null, 20, 'FaceColor',[0.7 0.7 0.7]); hold on;
xline(AUC, 'r-', 'LineWidth', 2);
xlabel('AUC (null)'); ylabel('Count');
title(sprintf('AUC: %.3f  p = %.3f', AUC, p_trial_auc));

sgtitle('SECONDARY: Choice ~ ICN state (trial level)');

%% Figure 5: Coefficient profiles — which ICNs drive predictions
figure('Name','Predictive ICN Coefficients','Position',[50 1250 1200 380]);

subplot(1,2,1);
bar(B_best_block(B_best_block ~= 0), 'FaceColor',[0.2 0.5 0.9]);
xlabel('Feature index (non-zero only)'); ylabel('\beta coefficient');
title('Primary: Block regression coefficients');

subplot(1,2,2);
bar(B_best_trial(B_best_trial ~= 0), 'FaceColor',[0.9 0.4 0.2]);
xlabel('Feature index (non-zero only)'); ylabel('\beta coefficient');
title('Secondary: Trial logistic coefficients');

sgtitle('Predictive ICN Features');

%% ========================================================================
%  STAGE 12: SAVE RESULTS
% =========================================================================

fprintf('=== STAGE 12: Saving results ===\n');

results = struct();

% Neural
results.coh_matrices    = coh_matrices;
results.seg_block       = seg_block;
results.seg_onset       = seg_onset;
results.icn_weights     = icn_weights;
results.icn_proj        = icn_proj;
results.icn_block_var   = icn_block_var;
results.n_ics           = n_ics;
results.freq_bands      = freq_bands;
results.band_names      = band_names;

% Primary analysis
results.primary.B            = B_best_block;
results.primary.intercept    = int_best_block;
results.primary.y_pred       = y_pred_block;
results.primary.R2           = R2_block;
results.primary.R2_null      = R2_null_block;
results.primary.p_val        = p_block;

% Secondary analysis
results.secondary.B          = B_best_trial;
results.secondary.intercept  = int_best_trial;
results.secondary.p_accept   = p_accept;
results.secondary.accuracy   = acc_trial;
results.secondary.AUC        = AUC;
results.secondary.acc_null   = acc_null;
results.secondary.AUC_null   = AUC_null;
results.secondary.p_acc      = p_trial_acc;
results.secondary.p_auc      = p_trial_auc;

save('LFP_pipeline_v2_results.mat', 'results', '-v7.3');
fprintf('Results saved to LFP_pipeline_v2_results.mat\n');
fprintf('Pipeline complete.\n');

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function lfp_clean = ica_artifact_removal(lfp)
%ICA_ARTIFACT_REMOVAL  FastICA-based artifact removal
%
%  Approximates the paper's approach:
%   - FastICA decomposition (symmetric, tanh nonlinearity)
%   - Artifact ICs identified by high kurtosis (proxy for large deflections)
%   - Top ~10% of ICs by kurtosis removed (matching paper's ~10% rate)
%   - Remaining ICs reconstructed to signal space
%
%  For production: replace kurtosis flagging with a trained logistic
%  classifier using power spectrum + amplitude distribution as features
%  (fitclinear or mnrfit), trained on manually labeled ICs

    fprintf('  Running FastICA...\n');
    lfp_norm = zscore(lfp, 0, 2);
    
    try
        [icasig, A, ~] = fastica(lfp_norm, ...
                                  'approach', 'symm', ...
                                  'g',        'tanh', ...
                                  'verbose',  'off');
    catch
        warning(['FastICA not found — skipping artifact removal.\n' ...
                 'Install from: https://research.ics.aalto.fi/ica/fastica/']);
        lfp_clean = lfp;
        return;
    end
    
    % Flag top 10% by kurtosis as artifact
    ic_kurt      = kurtosis(icasig, 0, 2);
    artifact_idx = ic_kurt > prctile(ic_kurt, 90);
    fprintf('    Artifact ICs removed: %d / %d\n', ...
            sum(artifact_idx), size(icasig,1));
    
    icasig_clean             = icasig;
    icasig_clean(artifact_idx,:) = 0;
    lfp_clean                = A * icasig_clean;
end


function x_surr = phase_randomize(x)
%PHASE_RANDOMIZE  Phase-randomized surrogate preserving power spectrum
%   Randomizes Fourier phases while maintaining conjugate symmetry
%   so that ifft returns a real-valued signal

    N  = length(x);
    X  = fft(x);
    
    if mod(N, 2) == 0
        n_r = N/2 - 1;
        rp  = 2*pi * rand(n_r, 1);
        ph  = [0; rp; 0; -flipud(rp)];
    else
        n_r = (N-1)/2;
        rp  = 2*pi * rand(n_r, 1);
        ph  = [0; rp; -flipud(rp)];
    end
    
    x_surr = real(ifft(abs(X) .* exp(1i * ph)));
end
