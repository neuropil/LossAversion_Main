%% step5_time_resolved_model.m
% Time-resolved sliding window: choice ~ gainZ + lossZ + neuralZ
% Runs per region × epoch × band, saves beta curves with 95% CI
% Requires: epoched data from step2

clear; clc;

project_dir  = 'Z:\LossAversion\LH_Data\JAT_ADAtest';
epoch_dir    = fullfile(project_dir, 'epochs');
results_dir  = fullfile(project_dir, 'results');

fs           = 500;
win_ms       = 200;
step_ms      = 20;
win_samp     = round(win_ms/1000 * fs);
step_samp    = round(step_ms/1000 * fs);

bands.theta  = [4   8];
bands.beta   = [13 30];
bands.hfa    = [70 150];
band_names   = fieldnames(bands);

band_min_ms.theta = 250;   % minimum window for reliable theta estimate
band_min_ms.beta  = 100;
band_min_ms.hfa   = 100;

regions      = {'LAMY','LMOF','RAMY','LLOF','RLOF','RMOF','RAH','LAH','RPH','LPH'};
epoch_types  = {'gamble','cue','outcome'};

for r = 1:numel(regions)
    rgn = regions{r};

    for e = 1:numel(epoch_types)
        ep = epoch_types{e};

        epoch_file = fullfile(epoch_dir, ...
            sprintf('*_%s_%s.mat', rgn, ep));
        files = dir(epoch_file);
        if isempty(files)
            fprintf('No files found: %s | %s\n', rgn, ep);
            continue
        end

        %% ── Stack neural data and behavioral tables across subjects ───────
        all_data = [];
        all_beh  = [];
        t_axis   = [];

        for f = 1:numel(files)
            fpath = fullfile(epoch_dir, files(f).name);
            load(fpath, 'out');

            % Validate expected fields
            if ~isfield(out, 'data') || ~isfield(out, 'beh') || ~isfield(out, 't_axis')
                fprintf('WARNING: missing fields in %s — skipping\n', files(f).name);
                continue
            end

            % Average over channels → [1 x time x trials]
            % Necessary because subjects have different numbers of bipolar pairs
            data_avg = mean(out.data, 1);

            if isempty(all_data)
                t_axis   = out.t_axis;
                all_data = data_avg;
            else
                % Verify time axis consistency across subjects
                if size(data_avg, 2) ~= size(all_data, 2)
                    fprintf('WARNING: %s time axis mismatch (%d vs %d samples) — skipping\n', ...
                        files(f).name, size(data_avg,2), size(all_data,2));
                    continue
                end
                all_data = cat(3, all_data, data_avg);
            end

            all_beh = stack_beh(all_beh, out.beh);
        end

        if isempty(all_data) || isempty(all_beh)
            fprintf('No data loaded for %s | %s — skipping\n', rgn, ep);
            continue
        end

        %% ── Validate trial count consistency ─────────────────────────────
        n_trials_neu = size(all_data, 3);
        n_trials_beh = height(all_beh);

        if n_trials_neu ~= n_trials_beh
            fprintf('WARNING: %s %s — neural trials (%d) ≠ beh trials (%d), skipping\n', ...
                rgn, ep, n_trials_neu, n_trials_beh);
            continue
        end
        n_trials = n_trials_neu;
        n_time   = size(all_data, 2);

        %% ── Build base behavioral table with within-subject z-scores ──────
        T_base         = table();
        T_base.choice  = double(all_beh.choice(:));
        T_base.subject = categorical(all_beh.subjectIndex(:));
        T_base.gainZ = zeros(n_trials, 1);
        T_base.lossZ = zeros(n_trials, 1);

        subjectSSs = unique(T_base.subject);
        n_subjects = numel(subjectSSs);

        for s = 1:numel(unique(T_base.subject))
            idx = T_base.subject == subjectSSs(s);

            gain_vals = double(all_beh.riskyGain(idx));
            loss_vals = double(all_beh.riskyLoss(idx));

            % Z-score within subject, guard against constant columns
            if std(gain_vals, 'omitnan') < 1e-6
                T_base.gainZ(idx) = 0;
            else
                T_base.gainZ(idx) = zscore(gain_vals);
            end

            if std(loss_vals, 'omitnan') < 1e-6
                T_base.lossZ(idx) = 0;
            else
                T_base.lossZ(idx) = zscore(loss_vals);
            end
        end

        % subjects   = unique(T_base.subject);
        % n_subjects = numel(subjects);
        % 
        % for s = 1:numel(subjects)
        %     idx = T_base.subject == subjects(s);
        %     T_base.gainZ(idx) = zscore(all_beh.gain(idx));
        %     T_base.lossZ(idx) = zscore(all_beh.loss(idx));
        % end

        %% ── Sliding window setup ─────────────────────────────────────────
        win_starts = 1 : step_samp : (n_time - win_samp + 1);
        n_wins     = numel(win_starts);
        t_centers  = t_axis(win_starts + floor(win_samp/2));

        % Preallocate [n_bands x n_wins]
        beta_t  = nan(numel(band_names), n_wins);
        pval_t  = nan(numel(band_names), n_wins);
        ci_lo_t = nan(numel(band_names), n_wins);
        ci_hi_t = nan(numel(band_names), n_wins);

        fprintf('Time-resolving: %s | %s  (%d windows, %d subjects, %d trials)...\n', ...
            rgn, ep, n_wins, n_subjects, n_trials);

        %% ── Sliding window loop ──────────────────────────────────────────
        for w = 1:n_wins
            win_idx = win_starts(w) : win_starts(w) + win_samp - 1;
            seg     = all_data(1, win_idx, :);   % [1 x win_samp x trials]

            for b = 1:numel(band_names)
                bname  = band_names{b};
                brange = bands.(bname);

                % Skip band if window too short for reliable estimate
                win_ms_actual = win_samp / fs * 1000;
                if win_ms_actual < band_min_ms.(bname)
                    continue
                end

                % Band power per trial — channel avg already done at load
                feat_vals = nan(n_trials, 1);
                for tr = 1:n_trials
                    seg_trial     = squeeze(seg(1, :, tr));   % [win_samp x 1]
                    feat_vals(tr) = mean_band_power(seg_trial, fs, brange);
                end

                % Within-subject z-score of neural feature
                feat_z = nan(n_trials, 1);
                for s = 1:numel(subjectSSs)
                    sidx = T_base.subject == subjectSSs(s);
                    vals = feat_vals(sidx);
                    if std(vals, 'omitnan') < 1e-6
                        feat_z(sidx) = 0;
                    else
                        feat_z(sidx) = zscore(vals);
                    end
                end

                % Skip degenerate windows
                if all(isnan(feat_z)) || std(feat_z, 'omitnan') < 1e-6
                    continue
                end

                % Assemble per-window table
                T         = T_base;
                T.neuralZ = feat_z;

                [formula, cov_pat] = build_formula('neuralZ', n_subjects, 'main');

                try
                    mdl = fitglme(T, formula, ...
                            'Distribution',     'Binomial', ...
                            'Link',             'logit', ...
                            'FitMethod',        'Laplace', ...
                            'CovariancePattern', cov_pat);

                    coef = mdl.Coefficients;
                    row  = strcmp(coef.Name, 'neuralZ');

                    beta_t(b, w)  = coef.Estimate(row);
                    pval_t(b, w)  = coef.pValue(row);
                    ci_lo_t(b, w) = coef.Estimate(row) - 1.96 * coef.SE(row);
                    ci_hi_t(b, w) = coef.Estimate(row) + 1.96 * coef.SE(row);

                catch ME
                    if w == 1
                        fprintf('  fitglme failed [%s | %s | %s | w=1]: %s\n', ...
                            rgn, ep, bname, ME.message);
                    end
                end
            end

            if mod(w, 50) == 0
                fprintf('  window %d / %d\n', w, n_wins);
            end
        end

        %% ── Save results ─────────────────────────────────────────────────
        tr_res.region     = rgn;
        tr_res.epoch      = ep;
        tr_res.t_centers  = t_centers;
        tr_res.band_names = band_names;
        tr_res.beta_t     = beta_t;
        tr_res.pval_t     = pval_t;
        tr_res.ci_lo_t    = ci_lo_t;
        tr_res.ci_hi_t    = ci_hi_t;
        tr_res.n_subjects = n_subjects;
        tr_res.n_trials   = n_trials;
        tr_res.n_wins     = n_wins;
        tr_res.win_ms     = win_ms;
        tr_res.step_ms    = step_ms;

        out_name = fullfile(results_dir, ...
            sprintf('time_resolved_%s_%s.mat', rgn, ep));
        save(out_name, 'tr_res', '-v7.3');
        fprintf('  Saved: %s\n', out_name);

        plot_beta_curve(tr_res , results_dir , rgn , ep);
    end
end

%% ═══════════════════════════════════════════════════════════════
%% HELPERS
%% ═══════════════════════════════════════════════════════════════

function [formula, cov_pat] = build_formula(feat, n_subjects, model_type)
% Mirrors step4 build_formula — keep in sync if RE threshold changes.
    if n_subjects >= 8
        re_str  = '(1 + lossZ | subject)';
        cov_pat = 'Diagonal';
    else
        re_str  = '(1 | subject)';
        cov_pat = 'FullCholesky';
    end
    switch model_type
        case 'main'
            formula = sprintf('choice ~ gainZ + lossZ + %s + %s', feat, re_str);
        otherwise
            error('build_formula: unknown model_type "%s"', model_type);
    end
end

function beh_out = stack_beh(beh_in, beh_new)
% Vertically concatenate behavioral tables across subjects.
    if isempty(beh_in)
        beh_out = beh_new;
        return
    end

    % Warn if schemas don't match
    vars_in  = beh_in.Properties.VariableNames;
    vars_new = beh_new.Properties.VariableNames;
    missing  = setdiff(vars_in,  vars_new);
    extra    = setdiff(vars_new, vars_in);
    if ~isempty(missing)
        fprintf('WARNING stack_beh: new table missing columns: %s\n', strjoin(missing, ', '));
    end
    if ~isempty(extra)
        fprintf('WARNING stack_beh: new table has extra columns: %s\n', strjoin(extra, ', '));
    end

    beh_out = vertcat(beh_in, beh_new);
end

function mp = mean_band_power(sig, fs, band)
    sig     = sig(:);                          % ensure column vector
    n_samp  = length(sig);
    
    % pwelch segment must not exceed signal length
    % target ~250ms for frequency resolution, cap at signal length
    win_len = min(round(fs * 0.25), n_samp);
    
    % Overlap: 50% of window length
    n_overlap = floor(win_len / 2);
    
    [pxx, f] = pwelch(sig, hamming(win_len), n_overlap, [], fs);
    f_idx    = f >= band(1) & f <= band(2);
    
    if ~any(f_idx)
        mp = NaN;
        return
    end
    mp = mean(pxx(f_idx));
end

function plot_beta_curve(tr_res , results_dir , rgn , ep)
    % Define colors as RGB [0-1] — hex strings not supported in all MATLAB versions
    band_colors = {[0.129 0.588 0.953], ...   % #2196F3 blue   (theta)
                   [1.000 0.341 0.133], ...   % #FF5722 orange (beta)
                   [0.298 0.686 0.314]};      % #4CAF50 green  (hfa)

    figure('Name', sprintf('Beta curve: %s | %s', tr_res.region, tr_res.epoch), ...
           'Position', [100 100 900 400]);
    hold on;

    for b = 1:numel(tr_res.band_names)
        t    = tr_res.t_centers;
        beta = tr_res.beta_t(b, :);
        lo   = tr_res.ci_lo_t(b, :);
        hi   = tr_res.ci_hi_t(b, :);
        col  = band_colors{b};

        % Skip if all NaN (e.g. theta skipped due to short window)
        if all(isnan(beta)), continue; end

        % Mask NaNs for clean plotting
        valid = ~isnan(beta) & ~isnan(lo) & ~isnan(hi);

        % Shaded CI
        t_v  = t(valid);
        lo_v = lo(valid);
        hi_v = hi(valid);
        fill([t_v fliplr(t_v)], [lo_v fliplr(hi_v)], col, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');

        % Beta curve
        plot(t_v, beta(valid), 'Color', col, 'LineWidth', 2, ...
            'DisplayName', tr_res.band_names{b});

        % Significant windows (uncorrected p<0.05 for display only)
        sig_idx = tr_res.pval_t(b, :) < 0.05 & valid;
        if any(sig_idx)
            scatter(t(sig_idx), beta(sig_idx), 30, col, 'filled', ...
                'HandleVisibility', 'off');
        end
    end

    yline(0, 'k--', 'LineWidth', 1);
    xline(0, 'k-',  'LineWidth', 1.5);
    xlabel('Time relative to epoch onset (s)');
    ylabel('\beta  (neuralZ \rightarrow choice)');
    title(sprintf('%s  |  %s epoch  |  n=%d subjects', ...
        tr_res.region, tr_res.epoch, tr_res.n_subjects));
    legend('Location', 'best');
    box off;
    hold off;

    fig_name = fullfile(results_dir, sprintf('time_resolved_%s_%s.png', rgn, ep));
    saveas(gcf, fig_name);
    close(gcf);

end