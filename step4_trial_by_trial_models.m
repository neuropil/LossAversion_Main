%% step4_trial_by_trial_models.m
% Logistic mixed-effects models: Choice ~ EV + neural features
% Runs per region × epoch, with interaction and time-resolved variants
% Requires: master_table.mat from step3

clear; clc;

project_dir  = 'Z:\LossAversion\LH_Data\JAT_ADAtest';
results_dir  = fullfile(project_dir, 'results');
load(fullfile(results_dir, 'master_table.mat'), 'master_table');

% regions      = unique(master_table.region);
% epochs       = unique(master_table.epoch);
% neural_feats = {'theta', 'alpha', 'beta', 'hfa'};
n_perm       = 1000;   % permutations for null distribution

subjects = unique(master_table.subject);
for s = 1:numel(subjects)
    idx = master_table.subject == subjects(s);
    master_table.gain_z(idx) = zscore(master_table.gain(idx));
    master_table.loss_z(idx) = zscore(master_table.loss(idx));
end

% Also z-score neural features within subject × region × epoch
regions      = unique(master_table.region);
epochs       = unique(master_table.epoch);
neural_feats = {'theta','alpha','beta','hfa'};

for s = 1:numel(subjects)
    for r = 1:numel(regions)
        for e = 1:numel(epochs)
            for f = 1:numel(neural_feats)
                feat = neural_feats{f};
                idx  = master_table.subject == subjects(s) & ...
                       master_table.region   == regions(r) & ...
                       master_table.epoch    == epochs(e);
                if sum(idx) > 1
                    master_table.([feat '_z'])(idx) = zscore(master_table.(feat)(idx));
                end
            end
        end
    end
end

rename_map = { ...
    'gain_z',  'gainZ';  ...
    'loss_z',  'lossZ';  ...
    'theta_z', 'thetaZ'; ...
    'alpha_z', 'alphaZ'; ...
    'beta_z',  'betaZ';  ...
    'hfa_z',   'hfaZ'};

for i = 1:size(rename_map, 1)
    old = rename_map{i,1};
    new = rename_map{i,2};
    if ismember(old, master_table.Properties.VariableNames)
        master_table.(new) = master_table.(old);
        master_table.(old) = [];
    end
end

% ── Ensure correct variable types ────────────────────────────────────────
master_table.subject = categorical(master_table.subject);
master_table.region  = categorical(master_table.region);
master_table.epoch   = categorical(master_table.epoch);
master_table.choice  = double(master_table.choice);



%% ═══════════════════════════════════════════════════════════════
%% MODEL 1 — Base model (behavioral only)
%% ═══════════════════════════════════════════════════════════════
fprintf('\n=== BASE MODEL (behavioral only) ===\n');

n_sub_full               = numel(unique(master_table.subject));
[base_formula, base_cov] = build_formula('', n_sub_full, 'base');

mdl_base = fitglme(master_table, base_formula, ...
             'Distribution',     'Binomial', ...
             'Link',             'logit', ...
             'FitMethod',        'Laplace', ...
             'CovariancePattern', base_cov)

results.base.model   = mdl_base;
results.base.AIC     = mdl_base.ModelCriterion.AIC;
results.base.formula = base_formula;
fprintf('Base AIC: %.2f\n', results.base.AIC);

%% ═══════════════════════════════════════════════════════════════
%% MODEL 2 — Neural main effect per region × epoch × feature
%% choice ~ gainZ + lossZ + neural_z + RE
%% ═══════════════════════════════════════════════════════════════
fprintf('\n=== NEURAL MAIN EFFECT MODELS ===\n');

for r = 1:numel(regions)
    rgn = char(regions(r));

    for e = 1:numel(epochs)
        ep   = char(epochs(e));
        mask = master_table.region == rgn & master_table.epoch == ep;
        T    = master_table(mask, :);

        if height(T) < 50
            fprintf('Skipping %s %s — too few rows (%d)\n', rgn, ep, height(T));
            continue
        end

        n_subjects = numel(unique(T.subject));

        for f = 1:numel(neural_feats_z)
            feat = neural_feats_z{f};

            if all(isnan(T.(feat))) || std(T.(feat), 'omitnan') < 1e-6
                fprintf('Skipping %s %s %s — no variance\n', rgn, ep, feat);
                continue
            end

            [formula, cov_pat] = build_formula(feat, n_subjects, 'main');

            try
                mdl = fitglme(T, formula, ...
                        'Distribution',     'Binomial', ...
                        'Link',             'logit', ...
                        'FitMethod',        'Laplace', ...
                        'CovariancePattern', cov_pat);

                coef_tbl = mdl.Coefficients;
                feat_row = strcmp(coef_tbl.Name, feat);

                res.model        = mdl;
                res.region       = rgn;
                res.epoch        = ep;
                res.feature      = feat;
                res.beta         = coef_tbl.Estimate(feat_row);
                res.SE           = coef_tbl.SE(feat_row);
                res.tstat        = coef_tbl.tStat(feat_row);
                res.pval         = coef_tbl.pValue(feat_row);
                res.AIC          = mdl.ModelCriterion.AIC;
                res.deltaAIC     = results.base.AIC - res.AIC;
                res.n_trials     = height(T);
                res.n_subjects   = n_subjects;
                res.re_structure = n_subjects >= 8;  % true = random slope

                key = sprintf('%s_%s_%s', rgn, ep, feat);
                results.neural.(key) = res;

                fprintf('[%s | %s | %s]  beta=%.3f  p=%.4f  dAIC=%.2f  RE=%s\n', ...
                    rgn, ep, feat, res.beta, res.pval, res.deltaAIC, ...
                    mat2str(res.re_structure));

            catch ME
                fprintf('Model failed: %s %s %s — %s\n', rgn, ep, feat, ME.message);
            end
        end
    end
end

%% ═══════════════════════════════════════════════════════════════
%% MODEL 3 — Interaction model (neural × loss)
%% choice ~ gainZ + lossZ + neural_z + neural_z:lossZ + RE
%% ═══════════════════════════════════════════════════════════════
fprintf('\n=== INTERACTION MODELS (neural_z × lossZ) ===\n');

for r = 1:numel(regions)
    rgn = char(regions(r));

    for e = 1:numel(epochs)
        ep   = char(epochs(e));
        mask = master_table.region == rgn & master_table.epoch == ep;
        T    = master_table(mask, :);

        if height(T) < 50, continue; end

        n_subjects = numel(unique(T.subject));

        for f = 1:numel(neural_feats_z)
            feat = neural_feats_z{f};

            if all(isnan(T.(feat))) || std(T.(feat), 'omitnan') < 1e-6
                continue
            end

            [formula, cov_pat] = build_formula(feat, n_subjects, 'interaction');
            int_name           = sprintf('%s:lossZ', feat);

            try
                mdl = fitglme(T, formula, ...
                        'Distribution',     'Binomial', ...
                        'Link',             'logit', ...
                        'FitMethod',        'Laplace', ...
                        'CovariancePattern', cov_pat);

                coef_tbl = mdl.Coefficients;
                feat_row = strcmp(coef_tbl.Name, feat);
                % MATLAB may flip interaction term order — check both
                int_name_a = sprintf('%s:lossZ', feat);
                int_name_b = sprintf('lossZ:%s', feat);
                int_row    = strcmp(coef_tbl.Name, int_name_a) | ...
                    strcmp(coef_tbl.Name, int_name_b);

                if ~any(int_row)
                    fprintf('WARNING: interaction term not found — %s %s %s\n', rgn, ep, feat);
                    fprintf('  Available: %s\n', strjoin(coef_tbl.Name, ', '));
                    continue
                end

                res.model        = mdl;
                res.region       = rgn;
                res.epoch        = ep;
                res.feature      = feat;
                res.beta_main    = coef_tbl.Estimate(feat_row);
                res.p_main       = coef_tbl.pValue(feat_row);
                res.beta_int     = coef_tbl.Estimate(int_row);
                res.SE_int       = coef_tbl.SE(int_row);
                res.tstat_int    = coef_tbl.tStat(int_row);
                res.p_int        = coef_tbl.pValue(int_row);
                res.AIC          = mdl.ModelCriterion.AIC;
                res.deltaAIC     = results.base.AIC - res.AIC;
                res.n_trials     = height(T);
                res.n_subjects   = n_subjects;
                res.re_structure = n_subjects >= 8;

                key = sprintf('%s_%s_%s_interaction', rgn, ep, feat);
                results.interaction.(key) = res;

                fprintf('[%s | %s | %s×lossZ]  beta_int=%.3f  p=%.4f  RE=%s\n', ...
                    rgn, ep, feat, res.beta_int, res.p_int, mat2str(res.re_structure));

            catch ME
                fprintf('Interaction model failed: %s %s %s — %s\n', ...
                    rgn, ep, feat, ME.message);
            end
        end
    end
end

%% ═══════════════════════════════════════════════════════════════
%% MODEL 4 — Permutation test
%% Shuffle choice within subject to build null distribution
%% ═══════════════════════════════════════════════════════════════
fprintf('\n=== PERMUTATION TESTING ===\n');

% Automatically target the top result from Model 2
top_row        = summary_tbl(1, :);   % already sorted by p_fdr
target_region  = top_row.region{1};
target_epoch   = top_row.epoch{1};
target_feature = top_row.feature{1};

fprintf('Running permutation on top result: %s | %s | %s\n', ...
    target_region, target_epoch, target_feature);

mask   = master_table.region == target_region & ...
         master_table.epoch  == target_epoch;
T_perm = master_table(mask, :);

n_subjects_perm              = numel(unique(T_perm.subject));
[formula_perm, cov_pat_perm] = build_formula(target_feature, n_subjects_perm, 'main');

% Observed beta
mdl_obs  = fitglme(T_perm, formula_perm, ...
             'Distribution',     'Binomial', ...
             'Link',             'logit', ...
             'FitMethod',        'Laplace', ...
             'CovariancePattern', cov_pat_perm);
coef_obs = mdl_obs.Coefficients;
beta_obs = coef_obs.Estimate(strcmp(coef_obs.Name, target_feature));

% Null distribution
null_betas = nan(n_perm, 1);
subjects   = unique(T_perm.subject);

fprintf('Running %d permutations for %s | %s | %s...\n', ...
    n_perm, target_region, target_epoch, target_feature);

for p = 1:n_perm
    T_shuf = T_perm;
    for s = 1:numel(subjects)
        sub_idx = find(T_shuf.subject == subjects(s));
        T_shuf.choice(sub_idx) = T_shuf.choice(sub_idx(randperm(numel(sub_idx))));
    end
    try
        mdl_null  = fitglme(T_shuf, formula_perm, ...
                      'Distribution',     'Binomial', ...
                      'Link',             'logit', ...
                      'FitMethod',        'Laplace', ...
                      'CovariancePattern', cov_pat_perm);
        coef_null = mdl_null.Coefficients;
        null_betas(p) = coef_null.Estimate(strcmp(coef_null.Name, target_feature));
    catch
        null_betas(p) = NaN;
    end
    if mod(p,100) == 0, fprintf('  %d/%d done\n', p, n_perm); end
end

perm_p = mean(abs(null_betas) >= abs(beta_obs), 'omitnan');

perm_key = sprintf('%s_%s_%s', target_region, target_epoch, target_feature);
results.permutation.(perm_key) = struct( ...
    'null_betas',    null_betas,     ...
    'beta_obs',      beta_obs,       ...
    'perm_pval',     perm_p,         ...
    'region',        target_region,  ...
    'epoch',         target_epoch,   ...
    'feature',       target_feature, ...
    'n_subjects',    n_subjects_perm, ...
    're_structure',  n_subjects_perm >= 8);

fprintf('Permutation p = %.4f  (observed beta = %.4f)\n', perm_p, beta_obs);

%% ═══════════════════════════════════════════════════════════════
%% SUMMARIZE — Results table with FDR correction
%% ═══════════════════════════════════════════════════════════════
fprintf('\n=== BUILDING SUMMARY TABLE ===\n');

keys    = fieldnames(results.neural);
n_tests = numel(keys);

sum_region   = cell(n_tests, 1);
sum_epoch    = cell(n_tests, 1);
sum_feature  = cell(n_tests, 1);
sum_beta     = nan(n_tests, 1);
sum_SE       = nan(n_tests, 1);
sum_tstat    = nan(n_tests, 1);
sum_pval     = nan(n_tests, 1);
sum_dAIC     = nan(n_tests, 1);
sum_ntrials  = nan(n_tests, 1);

for k = 1:n_tests
    res              = results.neural.(keys{k});
    sum_region{k}    = res.region;
    sum_epoch{k}     = res.epoch;
    sum_feature{k}   = res.feature;
    sum_beta(k)      = res.beta;
    sum_SE(k)        = res.SE;
    sum_tstat(k)     = res.tstat;
    sum_pval(k)      = res.pval;
    sum_dAIC(k)      = res.deltaAIC;
    sum_ntrials(k)   = res.n_trials;
end

sum_pval_fdr = fdr_bh(sum_pval);

summary_tbl = table(sum_region, sum_epoch, sum_feature, ...
                    sum_beta, sum_SE, sum_tstat, ...
                    sum_pval, sum_pval_fdr, sum_dAIC, sum_ntrials, ...
                    'VariableNames', ...
                    {'region','epoch','feature', ...
                     'beta','SE','tstat', ...
                     'p_raw','p_fdr','deltaAIC','n_trials'});

summary_tbl = sortrows(summary_tbl, 'p_fdr');

fprintf('\n--- Significant results (FDR q<0.05) ---\n');
disp(summary_tbl(summary_tbl.p_fdr < 0.05, :));

save(fullfile(results_dir, 'model_results.mat'), 'results', 'summary_tbl', '-v7.3');
writetable(summary_tbl, fullfile(results_dir, 'summary_table.csv'));
fprintf('Results saved.\n');

%% ═══════════════════════════════════════════════════════════════
%% HELPERS
%% ═══════════════════════════════════════════════════════════════

function [formula, cov_pat] = build_formula(feat, n_subjects, model_type)
% Build fitglme formula and CovariancePattern based on subject count.
% Uses random slope on lossZ when n_subjects >= 8, intercept-only otherwise.

    if n_subjects >= 8
        re_str  = '(1 + lossZ | subject)';
        cov_pat = 'Diagonal';
    else
        re_str  = '(1 | subject)';
        cov_pat = 'FullCholesky';
    end

    switch model_type
        case 'base'
            formula = sprintf('choice ~ gainZ + lossZ + %s', re_str);
        case 'main'
            formula = sprintf('choice ~ gainZ + lossZ + %s + %s', feat, re_str);
        case 'interaction'
            formula = sprintf('choice ~ gainZ + lossZ + %s + %s:lossZ + %s', ...
                              feat, feat, re_str);
        otherwise
            error('build_formula: unknown model_type "%s"', model_type);
    end
end

function padj = fdr_bh(pvals)
% Benjamini-Hochberg FDR correction.
    n     = numel(pvals);
    [~,I] = sort(pvals);
    padj  = pvals;
    for i = n:-1:1
        if i == n
            padj(I(i)) = pvals(I(i));
        else
            padj(I(i)) = min(pvals(I(i)) * n/i, padj(I(i+1)));
        end
    end
    padj = min(padj, 1);
end