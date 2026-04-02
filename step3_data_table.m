%% step3_build_model_table.m
% Assembles one tall table across all subjects, regions, epochs
% Ready to pass into fitglme for the logistic mixed-effects model

clear; clc;

project_dir = 'Z:\LossAversion\LH_Data\JAT_ADAtest';
feat_dir    = fullfile(project_dir, 'features');
results_dir = fullfile(project_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

feat_files = dir(fullfile(feat_dir, '*.mat'));
all_tables = {};

for f = 1:numel(feat_files)
    load(fullfile(feat_dir, feat_files(f).name), 'feat');

    n = feat.n_trials;
    b = feat.beh;   % b is the behavioral table (rows = trials)
    b.TrialID = transpose(1:135);

    T = table();
    T.subject = repmat({feat.subject}, n, 1);
    T.region  = repmat({feat.region},  n, 1);
    T.epoch   = repmat({feat.epoch},   n, 1);

    % ------------------------------------------------------------------ %
    %  Behavioral variables — mapped from b (the beh table)
    %  Column names match the spreadsheet headers exactly
    % ------------------------------------------------------------------ %

    % Core choice variable (1 = accepted risky gamble, 0 = rejected)
    T.choice   = b.choice(:);

    % Gamble parameters
    T.gain     = b.riskyGain(:);           % gain magnitude (col 2)
    T.loss     = abs(b.riskyLoss(:));      % loss magnitude, stored as negative → make positive
    T.alt      = b.alternative(:);         % safe alternative value (col 4)

    % Trial metadata
    T.triBlock      = b.triBlock(:);       % trial number within block
    T.ischecktrial  = logical(b.ischecktrial(:));  % attention/check trial flag
    T.TrialID = b.TrialID;

    % Outcome
    T.outcome  = b.outcome(:);             % realized outcome on that trial

    % Expected value of the risky option (assuming 50/50)
    T.EV = 0.5 * b.riskyGain(:) + 0.5 * b.riskyLoss(:);
    % Note: riskyLoss is already signed (negative), so this gives true EV

    % Loss-aversion parameters — stored per-subject but repeated across trials
    % T.lambda   = b.lambda(:);             % continuous λ estimate
    % T.LA_group = b.LA_group(:);           % categorical: 'low' / 'neutral' / 'high'

    % Optional response-side columns (useful for RT/motor confound checks)
    T.loc      = b.loc(:);               % response location (1 or 2)
    T.response = b.response(:);          % key pressed ('m' or 'z')

    % ------------------------------------------------------------------ %
    %  Neural features (channel-mean per trial, from step2)
    % ------------------------------------------------------------------ %
    T.theta = feat.theta_mean(:);
    T.alpha = feat.alpha_mean(:);
    T.beta  = feat.beta_mean(:);
    T.hfa   = feat.hfa_mean(:);

    % Z-score neural features across trials (within subject/region/epoch)
    for col = {'theta','alpha','beta','hfa'}
        T.(col{1}) = zscore(T.(col{1}));
    end

    % ------------------------------------------------------------------ %
    %  Exclude check trials from modeling (keep for QC if needed)
    % ------------------------------------------------------------------ %
    T = T(~T.ischecktrial, :);

    % ------------------------------------------------------------------ %
    %  Exclude trials with missing response
    % ------------------------------------------------------------------ %
    T = T(cellfun(@(x) ~isempty(x), T.response, 'UniformOutput', true), :);

    all_tables{end+1} = T;
end

%% Stack into one tall table
master_table = vertcat(all_tables{:});

%% Convert to categorical for fitglme
master_table.subject  = categorical(master_table.subject);
master_table.region   = categorical(master_table.region);
master_table.epoch    = categorical(master_table.epoch);
% master_table.LA_group = categorical(master_table.LA_group);
master_table.response = categorical(master_table.response);

%% Save
save(fullfile(results_dir, 'master_table.mat'), 'master_table', '-v7.3');

fprintf('Master table: %d rows, %d variables\n', height(master_table), width(master_table));
fprintf('Subjects: %d | Regions: %d | Epoch types: %d\n', ...
    numel(unique(master_table.subject)), ...
    numel(unique(master_table.region)),  ...
    numel(unique(master_table.epoch)));