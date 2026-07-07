function test_load_region_epochs(filepath, events)
%% test_load_region_epochs(filepath, events)
% Diagnostic for load_region_epochs on ONE real region file. Reports what the
% file actually contains, runs the loader, and checks the outputs are shaped
% the way the pipeline expects. Run this before the full pipeline.
%
% Usage:
%   test_load_region_epochs('CLASE001_L_AC_TrialDATA.mat')
%   test_load_region_epochs('CLASE001_L_AC_TrialDATA.mat', {'CHOICE','RESPONSEON','OUTCOME'})

if nargin < 2 || isempty(events), events = {'CHOICE','RESPONSEON','OUTCOME'}; end
fprintf('\n===== testing load_region_epochs on =====\n  %s\n\n', filepath);
pass = true;

%% 0. file exists + raw variable dump (so structural mismatches are visible)
if exist(filepath,'file') ~= 2
    error('File not found: %s', filepath);
end
info = whos('-file', filepath);
fprintf('Variables stored in the file:\n');
for i = 1:numel(info)
    fprintf('   %-18s %-12s size %s\n', info(i).name, info(i).class, mat2str(info(i).size));
end

S  = load(filepath);
fn = fieldnames(S);
main = [];
for i = 1:numel(fn)
    if safe_istable(S.(fn{i})) || isstruct(S.(fn{i})), main = S.(fn{i}); break; end
end
if isempty(main)
    fprintf('\n[!] No table/struct variable found — loader will error. Inspect the dump above.\n');
    return;
end
if safe_istable(main)
    fprintf('\nMain variable is a TABLE with columns: %s\n', strjoin(main.Properties.VariableNames, ', '));
    fprintf('   rows (trials): %d\n', size(main,1));
else
    fprintf('\nMain variable is a %s (size %s).\n', class(main), mat2str(size(main)));
end

%% 1. run the loader
try
    [labels, epochs, fs] = load_region_epochs(filepath); %#ok<ASGLU>
catch err
    fprintf('\n[FAIL] load_region_epochs threw: %s\n', err.message);
    fprintf('       -> edit field names in load_region_epochs.m to match the dump above.\n');
    return;
end
n = numel(labels);
fprintf('\nLoader returned %d epochs (labels) and %d matrices.\n', n, numel(epochs));
pass = pass && (n>0) && (numel(epochs)==n);

%% 2. events present?
u = unique(labels);
fprintf('\nUnique event labels found: %s\n', strjoin(cellstr(u), ', '));
for e = 1:numel(events)
    c = sum(strcmp(labels, events{e}));
    ok = c>0; pass = pass && ok;
    fprintf('   %-16s : %4d epochs   [%s]\n', events{e}, c, tern(ok,'OK','MISSING <-'));
end

%% 3. separate REAL epochs from placeholders, then check shapes
nch = cellfun(@(x) size(x,1), epochs);
nsm = cellfun(@(x) size(x,2), epochs);
fin = cellfun(@(x) all(isfinite(x(:))), epochs);
% a "real" epoch is full-length and finite; placeholders are short/NaN cells
isReal = (nsm >= 50) & fin;
nPlace = sum(~isReal);
fprintf('\nReal epochs: %d | placeholder/degenerate (short or NaN): %d\n', sum(isReal), nPlace);
if nPlace>0
    fprintf('   (this is normal if some trials lack an event; the pipeline skips them)\n');
end
if any(isReal)
    rc = nch(isReal); rs = nsm(isReal);
    fprintf('Among REAL epochs -> channels: min %d max %d (expect constant per region)\n', min(rc), max(rc));
    fprintf('                     samples : min %d max %d (1501 = %.2fs @500Hz)\n', min(rs), max(rs), max(rs)/500);
    if numel(unique(rc))>1
        fprintf('   [!] channel count varies among REAL epochs — unexpected.\n'); pass=false;
    end
    if max(rc) > min(rs)
        fprintf('   [!] more channels than samples — Volts may be TRANSPOSED.\n'); pass=false;
    end
else
    fprintf('   [FAIL] no real epochs found.\n'); pass=false;
end

%% 5. peek at one REAL epoch of the first event
cand = find(strcmp(labels, events{1}) & isReal, 1);
if ~isempty(cand)
    x = epochs{cand};
    gm = mean(x(:));
    fprintf('\nFirst real %s epoch: size %s | ch1 range [%.2f, %.2f] | grand mean %.3f\n', ...
        events{1}, mat2str(size(x)), min(x(1,:)), max(x(1,:)), gm);
    if abs(gm) > 5 || max(abs(x(:))) > 20
        fprintf('   -> looks like RAW voltage (not z-scored); keep pipeline ''zscore'',true.\n');
    end
end

%% 6. cross-check loader vs a manual extraction (best-effort)
try
    if safe_istable(main), col = main.TrialTablesZS; else, col = main; end
    if iscell(col), tr1 = col{1}; else, tr1 = col(1); end
    manual = double(tr1.(events{1}).Volts);
    loaded = epochs{find(strcmp(labels,events{1}),1)};   % first-trial extraction
    if isequal(manual, loaded)
        fprintf('Cross-check: loader output == manual extraction of first %s.  OK\n', events{1});
    else
        fprintf('Cross-check: MISMATCH (manual %s vs loader %s).\n', mat2str(size(manual)), mat2str(size(loaded)));
        pass = false;
    end
catch err
    fprintf('Cross-check skipped (%s)\n', err.message);
end

fprintf('\n===== %s =====\n\n', tern(pass,'PASS — loader looks good, safe to run the pipeline', ...
                                        'ISSUES ABOVE — fix load_region_epochs.m before running the pipeline'));
end

% ---- helpers ----
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
function tf = safe_istable(v)
if exist('istable','builtin') || exist('istable','file'), tf = istable(v); else, tf = false; end
end
