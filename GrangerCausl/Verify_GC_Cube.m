function out = Verify_GC_Cube(gcFile, dataDir, varargin)
%% out = Verify_GC_Cube(gcFile, dataDir, 'name',value, ...)
% Recomputes GC for a handful of trials and diffs the result against what is
% stored in <subject>_GC.mat, element by element.
%
% WHY
%   Diagnose_GC_Floor showed a contradiction that cannot both be true:
%     TEST 2/3 -- a FRESH computation on trial 1 gives mean off-diag GC
%                 0.00511 (1.53x the p/N floor), every order ratio >= 1.07,
%                 no negatives. The estimator is behaving correctly.
%     TEST 1   -- the STORED cube for the same subject, event, order and
%                 window has mean 0.00184 (0.55x floor) and 39% negatives.
%   GC = log(RSS_reduced / RSS_full) with nested models cannot be negative,
%   so the stored numbers are not what this estimator produces on this data.
%   Either the stored file came from a different run/settings, or something
%   between assemble_epoch and Compute_GC_V2 differs from the path
%   Get_GC_Epoch reproduces. This function tells you which, in one shot.
%
% INPUT
%   gcFile   path to <subject>_GC.mat
%   dataDir  folder of *_TrialDATA.mat
%
% PARAMETERS
%   'event'     ''   '' = first event in the file
%   'trials'    1:10 trial indices to check
%   'epochOpts' {}   MUST mirror the Run_GC_Pipeline options you used
%
% OUTPUT struct `out`
%   .stored .fresh   [nChan nChan nTrial] matched cubes
%   .maxAbsDiff .corr .agree
%   .verdict

ip = inputParser;
ip.addRequired('gcFile',@ischar);
ip.addRequired('dataDir',@ischar);
ip.addParameter('event','',@ischar);
ip.addParameter('trials',1:10,@isnumeric);
ip.addParameter('epochOpts',{},@iscell);
ip.parse(gcFile,dataDir,varargin{:});
o = ip.Results;

S = load(gcFile);
if isfield(S,'subjResult'), R = S.subjResult;
else, fn = fieldnames(S); R = S.(fn{1}); end
[~,bn] = fileparts(gcFile);
sub = regexprep(bn,'_GC$','');

events = fieldnames(R.events);
evt = o.event; if isempty(evt), evt = events{1}; end
p   = R.order;
nChan = numel(R.conn.nodeLabel);
offd  = ~eye(nChan);

fprintf('\n===============================================================\n');
fprintf(' Verify_GC_Cube   %s / %s   order %d, %d nodes\n', sub, evt, p, nChan);
fprintf('===============================================================\n');
if isfield(R,'mode'), fprintf('  file mode      : %s\n', R.mode); end
if isfield(R,'winSamples'), fprintf('  file N per fit : %g\n', R.winSamples); end
if isfield(R,'blockTrials'), fprintf('  file blockTrials: %d\n', R.blockTrials); end
fprintf('  epochOpts      : %s\n', tern_(isempty(o.epochOpts), ...
        '(none -- defaults)', 'supplied'));

stored = []; fresh = []; used = [];
args = [{'verbose',false}, o.epochOpts];

fprintf('\n  %6s %14s %14s %14s %10s\n','trial','stored mean','fresh mean','maxAbsDiff','stored<0');
for k = o.trials(:)'
    if k > size(R.events.(evt).GC,3), continue; end
    Gs = R.events.(evt).GC(:,:,k);
    if all(isnan(Gs(offd))), continue; end          % skipped trial

    M = Get_GC_Epoch(dataDir, sub, evt, k, args{:});
    if isempty(M), continue; end

    GC = Compute_GC_V2(M, 'order', p, 'nPerm', 0, 'zscore', false);
    Gf = GC.matrix;

    d = abs(Gs(offd) - Gf(offd));
    fprintf('  %6d %14.5f %14.5f %14.2e %9.1f%%\n', k, ...
            mean(Gs(offd),'omitnan'), mean(Gf(offd),'omitnan'), ...
            max(d), 100*mean(Gs(offd) < 0));

    stored = cat(3, stored, Gs);
    fresh  = cat(3, fresh,  Gf);
    used(end+1) = k; %#ok<AGROW>
end

if isempty(stored)
    error('Verify_GC_Cube:none','No comparable trials found.');
end

sv = stored(repmat(offd,[1 1 size(stored,3)]));
fv = fresh( repmat(offd,[1 1 size(fresh,3)]));
ok = isfinite(sv) & isfinite(fv);
sv = sv(ok); fv = fv(ok);

out.stored = stored; out.fresh = fresh; out.trials = used;
out.maxAbsDiff = max(abs(sv-fv));
cc = corrcoef(sv,fv); out.corr = cc(1,2);
out.agree = out.maxAbsDiff < 1e-9;

fprintf('\n  pooled over %d trials, %d off-diagonal values:\n', numel(used), numel(sv));
fprintf('    stored : mean %.5f, %.1f%% negative, min %.5f\n', ...
        mean(sv), 100*mean(sv<0), min(sv));
fprintf('    fresh  : mean %.5f, %.1f%% negative, min %.5f\n', ...
        mean(fv), 100*mean(fv<0), min(fv));
fprintf('    max |stored - fresh| = %.3e,  corr = %.4f\n', out.maxAbsDiff, out.corr);

fprintf('\n  VERDICT: ');
if out.agree
    fprintf(['stored and fresh are IDENTICAL. The negatives are real\n' ...
             '  output of this estimator on this data, which contradicts the\n' ...
             '  nested-model argument -- report this back, it means\n' ...
             '  Autoregressive_Process_V1 is not returning least-squares\n' ...
             '  residuals for every pair.\n']);
    out.verdict = 'identical - estimator itself produces negatives';
elseif out.corr > 0.95
    fprintf(['same shape, shifted/scaled. Stored and fresh correlate at\n' ...
             '  %.3f but differ by up to %.2e. Most likely the run used\n' ...
             '  different options than the epochOpts here -- check winSamples,\n' ...
             '  collapse, zscore, channelUnit.\n'], out.corr, out.maxAbsDiff);
    out.verdict = 'correlated but offset - options mismatch';
else
    fprintf(['stored and fresh DISAGREE (corr %.3f). The saved file was\n' ...
             '  not produced by this code path on this data. Most likely it\n' ...
             '  predates a change -- delete it and re-run the pipeline.\n'], out.corr);
    out.verdict = 'disagree - stale or mismatched file';
end
fprintf('===============================================================\n\n');
end

function s = tern_(c,a,b), if c, s=a; else, s=b; end, end
