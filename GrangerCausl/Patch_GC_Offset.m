function out = Patch_GC_Offset(gcFile, varargin)
%% out = Patch_GC_Offset(gcFile, 'name',value, ...)
% Corrects the dof-normalisation offset in a <subject>_GC.mat written by a
% Run_GC_Pipeline run that used the pre-fix Compute_Spectral_GC.
%
% THE BUG
%   local_fit_var normalised the residual covariance by (m - n*p). The reduced
%   model in pairwise GC is univariate (n=1) and the full model is bivariate
%   (n=2), so the two covariances were divided by different numbers and
%       td = log(RSS_r/RSS_f) + log((m-2p)/(m-p))
%   The trailing term is ~ -p/m, i.e. almost exactly the null bias floor p/N.
%   Every stored time-domain GC value is low by that constant.
%
% WHAT NEEDS FIXING, AND WHAT DOES NOT
%   .GC        affected ONLY if the run used the default 'specBand','' so that
%              .GC holds the time-domain td. If a band was selected, .GC holds
%              a spectral quantity and needs nothing.
%   .spec      NOT affected. Sigma appears in both halves of the Geweke ratio,
%   .bandGC    so any scalar multiple cancels (verified to 2.9e-14).
%   .sig/.pval NOT affected. Surrogates went through the identical code path,
%              so the offset cancels in the observed-vs-null comparison. Your
%              permutation results are valid as they stand -- this is why
%              patching beats re-running.
%   voltage-mode files are NOT affected at all (they use Compute_GC_V2).
%
% The patch is idempotent: it writes a .gcOffsetCorrected flag and refuses to
% apply twice.
%
% PARAMETERS
%   'specBand'  ''     set this if the run used a non-empty 'specBand', in
%                      which case nothing is done to .GC
%   'backup'    true   write <file>.prepatch.mat first
%   'dryRun'    false  report what would change without writing
%
% OUTPUT struct `out`: .correction .meanBefore .meanAfter .applied

ip = inputParser;
ip.addRequired('gcFile',@ischar);
ip.addParameter('specBand','',@ischar);
ip.addParameter('backup',true,@islogical);
ip.addParameter('dryRun',false,@islogical);
ip.parse(gcFile,varargin{:});
o = ip.Results;

S = load(gcFile);
fld = 'subjResult';
if ~isfield(S,fld), fn = fieldnames(S); fld = fn{1}; end
R = S.(fld);

fprintf('\n--- Patch_GC_Offset : %s ---\n', gcFile);

out.applied = false;
if isfield(R,'gcOffsetCorrected') && R.gcOffsetCorrected
    fprintf('  ALREADY CORRECTED (offset %.6f applied previously). Nothing to do.\n\n', ...
            R.gcOffsetValue);
    out.correction = R.gcOffsetValue; return;
end

md = 'voltage'; if isfield(R,'mode'), md = R.mode; end
if strcmpi(md,'voltage') && (~isfield(R,'blockTrials') || R.blockTrials == 1)
    fprintf(['  mode is ''voltage'' with per-trial fits -> GC came from\n' ...
             '  Compute_GC_V2, which never had this bug. Nothing to do.\n\n']);
    out.correction = 0; return;
end
if ~isempty(o.specBand)
    fprintf(['  run used specBand=''%s'', so .GC holds a spectral quantity,\n' ...
             '  which is unaffected. Nothing to do.\n\n'], o.specBand);
    out.correction = 0; return;
end

% ---- reconstruct m exactly as local_fit_var computed it ----
if ~isfield(R,'winSamples') || ~isfield(R,'order')
    error('Patch_GC_Offset:fields', ...
        ['File lacks .winSamples/.order, so m cannot be reconstructed. ' ...
         'Re-run the pipeline instead.']);
end
bt = 1; if isfield(R,'blockTrials'), bt = R.blockTrials; end
p  = R.order;
nSampEff = R.winSamples / bt;               % winSamples = nSampEff*blockTrials
m  = bt * (nSampEff - p);
if m <= 2*p
    error('Patch_GC_Offset:tiny','Reconstructed m=%g is too small.', m);
end
correction = -log((m - 2*p)/(m - p));       % ADD this to every stored td

fprintf('  mode %s | order %d | blockTrials %d | nSamp/trial %g -> m = %g\n', ...
        md, p, bt, nSampEff, m);
fprintf('  correction to ADD : %+.6f   (p/N floor = %.6f)\n', correction, p/R.winSamples);

events = fieldnames(R.events);
before = []; after = [];
for e = 1:numel(events)
    G = R.events.(events{e}).GC;
    if isempty(G), continue; end
    nChan = size(G,1); offd = ~eye(nChan);
    fl = reshape(G, nChan*nChan, []); v = fl(offd(:),:); v = v(isfinite(v));
    before = [before; v]; %#ok<AGROW>
    if ~o.dryRun
        R.events.(events{e}).GC = G + correction;
    end
    after = [after; v + correction]; %#ok<AGROW>
end

fprintf('  off-diagonal GC : mean %.5f -> %.5f  |  %%negative %.1f%% -> %.1f%%\n', ...
        mean(before), mean(after), 100*mean(before<0), 100*mean(after<0));
out.correction = correction;
out.meanBefore = mean(before); out.meanAfter = mean(after);

if o.dryRun
    fprintf('  DRY RUN -- nothing written.\n\n'); return;
end

if o.backup
    bk = regexprep(gcFile,'\.mat$','.prepatch.mat');
    if exist(bk,'file') ~= 2
        copyfile(gcFile, bk);
        fprintf('  backup    : %s\n', bk);
    else
        fprintf('  backup    : %s (already exists, kept)\n', bk);
    end
end

R.gcOffsetCorrected = true;
R.gcOffsetValue     = correction;
S.(fld) = R;
save(gcFile, '-struct', 'S', '-v7.3');
fprintf('  written   : %s\n\n', gcFile);
out.applied = true;
end
