function out = Diagnose_GC_Floor(source, varargin)
%% out = Diagnose_GC_Floor(source, 'name',value, ...)
% Resolves ONE question: why is the mean off-diagonal GC below the analytic
% null bias floor p/N, when that floor is supposed to be a lower bound?
%
% Background. For pairwise GC the reduced model (target's own p lags) and the
% full model (target + source, 2p lags) differ by exactly p free parameters.
% The likelihood ratio N*GC is asymptotically chi2(p), so E[GC | no coupling]
% = p/N. Simulation of THIS pipeline's exact estimator confirms it to within
% 5-7% across white, AR, 1/f, correlated, autocorrelated, and near-collinear
% pairs (rcond down to 7e-8). An observed mean below the floor is therefore
% not explicable by the null and means one of the assumptions is off.
%
% This function attacks it three ways:
%   TEST 1  distribution of the saved GC values -- how many are negative or
%           near zero, and where does the mean actually sit
%   TEST 2  empirical null ON YOUR OWN DATA: circularly shift the source far
%           enough to destroy any causal relation, recompute, and measure the
%           resulting mean. If that empirical null also lands near your
%           observed mean, the floor formula does not apply to your data and
%           the absolute levels are fine (only relative comparisons matter).
%           If it lands near p/N while your observed mean is half that, the
%           observed values are wrong.
%   TEST 3  order sweep -- GC and the floor should both scale ~linearly in p.
%           A flat or non-monotonic curve indicates a misspecified fit.
%
% INPUT
%   source : path to <subject>_GC.mat, or the subjResult struct
%
% PARAMETERS
%   'dataDir' ''    folder of *_TrialDATA.mat. Given this, the epoch is
%                   fetched automatically via Get_GC_Epoch -- you do NOT need
%                   to build a matrix yourself. This is the easy path:
%                     Diagnose_GC_Floor('CLASE001_GC.mat','dataDir',dataDir)
%   'subject' ''    '' = inferred from the .mat filename
%   'event'   ''    '' = first event stored in the file
%   'trial'   []    [] = first usable trial
%   'epochOpts' {}  cell of extra args for Get_GC_Epoch. MUST mirror whatever
%                   you passed to Run_GC_Pipeline (winSamples, collapse, mode,
%                   band, ...) or the matrix will not match what was fitted.
%   'epoch'   []    supply the nChan x nSamp matrix directly instead.
%   'order'   []    VAR order; [] = read from the struct
%   'nNull'   200   surrogates for TEST 2
%   'orders'  [2 5 10 20 30]  orders for TEST 3
%
% OUTPUT struct `out` with .observed .empiricalNull .analyticFloor .verdict

ip = inputParser;
ip.addRequired('source');
ip.addParameter('epoch',[],@(x)isempty(x)||isnumeric(x));
ip.addParameter('order',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('nNull',200,@isscalar);
ip.addParameter('orders',[2 5 10 20 30],@isnumeric);
ip.addParameter('dataDir','',@ischar);      % fetch the epoch automatically
ip.addParameter('subject','',@ischar);      % '' = infer from the filename
ip.addParameter('event','',@ischar);        % '' = first event in the file
ip.addParameter('trial',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('epochOpts',{},@iscell);    % passed through to Get_GC_Epoch
ip.parse(source,varargin{:});
o = ip.Results;

if ischar(source)
    S = load(source);
    if isfield(S,'subjResult'), R = S.subjResult;
    else, fn = fieldnames(S); R = S.(fn{1}); end
else
    R = source;
end
p  = o.order; if isempty(p), p = R.order; end
Nw = NaN;
if isfield(R,'winSamples'), Nw = R.winSamples; end
floorV = p / Nw;

fprintf('\n===============================================================\n');
fprintf(' Diagnose_GC_Floor    order p = %d,  N = %g,  p/N = %.5f\n', p, Nw, floorV);
fprintf('===============================================================\n');

%% ---------------- TEST 1: distribution of saved GC ----------------
fprintf('\nTEST 1  distribution of saved GC values\n');
events = fieldnames(R.events);
nChan  = numel(R.conn.nodeLabel);
offd   = ~eye(nChan);
out.observed = struct();
for e = 1:numel(events)
    evt = events{e};
    G = R.events.(evt).GC;
    v = logical(R.events.(evt).validEpoch(:));
    G = G(:,:,v);
    flat = reshape(G, nChan*nChan, []);
    x = flat(offd(:), :); x = x(:); x = x(isfinite(x));
    if isempty(x), continue; end
    fprintf(['  %-12s n=%5d  mean %.5f (%.2fx floor)  median %.5f\n' ...
             '                min %.5f  max %.5f  %%<0: %.1f%%  %%<floor: %.1f%%\n'], ...
        evt, numel(x), mean(x), mean(x)/floorV, median(x), min(x), max(x), ...
        100*mean(x<0), 100*mean(x<floorV));
    out.observed.(evt) = [mean(x) median(x) min(x) max(x) mean(x<0)];
end
fprintf(['\n  Reading: GC is mathematically >= 0 in-sample (nested models), so a\n' ...
         '  substantial %%<0 means the variance ratio is being computed on\n' ...
         '  demeaned residuals whose means differ between the two models.\n']);

if isfield(R,'rcondMin')
    fprintf('\n  min rcond over all fits: %.3e%s\n', R.rcondMin, ...
        tern_(R.rcondMin < 1e-6, '   <-- ill-conditioned', '   (healthy)'));
end

%% ---------------- fetch the epoch if we were given a dataDir ----------------
if isempty(o.epoch) && ~isempty(o.dataDir)
    sub = o.subject;
    if isempty(sub) && ischar(source)
        [~, bn] = fileparts(source);
        sub = regexprep(bn, '_GC$', '');
    end
    evt = o.event; if isempty(evt), evt = events{1}; end
    fprintf('\nFetching an epoch: subject %s, event %s\n', sub, evt);
    args = [{'verbose', false}, o.epochOpts];
    if isempty(o.trial)
        % walk forward to the first usable trial
        o.epoch = [];
        for k = 1:200
            try
                [Mk, ~] = Get_GC_Epoch(o.dataDir, sub, evt, k, args{:});
            catch err
                error('Diagnose_GC_Floor:fetch', ...
                      'Get_GC_Epoch failed: %s', err.message);
            end
            if ~isempty(Mk), o.epoch = Mk; o.trial = k; break; end
        end
        if isempty(o.epoch)
            error('Diagnose_GC_Floor:noTrial', ...
                  'No usable trial found in the first 200 of %s.', evt);
        end
    else
        o.epoch = Get_GC_Epoch(o.dataDir, sub, evt, o.trial, args{:});
    end
    fprintf('  using trial %d -> M is %d x %d\n', o.trial, ...
            size(o.epoch,1), size(o.epoch,2));
    if size(o.epoch,2) ~= Nw
        fprintf(['  *** WARNING: this epoch is %d samples but the saved file\n' ...
                 '      says N = %g. Your epochOpts do not match the run.\n'], ...
                 size(o.epoch,2), Nw);
    end
end

%% ---------------- TESTS 2 & 3 need the actual data ----------------
if isempty(o.epoch)
    fprintf(['\nTESTS 2 and 3 skipped -- they need the actual signal, which\n' ...
             'the saved file does not contain (it holds GC, not voltages).\n' ...
             'Re-run with the data folder and it will fetch an epoch itself:\n\n' ...
             '   Diagnose_GC_Floor(''%s'', ''dataDir'', dataDir)\n\n'], ...
             tern_(ischar(source), source, '<file>.mat'));
    out.verdict = 'incomplete - supply an epoch';
    fprintf('===============================================================\n\n');
    return;
end

M = o.epoch;
if size(M,1) ~= nChan
    warning('Diagnose_GC_Floor:chan','epoch has %d rows, file says %d nodes.', ...
            size(M,1), nChan);
end
nChan = size(M,1); nSamp = size(M,2); offd = ~eye(nChan);
M = (M - mean(M,2)) ./ max(std(M,0,2), eps);

fprintf('\nTEST 2  empirical null on YOUR data (%d x %d epoch, %d surrogates)\n', ...
        nChan, nSamp, o.nNull);
obs = pairwise_gc(M, p);
fprintf('  observed mean off-diag GC : %.5f  (%.2fx p/N)\n', ...
        mean(obs(offd)), mean(obs(offd))/(p/nSamp));

nullMeans = zeros(o.nNull,1);
for it = 1:o.nNull
    Ms = M;
    for c = 1:nChan
        sh = round(nSamp*0.25) + randi(round(nSamp*0.5));
        Ms(c,:) = circshift(M(c,:), [0 sh]);
    end
    g = pairwise_gc(Ms, p);
    nullMeans(it) = mean(g(offd));
end
empNull = mean(nullMeans);
fprintf('  empirical null mean       : %.5f  (%.2fx p/N)\n', empNull, empNull/(p/nSamp));
fprintf('  analytic floor p/N        : %.5f\n', p/nSamp);
out.empiricalNull  = empNull;
out.analyticFloor  = p/nSamp;
out.observedMean   = mean(obs(offd));

fprintf('\n  VERDICT: ');
if abs(empNull - out.observedMean) < 0.25*out.analyticFloor
    fprintf(['observed ~= empirical null. There is little or no net\n' ...
             '  directed coupling beyond bias, BUT the empirical null itself\n' ...
             '  sits below p/N on this data, so absolute levels are simply not\n' ...
             '  comparable to the analytic floor. Use paired contrasts only.\n']);
    out.verdict = 'observed matches empirical null; floor formula inapplicable';
elseif out.observedMean > empNull * 1.25
    fprintf(['observed clearly exceeds the empirical null -- there IS\n' ...
             '  directed structure. The analytic p/N floor is the wrong\n' ...
             '  reference for this data; use the empirical null instead.\n']);
    out.verdict = 'real signal above empirical null';
else
    fprintf(['observed is BELOW the empirical null. That points at the\n' ...
             '  estimator or the epoch, not the biology. Check TEST 3.\n']);
    out.verdict = 'observed below empirical null - investigate estimator';
end

%% ---------------- TEST 3: order sweep ----------------
fprintf('\nTEST 3  order sweep (GC and floor should both scale ~linearly in p)\n');
fprintf('  %4s %12s %12s %10s\n','p','mean GC','p/N','ratio');
sweep = nan(numel(o.orders),3);
for k = 1:numel(o.orders)
    pk = o.orders(k);
    if nSamp <= 3*pk*nChan, continue; end
    g = pairwise_gc(M, pk);
    mg = mean(g(offd)); fl = pk/nSamp;
    sweep(k,:) = [pk mg fl];
    fprintf('  %4d %12.5f %12.5f %9.2fx\n', pk, mg, fl, mg/fl);
end
out.sweep = sweep;
fprintf(['\n  Reading: the ratio column should be roughly CONSTANT and >= 1.\n' ...
         '  A ratio well below 1 at every order means the residual variances\n' ...
         '  are not behaving as nested least-squares fits should.\n']);
fprintf('===============================================================\n\n');
end

%% ---------------- helpers ----------------
function G = pairwise_gc(data, p)
% Mirrors Compute_GC_V2's statistic exactly, without the permutation loop.
n = size(data,1);
G = nan(n);
for tgt = 1:n
    for src = 1:n
        if tgt == src, continue; end
        [Rr,~,~] = Autoregressive_Process_V1(data(tgt,:), p);
        [Rf,~,~] = Autoregressive_Process_V1([data(tgt,:); data(src,:)], p);
        vr = var(Rr,0,2); vf = var(Rf,0,2);
        G(tgt,src) = log(vr ./ vf(1,1));
    end
end
end

function s = tern_(c,a,b), if c, s=a; else, s=b; end, end
