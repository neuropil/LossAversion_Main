function BS = Bootstrap_GC_TimeResolved(dataDir, subject, varargin)
%% BS = Bootstrap_GC_TimeResolved(dataDir, subject, 'name',value, ...)
% Trial bootstrap for time-resolved GC. Produces confidence bands and tests
% for the two questions the circular-shift null CANNOT answer:
%
%   (1) CHANGE   -- did coupling change from baseline after the event?
%   (2) CONTRAST -- do two events differ (e.g. CHOICE vs OUTCOME)?
%
% WHY THE CIRCULAR-SHIFT NULL IS NOT ENOUGH
%   Run_GC_TimeResolved's 'nPerm' shifts the source in time and asks "is
%   there ANY coupling in this window". With ~30k pooled rows per fit that
%   question is answered yes almost everywhere -- on real data it fired in
%   55-67% of PRE-event windows, i.e. it is detecting the coupled background
%   rather than the task response. Your hypotheses are about CHANGES on top
%   of that background, which needs a different null.
%
%   Pooling trials is what makes the time resolution possible, but it leaves
%   one GC value per window per event and therefore no trial-level variance.
%   The bootstrap puts that variance back: resample trials with replacement,
%   refit everything, and the spread across resamples is the sampling
%   distribution of the pooled estimate.
%
% INPUT / PARAMETERS
%   Takes the SAME arguments as Run_GC_TimeResolved (winLen, step, tZero,
%   order, epochOpts, mode, band, ...) plus:
%   'nBoot'     [500]   bootstrap resamples. 200 is enough to look at, 1000
%                       for publication. Runtime scales linearly.
%   'baseline'  [1 375] SAMPLES, as in Run_GC_TimeResolved. Required here --
%                       the change test is defined against it.
%   'contrast'  {'CHOICE','OUTCOME'}  the two events to difference (B - A).
%                       Set {} to skip and only do the change test.
%   'alpha'     [0.05]  two-sided level for the CIs.
%   'seed'      []      rng seed for reproducibility.
%   'peakWin'   [0 Inf] SECONDS relative to tZero. Window inside which the
%                       peak and its latency are located. Defaults to
%                       everything post-event.
%
% OUTPUT struct BS
%   .t                          [1 x nWin]
%   .change.(evt).mean          [nChan nChan nWin]  post-minus-baseline
%   .change.(evt).lo/.hi        CI bounds
%   .change.(evt).sig           CI excludes zero
%   .change.(evt).p             bootstrap two-sided p
%   .contrast.mean/.lo/.hi/.sig/.p    same, for event B minus event A
%   .contrastEvents .nBoot .alpha .order .nodeLabel
%
% MULTIPLE COMPARISONS
%   .sig is pointwise -- one test per edge per window, uncorrected. With
%   nChan^2 edges x nWin windows that is thousands of tests. BS.sigCluster
%   gives a cluster-corrected version: contiguous runs of pointwise-
%   significant windows are compared against the largest run arising in the
%   bootstrap null. Report that, not the pointwise flags.
%
% EXAMPLE
%   BS = Bootstrap_GC_TimeResolved(path2sub,'CLASE001', ...
%          'winLen',250,'step',25,'tZero',501,'baseline',[1 375], ...
%          'order',10,'nBoot',500,'contrast',{'CHOICE','OUTCOME'}, ...
%          'epochOpts',{'winSamples',1500});
%   Plot_GC_Bootstrap(BS);

ip = inputParser; ip.KeepUnmatched = true;
ip.addParameter('events',{'CHOICE','RESPONSEON','OUTCOME'},@iscell);
ip.addParameter('winLen',250,@isscalar);
ip.addParameter('step',25,@isscalar);
ip.addParameter('tZero',1,@isscalar);
ip.addParameter('order',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('maxOrder',30,@isscalar);
ip.addParameter('baseline',[1 375],@(x)numel(x)==2);
ip.addParameter('nBoot',500,@isscalar);
ip.addParameter('contrast',{'CHOICE','OUTCOME'},@iscell);
ip.addParameter('alpha',0.05,@isscalar);
ip.addParameter('seed',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('peakWin',[0 Inf],@(x)numel(x)==2);   % SECONDS, for latency
ip.addParameter('mode','voltage',@ischar);
ip.addParameter('band',[],@(x)isempty(x)||numel(x)==2);
ip.addParameter('epochOpts',{},@iscell);
ip.addParameter('fs',500,@isscalar);
ip.addParameter('verbose',true,@(x)islogical(x)||isnumeric(x));
ip.parse(varargin{:});
o = ip.Results;
V = o.verbose;
if ~isempty(o.seed), rng(o.seed); end

epArgs = [{'verbose',false,'events',o.events,'fs',o.fs, ...
           'mode',o.mode,'band',o.band}, o.epochOpts];

if V
    fprintf('\n=== Bootstrap_GC_TimeResolved : %s ===\n', subject);
    fprintf('  %d resamples, alpha %.3f\n', o.nBoot, o.alpha);
end

%% ---- load, restricted to trials valid in EVERY event ----
ALLe = cell(numel(o.events),1); VALe = cell(numel(o.events),1);
for e = 1:numel(o.events)
    [~, info] = Get_GC_Epoch(dataDir, subject, o.events{e}, [], epArgs{:});
    ALLe{e} = info.all; VALe{e} = logical(info.valid(:));
end
nT0 = min(cellfun(@numel, VALe));
keep = true(nT0,1);
for e = 1:numel(o.events), keep = keep & VALe{e}(1:nT0); end
idx = find(keep);
nTrial = numel(idx);
if nTrial < 20
    error('Bootstrap_GC_TimeResolved:few','Only %d common trials.', nTrial);
end
EP = cell(numel(o.events),1);
for e = 1:numel(o.events), EP{e} = ALLe{e}(idx); end

nodeLabel = info.nodeLabel; nChan = numel(nodeLabel);
nSamp = min(cellfun(@(c) min(cellfun(@(x) size(x,2), c)), EP));
fsA = info.fsAnalysis;

starts = 1 : o.step : (nSamp - o.winLen + 1);
nWin = numel(starts);
tC = (starts + o.winLen/2 - o.tZero) / fsA;
bMask = tC >= (o.baseline(1)-o.tZero)/fsA & tC <= (o.baseline(2)-o.tZero)/fsA;
if ~any(bMask)
    error('Bootstrap_GC_TimeResolved:baseline','Baseline matched no window.');
end

p = o.order;
if isempty(p)
    mid = starts(max(1,round(nWin/2)));
    pool = cellfun(@(x) x(:, mid:mid+o.winLen-1), EP{1}, 'UniformOutput', false);
    p = Select_Model_Order_MT(pool, o.maxOrder, 'BIC');
end
if V
    fprintf('  %d common trials, %d windows, order %d, %d baseline windows\n', ...
            nTrial, nWin, p, sum(bMask));
end

%% ---- pre-slice every window of every event (done once) ----
% W{e}{w} is [nChan x winLen x nTrial]
W = cell(numel(o.events),1);
for e = 1:numel(o.events)
    W{e} = cell(nWin,1);
    for w = 1:nWin
        s = starts(w);
        X3 = zeros(nChan, o.winLen, nTrial);
        for q = 1:nTrial, X3(:,:,q) = EP{e}{q}(:, s:s+o.winLen-1); end
        W{e}{w} = X3;
    end
end
clear EP ALLe;

cIdx = [];
if ~isempty(o.contrast)
    [tf, cIdx] = ismember(o.contrast, o.events);
    if ~all(tf)
        error('Bootstrap_GC_TimeResolved:contrast', ...
              'Contrast events {%s} not in the events list.', ...
              strjoin(o.contrast(~tf),', '));
    end
end

%% ---- point estimate + bootstrap ----
nE = numel(o.events);
chgHat = zeros(nChan,nChan,nWin,nE);
chgBoot = zeros(nChan,nChan,nWin,nE,o.nBoot);
if ~isempty(cIdx)
    conHat  = zeros(nChan,nChan,nWin);
    conBoot = zeros(nChan,nChan,nWin,o.nBoot);
end
% peak amplitude / latency, recomputed inside every resample so that any
% between-edge comparison is PAIRED on the same trial draw and their shared
% trial-level noise cancels in the difference
pkMask = tC >= o.peakWin(1) & tC <= o.peakWin(2);
if ~any(pkMask), pkMask = tC >= 0; end
ampHat  = nan(nChan,nChan,nE);      ampBoot  = nan(nChan,nChan,nE,o.nBoot);
latHat  = nan(nChan,nChan,nE);      latBoot  = nan(nChan,nChan,nE,o.nBoot);
riseHat = nan(nChan,nChan,nE);      riseBoot = nan(nChan,nChan,nE,o.nBoot);

tAll = tic;
for b = 0:o.nBoot
    if b == 0
        sel = 1:nTrial;                 % b == 0 is the point estimate
    else
        sel = randi(nTrial, [nTrial 1]);
    end
    gcB = zeros(nChan,nChan,nWin,nE);
    for e = 1:nE
        for w = 1:nWin
            S = Compute_Spectral_GC(W{e}{w}(:,:,sel), p, fsA, 'nFreq', 9, ...
                    'nPerm', 0, 'zscore', true, 'verbose', false);
            gcB(:,:,w,e) = S.td;
        end
    end
    % baseline-subtract within this resample, so baseline uncertainty is
    % propagated rather than treated as a fixed known quantity
    ch = gcB - repmat(mean(gcB(:,:,bMask,:),3), [1 1 nWin 1]);
    % peak statistics for this resample
    aTmp = nan(nChan,nChan,nE); lTmp = aTmp; rTmp = aTmp;
    for e = 1:nE
        for i = 1:nChan
            for j = 1:nChan
                if i == j, continue; end
                [aTmp(i,j,e), lTmp(i,j,e), rTmp(i,j,e)] = ...
                    peak_stats(squeeze(ch(i,j,:,e))', tC, pkMask);
            end
        end
    end

    if b == 0
        chgHat = ch;
        ampHat = aTmp; latHat = lTmp; riseHat = rTmp;
        if ~isempty(cIdx), conHat = ch(:,:,:,cIdx(2)) - ch(:,:,:,cIdx(1)); end
    else
        chgBoot(:,:,:,:,b) = ch;
        ampBoot(:,:,:,b) = aTmp; latBoot(:,:,:,b) = lTmp; riseBoot(:,:,:,b) = rTmp;
        if ~isempty(cIdx)
            conBoot(:,:,:,b) = ch(:,:,:,cIdx(2)) - ch(:,:,:,cIdx(1));
        end
    end
    if V && b > 0 && mod(b, max(1,floor(o.nBoot/10))) == 0
        el = toc(tAll);
        fprintf('    boot %4d/%4d  %6.1fs  ~%.0fs left\n', ...
                b, o.nBoot, el, el*(o.nBoot-b)/b);
    end
end

%% ---- assemble ----
BS.t = tC; BS.nodeLabel = nodeLabel; BS.order = p; BS.nBoot = o.nBoot;
BS.alpha = o.alpha; BS.events = o.events; BS.nTrial = nTrial;
BS.winLen = o.winLen; BS.step = o.step; BS.fs = fsA; BS.subject = subject;
BS.baselineWin = find(bMask);
% Peak statistics. The BOOT arrays are kept in full (not just summarised)
% because comparing two edges requires the paired per-resample values --
% see Compare_GC_Edges.
BS.peak.window   = o.peakWin;
BS.peak.amp      = ampHat;    BS.peak.ampBoot  = ampBoot;
BS.peak.latency  = latHat;    BS.peak.latBoot  = latBoot;
BS.peak.rise50   = riseHat;   BS.peak.riseBoot = riseBoot;

lo = 100*o.alpha/2; hi = 100*(1-o.alpha/2);
for e = 1:nE
    D = squeeze(chgBoot(:,:,:,e,:));                 % nChan nChan nWin nBoot
    evt = o.events{e};
    BS.change.(evt).mean = chgHat(:,:,:,e);
    BS.change.(evt).lo   = prctile_(D, lo, 4);
    BS.change.(evt).hi   = prctile_(D, hi, 4);
    BS.change.(evt).sig  = (BS.change.(evt).lo > 0) | (BS.change.(evt).hi < 0);
    BS.change.(evt).p    = boot_p(D, chgHat(:,:,:,e));
end
if ~isempty(cIdx)
    BS.contrast.mean = conHat;
    BS.contrast.lo   = prctile_(conBoot, lo, 4);
    BS.contrast.hi   = prctile_(conBoot, hi, 4);
    BS.contrast.sig  = (BS.contrast.lo > 0) | (BS.contrast.hi < 0);
    BS.contrast.p    = boot_p(conBoot, conHat);
    BS.contrastEvents = o.contrast;
    BS.sigCluster = cluster_correct(BS.contrast.sig, conBoot, o.alpha);
end

if V
    fprintf('  done in %.1f s\n', toc(tAll));
    if ~isempty(cIdx)
        od = ~eye(nChan);
        sp = BS.contrast.sig; spv = sp(repmat(od,[1 1 nWin]));
        sc = BS.sigCluster;   scv = sc(repmat(od,[1 1 nWin]));
        fprintf('  contrast %s - %s: %.1f%% of edge-windows pointwise sig, ', ...
                o.contrast{2}, o.contrast{1}, 100*mean(spv));
        fprintf('%.1f%% after cluster correction\n', 100*mean(scv));
        fprintf(['  (correction is family-wise over %d edges x %d windows; ' ...
                 'report .sigCluster)\n'], nChan*(nChan-1), nWin);
    end
    fprintf('\n');
end
end

%% ================= helpers =================
function v = prctile_(X, pct, dim)
% percentile along dim, toolbox-free
X = sort(X, dim);
n = size(X, dim);
k = max(1, min(n, ceil(pct/100*n)));
idx = repmat({':'}, 1, ndims(X)); idx{dim} = k;
v = X(idx{:});
end

function P = boot_p(D, hat)
% two-sided bootstrap p: how often does the resample distribution cross zero
nB = size(D,4);
frac = sum(D <= 0, 4) / nB;
frac = min(frac, 1-frac);
P = 2*frac;
P(~isfinite(hat)) = NaN;
P = min(P, 1);
P = max(P, 1/nB);
end

function SC = cluster_correct(sigPt, D, alpha)
% Cluster-extent correction over time AND across edges.
%
% Observed statistic : length of each contiguous run of pointwise-significant
%                      windows, per directed edge.
% Null               : the SAME statistic from centred bootstrap resamples
%                      (centring removes the true effect, leaving noise with
%                      the right covariance). For each resample we take the
%                      LONGEST run found across EVERY edge, so the threshold
%                      controls the family-wise error rate over all edges and
%                      all windows jointly -- not just over time within one
%                      edge. With 20 edges x 51 windows, per-edge correction
%                      alone would still leave ~50 false positives at 0.05.
[n1,n2,nW,nB] = size(D);
SC = false(n1,n2,nW);
Dc = D - repmat(mean(D,4), [1 1 1 nB]);          % centre -> null
sd = std(D, 0, 4); sd(sd < eps) = Inf;
thr = 1.96;

% ---- null: max run length over ALL edges, per resample ----
nullMax = zeros(nB,1);
for b = 1:nB
    mx = 0;
    for i = 1:n1
        for j = 1:n2
            if i == j, continue; end
            zz = abs(squeeze(Dc(i,j,:,b))) ./ squeeze(sd(i,j,:));
            r = max_run(zz > thr);
            if r > mx, mx = r; end
        end
    end
    nullMax(b) = mx;
end
cut = prctile_(reshape(nullMax,1,1,1,nB), 100*(1-alpha), 4);

% ---- keep observed runs that beat the global threshold ----
for i = 1:n1
    for j = 1:n2
        if i == j, continue; end
        [runs, s0, e0] = runs_of(squeeze(sigPt(i,j,:)));
        for r = 1:numel(runs)
            if runs(r) >= cut, SC(i,j,s0(r):e0(r)) = true; end
        end
    end
end
end

function [amp, tPeak, tRise] = peak_stats(y, t, mask)
% Extremum by MAGNITUDE inside mask (so a negative-going effect is handled),
% its latency, and the 50%-of-peak rise latency.
%
% tRise is usually the better latency measure for comparing edges: argmax is
% unstable when a curve is broad or flat, whereas the time at which the
% response first reaches half its peak is far less sensitive to noise near
% the top. Both are returned; use whichever the data supports.
amp = NaN; tPeak = NaN; tRise = NaN;
yy = y; yy(~mask) = NaN;
if all(isnan(yy)), return; end
[~, k] = max(abs(yy));
amp = yy(k); tPeak = t(k);
if ~isfinite(amp) || amp == 0, return; end
half = 0.5*abs(amp); s = sign(amp);
tRise = t(k);
for q = k:-1:2
    if s*yy(q-1) < half || ~isfinite(yy(q-1))
        y1 = s*yy(q-1); y2 = s*yy(q);
        if isfinite(y1) && y2 > y1
            f = (half - y1) / (y2 - y1);
            tRise = t(q-1) + f*(t(q) - t(q-1));
        else
            tRise = t(q);
        end
        return;
    end
end
tRise = t(find(mask,1,'first'));
end

function m = max_run(v)
m = 0; c = 0;
for k = 1:numel(v)
    if v(k), c = c+1; if c > m, m = c; end, else, c = 0; end
end
end

function [len, s0, e0] = runs_of(v)
len = []; s0 = []; e0 = []; k = 1; n = numel(v);
while k <= n
    if v(k)
        j = k; while j < n && v(j+1), j = j+1; end
        len(end+1) = j-k+1; s0(end+1) = k; e0(end+1) = j; %#ok<AGROW>
        k = j+1;
    else
        k = k+1;
    end
end
end
