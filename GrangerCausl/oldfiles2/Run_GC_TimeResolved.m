function TR = Run_GC_TimeResolved(dataDir, subject, varargin)
%% TR = Run_GC_TimeResolved(dataDir, subject, 'name',value, ...)
% Sliding-window Granger causality: GC(t) per directed edge, and GC(t,f) if
% you want the time-frequency version.
%
% HOW IT WORKS, AND WHY IT IS POOLED
%   A window short enough to resolve time (say 500 ms = 250 samples) is far
%   too short to fit a VAR per trial. So every window pools ALL valid trials
%   as separate realizations -- the lag matrix is built within each trial and
%   only then concatenated, so no design row crosses a trial boundary. With
%   ~128 trials a 250-sample window yields ~31k pooled rows, which is ample.
%   The null bias floor is p/m with m = nTrial*(winLen-p), i.e. tiny, so the
%   time course is dominated by signal rather than by window-length bias.
%
%   One consequence worth stating in a methods section: this is a
%   TRIAL-AVERAGED time course. There is one GC value per window, not one per
%   trial, so trial-level statistics are not available. Significance comes
%   from the circular-shift null instead (set 'nPerm').
%
% INPUT
%   dataDir  folder of *_TrialDATA.mat
%   subject  e.g. 'CLASE001'
%
% PARAMETERS
%   'events'    {'CHOICE','RESPONSEON','OUTCOME'}
%   'winLen'    [250]  window length in SAMPLES (500 ms at fs=500)
%   'step'      [25]   hop in SAMPLES (50 ms at fs=500)
%   'tZero'     [1]    sample index of event onset, for the time axis
%   'order'     []     VAR order; [] = select once from the middle window
%   'maxOrder'  [30]
%   'spectral'  [false] also keep GC(t,f). Costs memory: nChan^2*nFreq*nWin
%   'nFreq'     [65]
%   'bands'     []     band definitions when 'spectral' is on
%   'nPerm'     [0]    circular-shift surrogates PER WINDOW. Expensive:
%                      multiplies runtime by (1+nPerm). Start at 0.
%   'baseline'  []     [t1 t2] in SAMPLES relative to epoch start. If given,
%                      TR.gcBase holds GC minus the mean over that range.
%   'mode'      ['voltage'] | 'bandpower'
%   'band'      []     required for mode 'bandpower'
%   'epochOpts' {}     extra args passed to Get_GC_Epoch (winSamples etc.)
%   'verbose'   [true]
%
% OUTPUT struct TR
%   .t            [1 x nWin]  window CENTRE in seconds relative to tZero
%   .gc           [nChan x nChan x nWin x nEvent]   (i,j) = FROM j TO i
%   .gcBase       same, baseline-subtracted (if 'baseline' given)
%   .pval .sig    same size (if nPerm > 0)
%   .spec         [nChan x nChan x nFreq x nWin x nEvent] (if 'spectral')
%   .f .events .nodeLabel .order .winLen .step .nTrial .rcondMin

ip = inputParser;
ip.addParameter('events',{'CHOICE','RESPONSEON','OUTCOME'},@iscell);
ip.addParameter('winLen',250,@isscalar);
ip.addParameter('step',25,@isscalar);
ip.addParameter('tZero',1,@isscalar);
ip.addParameter('order',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('maxOrder',30,@isscalar);
ip.addParameter('spectral',false,@(x)islogical(x)||isnumeric(x));
ip.addParameter('nFreq',65,@isscalar);
ip.addParameter('bands',[],@(x)isempty(x)||isstruct(x));
ip.addParameter('nPerm',0,@isscalar);
ip.addParameter('baseline',[],@(x)isempty(x)||numel(x)==2);
ip.addParameter('mode','voltage',@ischar);
ip.addParameter('band',[],@(x)isempty(x)||numel(x)==2);
ip.addParameter('epochOpts',{},@iscell);
ip.addParameter('fs',500,@isscalar);
ip.addParameter('verbose',true,@(x)islogical(x)||isnumeric(x));
ip.parse(varargin{:});
o = ip.Results;
V = o.verbose;

if isempty(o.bands)
    o.bands = struct('theta',[4 8],'alpha',[8 12],'beta',[13 30], ...
                     'lowGamma',[30 70],'highGamma',[70 150]);
end
bfn = fieldnames(o.bands);
for b = 1:numel(bfn)
    if o.bands.(bfn{b})(2) >= o.fs/2, o.bands = rmfield(o.bands,bfn{b}); end
end

epArgs = [{'verbose',false,'events',o.events,'fs',o.fs, ...
           'mode',o.mode,'band',o.band}, o.epochOpts];

if V
    fprintf('\n=== Run_GC_TimeResolved : %s ===\n', subject);
    fprintf('  window %d samp (%.0f ms), step %d samp (%.0f ms)\n', ...
            o.winLen, 1000*o.winLen/o.fs, o.step, 1000*o.step/o.fs);
end

%% ---- load every event's epochs once ----
EP = cell(numel(o.events),1);
for e = 1:numel(o.events)
    [~, info] = Get_GC_Epoch(dataDir, subject, o.events{e}, [], epArgs{:});
    EP{e} = info.all(info.valid);
    if V, fprintf('  %-12s %d usable trials\n', o.events{e}, numel(EP{e})); end
end
nodeLabel = info.nodeLabel;
nChan = numel(nodeLabel);
nSamp = size(EP{1}{1},2);
fsA   = info.fsAnalysis;

%% ---- window grid ----
starts = 1 : o.step : (nSamp - o.winLen + 1);
nWin = numel(starts);
tCent = (starts + o.winLen/2 - o.tZero) / fsA;
if V
    fprintf('  %d windows spanning %.2f to %.2f s (fs %g)\n', ...
            nWin, tCent(1), tCent(end), fsA);
end

%% ---- model order: pick once, from the middle window of the first event ----
p = o.order;
if isempty(p)
    mid = starts(max(1,round(nWin/2)));
    pool = cellfun(@(x) x(:, mid:mid+o.winLen-1), EP{1}, 'UniformOutput', false);
    [p, oi] = Select_Model_Order_MT(pool, o.maxOrder, 'BIC');
    if V
        fprintf('  order %d (identifiable limit %d, from %d trials x %d samp)\n', ...
                p, oi.maxIdentifiable, numel(pool), o.winLen);
    end
elseif V
    fprintf('  order %d (fixed)\n', p);
end

mPool = numel(EP{1}) * (o.winLen - p);
if V
    fprintf('  pooled rows per window ~%d -> null floor p/m = %.6f\n', ...
            mPool, p/mPool);
    fprintf('  %d windows x %d events x %d pairs x %d = %.0fk VAR solves\n', ...
            nWin, numel(o.events), nChan*(nChan-1), 1+o.nPerm, ...
            nWin*numel(o.events)*nChan*(nChan-1)*(1+o.nPerm)/1000);
end

%% ---- sweep ----
TR.gc   = nan(nChan,nChan,nWin,numel(o.events));
TR.pval = nan(nChan,nChan,nWin,numel(o.events));
TR.sig  = false(nChan,nChan,nWin,numel(o.events));
if o.spectral
    TR.spec = nan(nChan,nChan,o.nFreq,nWin,numel(o.events));
end
rcondMin = Inf; fAxis = [];

for e = 1:numel(o.events)
    tE = tic;
    for w = 1:nWin
        s = starts(w);
        X3 = zeros(nChan, o.winLen, numel(EP{e}));
        for q = 1:numel(EP{e})
            X3(:,:,q) = EP{e}{q}(:, s:s+o.winLen-1);
        end
        S = Compute_Spectral_GC(X3, p, fsA, 'nFreq',o.nFreq, 'bands',o.bands, ...
                'nPerm',o.nPerm, 'zscore',true, 'verbose',false);
        TR.gc(:,:,w,e) = S.td;
        if isfield(S,'pval')
            TR.pval(:,:,w,e) = S.pval; TR.sig(:,:,w,e) = S.sig;
        end
        if o.spectral, TR.spec(:,:,:,w,e) = S.spec; fAxis = S.f; end
        rcondMin = min(rcondMin, S.rcondMin);

        if V && mod(w, max(1,floor(nWin/10))) == 0
            el = toc(tE);
            fprintf('    %-12s %3d/%3d windows  %5.1fs  ~%.0fs left\n', ...
                    o.events{e}, w, nWin, el, el*(nWin-w)/w);
        end
    end
end

%% ---- baseline ----
if ~isempty(o.baseline)
    bIdx = tCent >= (o.baseline(1)-o.tZero)/fsA & tCent <= (o.baseline(2)-o.tZero)/fsA;
    if ~any(bIdx)
        warning('Run_GC_TimeResolved:baseline','Baseline range matched no window.');
    else
        TR.gcBase = TR.gc - repmat(mean(TR.gc(:,:,bIdx,:),3), [1 1 nWin 1]);
        TR.baselineWin = find(bIdx);
    end
end

TR.t = tCent; TR.f = fAxis; TR.events = o.events; TR.nodeLabel = nodeLabel;
TR.order = p; TR.winLen = o.winLen; TR.step = o.step; TR.fs = fsA;
TR.nTrial = cellfun(@numel, EP); TR.rcondMin = rcondMin;
TR.nullFloor = p/mPool; TR.subject = subject; TR.mode = o.mode;
if strcmpi(o.mode,'bandpower'), TR.band = o.band; end

if V
    fprintf('  done. min rcond %.2e, null floor %.6f\n\n', rcondMin, TR.nullFloor);
end
end
