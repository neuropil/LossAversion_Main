function summary = Plot_GC_Summary(source, varargin)
%% summary = Plot_GC_Summary(source, 'name',value, ...)
% Summary + interpretation panel for ONE subject's GC results file
% (the <subject>_GC.mat written by Run_GC_Pipeline, containing `subjResult`).
%
% DIRECTION CONVENTION (verified on ground-truth simulation):
%   GC(i,j) = influence FROM j (source) TO i (target).
%   All heatmaps are therefore plotted rows = TARGET, cols = SOURCE.
%
% INPUT
%   source : path to <subject>_GC.mat, OR the subjResult struct itself.
%
% OPTIONS (defaults in brackets)
%   'events'   []        events to show; [] = all present in the file
%   'alpha'    [0.05]    per-epoch significance level used by the pipeline
%   'contrast' []        {eventA,eventB} for the change panel; [] = {first,last}
%   'nPerm'    [5000]    sign-flip permutations for the contrast test
%   'fdrQ'     [0.05]    FDR level for both the consistency and contrast tests
%   'topN'     [6]       # of edges listed in the text panel
%   'saveFig'  ''        if non-empty, print the figure to this path (.png)
%   'visible'  ['on']    'off' for headless batch use
%
% OUTPUT struct `summary` (all numbers behind the figure)
%   .nodeLabel, .order, .rcondMin
%   .event.(evt): .meanGC .semGC .sigRate .nValid .pConsistent .hConsistent
%   .contrast: .events .diff .p .h .nPaired
%   .byType.(evt): .intraL .intraR .inter  (mean +/- sem over pairs)
%   .topEdges, .topChanges
%
% STATISTICAL NOTES (these matter for interpretation)
%  - Mean GC is a biased-upward estimate; the bias depends on epoch length and
%    model order. Comparing events is only meaningful because Run_GC_Pipeline
%    fixes BOTH across events. Do not compare these values to another dataset
%    analysed with a different window or order.
%  - The pipeline's per-epoch .sig flags are UNCORRECTED. Here, each pair gets a
%    binomial test (is its significant-epoch rate above the alpha baseline?)
%    and those are FDR-corrected across pairs. Interpret .hConsistent, not the
%    raw sig rate.
%  - The contrast uses a paired sign-flip permutation test over trials valid in
%    BOTH events, so it is robust to the GC bias (which cancels in the paired
%    difference) and needs no distributional assumption.

%% ---------------- parse ----------------
ip = inputParser;
ip.addRequired('source');
ip.addParameter('events',[]);
ip.addParameter('alpha',0.05,@isscalar);
ip.addParameter('contrast',[]);
ip.addParameter('nPerm',5000,@isscalar);
ip.addParameter('fdrQ',0.05,@isscalar);
ip.addParameter('topN',6,@isscalar);
ip.addParameter('saveFig','',@ischar);
ip.addParameter('visible','on',@ischar);
ip.parse(source,varargin{:});
opt = ip.Results;

% ---- load ----
if ischar(source)
    S = load(source);
    if isfield(S,'subjResult'), R = S.subjResult;
    else
        fn = fieldnames(S); R = S.(fn{1});
    end
    [~,subName] = fileparts(source);
else
    R = source; subName = 'subject';
end
subName = strrep(subName,'_',' ');

if isempty(opt.events)
    events = fieldnames(R.events);
else
    events = opt.events;
    if ischar(events), events = {events}; end
end
events = events(:)';
nEvent = numel(events);

labels    = R.conn.nodeLabel(:);
dispLab   = strrep(labels,'_','-');       % '_' would render as subscript
nChan     = numel(labels);
offDiag   = ~eye(nChan);

summary.nodeLabel = labels;
summary.order     = R.order;
summary.rcondMin  = R.rcondMin;

%% ---------------- per-event stats ----------------
meanGC = cell(nEvent,1); sigRate = cell(nEvent,1);
for e = 1:nEvent
    evt = events{e};
    G   = R.events.(evt).GC;
    Sg  = R.events.(evt).sig;
    if isfield(R.events.(evt),'validEpoch')
        v = logical(R.events.(evt).validEpoch(:));
    else
        v = squeeze(any(any(~isnan(G),1),2));   % fall back: non-empty slices
    end
    Gv = G(:,:,v); Sv = Sg(:,:,v);
    n  = sum(v);

    m  = nanmean_(Gv,3);
    sd = nanstd_(Gv,3);
    se = sd ./ sqrt(max(n,1));
    sr = nanmean_(double(Sv),3);

    % per-pair consistency: is the sig-epoch rate above the alpha baseline?
    k  = round(sr * n);
    pC = nan(nChan);
    for i = 1:nChan
        for j = 1:nChan
            if i==j, continue; end
            pC(i,j) = binom_upper_tail(k(i,j), n, opt.alpha);
        end
    end
    hC = false(nChan);
    hC(offDiag) = bh_fdr(pC(offDiag), opt.fdrQ);

    meanGC{e} = m; sigRate{e} = sr;
    summary.event.(evt).meanGC      = m;
    summary.event.(evt).semGC       = se;
    summary.event.(evt).sigRate     = sr;
    summary.event.(evt).nValid      = n;
    summary.event.(evt).pConsistent = pC;
    summary.event.(evt).hConsistent = hC;
end

%% ---------------- contrast (the "change across task" question) ----------------
doContrast = nEvent >= 2;
if doContrast
    if isempty(opt.contrast), cEv = {events{1}, events{end}};
    else,                     cEv = opt.contrast; end
    A = R.events.(cEv{1}); B = R.events.(cEv{2});
    if isfield(A,'validEpoch'), vA = logical(A.validEpoch(:)); else, vA = true(size(A.GC,3),1); end
    if isfield(B,'validEpoch'), vB = logical(B.validEpoch(:)); else, vB = true(size(B.GC,3),1); end
    nT = min(numel(vA), numel(vB));
    both = vA(1:nT) & vB(1:nT);                 % trials usable in BOTH events
    D  = B.GC(:,:,both) - A.GC(:,:,both);       % paired difference, bias cancels
    nP = sum(both);

    dMean = nanmean_(D,3);
    pD = nan(nChan);
    for i = 1:nChan
        for j = 1:nChan
            if i==j, continue; end
            pD(i,j) = signflip_p(squeeze(D(i,j,:)), opt.nPerm);
        end
    end
    hD = false(nChan);
    hD(offDiag) = bh_fdr(pD(offDiag), opt.fdrQ);

    summary.contrast.events  = cEv;
    summary.contrast.diff    = dMean;
    summary.contrast.p       = pD;
    summary.contrast.h       = hD;
    summary.contrast.nPaired = nP;
end

%% ---------------- connection-type summary ----------------
mkT = {'intraL','intraR','inter'};
for e = 1:nEvent
    evt = events{e}; m = meanGC{e};
    for t = 1:numel(mkT)
        msk = R.conn.(mkT{t}) & offDiag;
        vals = m(msk); vals = vals(~isnan(vals));
        if isempty(vals)
            summary.byType.(evt).(mkT{t}) = [NaN NaN 0];
        else
            summary.byType.(evt).(mkT{t}) = [mean(vals), std(vals)/sqrt(numel(vals)), numel(vals)];
        end
    end
end

%% ---------------- top edges / top changes ----------------
e1 = events{1};
[summary.topEdges] = rank_edges(summary.event.(e1).meanGC, summary.event.(e1).hConsistent, ...
                                labels, opt.topN, 'desc');
if doContrast
    [summary.topChanges] = rank_edges(summary.contrast.diff, summary.contrast.h, ...
                                      labels, opt.topN, 'abs');
end

%% ---------------- figure ----------------
nCol = max(nEvent,3);
fig = figure('visible',opt.visible,'Color','w','Position',[80 80 460*nCol 1150]);

cmapSeq = seq_map(256);
cmapDiv = div_map(256);

% ---- Row 1: mean GC per event (shared colour scale) ----
allM = []; for e=1:nEvent, allM = [allM; meanGC{e}(offDiag)]; end %#ok<AGROW>
allM = allM(~isnan(allM));
cmax = prctile_(allM, 98); if ~isfinite(cmax) || cmax<=0, cmax = 1; end
for e = 1:nEvent
    subplot(3,nCol,e);
    M = meanGC{e}; M(1:nChan+1:end) = NaN;
    draw_mat(M, dispLab, [0 cmax], cmapSeq);
    title(sprintf('%s  mean GC  (n=%d)', events{e}, summary.event.(events{e}).nValid), ...
          'FontWeight','bold');
    if e==1, ylabel('TARGET  (to)','FontWeight','bold'); end
    xlabel('SOURCE  (from)');
end

% ---- Row 2: consistency (sig rate, FDR-surviving pairs outlined) ----
for e = 1:nEvent
    subplot(3,nCol,nCol+e);
    SR = sigRate{e}; SR(1:nChan+1:end) = NaN;
    draw_mat(SR, dispLab, [0 1], cmapSeq);
    hold on; outline_cells(summary.event.(events{e}).hConsistent);
    title(sprintf('%s  sig-epoch rate', events{e}));
    if e==1, ylabel('TARGET  (to)','FontWeight','bold'); end
    xlabel('SOURCE  (from)');
end

% ---- Row 3a: contrast ----
subplot(3,nCol,2*nCol+1);
if doContrast
    Dm = summary.contrast.diff; Dm(1:nChan+1:end) = NaN;
    lim = max(abs(Dm(offDiag & ~isnan(Dm)))); if isempty(lim)||lim<=0, lim=1; end
    draw_mat(Dm, dispLab, [-lim lim], cmapDiv);
    hold on; outline_cells(summary.contrast.h);
    title(sprintf('%s minus %s  (n=%d paired)', ...
        summary.contrast.events{2}, summary.contrast.events{1}, summary.contrast.nPaired), ...
        'FontWeight','bold');
    ylabel('TARGET  (to)','FontWeight','bold'); xlabel('SOURCE  (from)');
else
    axis off; text(0.1,0.5,'contrast needs >= 2 events','FontSize',11);
end

% ---- Row 3b: connection-type bars ----
subplot(3,nCol,2*nCol+2);
Mb = zeros(nEvent,3); Eb = zeros(nEvent,3);
for e=1:nEvent
    for t=1:3
        v = summary.byType.(events{e}).(mkT{t});
        Mb(e,t)=v(1); Eb(e,t)=v(2);
    end
end
hb = bar(Mb); hold on;
nb = size(Mb,2);
for t=1:nb
    if verLessThan_safe()
        xc = (1:nEvent) + (t-(nb+1)/2)*(0.8/nb);
    else
        try, xc = hb(t).XEndPoints; catch, xc = (1:nEvent) + (t-(nb+1)/2)*(0.8/nb); end
    end
    for e = 1:nEvent      % manual error bars (errorbar() arg forms differ by version)
        if isfinite(Eb(e,t)) && Eb(e,t)>0
            plot([xc(e) xc(e)], [Mb(e,t)-Eb(e,t), Mb(e,t)+Eb(e,t)], 'k-','LineWidth',1);
            w = 0.06;
            plot([xc(e)-w xc(e)+w], [Mb(e,t)-Eb(e,t) Mb(e,t)-Eb(e,t)], 'k-','LineWidth',1);
            plot([xc(e)-w xc(e)+w], [Mb(e,t)+Eb(e,t) Mb(e,t)+Eb(e,t)], 'k-','LineWidth',1);
        end
    end
end
set(gca,'XTick',1:nEvent,'XTickLabel',events);
ylabel('mean GC'); title('by connection type (+/- SEM over pairs)');
nL = sum(strcmp({R.channelInfo.hemi},'L')); nR = sum(strcmp({R.channelInfo.hemi},'R'));
legend({sprintf('intra-L (%d)',sum(R.conn.intraL(:))), ...
        sprintf('intra-R (%d)',sum(R.conn.intraR(:))), ...
        sprintf('inter (%d)',  sum(R.conn.inter(:)))}, 'Location','best');

% ---- Row 3c: diagnostics + interpretation text ----
subplot(3,nCol,2*nCol+3); axis off;
tx = {};
tx{end+1} = sprintf('%s', subName);
tx{end+1} = sprintf('nodes: %d   (L=%d, R=%d)', nChan, nL, nR);
tx{end+1} = sprintf('VAR order: %d', R.order);
tx{end+1} = sprintf('min rcond: %.1e%s', R.rcondMin, tern(R.rcondMin<1e-6,'   <-- ill-conditioned',''));
tx{end+1} = '';
tx{end+1} = 'valid epochs:';
for e=1:nEvent
    tx{end+1} = sprintf('   %-12s %d', events{e}, summary.event.(events{e}).nValid); %#ok<AGROW>
end
tx{end+1} = '';
tx{end+1} = sprintf('strongest edges (%s, FDR q=%.2f):', events{1}, opt.fdrQ);
for i=1:min(opt.topN,numel(summary.topEdges))
    tx{end+1} = sprintf('   %s', summary.topEdges{i}); %#ok<AGROW>
end
if doContrast
    tx{end+1} = '';
    tx{end+1} = sprintf('largest changes (%s-%s):', ...
        summary.contrast.events{2}, summary.contrast.events{1});
    for i=1:min(opt.topN,numel(summary.topChanges))
        tx{end+1} = sprintf('   %s', summary.topChanges{i}); %#ok<AGROW>
    end
end
tx{end+1} = '';
tx{end+1} = 'black outline = survives FDR';
tx{end+1} = 'matrix(i,j) = from j (col) to i (row)';
text(0,1,tx,'VerticalAlignment','top','FontName','FixedWidth','FontSize',9, ...
     'Interpreter','none');

if ~isempty(opt.saveFig)
    try
        print(fig, opt.saveFig, '-dpng', '-r120');
    catch err
        warning('Plot_GC_Summary:save','Could not save figure: %s', err.message);
    end
end
end

%% ======================= local helpers =======================

function draw_mat(M, labs, clim, cmap)
imagesc(M, clim); axis square; colormap(gca, cmap); colorbar;
n = size(M,1);
set(gca,'XTick',1:n,'XTickLabel',labs,'YTick',1:n,'YTickLabel',labs, ...
        'TickLength',[0 0],'FontSize',8);
try, xtickangle(45); catch, end
% grey out the diagonal so it reads as "not computed"
hold on;
for i=1:n, plot(i,i,'s','MarkerSize',10,'MarkerFaceColor',[.75 .75 .75], ...
                'MarkerEdgeColor',[.75 .75 .75]); end
end

function outline_cells(H)
[ii,jj] = find(H);
for k = 1:numel(ii)
    i = ii(k); j = jj(k);
    plot([j-.5 j+.5 j+.5 j-.5 j-.5], [i-.5 i-.5 i+.5 i+.5 i-.5], 'k-', 'LineWidth',1.6);
end
end

function m = nanmean_(X,dim)
N = sum(~isnan(X),dim); X(isnan(X)) = 0;
m = sum(X,dim) ./ max(N,1); m(N==0) = NaN;
end

function s = nanstd_(X,dim)
m = nanmean_(X,dim);
D = X - m; D(isnan(X)) = 0;
N = sum(~isnan(X),dim);
s = sqrt(sum(D.^2,dim) ./ max(N-1,1)); s(N<2) = NaN;
end

function p = binom_upper_tail(k, n, p0)
% P(X >= k), X ~ Binomial(n,p0); toolbox-free via gammaln
if n<=0, p = NaN; return; end
k = max(0,min(k,n));
if k==0, p = 1; return; end
ks = k:n;
lg = gammaln(n+1) - gammaln(ks+1) - gammaln(n-ks+1) ...
     + ks*log(max(p0,eps)) + (n-ks)*log(max(1-p0,eps));
p = min(1, sum(exp(lg)));
end

function h = bh_fdr(p, q)
% Benjamini-Hochberg; returns logical vector aligned with p
p = p(:); h = false(size(p));
ok = ~isnan(p); pv = p(ok);
if isempty(pv), return; end
[ps, idx] = sort(pv); m = numel(ps);
thr = (1:m)'/m * q;
below = find(ps <= thr, 1, 'last');
if ~isempty(below)
    keep = false(m,1); keep(1:below) = true;
    tmp = false(m,1); tmp(idx) = keep;
    h(ok) = tmp;
end
end

function p = signflip_p(d, nPerm)
% two-sided paired sign-flip permutation test on the mean of d
d = d(~isnan(d));
n = numel(d);
if n < 3, p = NaN; return; end
obs = abs(mean(d));
cnt = 0;
for it = 1:nPerm
    s = sign(randn(n,1)); s(s==0) = 1;
    if abs(mean(d.*s)) >= obs, cnt = cnt + 1; end
end
p = (1 + cnt) / (nPerm + 1);
end

function out = rank_edges(M, H, labels, topN, mode)
n = size(M,1); rows = {};
vals = []; ii = []; jj = [];
for i=1:n
    for j=1:n
        if i==j || isnan(M(i,j)), continue; end
        vals(end+1) = M(i,j); ii(end+1)=i; jj(end+1)=j; %#ok<AGROW>
    end
end
if isempty(vals), out = {}; return; end
switch mode
    case 'abs',  [~,ord] = sort(abs(vals),'descend');
    otherwise,   [~,ord] = sort(vals,'descend');
end
ord = ord(1:min(topN,numel(ord)));
for k = 1:numel(ord)
    a = ord(k); i = ii(a); j = jj(a);
    star = ''; if ~isempty(H) && H(i,j), star = ' *'; end
    rows{end+1,1} = sprintf('%s -> %s  %+.3f%s', labels{j}, labels{i}, vals(a), star); %#ok<AGROW>
end
out = rows;
end

function c = seq_map(n)
% white -> deep blue (sequential)
t = linspace(0,1,n)';
c = [1-0.85*t, 1-0.75*t, 1-0.35*t];
end

function c = div_map(n)
% blue <- white -> red (diverging)
h = floor(n/2); t = linspace(0,1,h)';
lower = [t, t, ones(h,1)];
upper = [ones(n-h,1), flipud(linspace(0,1,n-h)'), flipud(linspace(0,1,n-h)')];
c = [lower; upper];
end

function v = prctile_(x, pct)
x = sort(x(~isnan(x)));
if isempty(x), v = NaN; return; end
idx = max(1, min(numel(x), ceil(pct/100*numel(x))));
v = x(idx);
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end

function tf = verLessThan_safe()
tf = exist('OCTAVE_VERSION','builtin') > 0;
end