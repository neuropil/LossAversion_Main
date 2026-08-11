function out = Plot_Spectral_GC(source, varargin)
%% out = Plot_Spectral_GC(source, 'name',value, ...)
% Viewer for the frequency-resolved output written by Run_GC_Pipeline when
% run with 'mode','spectral' (or 'bandpower' with pooled fits). The existing
% Plot_GC_Summary / Plot_GC_Contrast2 only read the scalar .GC cube, so this
% is what you use to actually look at GC(f).
%
% Row 1  band-averaged GC heatmaps, one per band, shared colour scale.
%        Rows = TARGET, cols = SOURCE, matching the pipeline convention.
% Row 2  GC(f) spectra for the top-K edges, overlaid across events.
% Row 3  per-band paired contrast between two events + diagnostics.
%
% INPUT
%   source : path to <subject>_GC.mat, or the subjResult struct
%
% OPTIONS
%   'events'   []       events to show; [] = all present
%   'contrast' []       {evtA,evtB} for row 3; [] = {first,last}
%   'topK'     [6]      edges shown as spectra
%   'fdrQ'     [0.05]
%   'nPerm'    [5000]   sign-flip permutations for the band contrast
%   'fMax'     []       upper frequency limit for the spectra plots
%   'saveFig'  ''
%   'visible'  ['on']
%
% OUTPUT struct `out`
%   .band.(evt).(bandName)  mean band GC matrix
%   .contrast.(bandName)    [meanDiff ci p h] per band
%   .topEdges               labels of the plotted edges

ip = inputParser;
ip.addRequired('source');
ip.addParameter('events',[]);
ip.addParameter('contrast',[]);
ip.addParameter('topK',6,@isscalar);
ip.addParameter('fdrQ',0.05,@isscalar);
ip.addParameter('nPerm',5000,@isscalar);
ip.addParameter('fMax',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('saveFig','',@ischar);
ip.addParameter('visible','on',@ischar);
ip.parse(source,varargin{:});
opt = ip.Results;

if ischar(source)
    S = load(source);
    if isfield(S,'subjResult'), R = S.subjResult;
    else, fn = fieldnames(S); R = S.(fn{1}); end
    [~,subName] = fileparts(source);
else
    R = source; subName = 'subject';
end
subName = strrep(subName,'_',' ');

if isempty(opt.events), events = fieldnames(R.events);
else, events = opt.events; if ischar(events), events = {events}; end
end
events = events(:)';
nEvent = numel(events);

% require the spectral fields
if ~isfield(R.events.(events{1}),'spec')
    error('Plot_Spectral_GC:noSpec', ...
        ['No .spec field. Re-run Run_GC_Pipeline with ''mode'',''spectral'' ' ...
         '(this file looks like mode ''%s'').'], getfielddef(R,'mode','voltage'));
end

labels  = R.conn.nodeLabel(:);
dispLab = strrep(labels,'_','-');
nChan   = numel(labels);
offd    = ~eye(nChan);
f       = R.events.(events{1}).f;
fMax    = opt.fMax; if isempty(fMax), fMax = max(f); end
bandNames = fieldnames(R.events.(events{1}).bandGC);
nBand   = numel(bandNames);

% ---- band means per event ----
for e = 1:nEvent
    evt = events{e};
    v = logical(R.events.(evt).validEpoch(:));
    for b = 1:nBand
        C = R.events.(evt).bandGC.(bandNames{b});
        out.band.(evt).(bandNames{b}) = nanmean_(C(:,:,v),3);
    end
end

% ---- rank edges by band GC in the first event, summed over bands ----
tot = zeros(nChan);
for b = 1:nBand, tot = tot + out.band.(events{1}).(bandNames{b}); end
tot(1:nChan+1:end) = NaN;
[ii,jj] = deal([]);
for i=1:nChan, for j=1:nChan, if i~=j, ii(end+1)=i; jj(end+1)=j; end, end, end %#ok<AGROW>
lin  = sub2ind([nChan nChan], ii, jj);
[~,ord] = sort(tot(lin),'descend');
K = min(opt.topK, numel(ord)); sel = ord(1:K);
out.topEdges = arrayfun(@(k) sprintf('%s -> %s', dispLab{jj(sel(k))}, ...
                        dispLab{ii(sel(k))}), 1:K, 'UniformOutput',false)';

%% ================= figure =================
nCol = max([nBand, K, 3]);
fig = figure('visible',opt.visible,'Color','w','Position',[60 60 260*nCol 950]);
cmapSeq = seq_map(256);
evCols  = lines_(nEvent);

% ---- Row 1: band heatmaps (first event) ----
allv = [];
for b = 1:nBand
    M = out.band.(events{1}).(bandNames{b}); allv = [allv; M(offd)]; %#ok<AGROW>
end
allv = allv(isfinite(allv));
cmax = prctile_(allv,98); if ~isfinite(cmax)||cmax<=0, cmax = 1; end
for b = 1:nBand
    subplot(3,nCol,b);
    M = out.band.(events{1}).(bandNames{b}); M(1:nChan+1:end) = NaN;
    imagesc(M,[0 cmax]); axis square; colormap(gca,cmapSeq);
    set(gca,'XTick',1:nChan,'XTickLabel',dispLab,'YTick',1:nChan, ...
            'YTickLabel',dispLab,'TickLength',[0 0],'FontSize',7);
    try, xtickangle(45); catch, end
    title(sprintf('%s',bandNames{b}),'FontWeight','bold','FontSize',9);
    if b==1, ylabel('TARGET (to)','FontWeight','bold'); end
    if b==nBand, colorbar; end
end

% ---- Row 2: GC(f) for top edges, events overlaid ----
for k = 1:K
    subplot(3,nCol,nCol+k); hold on;
    for e = 1:nEvent
        sp = squeeze(R.events.(events{e}).spec(ii(sel(k)), jj(sel(k)), :));
        plot(f, sp, '-', 'Color', evCols(e,:), 'LineWidth', 1.4);
    end
    xlim([0 fMax]); box on; grid on;
    title(out.topEdges{k},'FontSize',8,'FontWeight','normal');
    xlabel('Hz','FontSize',8);
    if k==1
        ylabel('spectral GC','FontWeight','bold','FontSize',8);
        legend(events,'Location','best','FontSize',7); legend boxoff;
    end
end

% ---- Row 3a: per-band contrast ----
if nEvent >= 2
    if isempty(opt.contrast), cEv = {events{1}, events{end}};
    else, cEv = opt.contrast; end
    vA = logical(R.events.(cEv{1}).validEpoch(:));
    vB = logical(R.events.(cEv{2}).validEpoch(:));
    nT = min(numel(vA),numel(vB));
    both = vA(1:nT) & vB(1:nT); nP = sum(both);

    mu = zeros(nBand,1); ci = zeros(nBand,1); pv = nan(nBand,1);
    for b = 1:nBand
        CA = R.events.(cEv{1}).bandGC.(bandNames{b})(:,:,both);
        CB = R.events.(cEv{2}).bandGC.(bandNames{b})(:,:,both);
        D  = CB - CA;
        d  = zeros(nP,1);
        for q = 1:nP, sl = D(:,:,q); d(q) = nanmean_(sl(offd),1); end
        mu(b) = nanmean_(d,1);
        ci(b) = 1.96*nanstd_(d,1)/sqrt(max(nP,1));
        pv(b) = signflip_p(d, opt.nPerm);
    end
    hv = bh_fdr(pv, opt.fdrQ);
    for b = 1:nBand
        out.contrast.(bandNames{b}) = [mu(b) ci(b) pv(b) hv(b)];
    end

    subplot(3,nCol,2*nCol+1); hold on;
    bar(1:nBand, mu, 'FaceColor',[.6 .7 .85],'EdgeColor',[.25 .3 .45]);
    for b = 1:nBand
        plot([b b],[mu(b)-ci(b), mu(b)+ci(b)],'k-','LineWidth',1.2);
        if hv(b)
            text(b, mu(b)+ci(b), '*','HorizontalAlignment','center', ...
                 'FontSize',14,'FontWeight','bold');
        end
    end
    plot([0 nBand+1],[0 0],'k-');
    set(gca,'XTick',1:nBand,'XTickLabel',bandNames,'FontSize',8);
    try, xtickangle(45); catch, end
    xlim([0.4 nBand+0.6]); box on; grid on;
    ylabel('mean paired change');
    title(sprintf('%s - %s  (n=%d, * FDR q<%.2f)', cEv{2}, cEv{1}, nP, opt.fdrQ), ...
          'FontWeight','bold','FontSize',9);
end

% ---- Row 3b: diagnostics ----
subplot(3,nCol,2*nCol+2); axis off;
tx = {};
tx{end+1} = subName;
tx{end+1} = sprintf('mode: %s', getfielddef(R,'mode','?'));
if isfield(R,'band') && ~isempty(R.band)
    tx{end+1} = sprintf('band: %g-%g Hz', R.band(1), R.band(2));
end
tx{end+1} = sprintf('nodes: %d', nChan);
tx{end+1} = sprintf('VAR order: %d', R.order);
tx{end+1} = sprintf('fs (analysis): %g Hz', getfielddef(R,'fsAnalysis',NaN));
tx{end+1} = sprintf('trials per fit: %d', getfielddef(R,'blockTrials',1));
tx{end+1} = sprintf('N per fit: %d', getfielddef(R,'winSamples',NaN));
nf = getfielddef(R,'winSamples',NaN);
tx{end+1} = sprintf('null floor ~ p/N: %.4f', R.order/nf);
tx{end+1} = sprintf('min rcond: %.1e%s', R.rcondMin, ...
                    tern(R.rcondMin<1e-6,'  <-- ill-conditioned',''));
tx{end+1} = '';
tx{end+1} = 'fits per event:';
for e = 1:nEvent
    tx{end+1} = sprintf('   %-12s %d', events{e}, ...
        sum(logical(R.events.(events{e}).validEpoch(:)))); %#ok<AGROW>
end
tx{end+1} = '';
tx{end+1} = 'matrix(i,j) = from j (col) to i (row)';
text(0,1,tx,'VerticalAlignment','top','FontName','FixedWidth','FontSize',8, ...
     'Interpreter','none');

if ~isempty(opt.saveFig)
    try, print(fig, opt.saveFig, '-dpng','-r120');
    catch err, warning('Plot_Spectral_GC:save','%s',err.message); end
end
end

%% ---------------- helpers ----------------
function v = getfielddef(S,f,d)
if isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end
function m = nanmean_(X,dim)
N=sum(~isnan(X),dim); X(isnan(X))=0; m=sum(X,dim)./max(N,1); m(N==0)=NaN;
end
function s = nanstd_(X,dim)
m=nanmean_(X,dim); Dv=X-m; Dv(isnan(X))=0; N=sum(~isnan(X),dim);
s=sqrt(sum(Dv.^2,dim)./max(N-1,1)); s(N<2)=NaN;
end
function p = signflip_p(d,nPerm)
d=d(~isnan(d)); n=numel(d);
if n<3, p=NaN; return; end
obs=abs(mean(d)); cnt=0;
for it=1:nPerm
    s=sign(randn(n,1)); s(s==0)=1;
    if abs(mean(d.*s))>=obs, cnt=cnt+1; end
end
p=(1+cnt)/(nPerm+1);
end
function h = bh_fdr(p,q)
p=p(:); h=false(size(p)); ok=~isnan(p); pv=p(ok);
if isempty(pv), return; end
[ps,idx]=sort(pv); m=numel(ps); thr=(1:m)'/m*q;
b=find(ps<=thr,1,'last');
if ~isempty(b)
    keep=false(m,1); keep(1:b)=true; tmp=false(m,1); tmp(idx)=keep; h(ok)=tmp;
end
end
function c = seq_map(n)
t = linspace(0,1,n)';
c = [1-0.85*t, 1-0.75*t, 1-0.35*t];
end
function c = lines_(n)
base = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; ...
        0.49 0.18 0.56; 0.93 0.69 0.13; 0.30 0.75 0.93];
c = base(mod(0:n-1, size(base,1))+1, :);
end
function v = prctile_(x,pct)
x = sort(x(~isnan(x)));
if isempty(x), v = NaN; return; end
idx = max(1, min(numel(x), ceil(pct/100*numel(x))));
v = x(idx);
end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
