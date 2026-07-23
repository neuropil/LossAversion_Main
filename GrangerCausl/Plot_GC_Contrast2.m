function out = Plot_GC_Contrast2(source, evtA, evtB, varargin)
%% out = Plot_GC_Contrast(source, evtA, evtB, 'name',value, ...)
% Focused comparison of the GC network between TWO events (e.g. CHOICE vs
% OUTCOME) for one subject's <subject>_GC.mat.
%
% Four panels:
%   A  scatter of EVERY directed edge, evtA (x) vs evtB (y), with the y=x
%      line and the null-bias floor marked. No subsetting, so nothing is
%      hidden by selection.
%   B  dumbbell ("barbell") of the top-N edges by |paired change|, sorted,
%      with the null floor drawn as a vertical reference.
%   C  paired change +/- 95% CI for the SAME edges in the same row order.
%      This is the statistically trustworthy panel (see note below).
%   D  node-level in-strength / out-strength dumbbell (which regions become
%      more or less influential).
%
% WHY THE FLOOR LINE MATTERS
%   Pairwise GC is biased upward; for a truly unconnected pair the expected
%   value is approximately p/N (model order / window samples). At p=5,
%   N=1501 that is ~0.0033, which can be a large share of an observed value.
%   Absolute GC levels near the floor should not be read as connectivity.
%
% WHY PANEL C IS THE ONE TO TRUST
%   That bias is common to both events (same order, same window, enforced by
%   Run_GC_Pipeline), so it CANCELS in the within-trial paired difference.
%   Panel C is therefore essentially bias-free, whereas panels A/B still
%   carry the floor in their absolute positions.
%
% INPUT
%   source : path to <subject>_GC.mat, or the subjResult struct
%   evtA, evtB : event names, e.g. 'CHOICE','OUTCOME' (change = B - A)
%
% OPTIONS
%   'topN'       [15]    edges shown in panels B/C
%   'fdrQ'       [0.05]  FDR level for the paired test
%   'nPerm'      [5000]  sign-flip permutations
%   'winSamples' []      window length; [] = read from struct, else NaN floor
%   'nullFloor'  []      override the floor estimate (else p/winSamples)
%   'select'     ['change'] rank edges by 'change' (|B-A|) or 'level' (max)
%   'saveFig'    ''
%   'visible'    ['on']
%
% OUTPUT struct `out`
%   .edge      table-like struct arrays: label, a, b, diff, ci, p, h
%   .nodes     .label .inA .inB .outA .outB
%   .nullFloor .nPaired

ip = inputParser;
ip.addRequired('source');
ip.addRequired('evtA',@ischar);
ip.addRequired('evtB',@ischar);
ip.addParameter('topN',15,@isscalar);
ip.addParameter('fdrQ',0.05,@isscalar);
ip.addParameter('nPerm',5000,@isscalar);
ip.addParameter('winSamples',[]);
ip.addParameter('nullFloor',[]);
ip.addParameter('select','change',@ischar);
ip.addParameter('saveFig','',@ischar);
ip.addParameter('visible','on',@ischar);
ip.parse(source,evtA,evtB,varargin{:});
opt = ip.Results;

% ---- load ----
if ischar(source)
    S = load(source);
    if isfield(S,'subjResult'), R = S.subjResult; else, fn=fieldnames(S); R=S.(fn{1}); end
else
    R = source;
end
labels  = R.conn.nodeLabel(:);
dispLab = strrep(labels,'_','-');
nChan   = numel(labels);
offd    = ~eye(nChan);

% ---- null bias floor ----
if ~isempty(opt.nullFloor)
    floorVal = opt.nullFloor;
else
    if ~isempty(opt.winSamples), Nw = opt.winSamples;
    elseif isfield(R,'winSamples'), Nw = R.winSamples;
    else, Nw = NaN; end
    floorVal = R.order / Nw;     % E[GC | no coupling] ~ p/N
end

% ---- paired trials ----
A = R.events.(evtA); B = R.events.(evtB);
if isfield(A,'validEpoch'), vA=logical(A.validEpoch(:)); else, vA=true(size(A.GC,3),1); end
if isfield(B,'validEpoch'), vB=logical(B.validEpoch(:)); else, vB=true(size(B.GC,3),1); end
nT = min(numel(vA),numel(vB));
both = vA(1:nT) & vB(1:nT);
nPaired = sum(both);

GA = A.GC(:,:,both); GB = B.GC(:,:,both);
D  = GB - GA;                       % paired difference: bias cancels

mA = nanmean_(GA,3); mB = nanmean_(GB,3); mD = nanmean_(D,3);
sD = nanstd_(D,3);  seD = sD ./ sqrt(max(nPaired,1));
ciD = 1.96 * seD;

% ---- paired sign-flip test + FDR ----
pD = nan(nChan);
for i=1:nChan
    for j=1:nChan
        if i==j, continue; end
        pD(i,j) = signflip_p(squeeze(D(i,j,:)), opt.nPerm);
    end
end
hD = false(nChan);
hD(offd) = bh_fdr(pD(offd), opt.fdrQ);

% ---- OMNIBUS max-statistic permutation (FWER-controlled) ----
% Answers "did the network change AT ALL?" without counting edges, and gives
% each edge an FWER-corrected p from the same null. The max-statistic null
% automatically accounts for the 90 simultaneous tests.
obsAbs  = abs(mD); obsMax = max(obsAbs(offd));
nullMax = zeros(opt.nPerm,1);
for it = 1:opt.nPerm
    s = sign(randn(1,1,nPaired)); s(s==0)=1;
    Dp = nanmean_(bsxfun(@times, D, s), 3);
    tmp = abs(Dp); nullMax(it) = max(tmp(offd));
end
pOmni = (1 + sum(nullMax >= obsMax)) / (opt.nPerm + 1);
pFWER = nan(nChan);
for i=1:nChan
    for j=1:nChan
        if i==j, continue; end
        pFWER(i,j) = (1 + sum(nullMax >= obsAbs(i,j))) / (opt.nPerm + 1);
    end
end
out.omnibus.p = pOmni; out.omnibus.obsMax = obsMax;

% ---- NODE-LEVEL strength change (aggregates ~n-1 edges per node) ----
% in-strength  of node i = sum_j GC(i,j)  (row sum, incoming)
% out-strength of node j = sum_i GC(i,j)  (col sum, outgoing)
inD  = squeeze(nansum_(D,2));      % nChan x nTrial
outD = squeeze(nansum_(D,1));      % nChan x nTrial
if size(inD,1)~=nChan,  inD  = inD';  end
if size(outD,1)~=nChan, outD = outD'; end
nodeIn = struct('diff',[],'ci',[],'p',[]); nodeOut = nodeIn;
pIn = nan(nChan,1); pOut = nan(nChan,1);
dIn = nan(nChan,1); dOut = nan(nChan,1);
cIn = nan(nChan,1); cOut = nan(nChan,1);
for i=1:nChan
    di = inD(i,:)';
    dIn(i)  = nanmean_(di,1);  cIn(i)  = 1.96*nanstd_(di,1)/sqrt(nPaired);
    pIn(i)  = signflip_p(di, opt.nPerm);
    do_ = outD(i,:)';
    dOut(i) = nanmean_(do_,1); cOut(i) = 1.96*nanstd_(do_,1)/sqrt(nPaired);
    pOut(i) = signflip_p(do_, opt.nPerm);
end
hIn  = bh_fdr(pIn,  opt.fdrQ);
hOut = bh_fdr(pOut, opt.fdrQ);
out.node.label=labels; out.node.dIn=dIn; out.node.ciIn=cIn; out.node.pIn=pIn; out.node.hIn=hIn;
out.node.dOut=dOut; out.node.ciOut=cOut; out.node.pOut=pOut; out.node.hOut=hOut;

% ---- CONNECTION-TYPE change (intra-L / intra-R / inter) ----
tyName = {'intraL','intraR','inter'};
for t=1:3
    msk = R.conn.(tyName{t}) & offd;
    if ~any(msk(:)), out.type.(tyName{t}) = [NaN NaN NaN 0]; continue; end
    dt = zeros(nPaired,1);
    for k=1:nPaired
        sl = D(:,:,k); v = sl(msk); dt(k) = nanmean_(v(:),1);
    end
    out.type.(tyName{t}) = [nanmean_(dt,1), 1.96*nanstd_(dt,1)/sqrt(nPaired), ...
                            signflip_p(dt,opt.nPerm), sum(msk(:))];
end

% ---- flatten edges ----
ii=[]; jj=[];
for i=1:nChan, for j=1:nChan, if i~=j, ii(end+1)=i; jj(end+1)=j; end, end, end %#ok<AGROW>
lin = sub2ind([nChan nChan], ii, jj);
eLab = arrayfun(@(k) sprintf('%s -> %s', dispLab{jj(k)}, dispLab{ii(k)}), 1:numel(ii), 'UniformOutput',false);

eA = mA(lin)'; eB = mB(lin)'; eD = mD(lin)'; eCI = ciD(lin)';
eP = pD(lin)'; eH = hD(lin)';

% type per edge (for colour)
tIdx = zeros(numel(ii),1);
for k=1:numel(ii)
    if R.conn.intraL(ii(k),jj(k)), tIdx(k)=1;
    elseif R.conn.intraR(ii(k),jj(k)), tIdx(k)=2;
    else, tIdx(k)=3; end
end

out.edge.label = eLab(:); out.edge.a=eA; out.edge.b=eB; out.edge.diff=eD;
out.edge.ci=eCI; out.edge.p=eP; out.edge.h=eH; out.edge.type=tIdx;
out.edge.pFWER = pFWER(lin)';
out.edge.hFWER = out.edge.pFWER < 0.05;
out.nullFloor = floorVal; out.nPaired = nPaired;
out.events = {evtA, evtB};

% ---- node strengths (convention: matrix(i,j)=j->i) ----
inA  = nansum_(mA,2);  inB  = nansum_(mB,2);   % row sums  = incoming
outA = nansum_(mA,1)'; outB = nansum_(mB,1)';  % col sums  = outgoing
out.nodes.label=labels; out.nodes.inA=inA; out.nodes.inB=inB;
out.nodes.outA=outA; out.nodes.outB=outB;

% ---- rank edges ----
switch lower(opt.select)
    case 'level', [~,ord] = sort(max(eA,eB),'descend');
    otherwise,    [~,ord] = sort(abs(eD),'descend');
end
K = min(opt.topN, numel(ord)); sel = ord(1:K);
sel = flipud(sel(:));            % so the biggest ends up at the TOP of the axis

%% ================= figure =================
fig = figure('visible',opt.visible,'Color','w','Position',[60 60 1500 980]);
cIncr = [0.80 0.15 0.15]; cDecr = [0.15 0.30 0.80];
tCols = [0.20 0.55 0.85; 0.85 0.35 0.15; 0.95 0.72 0.20];

% ---- A: scatter, all edges ----
subplot(2,2,1); hold on;
lims = [0, max([eA;eB])*1.08];
plot(lims,lims,'k-','LineWidth',1);
if isfinite(floorVal)
    plot(lims,[floorVal floorVal],'--','Color',[.55 .55 .55]);
    plot([floorVal floorVal],lims,'--','Color',[.55 .55 .55]);
end
for t=1:3
    m = tIdx==t; if ~any(m), continue; end
    plot(eA(m), eB(m), 'o','MarkerSize',5,'MarkerEdgeColor',tCols(t,:), ...
        'MarkerFaceColor','none','LineWidth',1);
end
if any(eH)
    plot(eA(eH), eB(eH), 'o','MarkerSize',7,'MarkerEdgeColor','k', ...
        'MarkerFaceColor','k','LineWidth',1);
end
xlim(lims); ylim(lims); axis square; box on;
xlabel(sprintf('%s  mean GC',evtA)); ylabel(sprintf('%s  mean GC',evtB));
title(sprintf('all %d edges  (above line = stronger in %s)', numel(eA), evtB),'FontWeight','bold');
lg = {'y = x'};
if isfinite(floorVal), lg{end+1} = sprintf('null floor %.4f', floorVal); end
legend(lg,'Location','southeast'); legend boxoff;

% ---- B: dumbbell ----
subplot(2,2,3); hold on;
y = 1:K;
if isfinite(floorVal)
    plot([floorVal floorVal],[0 K+1],'--','Color',[.55 .55 .55],'LineWidth',1);
end
% thin profile lines linking each condition's points down the column
% plot(eA(sel), 1:K, '-', 'Color',[.45 .45 .45], 'LineWidth',0.75);
% plot(eB(sel), 1:K, '-', 'Color',[.15 .15 .15], 'LineWidth',0.75);
for k=1:K
    s = sel(k);
    if eD(s) >= 0, cc = cIncr; else, cc = cDecr; end
    plot([eA(s) eB(s)],[k k],'-','Color',cc,'LineWidth',2.5);
    plot(eA(s),k,'o','MarkerSize',7,'MarkerFaceColor','w','MarkerEdgeColor',[.35 .35 .35],'LineWidth',1.4);
    plot(eB(s),k,'o','MarkerSize',7,'MarkerFaceColor',cc,'MarkerEdgeColor',cc);
end
set(gca,'YTick',y,'YTickLabel',eLab(sel),'FontSize',8);
ylim([0 K+1]); box on; grid on;
xlabel('mean GC'); title(sprintf('top %d edges by |change|   (open=%s, filled=%s)',K,evtA,evtB),'FontWeight','bold');
% star FDR survivors
xr = xlim;
for k=1:K
    if eH(sel(k)), text(xr(2), k, ' *','FontSize',13,'FontWeight','bold','VerticalAlignment','middle'); end
end

% ---- C: paired change +/- CI (same rows) ----
subplot(2,2,4); hold on;
plot([0 0],[0 K+1],'k-','LineWidth',1);
for k=1:K
    s = sel(k);
    if eD(s) >= 0, cc = cIncr; else, cc = cDecr; end
    plot([eD(s)-eCI(s), eD(s)+eCI(s)],[k k],'-','Color',cc,'LineWidth',1.6);
    mk = 'o'; if eH(s), mk='s'; end
    plot(eD(s),k,mk,'MarkerSize',7,'MarkerFaceColor',cc,'MarkerEdgeColor',cc);
end
set(gca,'YTick',1:K,'YTickLabel',[],'FontSize',8);
ylim([0 K+1]); box on; grid on;
xlabel(sprintf('paired change  (%s - %s)',evtB,evtA));
title(sprintf('effect size +/- 95%% CI  (n=%d paired)  square = FDR q<%.2f',nPaired,opt.fdrQ),'FontWeight','bold');

% ---- D: node strength dumbbell ----
subplot(2,2,2); hold on;
yb = 1:nChan;
off = 0.18;
for i=1:nChan
    % outgoing (upper offset)
    cc = tern(outB(i)>=outA(i), cIncr, cDecr);
    plot([outA(i) outB(i)],[i+off i+off],'-','Color',cc,'LineWidth',2);
    plot(outA(i),i+off,'o','MarkerSize',6,'MarkerFaceColor','w','MarkerEdgeColor',[.35 .35 .35],'LineWidth',1.3);
    plot(outB(i),i+off,'o','MarkerSize',6,'MarkerFaceColor',cc,'MarkerEdgeColor',cc);
    % incoming (lower offset)
    cc2 = tern(inB(i)>=inA(i), cIncr, cDecr);
    plot([inA(i) inB(i)],[i-off i-off],'-','Color',cc2,'LineWidth',2);
    plot(inA(i),i-off,'s','MarkerSize',6,'MarkerFaceColor','w','MarkerEdgeColor',[.35 .35 .35],'LineWidth',1.3);
    plot(inB(i),i-off,'s','MarkerSize',6,'MarkerFaceColor',cc2,'MarkerEdgeColor',cc2);
end
set(gca,'YTick',1:nChan,'YTickLabel',dispLab,'FontSize',8);
ylim([0.4 nChan+0.6]); box on; grid on;
xlabel('summed GC'); title('node strength: o = OUT (row+), s = IN (row-)','FontWeight','bold');

if ~isempty(opt.saveFig)
    try, print(fig, opt.saveFig, '-dpng','-r120');
    catch err, warning('Plot_GC_Contrast:save','%s',err.message); end
end
end

%% ---------------- helpers ----------------
function m = nanmean_(X,dim)
N=sum(~isnan(X),dim); X(isnan(X))=0; m=sum(X,dim)./max(N,1); m(N==0)=NaN;
end
function s = nanstd_(X,dim)
m=nanmean_(X,dim); Dv=X-m; Dv(isnan(X))=0; N=sum(~isnan(X),dim);
s=sqrt(sum(Dv.^2,dim)./max(N-1,1)); s(N<2)=NaN;
end
function s = nansum_(X,dim)
X(isnan(X))=0; s=sum(X,dim);
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
function s = tern(c,a,b), if c, s=a; else, s=b; end, end