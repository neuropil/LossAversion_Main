function Plot_GC_TimeResolved(TR, varargin)
%% Plot_GC_TimeResolved(TR, 'name',value, ...)
% Views for the struct returned by Run_GC_TimeResolved.
%
% 'view' selects the layout:
%   'grid'  (default) nChan x nChan grid of GC(t), one panel per directed
%           edge, events overlaid. Rows = TARGET, cols = SOURCE, matching the
%           pipeline convention. The best overview of who drives whom, when.
%   'top'   GC(t) for the top-K edges only, larger panels with the null floor
%           and (if present) significance shading.
%   'tf'    time-frequency GC(t,f) for the top-K edges. Needs TR.spec, i.e.
%           Run_GC_TimeResolved(..., 'spectral', true).
%   'band'  band-averaged GC(t), one line per frequency band, for the top-K
%           edges. Also needs TR.spec. Use this to find WHICH rhythm carries
%           an effect seen in the broadband TR.gc trace.
%   'edge'  a single named edge, all events, full size.
%
% PARAMETERS
%   'view'    ['grid']
%   'topK'    [6]
%   'edge'    []      {'sourceLabel','targetLabel'} for view 'edge'
%   'events'  []      subset to plot; [] = all
%   'baseline' [false] plot TR.gcBase instead of TR.gc
%   'clim'    []      colour limits for 'tf'
%   'bands'   []      band definitions for view 'band'; [] = theta/alpha/
%                     beta/lowGamma/highGamma
%   'saveFig' ''
%   'visible' ['on']

ip = inputParser;
ip.addParameter('view','grid',@ischar);
ip.addParameter('topK',6,@isscalar);
ip.addParameter('edge',[],@(x)isempty(x)||iscell(x));
ip.addParameter('events',[],@(x)isempty(x)||iscell(x));
ip.addParameter('baseline',false,@(x)islogical(x)||isnumeric(x));
ip.addParameter('clim',[],@(x)isempty(x)||numel(x)==2);
ip.addParameter('bands',[],@(x)isempty(x)||isstruct(x));   % for view 'band'
ip.addParameter('saveFig','',@ischar);
ip.addParameter('visible','on',@ischar);
ip.parse(varargin{:});
o = ip.Results;

G = TR.gc;
if o.baseline
    if ~isfield(TR,'gcBase')
        error('Plot_GC_TimeResolved:nobase','TR has no .gcBase; pass ''baseline'' to Run_GC_TimeResolved.');
    end
    G = TR.gcBase;
end

evAll = TR.events;
if isempty(o.events), evIdx = 1:numel(evAll);
else, [~,evIdx] = ismember(o.events, evAll); evIdx = evIdx(evIdx>0); end
lab  = strrep(TR.nodeLabel,'_','-');
n    = numel(lab);
t    = TR.t;
cols = lines_(numel(evIdx));
floorV = TR.nullFloor;

switch lower(o.view)

case 'grid'
    fig = figure('visible',o.visible,'Color','w','Position',[40 40 240*n 200*n]);
    ylo = min(G(:)); yhi = max(G(:));
    if ~isfinite(ylo), ylo = 0; yhi = 1; end
    pad = 0.08*(yhi-ylo);
    for i = 1:n
        for j = 1:n
            subplot(n,n,(i-1)*n+j);
            if i == j
                axis off;
                text(0.5,0.5,lab{i},'HorizontalAlignment','center', ...
                     'FontWeight','bold','FontSize',9);
                continue;
            end
            hold on;
            if ~o.baseline && isfinite(floorV)
                plot([t(1) t(end)],[floorV floorV],'--','Color',[.6 .6 .6]);
            else
                plot([t(1) t(end)],[0 0],'-','Color',[.6 .6 .6]);
            end
            plot([0 0],[ylo-pad yhi+pad],'-','Color',[.8 .8 .8]);
            for k = 1:numel(evIdx)
                plot(t, squeeze(G(i,j,:,evIdx(k))), '-', ...
                     'Color', cols(k,:), 'LineWidth', 1.3);
            end
            xlim([t(1) t(end)]); ylim([ylo-pad yhi+pad]);
            box on; set(gca,'FontSize',7);
            if i == 1, title(sprintf('from %s', lab{j}),'FontSize',8); end
            if j == 1, ylabel(sprintf('to %s', lab{i}),'FontSize',8,'FontWeight','bold'); end
            if i == n, xlabel('time (s)','FontSize',7); end
            if i == 1 && j == 2
                legend(evAll(evIdx),'FontSize',6,'Location','best'); legend boxoff;
            end
        end
    end

case 'band'
    % Band-averaged GC(t): collapses TR.spec over frequency within each band
    % and plots one time course per band. Answers "which rhythm carries the
    % effect I see in the broadband trace".
    if ~isfield(TR,'spec') || isempty(TR.spec)
        error('Plot_GC_TimeResolved:nospec', ...
             ['This view needs TR.spec. Re-run with:\n' ...
              '   Run_GC_TimeResolved(..., ''spectral'', true)\n' ...
              'The broadband TR.gc alone cannot be split by frequency.']);
    end
    bands = o.bands;
    if isempty(bands)
        bands = struct('theta',[4 8],'alpha',[8 12],'beta',[13 30], ...
                       'lowGamma',[30 70],'highGamma',[70 150]);
    end
    bn = fieldnames(bands);
    bn = bn(cellfun(@(b) bands.(b)(1) < max(TR.f), bn));
    nB = numel(bn);

    m = mean(G(:,:,:,evIdx(1)),3); m(1:n+1:end) = NaN;
    [ii,jj] = deal([]);
    for i=1:n, for j=1:n, if i~=j, ii(end+1)=i; jj(end+1)=j; end, end, end %#ok<AGROW>
    v = arrayfun(@(k) m(ii(k),jj(k)), 1:numel(ii));
    [~,ord] = sort(v,'descend');
    K = min(o.topK, numel(ord)); sel = ord(1:K);
    nc = ceil(sqrt(K)); nr = ceil(K/nc);
    bCols = lines_(nB);

    ev = evIdx(1);
    fig = figure('visible',o.visible,'Color','w','Position',[50 50 380*nc 280*nr]);
    for k = 1:K
        subplot(nr,nc,k); hold on;
        for b = 1:nB
            rg = bands.(bn{b});
            fm = TR.f >= rg(1) & TR.f <= rg(2);
            if ~any(fm), continue; end
            y = squeeze(mean(TR.spec(ii(sel(k)), jj(sel(k)), fm, :, ev), 3));
            plot(t, y, '-', 'Color', bCols(b,:), 'LineWidth', 1.5);
        end
        yl = ylim; plot([0 0], yl, '-', 'Color', [.8 .8 .8]);
        xlim([t(1) t(end)]); box on; grid on; set(gca,'FontSize',8);
        xlabel('time (s)'); ylabel('spectral GC');
        title(sprintf('%s -> %s   (%s)', lab{jj(sel(k))}, lab{ii(sel(k))}, ...
              evAll{ev}), 'FontSize',9,'FontWeight','bold');
        if k == 1
            legend(bn,'FontSize',7,'Location','best'); legend boxoff;
        end
    end

case {'top','tf'}
    m = mean(G(:,:,:,evIdx(1)),3); m(1:n+1:end) = NaN;
    [ii,jj] = deal([]);
    for i=1:n, for j=1:n, if i~=j, ii(end+1)=i; jj(end+1)=j; end, end, end %#ok<AGROW>
    v = arrayfun(@(k) m(ii(k),jj(k)), 1:numel(ii));
    [~,ord] = sort(v,'descend');
    K = min(o.topK, numel(ord)); sel = ord(1:K);
    nc = ceil(sqrt(K)); nr = ceil(K/nc);

    if strcmpi(o.view,'tf')
        if ~isfield(TR,'spec') || isempty(TR.spec)
            error('Plot_GC_TimeResolved:nospec', ...
                  'TR has no .spec; re-run with ''spectral'',true.');
        end
        fig = figure('visible',o.visible,'Color','w','Position',[50 50 380*nc 300*nr]);
        allv = TR.spec(:); allv = allv(isfinite(allv));
        cl = o.clim; if isempty(cl), cl = [0 prctile_(allv,99)]; end
        for k = 1:K
            subplot(nr,nc,k);
            Z = squeeze(TR.spec(ii(sel(k)), jj(sel(k)), :, :, evIdx(1)));
            imagesc(t, TR.f, Z, cl); axis xy; colormap(gca, seq_map(256));
            hold on; plot([0 0],[TR.f(1) TR.f(end)],'w-','LineWidth',1);
            colorbar; set(gca,'FontSize',8);
            xlabel('time (s)'); ylabel('Hz');
            title(sprintf('%s -> %s   (%s)', lab{jj(sel(k))}, lab{ii(sel(k))}, ...
                  evAll{evIdx(1)}),'FontSize',9,'FontWeight','bold');
        end
    else
        fig = figure('visible',o.visible,'Color','w','Position',[50 50 360*nc 260*nr]);
        for k = 1:K
            subplot(nr,nc,k); hold on;
            for q = 1:numel(evIdx)
                y = squeeze(G(ii(sel(k)), jj(sel(k)), :, evIdx(q)));
                plot(t, y, '-', 'Color', cols(q,:), 'LineWidth', 1.6);
                if isfield(TR,'sig') && any(TR.sig(:))
                    sg = squeeze(TR.sig(ii(sel(k)), jj(sel(k)), :, evIdx(q)));
                    yy = y; yy(~sg) = NaN;
                    plot(t, yy, '-', 'Color', cols(q,:), 'LineWidth', 4);
                end
            end
            if ~o.baseline && isfinite(floorV)
                plot([t(1) t(end)],[floorV floorV],'--','Color',[.55 .55 .55]);
            else
                plot([t(1) t(end)],[0 0],'-','Color',[.55 .55 .55]);
            end
            yl = ylim; plot([0 0],yl,'-','Color',[.8 .8 .8]);
            xlim([t(1) t(end)]); box on; grid on; set(gca,'FontSize',8);
            xlabel('time (s)'); ylabel('GC');
            title(sprintf('%s -> %s', lab{jj(sel(k))}, lab{ii(sel(k))}), ...
                  'FontSize',9,'FontWeight','bold');
            if k == 1
                legend(evAll(evIdx),'FontSize',7,'Location','best'); legend boxoff;
            end
        end
    end

case 'edge'
    if isempty(o.edge) || numel(o.edge) ~= 2
        error('Plot_GC_TimeResolved:edge','Pass ''edge'',{''source'',''target''}.');
    end
    j = find(strcmp(TR.nodeLabel, o.edge{1}));
    i = find(strcmp(TR.nodeLabel, o.edge{2}));
    if isempty(i) || isempty(j)
        error('Plot_GC_TimeResolved:edge','Labels not found. Available: %s', ...
              strjoin(TR.nodeLabel', ', '));
    end
    fig = figure('visible',o.visible,'Color','w','Position',[80 80 800 420]); hold on;
    for q = 1:numel(evIdx)
        plot(t, squeeze(G(i,j,:,evIdx(q))), '-','Color',cols(q,:),'LineWidth',2);
    end
    if ~o.baseline && isfinite(floorV)
        plot([t(1) t(end)],[floorV floorV],'--','Color',[.55 .55 .55]);
    end
    yl = ylim; plot([0 0],yl,'k-');
    xlim([t(1) t(end)]); box on; grid on;
    xlabel('time (s)'); ylabel('GC');
    title(sprintf('%s: %s -> %s', TR.subject, lab{j}, lab{i}),'FontWeight','bold');
    legend(evAll(evIdx),'Location','best'); legend boxoff;

otherwise
    error('Plot_GC_TimeResolved:view','Unknown view ''%s''.', o.view);
end

if ~isempty(o.saveFig)
    try, print(fig, o.saveFig, '-dpng','-r120');
    catch err, warning('Plot_GC_TimeResolved:save','%s',err.message); end
end
end

%% ---- helpers ----
function c = lines_(n)
base = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; ...
        0.49 0.18 0.56; 0.93 0.69 0.13; 0.30 0.75 0.93];
c = base(mod(0:n-1, size(base,1))+1, :);
end
function c = seq_map(n)
t = linspace(0,1,n)';
c = [1-0.85*t, 1-0.75*t, 1-0.35*t];
end
function v = prctile_(x,pct)
x = sort(x(~isnan(x)));
if isempty(x), v = NaN; return; end
v = x(max(1,min(numel(x),ceil(pct/100*numel(x)))));
end
