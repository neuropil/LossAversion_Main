function Plot_GC_Bootstrap(BS, varargin)
%% Plot_GC_Bootstrap(BS, 'name',value, ...)
% Plots the output of Bootstrap_GC_TimeResolved as time courses with
% bootstrap confidence bands, with significant stretches shaded.
%
% 'view':
%   'contrast' (default) CHOICE-minus-OUTCOME (or whichever pair was run),
%              with its CI band. Shading marks cluster-corrected windows.
%   'change'   per-event change from baseline, events overlaid, CI bands.
%   'both'     two rows: change on top, contrast below.
%
% PARAMETERS
%   'view'      ['contrast']
%   'topK'      [6]      edges shown, ranked by peak |effect|
%   'cluster'   [true]   shade cluster-corrected windows; false = pointwise
%   'shareY'    [true]
%   'ylim'      []
%   'saveFig'   ''
%   'visible'   ['on']
%
% Reading the figure: a band that excludes zero over a contiguous stretch is
% the result. A point estimate far from zero whose band straddles it is not.

ip = inputParser;
ip.addParameter('view','contrast',@ischar);
ip.addParameter('topK',6,@isscalar);
ip.addParameter('cluster',true,@(x)islogical(x)||isnumeric(x));
ip.addParameter('shareY',true,@(x)islogical(x)||isnumeric(x));
ip.addParameter('ylim',[],@(x)isempty(x)||numel(x)==2);
ip.addParameter('saveFig','',@ischar);
ip.addParameter('visible','on',@ischar);
ip.parse(varargin{:});
o = ip.Results;

t   = BS.t;
lab = strrep(BS.nodeLabel,'_','-');
n   = numel(lab);
evs = BS.events;

hasCon = isfield(BS,'contrast');
if strcmpi(o.view,'contrast') && ~hasCon
    error('Plot_GC_Bootstrap:nocontrast', ...
        'BS has no contrast. Re-run with ''contrast'',{''A'',''B''}.');
end

% ---- rank edges by peak |effect| of whatever is being shown ----
if hasCon && ~strcmpi(o.view,'change')
    R = BS.contrast.mean;
else
    R = BS.change.(evs{1}).mean;
end
[ii,jj] = deal([]);
for i=1:n, for j=1:n, if i~=j, ii(end+1)=i; jj(end+1)=j; end, end, end %#ok<AGROW>
pk = arrayfun(@(k) max(abs(squeeze(R(ii(k),jj(k),:)))), 1:numel(ii));
[~,ord] = sort(pk,'descend');
K = min(o.topK, numel(ord)); sel = ord(1:K);
nc = ceil(sqrt(K)); nr = ceil(K/nc);

cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; ...
        0.49 0.18 0.56; 0.93 0.69 0.13];

switch lower(o.view)
    case 'contrast', rows = {'contrast'};
    case 'change',   rows = {'change'};
    case 'both',     rows = {'change','contrast'};
    otherwise, error('Plot_GC_Bootstrap:view','Unknown view ''%s''.', o.view);
end
nR = numel(rows);

% ---- shared y ----
YL = o.ylim;
if isempty(YL) && o.shareY
    v = [];
    for r = 1:nR
        for k = 1:K
            [~,L,H] = get_series(BS, rows{r}, evs, ii(sel(k)), jj(sel(k)));
            v = [v; L(:); H(:)]; %#ok<AGROW>
        end
    end
    v = v(isfinite(v));
    if ~isempty(v)
        pad = 0.08*(max(v)-min(v));
        YL = [min(v)-pad, max(v)+pad];
    end
end

fig = figure('visible',o.visible,'Color','w', ...
             'Position',[40 40 400*nc 300*nr*nR]);

for r = 1:nR
    for k = 1:K
        subplot(nR*nr, nc, (r-1)*nr*nc + k); hold on;
        i = ii(sel(k)); j = jj(sel(k));

        if strcmp(rows{r},'contrast')
            m = squeeze(BS.contrast.mean(i,j,:))';
            lo = squeeze(BS.contrast.lo(i,j,:))';
            hi = squeeze(BS.contrast.hi(i,j,:))';
            if o.cluster && isfield(BS,'sigCluster')
                sg = squeeze(BS.sigCluster(i,j,:))';
            else
                sg = squeeze(BS.contrast.sig(i,j,:))';
            end
            band(t, lo, hi, cols(1,:));
            % shade significant stretches
            yl = YL; if isempty(yl), yl = [min(lo) max(hi)]; end
            shade_runs(t, sg, yl, [0.95 0.90 0.55]);
            plot(t, m, '-', 'Color', cols(1,:), 'LineWidth', 2);
            ttl = sprintf('%s - %s', BS.contrastEvents{2}, BS.contrastEvents{1});
        else
            hE = zeros(1,numel(evs));
            for e = 1:numel(evs)
                m  = squeeze(BS.change.(evs{e}).mean(i,j,:))';
                lo = squeeze(BS.change.(evs{e}).lo(i,j,:))';
                hi = squeeze(BS.change.(evs{e}).hi(i,j,:))';
                band(t, lo, hi, cols(min(e,size(cols,1)),:));
                hE(e) = plot(t, m, '-', 'Color', cols(min(e,size(cols,1)),:), ...
                             'LineWidth',1.8);
            end
            % legend bound to the trace handles -- the CI patches are drawn
            % first, so an unbound legend() would label those instead
            if k == 1
                legend(hE, evs, 'FontSize',7, 'Location','best'); legend boxoff;
            end
            ttl = 'change vs baseline';
        end

        plot([t(1) t(end)],[0 0],'-','Color',[.45 .45 .45]);
        if ~isempty(YL), ylim(YL); end
        yl = ylim; plot([0 0],yl,'-','Color',[.8 .8 .8]);
        xlim([t(1) t(end)]); box on; grid on; set(gca,'FontSize',8);
        xlabel('time (s)'); ylabel('\DeltaGC');
        title(sprintf('%s -> %s   (%s)', lab{j}, lab{i}, ttl), ...
              'FontSize',9,'FontWeight','bold');
    end
end

if ~isempty(o.saveFig)
    try, print(fig, o.saveFig, '-dpng','-r120');
    catch err, warning('Plot_GC_Bootstrap:save','%s',err.message); end
end
end

%% ---- helpers ----
function [M,L,H] = get_series(BS, which_, evs, i, j)
if strcmp(which_,'contrast')
    M = squeeze(BS.contrast.mean(i,j,:));
    L = squeeze(BS.contrast.lo(i,j,:));
    H = squeeze(BS.contrast.hi(i,j,:));
else
    M = []; L = []; H = [];
    for e = 1:numel(evs)
        M = [M; squeeze(BS.change.(evs{e}).mean(i,j,:))]; %#ok<AGROW>
        L = [L; squeeze(BS.change.(evs{e}).lo(i,j,:))];   %#ok<AGROW>
        H = [H; squeeze(BS.change.(evs{e}).hi(i,j,:))];   %#ok<AGROW>
    end
end
end

function band(t, lo, hi, c)
ok = isfinite(lo) & isfinite(hi);
if ~any(ok), return; end
t = t(ok); lo = lo(ok); hi = hi(ok);
px = [t, fliplr(t)]; py = [lo, fliplr(hi)];
h = fill(px, py, c, 'EdgeColor','none');
try, set(h,'FaceAlpha',0.20); catch, set(h,'FaceColor', 1-0.25*(1-c)); end
end

function shade_runs(t, sg, yl, c)
k = 1; nW = numel(sg);
while k <= nW
    if sg(k)
        j = k; while j < nW && sg(j+1), j = j+1; end
        x1 = t(max(k-1,1)); x2 = t(min(j+1,nW));
        h = fill([x1 x2 x2 x1],[yl(1) yl(1) yl(2) yl(2)], c, 'EdgeColor','none');
        try, set(h,'FaceAlpha',0.30); catch, end
        uistack(h,'bottom');
        k = j+1;
    else
        k = k+1;
    end
end
end
