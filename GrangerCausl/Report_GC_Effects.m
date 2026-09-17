function T = Report_GC_Effects(BS, varargin)
%% T = Report_GC_Effects(BS, 'name',value, ...)
% Turns the bootstrap result into a table of effects: one row per
% significant cluster, with onset, offset, duration, peak, latency to peak,
% and mean change over the cluster.
%
% Reports the CONTRAST by default (what BS.contrast holds) and optionally
% each event's change from baseline.
%
% All quantities are measured WITHIN the significant cluster, not over the
% whole epoch. A peak located outside the window where the effect is
% detectable is not a property of the effect, and averaging across
% non-significant stretches dilutes it toward zero.
%
% PARAMETERS
%   'what'     ['contrast'] 'contrast' | 'change' | 'both'
%   'cluster'  [true]    use cluster-corrected significance (recommended).
%                        false = pointwise, which at 20 edges x 51 windows
%                        will include roughly 50 false positives at 0.05.
%   'minWin'   [1]       ignore clusters shorter than this many windows
%   'metric'   ['rise50'] latency measure reported alongside the peak:
%                        'rise50' = time to 50% of peak, 'latency' = argmax.
%                        rise50 is the more stable of the two on broad peaks.
%   'sortBy'   ['peak']  'peak' | 'onset' | 'duration' | 'edge'
%   'print'    [true]
%   'csv'      ''        if non-empty, write the table to this path
%
% OUTPUT
%   T : struct array, one element per cluster, with fields
%       .source .target .edge .kind .event
%       .onset .offset .duration      seconds (cluster extent)
%       .peak .tPeak .tRise           extremum inside the cluster and its
%                                     latency (argmax and 50%-rise)
%       .meanChange                   mean over the cluster
%       .peakLo .peakHi               bootstrap CI at the peak window
%       .nWin
%
% NOTE ON WHAT THE NUMBERS MEAN
%   BS is a within-subject, trial-level bootstrap. Every interval here
%   describes THIS subject's data. They are not group statistics and do not
%   generalise across patients without a second-level model.

ip = inputParser;
ip.addParameter('what','contrast',@ischar);
ip.addParameter('cluster',true,@(x)islogical(x)||isnumeric(x));
ip.addParameter('minWin',1,@isscalar);
ip.addParameter('metric','rise50',@ischar);
ip.addParameter('sortBy','peak',@ischar);
ip.addParameter('print',true,@(x)islogical(x)||isnumeric(x));
ip.addParameter('csv','',@ischar);
ip.parse(varargin{:});
o = ip.Results;

t   = BS.t(:)';
lab = BS.nodeLabel(:);
n   = numel(lab);

blocks = {};
if any(strcmpi(o.what,{'contrast','both'}))
    if ~isfield(BS,'contrast')
        error('Report_GC_Effects:nocontrast','BS has no contrast.');
    end
    if o.cluster && isfield(BS,'sigCluster'), sg = BS.sigCluster;
    else,                                     sg = BS.contrast.sig; end
    blocks{end+1} = struct('kind','contrast', ...
        'event', sprintf('%s - %s', BS.contrastEvents{2}, BS.contrastEvents{1}), ...
        'M', BS.contrast.mean, 'LO', BS.contrast.lo, 'HI', BS.contrast.hi, 'SG', sg);
end
if any(strcmpi(o.what,{'change','both'}))
    for e = 1:numel(BS.events)
        ev = BS.events{e};
        blocks{end+1} = struct('kind','change','event',ev, ...
            'M', BS.change.(ev).mean, 'LO', BS.change.(ev).lo, ...
            'HI', BS.change.(ev).hi, 'SG', BS.change.(ev).sig); %#ok<AGROW>
    end
end

T = struct('source',{},'target',{},'edge',{},'kind',{},'event',{}, ...
           'onset',{},'offset',{},'duration',{},'peak',{},'tPeak',{}, ...
           'tRise',{},'meanChange',{},'peakLo',{},'peakHi',{},'nWin',{});

for b = 1:numel(blocks)
    B = blocks{b};
    for i = 1:n
        for j = 1:n
            if i == j, continue; end
            s = squeeze(B.SG(i,j,:))';
            y = squeeze(B.M(i,j,:))';
            [len, s0, e0] = runs_of(s);
            for r = 1:numel(len)
                if len(r) < o.minWin, continue; end
                k1 = s0(r); k2 = e0(r);
                seg = y(k1:k2); tseg = t(k1:k2);
                [~, kx] = max(abs(seg));
                pk = seg(kx); tp = tseg(kx);
                tr = rise50(seg, tseg, pk);

                rec.source = lab{j};        % matrix is (target, source)
                rec.target = lab{i};
                rec.edge   = sprintf('%s -> %s', lab{j}, lab{i});
                rec.kind   = B.kind;
                rec.event  = B.event;
                rec.onset  = t(k1);
                rec.offset = t(k2);
                rec.duration = t(k2) - t(k1);
                rec.peak   = pk;
                rec.tPeak  = tp;
                rec.tRise  = tr;
                rec.meanChange = mean(seg);
                rec.peakLo = B.LO(i,j,k1+kx-1);
                rec.peakHi = B.HI(i,j,k1+kx-1);
                rec.nWin   = len(r);
                T(end+1) = rec; %#ok<AGROW>
            end
        end
    end
end

if isempty(T)
    if o.print
        fprintf('\n  No significant clusters');
        if o.cluster, fprintf(' (cluster-corrected)'); end
        fprintf('.\n');
        if o.cluster
            fprintf(['  Try ''cluster'',false to see pointwise results, but treat\n' ...
                     '  those as exploratory -- they are uncorrected.\n']);
        end
        fprintf('\n');
    end
    return;
end

% ---- sort ----
switch lower(o.sortBy)
    case 'onset',    [~,ord] = sort([T.onset]);
    case 'duration', [~,ord] = sort([T.duration],'descend');
    case 'edge',     [~,ord] = sort({T.edge});
    otherwise,       [~,ord] = sort(abs([T.peak]),'descend');
end
T = T(ord);

% ---- print ----
if o.print
    latName = tern_(strcmpi(o.metric,'latency'),'t_peak','t_rise50');
    fprintf('\n===============================================================================\n');
    fprintf(' GC effects: %s   (%s significance, %d bootstrap resamples)\n', ...
            BS.subject, tern_(o.cluster,'cluster-corrected','POINTWISE/uncorrected'), ...
            BS.nBoot);
    fprintf('===============================================================================\n');
    fprintf('%-22s %-18s %7s %7s %7s %10s %8s %10s\n', ...
            'edge','event','onset','offset','dur','peak',latName,'mean');
    for k = 1:numel(T)
        lat = tern_(strcmpi(o.metric,'latency'), T(k).tPeak, T(k).tRise);
        fprintf('%-22s %-18s %7.2f %7.2f %7.2f %10.5f %8.2f %10.5f\n', ...
                T(k).edge, T(k).event, T(k).onset, T(k).offset, ...
                T(k).duration, T(k).peak, lat, T(k).meanChange);
    end
    fprintf('-------------------------------------------------------------------------------\n');
    fprintf(' times in seconds relative to event onset; peak/mean measured INSIDE the cluster\n');
    if BS.nBoot < 500
        fprintf(' NOTE: %d resamples -- CI edges are coarse; use >=500 before reporting\n', ...
                BS.nBoot);
    end
    fprintf('===============================================================================\n\n');
end

% ---- csv ----
if ~isempty(o.csv)
    fid = fopen(o.csv,'w');
    if fid < 0
        warning('Report_GC_Effects:csv','Could not open %s', o.csv);
    else
        fprintf(fid,['subject,source,target,kind,event,onset_s,offset_s,' ...
                     'duration_s,peak,t_peak_s,t_rise50_s,mean_change,' ...
                     'peak_lo,peak_hi,n_windows\n']);
        for k = 1:numel(T)
            fprintf(fid,'%s,%s,%s,%s,%s,%.4f,%.4f,%.4f,%.6f,%.4f,%.4f,%.6f,%.6f,%.6f,%d\n', ...
                BS.subject, T(k).source, T(k).target, T(k).kind, T(k).event, ...
                T(k).onset, T(k).offset, T(k).duration, T(k).peak, ...
                T(k).tPeak, T(k).tRise, T(k).meanChange, ...
                T(k).peakLo, T(k).peakHi, T(k).nWin);
        end
        fclose(fid);
        if o.print, fprintf('  wrote %s\n\n', o.csv); end
    end
end
end

%% ---- helpers ----
function tr = rise50(seg, tseg, pk)
tr = NaN;
if ~isfinite(pk) || pk == 0, return; end
[~,k] = max(abs(seg));
half = 0.5*abs(pk); s = sign(pk);
tr = tseg(k);
for q = k:-1:2
    if s*seg(q-1) < half
        y1 = s*seg(q-1); y2 = s*seg(q);
        if y2 > y1
            f = (half - y1)/(y2 - y1);
            tr = tseg(q-1) + f*(tseg(q) - tseg(q-1));
        else
            tr = tseg(q);
        end
        return;
    end
end
tr = tseg(1);
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

function s = tern_(c,a,b), if c, s=a; else, s=b; end, end
