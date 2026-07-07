function results = Run_GC_Pipeline(dataDir, outDir, varargin)
%% results = Run_GC_Pipeline(dataDir, outDir, 'name',value, ...)
% Drives Granger-causality connectivity across brain regions, per subject,
% from a folder of region-split files named {SUBJECT}_{HEMI}_{REGION}_TrialDATA.mat
% (e.g. CLASE001_L_AC_TrialDATA.mat). Regions available PER SUBJECT are
% discovered from the filenames, so subjects may have different region sets
% and one or two hemispheres. GC is computed within subject across all of
% that subject's regions, and every channel pair is tagged intra-L / intra-R
% / inter-hemisphere.
%
% Depends on: load_region_epochs.m (EDIT to your file internals),
%             Compute_GC_V2.m, Select_Model_Order.m, Autoregressive_Process_V1.m
%
% OPTIONS (defaults in brackets)
%   'events'      [{'CHOICE','RESPONSEON','OUTCOME'}] event-field names to
%                 analyse (must match the struct fields in TrialTablesZS;
%                 note RESPONSEON vs RESPONSEBUTTON is a choice you must make)
%   'fs'          [500]
%   'channelUnit' ['region'] 'region' = one node per region (contacts collapsed);
%                            'contact' = every contact is its own node
%   'collapse'    ['mean']  region-unit collapse: 'mean' | 'first' | 'pca1'
%   'order'       []        fixed VAR order; [] = auto-select ONCE per subject
%   'maxOrder'    [30]      max order for auto-selection
%   'criterion'   ['BIC']   'AIC' | 'BIC'
%   'winStart'    [1]       first sample of the common analysis window
%   'winSamples'  []        window length; [] = shortest epoch in the subject
%                           (matches epoch length across events -> unbiased contrast)
%   'minEpochSamp'[50]      skip epochs shorter than this
%   'orderEpochs' [10]      # epochs concatenated for order selection
%   'nPerm'       [100]     circular-shift surrogates (0 = skip significance)
%   'alpha'       [0.05]
%   'zscore'      [true]    z-score each channel within epoch. ON by default:
%                           Volts are raw voltage and the VAR has no intercept,
%                           so demeaning is required; the unit-variance step is
%                           GC-neutral (scale-invariant) but improves the rcond
%                           of the fit when regions differ in amplitude.
%   'subjects'    {}        restrict to these subject IDs ({} = all found)
%   'save'        [true]    write <outDir>/<subject>_GC.mat per subject
%
% OUTPUT struct `results`, with results.(subject) containing:
%   .channelInfo  struct array (node -> hemi, region, contact, label)
%   .order        VAR order used
%   .events.(evt).GC   nChan x nChan x nEpoch, entry (i,j)=j->i (src->tgt)
%   .events.(evt).sig  logical, same size
%   .events.(evt).pval same size
%   .conn         masks (.intraL/.intraR/.inter, nChan x nChan) + .typeLabel + .nodeLabel
%   .rcondMin     smallest rcond seen (watch for tiny values)

ip = inputParser;
ip.addRequired('dataDir',@ischar);
ip.addRequired('outDir',@ischar);
ip.addParameter('events',{'CHOICE','RESPONSEON','OUTCOME'},@iscell);
ip.addParameter('fs',500,@isscalar);
ip.addParameter('channelUnit','region',@ischar);
ip.addParameter('collapse','mean',@ischar);
ip.addParameter('order',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('maxOrder',30,@isscalar);
ip.addParameter('criterion','BIC',@ischar);
ip.addParameter('winStart',1,@isscalar);
ip.addParameter('winSamples',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('minEpochSamp',50,@isscalar);
ip.addParameter('orderEpochs',10,@isscalar);
ip.addParameter('nPerm',100,@isscalar);
ip.addParameter('alpha',0.05,@isscalar);
ip.addParameter('zscore',true,@islogical);   % Volts are raw -> z-score per epoch/channel
ip.addParameter('subjects',{},@iscell);
ip.addParameter('save',true,@islogical);
ip.parse(dataDir,outDir,varargin{:});
opt = ip.Results;
regionUnit = strcmpi(opt.channelUnit,'region');

if ~isfolder(outDir), mkdir(outDir); end

% ---------- discover & parse files ----------
D = dir(fullfile(dataDir,'*_TrialDATA.mat'));
if isempty(D), error('No *_TrialDATA.mat files in %s', dataDir); end
nF = numel(D);
[subj,hemi,region,fpath] = deal(cell(nF,1));
for i = 1:nF
    [s,h,r] = parse_region_filename(D(i).name);
    subj{i}=s; hemi{i}=h; region{i}=r;
    fpath{i}=fullfile(D(i).folder,D(i).name);
end

subjects = unique(subj);
if ~isempty(opt.subjects), subjects = intersect(subjects, opt.subjects); end

results = struct();
nEvent  = numel(opt.events);

% ---------- per subject ----------
for si = 1:numel(subjects)
    sub = subjects{si};
    idx = find(strcmp(subj,sub));
    fprintf('\n=== %s : %d region files ===\n', sub, numel(idx));

    % load every region for this subject
    nReg = numel(idx);
    R = struct('hemi',{},'region',{},'nContact',{},'byEvent',{});
    for r = 1:nReg
        f = idx(r);
        [labels, epochs] = load_region_epochs(fpath{f});
        R(r).hemi   = hemi{f};
        R(r).region = region{f};
        % channel count from VALID (full-length, finite) epochs only -- some
        % trials store a degenerate placeholder (e.g. 1x1 NaN) when an event
        % is absent; those must not set the channel count.
        chCounts = [];
        for kk = 1:numel(epochs)
            x = epochs{kk};
            if ~isempty(x) && size(x,2) >= opt.minEpochSamp && all(isfinite(x(:)))
                chCounts(end+1) = size(x,1); %#ok<AGROW>
            end
        end
        if isempty(chCounts)
            warning('%s %s_%s: no valid epochs found.', sub, hemi{f}, region{f});
            R(r).nContact = 0;
        else
            R(r).nContact = mode(chCounts);
        end
        R(r).byEvent  = cell(nEvent,1);
        for e = 1:nEvent
            sel = strcmp(labels, opt.events{e});
            R(r).byEvent{e} = epochs(sel);
        end
        fprintf('   %s_%s : %d contacts | epochs/event = [%s]\n', ...
            R(r).hemi, R(r).region, R(r).nContact, ...
            strtrim(sprintf('%d ', cellfun(@numel, R(r).byEvent))));
    end

    % epochs per event (aligned across regions -> use the common minimum)
    nEpoch = zeros(nEvent,1);
    for e = 1:nEvent
        counts = arrayfun(@(x) numel(x.byEvent{e}), R);
        nEpoch(e) = min(counts);
        if numel(unique(counts))>1
            warning('%s event %s: region epoch counts differ (%s); using %d.', ...
                sub, opt.events{e}, mat2str(counts), nEpoch(e));
        end
    end

    % common window length across ALL epochs/events (matches length -> unbiased)
    if isempty(opt.winSamples)
        allLen = [];
        for e = 1:nEvent
            for k = 1:nEpoch(e)
                for r = 1:nReg
                    x = R(r).byEvent{e}{k};
                    if ~isempty(x) && size(x,1)==R(r).nContact && ...
                       size(x,2) >= opt.minEpochSamp && all(isfinite(x(:)))
                        allLen(end+1) = size(x,2); %#ok<AGROW>
                    end
                end
            end
        end
        if isempty(allLen)
            warning('%s: no valid epochs across events; skipping subject.', sub);
            continue;
        end
        commonLen = min(allLen);
    else
        commonLen = opt.winSamples;
    end
    win = opt.winStart : (opt.winStart + commonLen - 1);
    fprintf('   common window: %d samples (%.0f ms)\n', commonLen, 1000*commonLen/opt.fs);

    % ---------- channel bookkeeping (fixed across epochs) ----------
    channelInfo = struct('hemi',{},'region',{},'contact',{},'label',{});
    node = 0;
    for r = 1:nReg
        if regionUnit
            node = node+1;
            channelInfo(node) = mk_node(R(r).hemi,R(r).region,0);
        else
            for c = 1:R(r).nContact
                node = node+1;
                channelInfo(node) = mk_node(R(r).hemi,R(r).region,c);
            end
        end
    end
    nChan = node;

    % ---------- model order: select ONCE, then fix ----------
    if isempty(opt.order)
        stack = []; got = 0; eSel = 1;              % scan first event for VALID epochs
        for k = 1:nEpoch(eSel)
            if got >= opt.orderEpochs, break; end
            M = assemble_epoch(R,eSel,k,win,regionUnit,opt.collapse,nChan);
            if isempty(M), continue; end
            if opt.zscore, M = normalize(M,2,'zscore'); end
            stack = [stack, M]; got = got+1; %#ok<AGROW>
        end
        if isempty(stack)
            warning('%s: no valid epochs for order selection; skipping subject.', sub);
            continue;
        end
        p = Select_Model_Order(stack, opt.maxOrder, opt.criterion);
    else
        p = opt.order;
    end
    fprintf('   VAR order: %d\n', p);

    % ---------- GC per event x epoch ----------
    rcondMin = Inf;
    for e = 1:nEvent
        evt = opt.events{e};
        GCcube = nan(nChan,nChan,nEpoch(e));
        SIGcube= false(nChan,nChan,nEpoch(e));
        Pcube  = nan(nChan,nChan,nEpoch(e));
        valid = false(nEpoch(e),1);
        for k = 1:nEpoch(e)
            M = assemble_epoch(R,e,k,win,regionUnit,opt.collapse,nChan);
            if isempty(M), continue; end
            GC = Compute_GC_V2(M, 'order',p, 'nPerm',opt.nPerm, ...
                               'alpha',opt.alpha, 'zscore',opt.zscore);
            GCcube(:,:,k) = GC.matrix;
            SIGcube(:,:,k)= GC.sig;
            Pcube(:,:,k)  = GC.pval;
            rcondMin = min(rcondMin, GC.rcondMin);
            valid(k) = true;
        end
        results.(sub).events.(evt).GC       = GCcube;
        results.(sub).events.(evt).sig      = SIGcube;
        results.(sub).events.(evt).pval     = Pcube;
        results.(sub).events.(evt).validEpoch = valid;   % which trials were usable
        fprintf('   %-14s %d/%d epochs used (%d skipped: placeholder/short/NaN), %dx%d\n', ...
            evt, sum(valid), nEpoch(e), nEpoch(e)-sum(valid), nChan, nChan);
    end

    % ---------- connection classification ----------
    results.(sub).channelInfo = channelInfo;
    results.(sub).conn        = classify_connections(channelInfo);
    results.(sub).order       = p;
    results.(sub).rcondMin    = rcondMin;
    results.(sub).events_list = opt.events;

    if opt.save
        outFile = fullfile(outDir, sprintf('%s_GC.mat', sub));
        subjResult = results.(sub); %#ok<NASGU>
        save(outFile, 'subjResult', '-v7.3');
        fprintf('   saved -> %s\n', outFile);
    end
end
end

% ======================= local functions =======================

function [s,h,r] = parse_region_filename(name)
% CLASE001_L_AC_TrialDATA.mat -> s=CLASE001 h=L r=AC
base = regexprep(name, '_TrialDATA\.mat$', '');
tok  = strsplit(base, '_');
if numel(tok) < 3
    error('parse_region_filename:bad','Cannot parse "%s".', name);
end
s = tok{1};
h = tok{2};
r = strjoin(tok(3:end), '_');   % region may (rarely) contain '_'
if ~any(strcmp(h,{'L','R'}))
    warning('Unexpected hemisphere "%s" in %s', h, name);
end
end

function nd = mk_node(hemi,region,contact)
nd.hemi=hemi; nd.region=region; nd.contact=contact;
if contact>0, nd.label=sprintf('%s_%s_%d',hemi,region,contact);
else,         nd.label=sprintf('%s_%s',hemi,region); end
end

function M = assemble_epoch(R,e,k,win,regionUnit,collapse,nChan)
% stack all regions' nodes for event e, epoch k into nChan x numel(win)
M = zeros(nChan, numel(win));
row = 0; ok = true;
for r = 1:numel(R)
    x = R(r).byEvent{e}{k};                 % nContact x nSamp
    if isempty(x) || size(x,1) ~= R(r).nContact || size(x,2) < win(end)
        ok=false; break;                    % empty / placeholder / too short
    end
    x = x(:, win);
    if ~all(isfinite(x(:))), ok=false; break; end   % drop NaN/Inf epochs
    if regionUnit
        x = collapse_contacts(x, collapse);  % 1 x nSamp
    end
    nr = size(x,1);
    M(row+1:row+nr, :) = x;
    row = row + nr;
end
if ~ok || row ~= nChan, M = []; end
end

function y = collapse_contacts(x, method)
if size(x,1)==1, y=x; return; end
switch lower(method)
    case 'mean',  y = mean(x,1);
    case 'first', y = x(1,:);
    case 'pca1'
        xc = x - mean(x,2);
        [U,~,~] = svd(xc,'econ');
        y = U(:,1)' * xc;
    otherwise, error('collapse_contacts:method','unknown method %s',method);
end
end

function conn = classify_connections(ci)
n = numel(ci);
conn.nodeLabel = {ci.label}';
conn.intraL = false(n); conn.intraR = false(n); conn.inter = false(n);
conn.typeLabel = repmat({''}, n, n);
for i = 1:n            % i = target
    for j = 1:n        % j = source ; matrix(i,j) = j->i
        if i==j, continue; end
        hi = ci(i).hemi; hj = ci(j).hemi;
        if strcmp(hi,'L') && strcmp(hj,'L')
            conn.intraL(i,j)=true; conn.typeLabel{i,j}='intra-L';
        elseif strcmp(hi,'R') && strcmp(hj,'R')
            conn.intraR(i,j)=true; conn.typeLabel{i,j}='intra-R';
        else
            conn.inter(i,j)=true;  conn.typeLabel{i,j}='inter';
        end
    end
end
end
