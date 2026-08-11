function [M, info] = Get_GC_Epoch(dataDir, subject, event, trialIdx, varargin)
%% [M, info] = Get_GC_Epoch(dataDir, subject, event, trialIdx, 'name',value, ...)
% Rebuilds the SAME nChan x nSamp node matrix that Run_GC_Pipeline feeds to
% Compute_GC_V2, for one subject / event / trial. Exists because the driver's
% assemble_epoch is a local function, so there was previously no way to get
% at the matrix the pipeline actually fits -- which you need for any
% diagnostic that has to re-run the estimator.
%
% Pass the SAME options you passed to Run_GC_Pipeline, or the matrix will not
% match what the pipeline used.
%
% INPUT
%   dataDir  folder of {SUBJECT}_{HEMI}_{REGION}_TrialDATA.mat files
%   subject  e.g. 'CLASE001'
%   event    e.g. 'CHOICE'
%   trialIdx trial number, or [] to return every usable trial in info.all
%
% PARAMETERS (defaults match Run_GC_Pipeline)
%   'events' {'CHOICE','RESPONSEON','OUTCOME'}   needed to reproduce the
%            common window length, which is the min over ALL events
%   'fs' [500]  'channelUnit' ['region']  'collapse' ['mean']
%   'winStart' [1]  'winSamples' []  'minEpochSamp' [50]  'zscore' [true]
%   'mode' ['voltage']  'band' []  'envLowpass' []
%   'verbose' [true]
%
% OUTPUT
%   M     nChan x nSamp, z-scored, windowed, collapsed -- ready for GC
%   info  .channelInfo .nodeLabel .commonLen .win .nEpoch .valid .all
%         .fsAnalysis .why (skip reason if M is empty)
%
% EXAMPLE
%   M = Get_GC_Epoch('/data/regions','CLASE001','CHOICE',1);
%   out = Diagnose_GC_Floor('CLASE001_GC.mat','epoch',M);

ip = inputParser;
ip.addParameter('events',{'CHOICE','RESPONSEON','OUTCOME'},@iscell);
ip.addParameter('fs',500,@isscalar);
ip.addParameter('channelUnit','region',@ischar);
ip.addParameter('collapse','mean',@ischar);
ip.addParameter('winStart',1,@isscalar);
ip.addParameter('winSamples',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('minEpochSamp',50,@isscalar);
ip.addParameter('zscore',true,@islogical);
ip.addParameter('mode','voltage',@ischar);
ip.addParameter('band',[],@(x)isempty(x)||numel(x)==2);
ip.addParameter('envLowpass',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('verbose',true,@(x)islogical(x)||isnumeric(x));
ip.parse(varargin{:});
opt = ip.Results;
regionUnit = strcmpi(opt.channelUnit,'region');
V = opt.verbose;

bp = [];
if strcmpi(opt.mode,'bandpower')
    if isempty(opt.band)
        error('Get_GC_Epoch:band','mode ''bandpower'' requires ''band'',[f1 f2].');
    end
    bp = struct('fs',opt.fs,'band',opt.band,'envLowpass',opt.envLowpass);
end

% ---------- discover this subject's region files ----------
D = dir(fullfile(dataDir,'*_TrialDATA.mat'));
if isempty(D), error('Get_GC_Epoch:noFiles','No *_TrialDATA.mat in %s', dataDir); end
keep = {}; hemis = {}; regs = {};
for i = 1:numel(D)
    base = regexprep(D(i).name, '_TrialDATA\.mat$', '');
    tok  = strsplit(base,'_');
    if numel(tok) < 3, continue; end
    if ~strcmp(tok{1}, subject), continue; end
    keep{end+1}  = fullfile(D(i).folder, D(i).name); %#ok<AGROW>
    hemis{end+1} = tok{2};                            %#ok<AGROW>
    regs{end+1}  = strjoin(tok(3:end),'_');           %#ok<AGROW>
end
if isempty(keep)
    error('Get_GC_Epoch:noSubject','No files for subject %s in %s', subject, dataDir);
end
nReg   = numel(keep);
nEvent = numel(opt.events);
eSel   = find(strcmp(opt.events, event));
if isempty(eSel)
    error('Get_GC_Epoch:event','Event ''%s'' not in the events list.', event);
end

% ---------- load, exactly as the driver does ----------
R = struct('hemi',{},'region',{},'nContact',{},'byEvent',{});
for r = 1:nReg
    [labels, epochs] = load_region_epochs(keep{r});
    R(r).hemi = hemis{r}; R(r).region = regs{r};
    chCounts = [];
    for kk = 1:numel(epochs)
        x = epochs{kk};
        if ~isempty(x) && size(x,2) >= opt.minEpochSamp && all(isfinite(x(:)))
            chCounts(end+1) = size(x,1); %#ok<AGROW>
        end
    end
    if isempty(chCounts), R(r).nContact = 0; else, R(r).nContact = mode(chCounts); end
    R(r).byEvent = cell(nEvent,1);
    for e = 1:nEvent
        R(r).byEvent{e} = epochs(strcmp(labels, opt.events{e}));
    end
    if V
        fprintf('   %s_%s : %d contacts\n', R(r).hemi, R(r).region, R(r).nContact);
    end
end

nEpoch = zeros(nEvent,1);
for e = 1:nEvent
    nEpoch(e) = min(arrayfun(@(x) numel(x.byEvent{e}), R));
end

% ---------- common window: min over ALL events, as the driver does ----------
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
    commonLen = min(allLen);
else
    commonLen = opt.winSamples;
end
win = opt.winStart : (opt.winStart + commonLen - 1);

% ---------- node bookkeeping ----------
channelInfo = struct('hemi',{},'region',{},'contact',{},'label',{});
node = 0;
for r = 1:nReg
    if regionUnit
        node = node+1;
        channelInfo(node).hemi=R(r).hemi; channelInfo(node).region=R(r).region;
        channelInfo(node).contact=0;
        channelInfo(node).label=sprintf('%s_%s',R(r).hemi,R(r).region);
    else
        for c = 1:R(r).nContact
            node = node+1;
            channelInfo(node).hemi=R(r).hemi; channelInfo(node).region=R(r).region;
            channelInfo(node).contact=c;
            channelInfo(node).label=sprintf('%s_%s_%d',R(r).hemi,R(r).region,c);
        end
    end
end
nChan = node;

fsAnalysis = opt.fs;
if ~isempty(bp)
    [~, fsAnalysis] = Extract_Band_Power(zeros(1,commonLen), opt.fs, opt.band, ...
        'envLowpass', opt.envLowpass, 'verbose', false);
end

info = struct('channelInfo',channelInfo,'nodeLabel',{{channelInfo.label}'}, ...
              'commonLen',commonLen,'win',win,'nEpoch',nEpoch(eSel), ...
              'fsAnalysis',fsAnalysis,'why','');

% ---------- assemble ----------
if isempty(trialIdx)
    allM = cell(nEpoch(eSel),1); valid = false(nEpoch(eSel),1);
    for k = 1:nEpoch(eSel)
        [m, w] = build_one(R,eSel,k,win,regionUnit,opt,nChan,bp);
        if ~isempty(m), allM{k} = m; valid(k) = true; end
        info.why = w;
    end
    info.all = allM; info.valid = valid;
    M = allM(valid);
    if V, fprintf('   %s: %d/%d trials usable\n', event, sum(valid), nEpoch(eSel)); end
    return;
end

[M, why] = build_one(R,eSel,trialIdx,win,regionUnit,opt,nChan,bp);
info.why = why; info.valid = ~isempty(M);
if isempty(M)
    warning('Get_GC_Epoch:skip','Trial %d of %s unusable (%s).', trialIdx, event, why);
elseif V
    fprintf('   built M: %d nodes x %d samples (%s)\n', size(M,1), size(M,2), why);
end
end

% ================================================================
function [M, why] = build_one(R,e,k,win,regionUnit,opt,nChan,bp)
% Byte-identical to Run_GC_Pipeline's assemble_epoch + its z-score step.
M = []; why = 'ok';
rows = {}; nr = 0;
if k < 1 || k > numel(R(1).byEvent{e}), why = 'out-of-range'; return; end
for r = 1:numel(R)
    x = R(r).byEvent{e}{k};
    if isempty(x),                 why = 'empty';       return; end
    if size(x,1) ~= R(r).nContact, why = 'placeholder'; return; end
    if size(x,2) < win(end),       why = 'short';       return; end
    x = x(:, win);
    if ~all(isfinite(x(:))),       why = 'NaN/Inf';     return; end
    if ~isempty(bp)
        x = Extract_Band_Power(x, bp.fs, bp.band, 'output','logpower', ...
                'envLowpass', bp.envLowpass, 'verbose', false);
        if isempty(x) || ~all(isfinite(x(:))), why = 'bandpower-fail'; return; end
    end
    if regionUnit, x = collapse_one(x, opt.collapse); end
    rows{end+1} = x; %#ok<AGROW>
    nr = nr + size(x,1);
end
if nr ~= nChan, why = 'node-count'; return; end
L = min(cellfun(@(a) size(a,2), rows));
for q = 1:numel(rows), rows{q} = rows{q}(:,1:L); end
M = vertcat(rows{:});
if opt.zscore
    mu = mean(M,2); sd = std(M,0,2); sd(sd < eps) = 1;
    M = (M - mu(:,ones(1,size(M,2)))) ./ sd(:,ones(1,size(M,2)));
end
end

function y = collapse_one(x, method)
if size(x,1)==1, y=x; return; end
switch lower(method)
    case 'mean',  y = mean(x,1);
    case 'first', y = x(1,:);
    case 'pca1'
        xc = x - mean(x,2);
        [U,~,~] = svd(xc,'econ');
        y = U(:,1)' * xc;
    otherwise, error('collapse_one:method','unknown method %s',method);
end
end
