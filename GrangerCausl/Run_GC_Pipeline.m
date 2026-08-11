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
%             Compute_GC_V2.m, Autoregressive_Process_V1.m,
%             Select_Model_Order_MT.m  (multi-trial, seam-free; replaces the
%                 single-segment Select_Model_Order for order selection),
%             Compute_Spectral_GC.m, Extract_Band_Power.m  (modes below)
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
%   'mode'        ['voltage'] 'voltage'   = time-domain GC on raw voltage (original)
%                             'bandpower' = GC on log band power (amplitude
%                                           envelope coupling; needs 'band')
%                             'spectral'  = Geweke frequency-resolved GC from a
%                                           broadband VAR; NO filtering applied
%   'band'        []        [f1 f2] Hz, required for mode 'bandpower'
%   'envLowpass'  []        Hz; envelope timescale. [] = min(0.5*bandwidth,25).
%                           This is a scientific choice, not a tuning knob.
%   'bands'       []        struct of band definitions for mode 'spectral'.
%                           [] = delta/theta/alpha/beta/lowGamma/highGamma
%   'specBand'    ''        which band's GC goes into .GC so the existing
%                           plotters work unchanged. '' = time-domain GC.
%   'blockTrials' []        trials pooled per VAR fit when per-trial fitting is
%                           not viable (short decimated envelopes). [] = auto.
%   'minFitSamples' [600]   target samples per VAR fit, drives the auto block size
%   'nFreq'       [257]     frequency bins for mode 'spectral'
%   'verbose'     [1]       0 = silent, 1 = normal progress,
%                           2 = detailed (per-fit progress, ETA, diagnostics)
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
ip.addParameter('mode','voltage',@ischar);
ip.addParameter('band',[],@(x)isempty(x)||numel(x)==2);
ip.addParameter('envLowpass',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('bands',[],@(x)isempty(x)||isstruct(x));
ip.addParameter('specBand','',@ischar);
ip.addParameter('blockTrials',[],@(x)isempty(x)||isscalar(x));
ip.addParameter('minFitSamples',600,@isscalar);
ip.addParameter('nFreq',257,@isscalar);
ip.addParameter('verbose',1,@isscalar);   % 0 silent | 1 normal | 2 detailed
ip.addParameter('zscore',true,@islogical);   % Volts are raw -> z-score per epoch/channel
ip.addParameter('subjects',{},@iscell);
ip.addParameter('save',true,@islogical);
ip.parse(dataDir,outDir,varargin{:});
opt = ip.Results;
regionUnit = strcmpi(opt.channelUnit,'region');

% ---------- mode validation ----------
opt.mode = lower(opt.mode);
if ~any(strcmp(opt.mode,{'voltage','bandpower','spectral'}))
    error('Run_GC_Pipeline:mode','mode must be voltage | bandpower | spectral.');
end
if strcmp(opt.mode,'bandpower') && isempty(opt.band)
    error('Run_GC_Pipeline:band','mode ''bandpower'' requires ''band'',[f1 f2].');
end
if strcmp(opt.mode,'bandpower') && opt.band(2) >= opt.fs/2
    error('Run_GC_Pipeline:nyquist','band upper edge %.1f >= Nyquist %.1f.', ...
          opt.band(2), opt.fs/2);
end
if isempty(opt.bands)
    opt.bands = struct('delta',[1 4],'theta',[4 8],'alpha',[8 12], ...
                       'beta',[13 30],'lowGamma',[30 70],'highGamma',[70 150]);
end
% drop spectral bands that exceed Nyquist
bfn = fieldnames(opt.bands);
for b = 1:numel(bfn)
    if opt.bands.(bfn{b})(2) >= opt.fs/2
        fprintf('   note: dropping band %s (above Nyquist %.0f Hz)\n', bfn{b}, opt.fs/2);
        opt.bands = rmfield(opt.bands, bfn{b});
    end
end
% band-power transform spec, threaded into assemble_epoch
bp = [];
if strcmp(opt.mode,'bandpower')
    bp = struct('fs',opt.fs,'band',opt.band,'envLowpass',opt.envLowpass);
end
V = opt.verbose;
tPipe = tic;
if V >= 1
    fprintf('\n');
    fprintf('===============================================================\n');
    fprintf(' Run_GC_Pipeline\n');
    fprintf('===============================================================\n');
    fprintf('  mode          : %s\n', opt.mode);
    if strcmp(opt.mode,'bandpower')
        fprintf('  band          : %g - %g Hz  (log power)\n', opt.band(1), opt.band(2));
        if isempty(opt.envLowpass)
            fprintf('  envLowpass    : auto = %g Hz\n', ...
                    min(0.5*(opt.band(2)-opt.band(1)), 25));
        else
            fprintf('  envLowpass    : %g Hz\n', opt.envLowpass);
        end
    end
    fprintf('  events        : %s\n', strjoin(opt.events, ', '));
    fprintf('  fs            : %g Hz\n', opt.fs);
    fprintf('  node unit     : %s', opt.channelUnit);
    if regionUnit, fprintf(' (collapse = %s)', opt.collapse); end
    fprintf('\n');
    if isempty(opt.order)
        fprintf('  order         : auto (%s, maxOrder %d)\n', opt.criterion, opt.maxOrder);
    else
        fprintf('  order         : %d (fixed)\n', opt.order);
    end
    fprintf('  nPerm         : %d%s\n', opt.nPerm, ...
            tern_(opt.nPerm==0, '   (significance skipped)', ''));
    fprintf('  dataDir       : %s\n', dataDir);
    fprintf('  outDir        : %s\n', outDir);
    fprintf('===============================================================\n');
end

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
if V >= 1
    fprintf('\nFound %d region files across %d subjects; analysing %d.\n', ...
            nF, numel(unique(subj)), numel(subjects));
end

results = struct();
nEvent  = numel(opt.events);

% ---------- per subject ----------
failed = {};
for si = 1:numel(subjects)
  try
    sub = subjects{si};
    idx = find(strcmp(subj,sub));
    tSub = tic;
    if V >= 1
        fprintf('\n---------------------------------------------------------------\n');
        fprintf('[%d/%d] %s : %d region files\n', si, numel(subjects), sub, numel(idx));
        fprintf('---------------------------------------------------------------\n');
    end

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
        if V >= 1
            fprintf('   load %-10s %2d contacts | epochs/event = [%s]', ...
                [R(r).hemi '_' R(r).region], R(r).nContact, ...
                strtrim(sprintf('%d ', cellfun(@numel, R(r).byEvent))));
            if V >= 2
                fprintf(' | %d raw epochs, %d valid', numel(epochs), numel(chCounts));
                if numel(unique(chCounts)) > 1
                    fprintf(' | WARNING contact count varies: %s', ...
                            mat2str(unique(chCounts)));
                end
            end
            fprintf('\n');
        end
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
    if V >= 1
        fprintf('   window     : %d samples (%.0f ms) starting at sample %d\n', ...
                commonLen, 1000*commonLen/opt.fs, opt.winStart);
    end

    % ---------- band-power geometry (same for every region -> compute once) ----------
    fsAnalysis = opt.fs; nSampEff = commonLen; bpInfo = [];
    if ~isempty(bp)
        [~, fsAnalysis, bpInfo] = Extract_Band_Power(zeros(1,commonLen), opt.fs, ...
            opt.band, 'envLowpass', opt.envLowpass, 'verbose', false);
        nSampEff = bpInfo.nSamp;
        if V >= 1
            fprintf(['   bandpower  : trim %d/side, decimate x%d, ' ...
                     'fs %g -> %g Hz, %d samp/epoch\n'], bpInfo.nTrim, ...
                     bpInfo.step, opt.fs, fsAnalysis, nSampEff);
        end
    end

    % ---------- per-trial vs pooled fitting ----------
    % A decimated envelope can leave far too few samples per trial to fit a VAR.
    % When that happens we pool trials into blocks: the lag matrix is still built
    % WITHIN each trial and only then concatenated, so no lags cross trial seams.
    perTrial = ~strcmp(opt.mode,'spectral') && nSampEff >= opt.minFitSamples;
    if ~isempty(opt.blockTrials)
        blockTrials = max(1, round(opt.blockTrials));
        perTrial = perTrial && blockTrials == 1;
    elseif perTrial
        blockTrials = 1;
    else
        blockTrials = max(1, ceil(opt.minFitSamples / nSampEff));
    end
    if blockTrials > 1
        if V >= 1
            fprintf(['   pooling    : %d trials per VAR fit ' ...
                     '(%d samp/trial < minFitSamples %d)\n'], ...
                     blockTrials, nSampEff, opt.minFitSamples);
        end
    end

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

    % ---------- PASS 1: assemble every epoch of every event ----------
    % Done up front so that (a) order selection sees real epochs rather than a
    % concatenated stack, and (b) pooled blocks can be built from the trials
    % valid in ALL events, which is what makes the paired contrast genuine.
    EPall = cell(nEvent,1); okAll = cell(nEvent,1);
    if V >= 1, fprintf('   assembling epochs...\n'); end
    for e = 1:nEvent
        EP = cell(nEpoch(e),1); okEp = false(nEpoch(e),1);
        why = cell(nEpoch(e),1);
        for k = 1:nEpoch(e)
            [M, w] = assemble_epoch(R,e,k,win,regionUnit,opt.collapse,nChan,bp);
            why{k} = w;
            if isempty(M), continue; end
            if opt.zscore, M = zsc_(M); end   % after the log transform, not before
            EP{k} = M; okEp(k) = true;
        end
        EPall{e} = EP; okAll{e} = okEp;
        if V >= 1
            fprintf('      %-14s %3d/%3d usable', opt.events{e}, sum(okEp), nEpoch(e));
            bad = why(~okEp);
            if ~isempty(bad)
                u = unique(bad);
                parts = {};
                for q = 1:numel(u)
                    parts{end+1} = sprintf('%s x%d', u{q}, sum(strcmp(bad,u{q}))); %#ok<AGROW>
                end
                fprintf('   skipped: %s', strjoin(parts, ', '));
            end
            fprintf('\n');
        end
    end

    % trials usable in EVERY event -> the only ones that can be truly paired
    nTrCommon = min(cellfun(@numel, okAll));
    commonOK  = true(nTrCommon,1);
    for e = 1:nEvent, commonOK = commonOK & okAll{e}(1:nTrCommon); end
    idxCommon = find(commonOK);
    if isempty(idxCommon)
        warning('%s: no trial is valid across all events; skipping subject.', sub);
        continue;
    end
    if V >= 1
        fprintf('   common     : %d trials valid in all %d events\n', ...
                numel(idxCommon), nEvent);
    end

    % ---------- model order: select ONCE, then fix ----------
    if isempty(opt.order)
        % NOTE: order selection runs on the TRANSFORMED data. The optimal order
        % for a decimated log envelope is a completely different number from the
        % optimal order for voltage. Epochs are passed as separate realizations
        % so no design row straddles an epoch boundary.
        pool = {}; got = 0;
        nWant = opt.orderEpochs;
        if nSampEff < opt.minFitSamples          % short epochs -> need more of them
            nWant = max(nWant, ceil(4*opt.minFitSamples / max(nSampEff,1)));
        end
        for k = 1:numel(idxCommon)
            if got >= nWant, break; end
            pool{end+1} = EPall{1}{idxCommon(k)}; got = got+1; %#ok<AGROW>
        end
        if isempty(pool)
            warning('%s: no valid epochs for order selection; skipping subject.', sub);
            continue;
        end
        [p, orderInfo] = Select_Model_Order_MT(pool, opt.maxOrder, opt.criterion);
        if V >= 1
            fprintf(['   order sel  : %d epochs x %d samp | limits: rows %d, ' ...
                     'len %d -> tested 1..%d\n'], numel(pool), nSampEff, ...
                     orderInfo.maxByRows, orderInfo.maxByLen, ...
                     numel(orderInfo.orders));
        end
        if V >= 2
            fprintf('              AIC picks %d, BIC picks %d\n', ...
                    orderInfo.aicOrder, orderInfo.bicOrder);
            oo = orderInfo.orders; cc = orderInfo.(upper(opt.criterion));
            kshow = unique(round(linspace(1, numel(oo), min(8,numel(oo)))));
            fprintf('              %s curve:', upper(opt.criterion));
            for q = kshow, fprintf('  p%d=%.4f', oo(q), cc(q)); end
            fprintf('\n');
        end
        if orderInfo.railed
            fprintf(['   *** order %d sits AT the tested ceiling -- the true ' ...
                     'order is probably higher ***\n'], p);
        end
    else
        p = opt.order; orderInfo = [];
    end
    if V >= 1
        fprintf('   VAR order  : %d   (null GC bias floor ~ p/N = %.5f)\n', ...
                p, p/(nSampEff*blockTrials));
    end

    % ---------- PASS 2: GC per event ----------
    rcondMin = Inf;
    bandNames = fieldnames(opt.bands);
    if V >= 1
        if blockTrials == 1
            nFitTot = 0;
            for e = 1:nEvent, nFitTot = nFitTot + sum(okAll{e}); end
        else
            nFitTot = nEvent * floor(numel(idxCommon)/blockTrials);
        end
        nPairs  = nChan*(nChan-1);
        nVarFit = nFitTot * nPairs * (1 + opt.nPerm);
        fprintf(['   computing GC: %d fits x %d ordered pairs x %d ' ...
                 '(1 obs + %d perm) = %s VAR solves\n'], ...
                 nFitTot, nPairs, 1+opt.nPerm, opt.nPerm, human_(nVarFit));
    end
    for e = 1:nEvent
        evt = opt.events{e};

        if blockTrials == 1
            % ---- one fit per trial, indexed BY TRIAL so slice k is trial k.
            % This preserves the original layout exactly: Plot_GC_Summary and
            % Plot_GC_Contrast2 pair events by trial index, so the cube must
            % stay trial-indexed with NaN in the unusable slots.
            nSlice = nEpoch(e);
            blocks = cell(nSlice,1);
            for k = 1:nSlice
                if okAll{e}(k), blocks{k} = k; end
            end
            sliceValid = okAll{e};
        else
            % ---- pooled fits over the COMMON trial set, so block b holds the
            % same trials in every event and the paired contrast stays paired.
            nSlice = floor(numel(idxCommon) / blockTrials);
            blocks = cell(nSlice,1);
            for b = 1:nSlice
                blocks{b} = idxCommon((b-1)*blockTrials+1 : b*blockTrials);
            end
            sliceValid = true(nSlice,1);
        end

        GCcube  = nan(nChan,nChan,nSlice);
        SIGcube = false(nChan,nChan,nSlice);
        Pcube   = nan(nChan,nChan,nSlice);
        specSum = zeros(nChan,nChan,opt.nFreq); specN = 0;
        bandCube = struct();
        for b = 1:numel(bandNames), bandCube.(bandNames{b}) = nan(nChan,nChan,nSlice); end
        fAxis = [];

        tEvt = tic; nDone = 0; nRep = 0;
        nToDo = sum(~cellfun(@isempty, blocks));
        for b = 1:nSlice
            selK = blocks{b};
            if isempty(selK), continue; end
            L  = min(cellfun(@(a) size(a,2), EPall{e}(selK)));
            X3 = zeros(nChan, L, numel(selK));
            for q = 1:numel(selK), X3(:,:,q) = EPall{e}{selK(q)}(:,1:L); end

            if blockTrials == 1 && ~strcmp(opt.mode,'spectral')
                % ---- original path: time-domain GC, one trial ----
                GC = Compute_GC_V2(X3(:,:,1), 'order',p, 'nPerm',opt.nPerm, ...
                                   'alpha',opt.alpha, 'zscore',false);
                GCcube(:,:,b)  = GC.matrix;
                SIGcube(:,:,b) = GC.sig;
                Pcube(:,:,b)   = GC.pval;
                rcondMin = min(rcondMin, GC.rcondMin);
            else
                % ---- pooled VAR: time-domain + Geweke spectral decomposition ----
                Sg = Compute_Spectral_GC(X3, p, fsAnalysis, 'nFreq',opt.nFreq, ...
                        'bands',opt.bands, 'nPerm',opt.nPerm, 'alpha',opt.alpha, ...
                        'zscore',false, 'verbose',false);
                fAxis   = Sg.f;
                specSum = specSum + fill_nan(Sg.spec); specN = specN + 1;
                for bb = 1:numel(bandNames)
                    bandCube.(bandNames{bb})(:,:,b) = Sg.band.(bandNames{bb});
                end
                if isempty(opt.specBand)
                    GCcube(:,:,b) = Sg.td;
                    if isfield(Sg,'sig')
                        SIGcube(:,:,b) = Sg.sig; Pcube(:,:,b) = Sg.pval;
                    end
                else
                    GCcube(:,:,b) = Sg.band.(opt.specBand);
                    if isfield(Sg,'bandSig')
                        SIGcube(:,:,b) = Sg.bandSig.(opt.specBand);
                        Pcube(:,:,b)   = Sg.bandPval.(opt.specBand);
                    end
                end
                rcondMin = min(rcondMin, Sg.rcondMin);
            end

            % ---- progress with ETA (report ~10 times, never more) ----
            nDone = nDone + 1;
            if V >= 1 && nToDo > 0
                step = max(1, floor(nToDo/10));
                if mod(nDone, step) == 0 || nDone == nToDo
                    el  = toc(tEvt);
                    eta = el * (nToDo - nDone) / max(nDone,1);
                    nRep = nRep + 1;
                    fprintf('      %-12s [%s] %3d/%3d  %5.1fs elapsed', ...
                            evt, bar_(nDone/nToDo, 20), nDone, nToDo, el);
                    if nDone < nToDo
                        fprintf('  ~%.0fs left', eta);
                    end
                    fprintf('\n');
                end
            end
        end

        % ---- what actually came out ----
        if V >= 1
            offd = ~eye(nChan);
            flat  = reshape(GCcube, nChan*nChan, []);
            % count NaN only among slices that actually held a fit -- the
            % trial-indexed cube is NaN by construction wherever an epoch was
            % skipped, and calling that "ill-conditioned" is simply wrong.
            offdF = flat(offd(:), logical(sliceValid(:))');
            nBad  = sum(~isfinite(offdF(:)));
            mG    = mean(offdF(isfinite(offdF)));
            floorV = p / (nSampEff*blockTrials);
            fprintf('      %-12s mean off-diag GC %.5f (floor %.5f, %.2fx)', ...
                    evt, mG, floorV, mG/floorV);
            if opt.nPerm > 0
                sflat = reshape(SIGcube, nChan*nChan, []);
                sv = sflat(offd(:),:);
                fprintf('  |  sig rate %.1f%%', 100*mean(double(sv(:))));
            end
            if nBad > 0
                fprintf('  |  %d/%d NaN among FITTED slices (ill-conditioned)', ...
                        nBad, numel(offdF));
            end
            if isfinite(mG) && mG < floorV
                fprintf('\n      %-12s *** mean is BELOW the null bias floor -- ', evt);
                fprintf('run Diagnose_GC_Floor ***');
            end
            fprintf('\n');
        end

        results.(sub).events.(evt).GC         = GCcube;
        results.(sub).events.(evt).sig        = SIGcube;
        results.(sub).events.(evt).pval       = Pcube;
        results.(sub).events.(evt).validEpoch = sliceValid;
        results.(sub).events.(evt).blockTrial = blocks;   % provenance per slice
        if specN > 0
            results.(sub).events.(evt).spec   = specSum / specN;
            results.(sub).events.(evt).f      = fAxis;
            results.(sub).events.(evt).bandGC = bandCube;
        end

    end

    % ---------- connection classification ----------
    results.(sub).channelInfo = channelInfo;
    results.(sub).conn        = classify_connections(channelInfo);
    results.(sub).order       = p;
    results.(sub).rcondMin    = rcondMin;
    results.(sub).events_list = opt.events;
    % winSamples = effective N behind ONE VAR fit. Plot_GC_Contrast2 reads this
    % to draw the null bias floor (~order/N); it was never written before, so
    % the floor silently evaluated to NaN.
    results.(sub).winSamples  = nSampEff * blockTrials;
    results.(sub).mode        = opt.mode;
    results.(sub).fsAnalysis  = fsAnalysis;
    results.(sub).blockTrials = blockTrials;
    results.(sub).bandInfo    = bpInfo;
    results.(sub).orderInfo   = orderInfo;
    results.(sub).nTrialCommon= numel(idxCommon);
    if strcmp(opt.mode,'bandpower'), results.(sub).band = opt.band; end
    if ~strcmp(opt.mode,'voltage'), results.(sub).bands = opt.bands; end

    if V >= 1
        fprintf('   min rcond  : %.2e%s\n', rcondMin, ...
                tern_(rcondMin < 1e-6, '   *** ILL-CONDITIONED ***', ''));
        nL = sum(strcmp({channelInfo.hemi},'L'));
        nR = sum(strcmp({channelInfo.hemi},'R'));
        cn = results.(sub).conn;
        fprintf('   nodes      : %d  (L=%d, R=%d)  edges: intra-L %d, intra-R %d, inter %d\n', ...
                nChan, nL, nR, sum(cn.intraL(:)), sum(cn.intraR(:)), sum(cn.inter(:)));
        fprintf('   subject done in %.1f s\n', toc(tSub));
    end

    if opt.save
        outFile = fullfile(outDir, sprintf('%s_GC.mat', sub));
        subjResult = results.(sub); %#ok<NASGU>
        save(outFile, 'subjResult', '-v7.3');
        if V >= 1, fprintf('   saved      : %s\n', outFile); end
    end

  catch err
    % A batch run is hours long. Losing all of it to one malformed file is
    % not acceptable, so record and move on rather than aborting.
    failed{end+1} = subjects{si}; %#ok<AGROW>
    fprintf(2, '\n   *** %s FAILED: %s\n', subjects{si}, err.message);
    if ~isempty(err.stack)
        fprintf(2, '       at %s line %d\n', err.stack(1).name, err.stack(1).line);
    end
  end
end
results_failed = failed;
if V >= 1
    fprintf('\n===============================================================\n');
    fprintf(' Done: %d/%d subject(s) in %.1f s\n', ...
            numel(fieldnames(results)), numel(subjects), toc(tPipe));
    if ~isempty(failed)
        fprintf(' FAILED (%d): %s\n', numel(failed), strjoin(failed, ', '));
    end
    % Cross-subject comparability check -- absolute GC is only comparable
    % across subjects when BOTH the model order and the window match, since
    % the bias floor is p/N.
    fn = fieldnames(results);
    if numel(fn) > 1
        ords = zeros(numel(fn),1); wins = zeros(numel(fn),1);
        for q = 1:numel(fn)
            ords(q) = results.(fn{q}).order;
            wins(q) = results.(fn{q}).winSamples;
        end
        if numel(unique(ords)) > 1 || numel(unique(wins)) > 1
            fprintf([' *** WARNING: order varies %s, N varies %s across subjects.\n' ...
                     '     The null bias floor is p/N, so ABSOLUTE GC is not\n' ...
                     '     comparable between them. For group analysis re-run with\n' ...
                     '     fixed ''order'' and ''winSamples''. Within-subject event\n' ...
                     '     contrasts are unaffected. ***\n'], ...
                     mat2str(unique(ords)'), mat2str(unique(wins)'));
        end
    end
    fprintf('===============================================================\n\n');
end
end

% ======================= local functions =======================

function s = bar_(frac, width)
n = max(0, min(width, round(frac*width)));
s = [repmat('=',1,n), repmat('.',1,width-n)];
end

function s = human_(x)
if x >= 1e6, s = sprintf('%.1fM', x/1e6);
elseif x >= 1e3, s = sprintf('%.1fk', x/1e3);
else, s = sprintf('%d', round(x)); end
end

function s = tern_(c,a,b), if c, s=a; else, s=b; end, end

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

function [M, why] = assemble_epoch(R,e,k,win,regionUnit,collapse,nChan,bp)
% Stack all regions' nodes for event e, epoch k.
% If bp is non-empty, band power is extracted PER CONTACT and only then are
% contacts collapsed. Order matters: averaging raw voltage across bipolar
% derivations first lets their band activity phase-cancel, which suppresses
% the very amplitude signal we are trying to measure. Envelope-then-average
% is the standard and is what is done here.
if nargin < 8, bp = []; end
M = []; why = 'ok';
rows = {}; nr = 0;
for r = 1:numel(R)
    x = R(r).byEvent{e}{k};                 % nContact x nSamp
    if isempty(x),                    why = 'empty';       return; end
    if size(x,1) ~= R(r).nContact,    why = 'placeholder'; return; end
    if size(x,2) < win(end),          why = 'short';       return; end
    x = x(:, win);
    if ~all(isfinite(x(:))),          why = 'NaN/Inf';     return; end
    if ~isempty(bp)
        x = Extract_Band_Power(x, bp.fs, bp.band, 'output','logpower', ...
                'envLowpass', bp.envLowpass, 'verbose', false);
        if isempty(x) || ~all(isfinite(x(:))), why = 'bandpower-fail'; return; end
    end
    if regionUnit
        x = collapse_contacts(x, collapse);  % 1 x nSamp
    end
    rows{end+1} = x; %#ok<AGROW>
    nr = nr + size(x,1);
end
if nr ~= nChan, why = 'node-count'; return; end
L = min(cellfun(@(a) size(a,2), rows));
for q = 1:numel(rows), rows{q} = rows{q}(:,1:L); end
M = vertcat(rows{:});
end

function Z = zsc_(M)
% per-row z-score, toolbox-free (normalize() is not in every Octave build)
mu = mean(M,2); sd = std(M,0,2); sd(sd < eps) = 1;
Z = (M - mu(:,ones(1,size(M,2)))) ./ sd(:,ones(1,size(M,2)));
end

function Y = fill_nan(X)
Y = X; Y(~isfinite(Y)) = 0;
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
