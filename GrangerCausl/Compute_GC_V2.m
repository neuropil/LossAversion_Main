function GC = Compute_GC_V2(data,varargin)
%% GC = Compute_GC_V2(data, 'name', value, ...)
% Pairwise time-domain Granger causality with the pipeline's VERIFIED index
% convention, model-order handling, a conditioning guard, and a
% circular-shift null that preserves the source's autocorrelation.
%
% ---------------------------------------------------------------------------
% CONVENTION (confirmed on a 2-node ground-truth simulation):
%   GC.matrix(i,j) = Granger causality FROM channel j TO channel i.
%   The entry at (target, source) is the influence of SOURCE on TARGET.
% ---------------------------------------------------------------------------
%
% INPUT
%   data : nChan x nSamp, ONE epoch. Trim to a length common to all
%          conditions before calling (GC bias is length-dependent, so
%          unequal lengths bias across-event comparisons). Pipeline already
%          z-scores; set 'zscore',true to do it here instead.
%
% OPTIONS (defaults in brackets)
%   'order'        [] fixed VAR order; if empty, auto-select via BIC. Passing
%                     a single fixed order across all epochs/conditions is the
%                     recommended, bias-matched usage.
%   'maxOrder'     [30]   max order for auto-selection
%   'criterion'    ['BIC'] 'AIC' or 'BIC' for auto-selection
%   'nPerm'        [200]  surrogates for the null; 0 skips significance
%   'alpha'        [0.05] one-sided significance level
%   'minShiftFrac' [0.05] min circular shift as a fraction of nSamp
%   'rcondTol'     [1e-10] below this the pair is ill-conditioned -> NaN
%   'zscore'       [false] z-score each channel (row) before fitting
%
% OUTPUT (struct GC)
%   .matrix   nChan x nChan, (i,j)=j->i, NaN on the diagonal / bad pairs
%   .sig      logical, obs >= (1-alpha) percentile of the null
%   .pval     one-sided permutation p-value = (1+#{surr>=obs})/(nPerm+1)
%   .order    VAR order actually used
%   .rcondMin smallest rcond seen (watch for tiny values at high order)
%   .chanInfo cell of 'FROM src | TO tgt' labels matching .matrix
%
% NOTE: significance here is per-pair, per-epoch and UNCORRECTED. Pool
% .pval across epochs / pairs and apply your multiple-comparisons control
% (e.g. FDR over the ~nChan^2 x nEvents tests) downstream.

ip = inputParser;
ip.addRequired('data',@isnumeric);
ip.addParameter('order',[],@(x)isempty(x)||(isscalar(x)&&x>0));
ip.addParameter('maxOrder',30,@isscalar);
ip.addParameter('criterion','BIC',@ischar);
ip.addParameter('nPerm',200,@isscalar);
ip.addParameter('alpha',0.05,@isscalar);
ip.addParameter('minShiftFrac',0.05,@isscalar);
ip.addParameter('rcondTol',1e-10,@isscalar);
ip.addParameter('zscore',false,@islogical);
ip.parse(data,varargin{:});
opt = ip.Results;

if opt.zscore
    data = normalize(data,2,'zscore');
end
[nChan,nSamp] = size(data);

% ---- model order ----
if isempty(opt.order)
    p = Select_Model_Order(data,opt.maxOrder,opt.criterion);
else
    p = opt.order;
end

% ---- safe circular-shift range (avoid trivial near-zero / near-full shifts) ----
minShift = max(p, round(opt.minShiftFrac*nSamp));
maxRange = nSamp - 2*minShift;
if maxRange < 1
    minShift = p;  maxRange = nSamp - 2*minShift;
    if maxRange < 1, minShift = 1; maxRange = nSamp - 2; end
end

GC.matrix   = nan(nChan);
GC.sig      = false(nChan);
GC.pval     = nan(nChan);
GC.order    = p;
GC.rcondMin = Inf;
GC.chanInfo = cell(nChan);

for tgt = 1:nChan
    for src = 1:nChan
        if tgt == src, continue; end
        GC.chanInfo{tgt,src} = sprintf('FROM %d | TO %d',src,tgt);

        % reduced: target from its own past. full: target + source.
        reduced = data(tgt,:);
        full    = [data(tgt,:); data(src,:)];

        [Rr,~,rcR] = Autoregressive_Process_V1(reduced,p);
        [Rf,~,rcF] = Autoregressive_Process_V1(full,p);
        GC.rcondMin = min([GC.rcondMin, rcR, rcF]);

        if rcF < opt.rcondTol || rcR < opt.rcondTol
            warning('Compute_GC_V2:illConditioned',...
                'rcond below tol at target=%d src=%d (rcondF=%.2e). GC left NaN.',...
                tgt,src,rcF);
            continue;
        end

        vr  = var(Rr,0,2);
        vf  = var(Rf,0,2);
        obs = log(vr./vf(1,1));            % (tgt,src) = src -> tgt
        GC.matrix(tgt,src) = obs;

        % ---- circular-shift null on the SOURCE (preserves its autocorr) ----
        if opt.nPerm > 0
            surr = zeros(opt.nPerm,1);
            for it = 1:opt.nPerm
                s   = minShift + randi(maxRange);
                idx = mod((0:nSamp-1) + s, nSamp) + 1;
                fullS = [data(tgt,:); data(src,idx)];
                [RfS,~,~] = Autoregressive_Process_V1(fullS,p);
                vfS = var(RfS,0,2);
                surr(it) = log(vr./vfS(1,1));
            end
            thr = prctile(surr,100*(1-opt.alpha));
            GC.sig(tgt,src)  = obs >= thr;
            GC.pval(tgt,src) = (1 + sum(surr >= obs)) / (opt.nPerm + 1);
        end
    end
end
end
