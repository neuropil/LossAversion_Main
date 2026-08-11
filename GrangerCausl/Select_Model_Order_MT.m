function [pSel, info] = Select_Model_Order_MT(X, maxOrder, criterion)
%% [pSel, info] = Select_Model_Order_MT(X, maxOrder, criterion)
% Multi-trial VAR order selection by AIC/BIC. Drop-in replacement for
% Select_Model_Order when the data are many short epochs rather than one
% long segment.
%
% WHY THIS EXISTS
%   Select_Model_Order takes ONE nChan x nSamp matrix, so the driver used to
%   concatenate epochs horizontally. That fabricates lags across every epoch
%   boundary: with E epochs of length L and order p, (E-1)*p of the design
%   rows straddle a seam and pair unrelated samples. For 3 s voltage epochs
%   that is ~0.3% of rows and harmless; for a decimated 8-12 Hz envelope
%   (19 samples/epoch) it is over 20% of rows, which is not.
%   Here the lag matrix is built WITHIN each epoch and only then concatenated,
%   so no row ever crosses a boundary.
%
% INPUT
%   X         [nChan x nSamp x nTrial] numeric, or cell of [nChan x nSamp].
%             Epochs may differ in length in the cell form.
%   maxOrder  largest order to test (default 30)
%   criterion 'AIC' or 'BIC' (default 'BIC')
%
% OUTPUT
%   pSel  selected order
%   info  .orders .AIC .BIC .aicOrder .bicOrder .nEff .nChan
%         .maxIdentifiable  order ceiling from the sample budget AND from the
%                           length of a single epoch (see the cap below)
%         .maxByRows .maxByLen  the two constraints separately
%         .railed           true if the selection sat at the ceiling
%
% Criterion form matches Select_Model_Order:
%   AIC = log det(Sigma) + 2k/m ;  BIC = log det(Sigma) + k log(m)/m
% with k = nChan^2 * p and m the pooled effective sample count. log det is
% taken via Cholesky rather than det() so it does not underflow to -Inf on
% the near-singular fits that show up at high order.

if nargin < 3 || isempty(criterion), criterion = 'BIC'; end
if nargin < 2 || isempty(maxOrder),  maxOrder  = 30;    end

% ---- normalize input to a cell of [nChan x nSamp] ----
if ~iscell(X)
    if ndims(X) == 2, X = {X};
    else
        C = cell(1, size(X,3));
        for t = 1:size(X,3), C{t} = X(:,:,t); end
        X = C;
    end
end
good = false(1,numel(X));
for t = 1:numel(X)
    good(t) = ~isempty(X{t}) && ismatrix(X{t}) && all(isfinite(X{t}(:)));
end
X = X(good);
if isempty(X), error('Select_Model_Order_MT:noData','No valid epochs.'); end

nChan = size(X{1},1);
lens  = cellfun(@(a) size(a,2), X);
for t = 1:numel(X)
    if size(X{t},1) ~= nChan
        error('Select_Model_Order_MT:chan','Epoch %d has %d channels, expected %d.', ...
              t, size(X{t},1), nChan);
    end
end

% ---- cap maxOrder at what is identifiable ----
% TWO constraints, and the second is the one that bites on short epochs:
%   (a) pooled rows at order p = sum(L_t) - T*p, and each equation estimates
%       nChan*p coefficients; require >= 10 rows per parameter.
%   (b) p must be small relative to the length of a SINGLE epoch. Pooling more
%       epochs raises the total sample count but cannot buy lag depth inside a
%       short one: at p near L the surviving rows are dominated by each epoch's
%       initial conditions, the residual covariance shrinks spuriously, and BIC
%       rails at whatever ceiling it is given. Verified on a known VAR(4): with
%       19-sample epochs and constraint (a) alone, BIC picked 16; adding
%       p <= L/4 recovered 4 exactly.
maxByRows = max(1, floor(sum(lens) / (numel(lens) + 10*nChan)));
maxByLen  = max(1, floor(min(lens) / 4));
maxIdent  = min(maxByRows, maxByLen);
maxOrderUse = min(maxOrder, maxIdent);

orders = 1:maxOrderUse;
AIC = nan(numel(orders),1);
BIC = nan(numel(orders),1);
mEff = nan(numel(orders),1);

for oi = 1:numel(orders)
    p = orders(oi);
    use = lens > p + 1;
    if ~any(use), break; end

    Ys = {}; Zs = {};
    for t = find(use)
        x = X{t}; N = size(x,2); m = N - p;
        Ys{end+1} = x(:, p+1:N); %#ok<AGROW>
        Z = zeros(nChan*p, m);
        for k = 1:p
            Z((k-1)*nChan+1 : k*nChan, :) = x(:, p+1-k : N-k);
        end
        Zs{end+1} = Z; %#ok<AGROW>
    end
    Y = [Ys{:}]; Z = [Zs{:}];
    m = size(Y,2);
    if m <= nChan*p, break; end

    G = Z*Z.';
    if ~isfinite(rcond(G)) || rcond(G) < eps, break; end
    A = (Y*Z.') / G;
    R = Y - A*Z;
    Sigma = (R*R.') / m;

    ld = logdet_(Sigma);
    if ~isfinite(ld), break; end
    k = nChan*nChan*p;
    AIC(oi) = ld + 2*k/m;
    BIC(oi) = ld + k*log(m)/m;
    mEff(oi) = m;
end

okA = ~isnan(AIC); okB = ~isnan(BIC);
if ~any(okB) && ~any(okA)
    error('Select_Model_Order_MT:noFit','No order could be fit; data too short.');
end
ia = pick_min(AIC); ib = pick_min(BIC);

info.orders   = orders;
info.AIC      = AIC;
info.BIC      = BIC;
info.aicOrder = orders(ia);
info.bicOrder = orders(ib);
info.nEff     = mEff;
info.nChan    = nChan;
info.nTrial   = numel(X);
info.maxIdentifiable = maxIdent;
info.maxByRows       = maxByRows;
info.maxByLen        = maxByLen;

if strcmpi(criterion,'AIC'), pSel = orders(ia); else, pSel = orders(ib); end
info.railed = (pSel >= maxOrderUse) && (maxOrderUse >= 2);
if info.railed
    warning('Select_Model_Order_MT:rail', ...
        ['Selected order %d sits at the tested ceiling %d (maxOrder=%d, ' ...
         'identifiable limit=%d). The true order is probably higher; ' ...
         'lengthen the window or pool more trials.'], ...
        pSel, maxOrderUse, maxOrder, maxIdent);
end
end

% ---------------------------------------------------------------
function i = pick_min(v)
ok = find(~isnan(v));
if isempty(ok), i = 1; return; end
[~, j] = min(v(ok)); i = ok(j);
end

function ld = logdet_(S)
% log(det(S)) via Cholesky; det() underflows to 0 for near-singular Sigma,
% which would make log det -Inf and hand every high order a spurious win.
S = (S + S.')/2;
[L, flag] = chol(S);
if flag == 0
    ld = 2*sum(log(diag(L)));
else
    ev = eig(S); ev = ev(ev > 0);
    if isempty(ev), ld = NaN; else, ld = sum(log(ev)); end
end
end
