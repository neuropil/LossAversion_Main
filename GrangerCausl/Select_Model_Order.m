function [pSel,info] = Select_Model_Order(data,maxOrder,criterion)
%% [pSel,info] = Select_Model_Order(data,maxOrder,criterion)
% Selects VAR model order by AIC or BIC on the FULL multivariate data.
%   data      : nChan x nSamp
%   maxOrder  : largest order to test (default 30)
%   criterion : 'AIC' or 'BIC' (default 'BIC')
% Returns the selected order and an info struct with the full curves.
%
% Recommended: run this ONCE on representative/concatenated data, then pass
% the chosen order as a FIXED order to every epoch and condition, so
% finite-sample GC bias is matched across the comparison.
if nargin < 3 || isempty(criterion), criterion = 'BIC'; end
if nargin < 2 || isempty(maxOrder),  maxOrder  = 30;    end
nChan  = size(data,1);
orders = 1:maxOrder;
AIC = zeros(numel(orders),1);
BIC = zeros(numel(orders),1);
for oi = 1:numel(orders)
    p = orders(oi);
    [R,~,~] = Autoregressive_Process_V1(data,p);
    m     = size(R,2);              % effective sample size
    Sigma = (R*R')/m;              % residual covariance
    k     = nChan*nChan*p;          % number of AR parameters
    ld    = log(det(Sigma));
    AIC(oi) = ld + 2*k/m;
    BIC(oi) = ld + k*log(m)/m;
end
[~,ia] = min(AIC);
[~,ib] = min(BIC);
info.orders   = orders;
info.AIC      = AIC;
info.BIC      = BIC;
info.aicOrder = orders(ia);
info.bicOrder = orders(ib);
if strcmpi(criterion,'AIC'), pSel = orders(ia); else, pSel = orders(ib); end
end
