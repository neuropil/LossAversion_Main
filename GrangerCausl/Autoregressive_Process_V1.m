function [R,err,rc] = Autoregressive_Process_V1(data,history)
%% [R,err,rc] = Autoregressive_Process_V1(data,history)
% Fits a VAR(history) to DATA (nChan x nSamp) and returns the residuals.
%   R   : residual matrix (nChan x (nSamp-history))
%   err : rms of residuals (toolbox-free; identical to rms(R(:)))
%   rc  : rcond of the lag covariance C_ = X_*X_'  (conditioning check)
% NOTE: 'history' is the VAR model ORDER (number of lags), not a duration.
[nChan,nSamp] = size(data);
nX2 = nSamp - history;
X_ = zeros(nChan*history,nX2);
k = 1;
for i=1:nChan
    for j=1:history
        X_(k,:) = data(i,j-1+(1:nX2));
        k = k + 1;
    end
end
X  = data(:,history+1:end);
C  = X*X_';
C_ = X_*X_';
rc = rcond(C_);            % <-- previously commented out; now returned
A  = C / C_;
R  = X - A*X_;
err = sqrt(mean(R(:).^2)); % <-- was rms(R(:)); identical, no Signal toolbox dep
end
