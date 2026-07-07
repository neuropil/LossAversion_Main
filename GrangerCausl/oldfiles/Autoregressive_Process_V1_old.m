function [R,err] = Autoregressive_Process_V1(data,history)

%% function [R,err] = Autoregressive_Process_V1(data,history)
% Notes for later.

%%
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

X = data(:,history+1:end);

C = X*X_';
C_ = X_*X_';

A = C / C_;
% zz = rcond(C_);

R = X - A*X_;
err = rms(R(:));

end