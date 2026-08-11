function test_GC_Extensions()
% TEST_GC_EXTENSIONS  Ground-truth checks for Compute_Spectral_GC and
%                     Extract_Band_Power. Run this before touching real data.
%
%   test_GC_Extensions()
%
% Replicates the two validation cases used during development. Expected:
%   TEST 1  spectral GC peaks at 40 Hz, correct direction, Geweke identity holds
%   TEST 2  raw-voltage GC ~0 both ways; envelope GC recovers the true direction

fprintf('\n================ TEST 1: spectral GC ================\n');
fprintf('Truth: node2 -> node1, node2 resonant at 40 Hz.\n');
rand('seed',7); randn('seed',7);

fs = 500; N = 1501; T = 20;
f0 = 40; r = 0.96;
a1 = 2*r*cos(2*pi*f0/fs); a2 = -r^2;

X = zeros(2,N,T);
for t = 1:T
    e = randn(2, N+200);
    x1 = zeros(1,N+200); x2 = zeros(1,N+200);
    for k = 4:N+200
        x2(k) = a1*x2(k-1) + a2*x2(k-2) + e(2,k);
        x1(k) = 0.30*x1(k-1) + 0.60*x2(k-2) + e(1,k);   % 2 -> 1
    end
    X(1,:,t) = x1(201:end);
    X(2,:,t) = x2(201:end);
end

S = Compute_Spectral_GC(X, 8, fs, 'nFreq', 257, 'verbose', true);

fprintf('\n  time-domain  2->1 = %.4f   <- should be LARGE\n', S.td(1,2));
fprintf('  time-domain  1->2 = %.4f   <- should be ~0\n',    S.td(2,1));
[pk, pi_] = max(squeeze(S.spec(1,2,:)));
fprintf('  spectral 2->1 peaks at %.1f Hz (truth 40), value %.3f\n', S.f(pi_), pk);
fprintf('  spectral 1->2 max     = %.4f\n', max(squeeze(S.spec(2,1,:))));
fprintf('\n  Geweke identity (mean over f vs time-domain):\n');
fprintf('    2->1  spec-mean %.4f  vs  td %.4f\n', S.specCheck(1,2), S.td(1,2));
fprintf('    1->2  spec-mean %.4f  vs  td %.4f\n', S.specCheck(2,1), S.td(2,1));

bn = fieldnames(S.band);
fprintf('\n  band averages:\n');
for b = 1:numel(bn)
    fprintf('    %-10s 2->1 = %7.4f    1->2 = %7.4f\n', ...
        bn{b}, S.band.(bn{b})(1,2), S.band.(bn{b})(2,1));
end

ok1 = S.td(1,2) > 10*max(S.td(2,1),1e-6) && abs(S.f(pi_)-40) < 5 && ...
      abs(S.specCheck(1,2)-S.td(1,2)) < 0.05;
fprintf('\n  TEST 1: %s\n', tf(ok1));

fprintf('\n================ TEST 2: envelope GC ================\n');
fprintf(['Truth: envelope of node2 drives envelope of node1 at 40 ms lag,\n' ...
         'with INDEPENDENT high-gamma carriers (no voltage coupling at all).\n']);
rand('seed',11); randn('seed',11);

T2 = 30; band = [70 150]; LAG = 20;
X2 = zeros(2,N,T2);
for t = 1:T2
    e = randn(2, N+400);
    s2 = zeros(1,N+400); s1 = zeros(1,N+400);
    for k = LAG+2 : N+400
        s2(k) = 0.995*s2(k-1) + 0.06*e(2,k);
        s1(k) = 0.995*s1(k-1) + 0.06*e(1,k) + 0.45*s2(k-LAG);
    end
    c  = randn(2,N);
    cb = Extract_Band_Power(c, fs, band, 'output','filtered', ...
                            'downsample',false, 'trimCycles',0, 'verbose',false);
    X2(1,:,t) = exp(s1(401:end)) .* cb(1,:);
    X2(2,:,t) = exp(s2(401:end)) .* cb(2,:);
end

Sraw = Compute_Spectral_GC(X2, 10, fs, 'nFreq', 65, 'verbose', false);
fprintf('\n  RAW VOLTAGE GC (expect ~0 both ways):\n');
fprintf('    2->1 = %.4f     1->2 = %.4f\n', Sraw.td(1,2), Sraw.td(2,1));

fprintf('\n  BAND-POWER TRANSFORM:\n');
[Y, fs2, info] = Extract_Band_Power(X2, fs, band, 'output','logpower');

fprintf('\n  ENVELOPE GC (true lag = %.1f samples after decimation):\n', LAG/info.step);
for p = [3 5 8 12]
    Se = Compute_Spectral_GC(Y, p, fs2, 'nFreq', 65, 'verbose', false);
    fprintf('    order %2d: 2->1 = %.4f   1->2 = %.4f   ratio %5.1fx   rcond %.1e\n', ...
        p, Se.td(1,2), Se.td(2,1), Se.td(1,2)/max(Se.td(2,1),1e-9), Se.rcondMin);
end

Se = Compute_Spectral_GC(Y, 8, fs2, 'nFreq', 65, 'verbose', false);
ok2 = abs(Sraw.td(1,2)) < 0.01 && abs(Sraw.td(2,1)) < 0.01 && ...
      Se.td(1,2) > 3*Se.td(2,1);
fprintf('\n  TEST 2: %s\n', tf(ok2));

fprintf('\n  Sample budget per band (3 s epochs at %d Hz):\n', fs);
for b = {[4 8],[8 12],[13 30],[30 70],[70 150]}
    [~, fo, inf2] = Extract_Band_Power(X2(:,:,1), fs, b{1}, 'verbose', false);
    fprintf('    %3d-%3d Hz -> fs %6.1f Hz, %4d samp/trial (trim %d/side)\n', ...
        b{1}(1), b{1}(2), fo, inf2.nSamp, inf2.nTrim);
end

fprintf('\n================ OVERALL: %s ================\n\n', tf(ok1 && ok2));
end

function s = tf(x)
if x, s = 'PASS'; else, s = 'FAIL'; end
end
