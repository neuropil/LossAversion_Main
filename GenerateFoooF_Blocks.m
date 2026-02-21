function [] = GenerateFoooF_Blocks()

cd('Z:\LossAversion\LH_Data\JAT_BlockData')

matD1 = dir('*.mat');
matD2 = struct2table(matD1);
matD3 = matD2.name;

for mmii = 1:length(matD3)

    tmpFile = matD3{mmii};
    pBlklist = {'PreBlock1','PreBlock2','PreBlock3','PreBlock4','PreBlock5'};
    load(tmpFile,'baselineBlock')

    startsS = 1:500:2500;
    stopsS = [startsS(2:end)-1,2500];

    for pii = 1:length(pBlklist)
        tmpBlk = baselineBlock.(pBlklist{pii}).LFP;
        fooofblock = cell(5,1);
        specPblock = cell(5,1);
        for secI = 1:5

            tmpSec = tmpBlk(:,startsS(secI):stopsS(secI));

            [fooof_results] = testPSD_fooof(tmpSec);
            [specParm_results] = jatSpecPARAM(tmpSec);

            fooofblock{secI} = fooof_results;
            specPblock{secI} = specParm_results;

        end
        preBlockFooof.FoooF.(['B_',num2str(pii)]) = fooofblock;
        preBlockFooof.JATSpecP.(['B_',num2str(pii)]) = specPblock;
    end
    % SAVE file append
    save(tmpFile,"preBlockFooof",'-append');
end









end















function [fooof_results] = testPSD_fooof(inputLFP)

% inputLFP = detrend(inputLFP,'linear');   % remove slow trends
% inputLFP = inputLFP - mean(inputLFP);           % remove DC offset
% hpFilt = designfilt('highpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency',1,'SampleRate',500);
% inputLFP = filtfilt(hpFilt,inputLFP);
% [pxx,f] = pmtm(inputLFP,2.5,1024,500);
% plot(f,log10(pxx))
% xlim([1 100])

% Inputs you need:
%   x        : 1-second signal vector
%   Fs       : 500
%   f_range  : e.g. [3 45] for FOOOF
%   settings : MATLAB struct (your settings)

Fs = 500;
f_range = [3 55]; % THIS IS NEW ----  WAS [3 45]

allPSD = zeros(1001,height(inputLFP));
for ii = 1:height(inputLFP)

    % --- preprocess (matches what improved your PSD) ---
    inputLFPi = detrend(inputLFP(ii,:),'linear');
    inputLFPi = inputLFPi - mean(inputLFPi);
    hp = designfilt('highpassiir','FilterOrder',4,...
        'HalfPowerFrequency',2,'SampleRate',Fs);
    inputLFPi = filtfilt(hp,inputLFPi);

    % --- PSD (multitaper recommended for 1-s) ---
    % nfft = 1024;
    % [psd, freqs] = pmtm(inputLFP, 3, nfft, Fs);

    [allPSD(:,ii), freqs] = pspectrum(inputLFPi,Fs,"power","FrequencyResolution",2,"Leakage",0.8,'FrequencyLimits',[2 55]);

    % --- restrict to f_range for plotting ---
    % idx = freqs >= f_range(1) & freqs <= f_range(2);
    % fplt = freqs(idx);
    % pplt = allPSD(idx);

end

meanPSD = mean(allPSD,2);

settings = struct( ...
    'peak_width_limits', [2, 10], ...   % narrower; wide peaks in 1 s = noise
    'max_n_peaks', 5, ...               % cap aggressively
    'min_peak_height', 0.15, ...        % raise floor
    'peak_threshold', 2, ...          % stricter than default 2
    'aperiodic_mode', 'fixed', ...      % knee not reliable here
    'verbose', true);

% --- run FOOOF (your wrapper expects full vectors + range) ---
fooof_results = fooof(freqs, meanPSD, f_range, settings, true);

% % --- overlay plot (log10 power is standard) ---
% figure; hold on;
% plot(fplt, log10(pplt), 'k', 'LineWidth', 1.5);                     % raw PSD
% plot(fooof_results.freqs, log10(fooof_results.ap_fit), 'b', 'LineWidth', 1.5);   % aperiodic fit
% plot(fooof_results.freqs, log10(fooof_results.fooofed_spectrum), 'r', 'LineWidth', 1.5); % full model
% plot(fooof_results.freqs, log10(fooof_results.power_spectrum - fooof_results.ap_fit), ...
%      'Color', [0.2 0.6 0.2], 'LineWidth', 1);                       % residual-ish (careful)
% 
% xlabel('Frequency (Hz)');
% ylabel('log_{10} Power');
% title('FOOOF fit overlay');
% legend({'PSD','Aperiodic fit','FOOOFed fit','(PSD - AP fit)'}, 'Location','northeast');
% xlim([f_range(1) f_range(2)]);
% grid on;

% --- print summary + peaks ---
disp('--- SETTINGS ---');
disp(settings);

disp('--- FOOOF RESULTS ---');
fprintf('Exponent: %.3f\n', fooof_results.aperiodic_params(2));
fprintf('Offset  : %.3f\n', fooof_results.aperiodic_params(1));
fprintf('R^2     : %.3f\n', fooof_results.r_squared);
fprintf('Error   : %.3f\n', fooof_results.error);

if isfield(fooof_results,'peak_params') && ~isempty(fooof_results.peak_params)
    peaks = fooof_results.peak_params;
    % columns: [CF, PW, BW] for fooof typically
    fprintf('\n--- PEAKS [CF  Amp  BW] ---\n');
    disp(peaks);
else
    disp('No peaks detected.');
end


end







function [specParm_results] = jatSpecPARAM(inputLFP)

Fs = 500;

allPSD = zeros(1001,height(inputLFP));
for ii = 1:height(inputLFP)

    % --- preprocess (matches what improved your PSD) ---
    inputLFPi = detrend(inputLFP(ii,:),'linear');
    inputLFPi = inputLFPi - mean(inputLFPi);
    hp = designfilt('highpassiir','FilterOrder',4,...
        'HalfPowerFrequency',2,'SampleRate',Fs);
    inputLFPi = filtfilt(hp,inputLFPi);

    [allPSD(:,ii), fxx] = pspectrum(inputLFPi,Fs,"power","FrequencyResolution",2,"Leakage",0.8,'FrequencyLimits',[2 55]);

end

Pxx = mean(allPSD,2);

% [pxx, f] = pwelch(signal, fs*2, fs, [], fs);
logP = log10(Pxx);
logF = log10(fxx);

% freq_mask = (fxx >= 2) & (fxx <= 90);
% maxIter = 5;
% exclude = true(size(fxx));
% for i = 1:maxIter
%     fit_idx = freq_mask & exclude;
%     p = polyfit(logF(fit_idx), logP(fit_idx), 1);
%     aperiodic_fit = polyval(p, logF);
%     residuals = logP - aperiodic_fit;
%     exclude = ~(residuals > 2*std(residuals(fit_idx)));
% end
% final_idx = freq_mask & exclude;

% ----- basic masks -----
valid = (fxx > 0) & isfinite(logP) & isfinite(logF);
band  = valid & (fxx >= 3) & (fxx <= 55);     % match your pspectrum limits

% =========================
% FIXED: keep your robust exclusion (optional)
% =========================
maxIter = 5;
exclude = true(size(fxx));
for i = 1:maxIter
    fit_idx = band & exclude;
    p = polyfit(logF(fit_idx), logP(fit_idx), 1);
    base_fit = polyval(p, logF);
    resid = logP - base_fit;
    exclude = ~(resid > 2*std(resid(fit_idx)));
end
final_idx_fixed = band & exclude;

% =========================
% KNEE: use a knee-appropriate mask (do NOT reuse fixed exclusion)
% =========================
final_idx_knee = band;   % start simple; you can add knee-robustness later

% % --- models ---
% % knee_fun  = @(b,f) b(1) - log10(b(2) + f.^b(3));   % b = [offset, knee, exponent]
% knee_fun_stable = @(b,f) b(1) - log10( 10.^b(2) + f.^b(3) );
% fixed_fun = @(b,f) b(1) - b(2).*log10(f);          % b = [offset, exponent]

% ----- models -----
knee_fun = @(b,f) b(1) - log10( 10.^b(2) + f.^b(3) );   % b=[offset, log10(knee), exponent]
fixed_fun = @(b,f) b(1) - b(2).*log10(f);              % b=[offset, exponent]

% ----- initial guesses -----
% offset guess: median power in band (more robust than mean)
off0 = median(logP(final_idx_knee));

% exponent guess: from a quick fixed slope over higher freqs (where knee matters less)
hi = valid & (fxx >= 20) & (fxx <= 55);
p_hi = polyfit(logF(hi), logP(hi), 1);
exp0 = max(0, min(5, -p_hi(1)));

% knee guess: pick a knee frequency fk ~ 5–10 Hz and back out k ~= fk^exp
fk = 8;  % Hz
log10k0 = log10(fk^exp0);

b0_knee  = [off0, log10k0, exp0];
b0_fixed = [median(logP(final_idx_fixed)), max(0,min(5,-p(1)))];  % p from last fixed iter

% --- initial guesses ---
% b0_knee  = [mean(logP(final_idx)), 1, 1];
% b0_knee = [mean(logP(fit_mask)), 0, 1];     % log10(knee)=0 -> knee=1
% b0_fixed = [mean(logP(final_idx)), 1];

% --- bounds ---
% lb_knee = [-inf, 0, 0];
% ub_knee = [ inf, inf, 5];
lb_knee = [-inf, -12, 0];
ub_knee = [ inf,  12, 5];

lb_fixed = [-inf, 0];
ub_fixed = [ inf, 5];

% --- fits ---
% b_fit_knee  = lsqcurvefit(knee_fun,  b0_knee,  fxx(final_idx), logP(final_idx), lb_knee,  ub_knee);
% b_fit_knee = lsqcurvefit(knee_fun_stable, b0_knee, fxx(fit_mask), logP(fit_mask), lb_knee, ub_knee);
% b_fit_fixed = lsqcurvefit(fixed_fun, b0_fixed, fxx(final_idx), logP(final_idx), lb_fixed, ub_fixed);

% ----- fits -----
b_fit_knee  = lsqcurvefit(knee_fun,  b0_knee,  fxx(final_idx_knee),  logP(final_idx_knee),  lb_knee,  ub_knee);
b_fit_fixed = lsqcurvefit(fixed_fun, b0_fixed, fxx(final_idx_fixed), logP(final_idx_fixed), lb_fixed, ub_fixed);

% --- full-band fits for plotting/residuals ---
% fit_knee_full = knee_fun_stable(b_fit_knee, fxx);
% fit_fixed_full = fixed_fun(b_fit_fixed, fxx);
% ----- full-band fits -----
fit_knee_full  = knee_fun(b_fit_knee,  fxx);
fit_fixed_full = fixed_fun(b_fit_fixed, fxx);

residuals_knee  = logP - fit_knee_full;
residuals_fixed = logP - fit_fixed_full;

% --- parameters ---
kneeParams.offset   = b_fit_knee(1);
kneeParams.knee     = 10.^b_fit_knee(2);
kneeParams.exponent = b_fit_knee(3);

fixedParams.offset   = b_fit_fixed(1);
fixedParams.exponent = b_fit_fixed(2);

% ----- fit quality (use model-specific masks) -----
y = logP(final_idx_fixed);
yhat = fit_fixed_full(final_idx_fixed);
res = y - yhat;
fixedParams.Error = sqrt(mean(res.^2));
fixedParams.Rsq   = 1 - sum(res.^2) / sum((y - mean(y)).^2);

y = logP(final_idx_knee);
yhat = fit_knee_full(final_idx_knee);
res = y - yhat;
kneeParams.Error = sqrt(mean(res.^2));
kneeParams.Rsq   = 1 - sum(res.^2) / sum((y - mean(y)).^2);

% ------ Find peaks ---------------------------------

search_range = (fxx >= 6) & (fxx <= 95);
% [peak_vals, locs, widths] = findpeaks(residuals_knee(search_range), fxx(search_range), ...
%     'MinPeakProminence', 0.1, 'MinPeakDistance', 4)

knee_x = fxx(search_range);
knee_y = residuals_knee(search_range);

% Robust noise estimate of residuals (avoid hard-coding 0.1)
knee_sigma = 1.4826 * mad(knee_y, 1);          % ~robust std
% minProm = max(0.1, 3*sigma);         % keep your 0.1 as a floor, or drop it
% knee_minProm = max(0.1, 3*knee_sigma)/2;
knee_minProm = min(0.15, max(0.08, 2*knee_sigma));

[knee_peak_vals, knee_locs, knee_widths, knee_proms] = findpeaks(knee_y, knee_x, ...
    'MinPeakProminence', knee_minProm, ...
    'MinPeakHeight',     0, ...
    'MinPeakDistance',   2, ...        % Hz (see note below)
    'WidthReference',    'halfheight', ...
    'MinPeakWidth',      2, ...        % Hz (match FOOOF lower BW limit)
    'MaxPeakWidth',      10);          % Hz (match FOOOF upper BW limit)

kpeaks.params.MinPeakProm = knee_minProm;
kpeaks.params.MinPeakHeight = 0;
kpeaks.params.MinPeakDistance = 2;
kpeaks.params.MinPeakWidth = 2;
kpeaks.params.MaxPeakWidth = 10;

kAmpBC = knee_proms;
kHz = knee_locs;
kBW = knee_widths;
kAmpRaw = knee_peak_vals;
kpeaks.Summary = table(kAmpBC,kHz,kBW,kAmpRaw,'VariableNames',...
    {'AmpBC','Hz','BW','AmpRAW'});

max_n_peaks = 5;
% ------- BUILD Fooofed_spectrum ------------------------------------------
if ~isempty(knee_locs)
    [~, ord] = sort(knee_proms, 'descend');           % rank by prominence
    ord = ord(1:min(max_n_peaks, numel(ord)));
    knee_locs   = knee_locs(ord);
    knee_widths = knee_widths(ord);
    knee_proms = knee_proms(ord);
end
% --- build Gaussian peak fit on the FULL frequency vector (same space as logP) ---
f = fxx(:);
peak_fit_full_knee = zeros(size(f));

for knumpks = 1:numel(knee_locs)
    sig = knee_widths(knumpks) / 2.3548;                     % FWHM -> sigma
    A   = knee_proms(knumpks);                              % use prominence as "Amp"
    peak_fit_full_knee = peak_fit_full_knee + A .* exp(-(f - knee_locs(knumpks)).^2 ./ (2*sig^2));
end

% --- reconstructed model spectrum (FOOOF-like) ---
fooofed_spectrum_knee = fit_knee_full(:) + peak_fit_full_knee;

% Optional: store peak table comparable to FOOOF [CF Amp BW]
% peaks_table_knee = [cf(:), prom(:), fwhm(:)];

fixed_x = fxx(search_range);
fixed_y = residuals_fixed(search_range);

% Robust noise estimate of residuals (avoid hard-coding 0.1)
fixed_sigma = 1.4826 * mad(fixed_y, 1);          % ~robust std
% minProm = max(0.1, 3*sigma);         % keep your 0.1 as a floor, or drop it
fixed_minProm = min(0.15, max(0.08, 2*fixed_sigma));

[fixed_peak_vals, fixed_locs, fixed_widths, fixed_proms] = findpeaks(fixed_y, fixed_x, ...
    'MinPeakProminence', fixed_minProm, ...
    'MinPeakHeight',     0, ...
    'MinPeakDistance',   2, ...        % Hz (see note below)
    'WidthReference',    'halfheight', ...
    'MinPeakWidth',      2, ...        % Hz (match FOOOF lower BW limit)
    'MaxPeakWidth',      10);          % Hz (match FOOOF upper BW limit)

fpeaks.params.MinPeakProm = fixed_minProm;
fpeaks.params.MinPeakHeight = 0;
fpeaks.params.MinPeakDistance = 2;
fpeaks.params.MinPeakWidth = 2;
fpeaks.params.MaxPeakWidth = 10;

fAmpBC = fixed_proms;
fHz = fixed_locs;
fBW = fixed_widths;
fAmpRaw = fixed_peak_vals;
fpeaks.Summary = table(fAmpBC,fHz,fBW,fAmpRaw,'VariableNames',...
    {'AmpBC','Hz','BW','AmpRAW'});

% ------- BUILD Fooofed_spectrum ------------------------------------------
if ~isempty(fixed_locs)
    [~, ord] = sort(fixed_proms, 'descend');           % rank by prominence
    ord = ord(1:min(max_n_peaks, numel(ord)));
    fixed_locs   = fixed_locs(ord);
    fixed_widths = fixed_widths(ord);
    fixed_proms = fixed_proms(ord);
end
% --- build Gaussian peak fit on the FULL frequency vector (same space as logP) ---
f = fxx(:);
peak_fit_full_fixed = zeros(size(f));

for fnumpks = 1:numel(fixed_locs)
    sig = fixed_widths(fnumpks) / 2.3548;                     % FWHM -> sigma
    A   = fixed_proms(fnumpks);                              % use prominence as "Amp"
    peak_fit_full_fixed = peak_fit_full_fixed + A .* exp(-(f - fixed_locs(fnumpks)).^2 ./ (2*sig^2));
end

% --- reconstructed model spectrum (FOOOF-like) ---
fooofed_spectrum_fixed = fit_fixed_full(:) + peak_fit_full_fixed;

% ------ Outputs ------------------------------------

specParm_results.FIXED.residuals = residuals_fixed;
specParm_results.FIXED.PeakSpec = fooofed_spectrum_fixed;
specParm_results.FIXED.params = fixedParams;
specParm_results.KNEE.residuals = residuals_knee;
specParm_results.KNEE.PeakSpec = fooofed_spectrum_knee;
specParm_results.KNEE.params = kneeParams;
specParm_results.FIXED.fit = fit_fixed_full;
specParm_results.KNEE.fit = fit_knee_full;
specParm_results.General.PSD = logP;
specParm_results.General.Freqs = fxx;
specParm_results.FIXED.peaks = fpeaks;
specParm_results.KNEE.peaks = kpeaks;

end