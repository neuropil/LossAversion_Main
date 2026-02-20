

%% 1. Simulate Data
fs = 500; t = 0:1/fs:10; n = length(t);
f_fft = (0:n-1)*(fs/n); exponent = 1.1;
aperiodic_spectrum = 1 ./ ((f_fft + 1).^exponent);
rand_phases = exp(1i*2*pi*rand(1,n));
spectrum = aperiodic_spectrum .* rand_phases;
background = real(ifft(spectrum)); background = background/std(background);

f_alpha = 10; f_beta = 22;
signal = background + 0.5*sin(2*pi*f_alpha*t) + 0.25*sin(2*pi*f_beta*t);
signal = signal + 0.2*randn(size(t));

%% 2. Compute PSD
[pxx, f] = pwelch(signal, fs*2, fs, [], fs);
logP = log10(pxx); logF = log10(f);

%% 3. Iterative Peak Exclusion (linear log-log fit)
freq_mask = (f >= 2) & (f <= 50);
maxIter = 5; exclude = true(size(f));
for i = 1:maxIter
    fit_idx = freq_mask & exclude;
    p = polyfit(logF(fit_idx), logP(fit_idx), 1);
    aperiodic_fit = polyval(p, logF);
    residuals = logP - aperiodic_fit;
    exclude = ~(residuals > 2*std(residuals(fit_idx)));
end
final_idx = freq_mask & exclude;

%% 4. Nonlinear 'Knee' Fit
fooof_fun = @(b, x) b(1) - log10(b(2) + x.^b(3)); % [offset, knee, exponent]
b0 = [mean(logP(final_idx)), 1, 1];
b_fit = lsqcurvefit(fooof_fun, b0, f(final_idx), logP(final_idx), [0 0 0], [inf inf 5]);
fooof_fit = fooof_fun(b_fit, f);

% Aperiodic parameters
offset = b_fit(1);
knee   = b_fit(2);
exponent = b_fit(3);

%% 5. Residuals (Periodic Component)
residuals = logP - fooof_fit;

%% 6. Periodic Peaks (using 6–40 Hz)
search_range = (f >= 6) & (f <= 40);
[peak_vals, locs, widths] = findpeaks(residuals(search_range), f(search_range), ...
    'MinPeakProminence', 0.05, 'MinPeakDistance', 2);

%% 7. Tidy Output Tables

% Power stats for raw, fit, and residuals
band_idx = search_range; % 6–40 Hz analysis band
tbl_power = table;
tbl_power.Band = "6–40 Hz";
tbl_power.MeanRaw   = mean(logP(band_idx));
tbl_power.SDRaw     = std(logP(band_idx));
tbl_power.MeanFit   = mean(fooof_fit(band_idx));
tbl_power.SDFit     = std(fooof_fit(band_idx));
tbl_power.MeanResid = mean(residuals(band_idx));
tbl_power.SDResid   = std(residuals(band_idx));

% Aperiodic parameters
tbl_aperiodic = table;
tbl_aperiodic.Offset  = offset;
tbl_aperiodic.Knee    = knee;
tbl_aperiodic.Exponent = exponent;

% Periodic (peak) summary
tbl_periodic = table;
tbl_periodic.Frequency_Hz = locs';
tbl_periodic.Peak_Amplitude = peak_vals';
tbl_periodic.Bandwidth_Hz = widths';

% Overall summary
disp('--- Aperiodic Fit Parameters ---');
disp(tbl_aperiodic);
disp(' ');
disp('--- Power Summary (6–40 Hz) ---');
disp(tbl_power);
disp(' ');
disp('--- Periodic Peaks ---');
disp(tbl_periodic);

%% Optional: Plot
figure;
subplot(3,1,1);
plot(f, logP, 'k'); hold on;
plot(f, fooof_fit, 'r-');
xlim([2 50]); ylabel('log_{10} Power');
legend('PSD','Aperiodic (knee)');
title('Power Spectrum Decomposition');

subplot(3,1,2);
plot(f, residuals, 'k'); hold on;
findpeaks(residuals(search_range), f(search_range));
xlim([2 50]); ylabel('Residuals');
title('Residuals (Oscillatory Peaks)');

subplot(3,1,3);
plot(f, pxx, 'k');
xlim([2 50]); xlabel('Frequency (Hz)'); ylabel('Power');
title('Raw Power Spectrum');
