function [] = SPRiNT_SpecParmJ_V1(signal,fs ,plotUSE)



T       = length(signal)/fs;    % Duration (seconds)
N       = T * fs;               % Number of samples
t       = (0:N-1) / fs; % Time vector


%% =========================================================
%  SECTION 2: SPRINT PARAMETERS
% =========================================================

sprint.win_len   = 1.5;       % Window length (seconds) % was 2
sprint.win_step  = 0.25;     % Step size between windows (seconds) % was 0.5
sprint.freq_range = [3 40]; % Frequency range for fitting (Hz)
sprint.min_peak_height = 0.15;  % Minimum peak height (log power)
sprint.peak_threshold  = 2.0;  % SD threshold above aperiodic for peak detection
sprint.max_n_peaks     = 7;    % Maximum number of peaks to fit
sprint.peak_width_limits = [1 6]; % [min max] peak width (Hz, BW = 2*sigma*sqrt(2))
sprint.aperiodic_mode  = 'knee'; % 'fixed' (no knee) or 'knee'

% Derived
sprint.win_samps  = round(sprint.win_len  * fs);
sprint.step_samps = round(sprint.win_step * fs);

% Window centers
% Pad signal symmetrically so windows cover the full time range
pad_samps = round(sprint.win_samps / 2);
signal_padded = [zeros(1, pad_samps), signal, zeros(1, pad_samps)];

% Shift window start indices to account for padding
win_starts = 1 : sprint.step_samps : (length(signal_padded) - sprint.win_samps + 1);
win_centers = ((win_starts + sprint.win_samps/2 - 1) - pad_samps) / fs;  % centres map back to original time
n_wins = length(win_starts);
signal_padded = [fliplr(signal(1:pad_samps)), signal, fliplr(signal(end-pad_samps+1:end))];

fprintf('=== SPRiNT Configuration ===\n');
fprintf('  Window: %.1f s, Step: %.1f s\n', sprint.win_len, sprint.win_step);
fprintf('  Freq range: %.0f–%.0f Hz\n', sprint.freq_range(1), sprint.freq_range(2));
fprintf('  Number of windows: %d\n\n', n_wins);


%% =========================================================
%  SECTION 3: COMPUTE SHORT-TIME POWER SPECTRA (STFT-based)
% =========================================================

fprintf('=== Computing short-time power spectra ===\n');

% Use Hann window + zero-padding for frequency resolution
nfft   = 2^nextpow2(sprint.win_samps * 4); % zero-pad 4x for resolution
hann_w = hann(sprint.win_samps)';

% Frequency vector
freqs_all = (0 : nfft/2) * fs / nfft;
freq_mask = freqs_all >= sprint.freq_range(1) & freqs_all <= sprint.freq_range(2);
freqs     = freqs_all(freq_mask);
n_freqs   = length(freqs);



% Pre-allocate spectrogram
psd_matrix = zeros(n_freqs, n_wins);

% for wi = 1:n_wins
%     idx = win_starts(wi) : win_starts(wi) + sprint.win_samps - 1;
%     seg = signal(idx) .* hann_w;
%     X   = fft(seg, nfft);
%     psd = (2 / (fs * sum(hann_w.^2))) * abs(X(1:nfft/2+1)).^2;
%     psd(2:end-1) = psd(2:end-1); % one-sided
%     psd_matrix(:, wi) = psd(freq_mask);
% end

% Use Welch's method within each window for smoother PSDs
for wi = 1:n_wins
    idx = win_starts(wi) : win_starts(wi) + sprint.win_samps - 1;
    seg = signal_padded(idx);

    % Welch: 50% overlapping sub-windows within the sliding window
    sub_win_len = round(sprint.win_samps / 2);  % sub-window = half of main window
    [pxx, f_welch] = pwelch(seg, hann(sub_win_len), round(sub_win_len*0.5), nfft, fs);

    psd_matrix(:, wi) = pxx(freq_mask);
end

fprintf('  Frequency resolution: %.3f Hz\n', freqs(2)-freqs(1));
fprintf('  Spectrogram size: %d freqs x %d windows\n\n', n_freqs, n_wins);

%% =========================================================
%  SECTION 4: FIT FOOOF MODEL TO EACH TIME WINDOW
% =========================================================

fprintf('=== Fitting FOOOF model to each window ===\n');

% Pre-allocate results
results.time        = win_centers;
results.offset      = NaN(1, n_wins);
results.exponent    = NaN(1, n_wins);
results.r_squared   = NaN(1, n_wins);
results.error       = NaN(1, n_wins);
results.peaks       = cell(1, n_wins);  % each cell: [cf, amplitude, bw]
results.knee = NaN(1, n_wins);  % initialise before the loop

log_freqs = log10(freqs);

for wi = 1:n_wins
    psd      = psd_matrix(:, wi)';
    log_psd  = log10(psd);

    % --- Step 1: Fit aperiodic component (initial) ---
    % ap_params = fit_aperiodic(log_freqs, log_psd, sprint.aperiodic_mode);

    % Exclude canonical oscillatory bands from initial aperiodic fit
    excl = (freqs >= 7 & freqs <= 14) | (freqs >= 15 & freqs <= 30);
    incl = ~excl;
    ap_params = fit_aperiodic(log_freqs(incl), log_psd(incl), sprint.aperiodic_mode);

    % --- Step 2: Iteratively find and remove peaks ---
    peaks_found = zeros(0, 3); % [cf, amplitude, bw]
    log_psd_flat = log_psd - aperiodic_model(log_freqs, ap_params, sprint.aperiodic_mode);

    for pk = 1:sprint.max_n_peaks
        % Find highest point above threshold
        [max_val, max_idx] = max(log_psd_flat);

        if max_val < sprint.min_peak_height
            break;
        end

        if freqs(max_idx) < sprint.freq_range(1)
            log_psd_flat(max_idx) = 0;  % suppress and continue
            continue;
        end

        % Check if above threshold relative to local noise
        % local_std = std(log_psd_flat);
        % if max_val < sprint.peak_threshold * local_std
        %     break;
        % end

        % Estimate initial peak params: [cf, amp, bw_sigma]
        cf_init = freqs(max_idx);
        amp_init = max_val;
        bw_init  = 2.0; % initial sigma guess (Hz)

        % Fit Gaussian to the peak
        try
            peak_p = fit_gaussian_peak(freqs, log_psd_flat, cf_init, amp_init, bw_init, sprint);
            if ~isempty(peak_p)
                peaks_found(end+1, :) = peak_p; %#ok<AGROW>
                % Subtract fitted peak from flattened spectrum
                log_psd_flat = log_psd_flat - gaussian_peak(freqs, peak_p);
            end
        catch
            % Skip bad fits
        end

        % Refit aperiodic with peaks removed
        log_psd_no_peaks = log_psd - sum(gaussian_peak_matrix(freqs, peaks_found), 1);
        ap_params = fit_aperiodic(log_freqs, log_psd_no_peaks, sprint.aperiodic_mode);
        log_psd_flat = log_psd - aperiodic_model(log_freqs, ap_params, sprint.aperiodic_mode) ...
            - sum(gaussian_peak_matrix(freqs, peaks_found), 1);
    end

    % --- Step 3: Final refit of aperiodic with all peaks removed ---
    if ~isempty(peaks_found)
        log_psd_no_peaks = log_psd - sum(gaussian_peak_matrix(freqs, peaks_found), 1);
        ap_params = fit_aperiodic(log_freqs, log_psd_no_peaks, sprint.aperiodic_mode);
    end

    % --- Step 4: Compute goodness of fit ---
    log_psd_hat = aperiodic_model(log_freqs, ap_params, sprint.aperiodic_mode) ...
        + sum(gaussian_peak_matrix(freqs, peaks_found), 1);
    ss_res = sum((log_psd - log_psd_hat).^2);
    ss_tot = sum((log_psd - mean(log_psd)).^2);
    r2     = 1 - ss_res / ss_tot;
    rmse   = sqrt(mean((log_psd - log_psd_hat).^2));

    % Store
    results.r_squared(wi) = r2;
    results.error(wi)     = rmse;
    results.peaks{wi}     = peaks_found;

    % Inside the loop, after fitting:
    if strcmp(sprint.aperiodic_mode, 'knee')
        results.offset(wi)   = ap_params(1);
        results.knee(wi)     = ap_params(2);
        results.exponent(wi) = ap_params(3);
    else
        results.offset(wi)   = ap_params(1);
        results.exponent(wi) = ap_params(2);
    end


    if mod(wi, 20) == 0
        fprintf('  Window %d / %d (t=%.1fs): exp=%.2f, R²=%.3f, %d peaks\n', ...
            wi, n_wins, win_centers(wi), ap_params(end), r2, size(peaks_found,1));
    end
end

fprintf('\n  Done. Mean R² = %.3f\n\n', nanmean(results.r_squared));

%% =========================================================
%  SECTION 5: POST-PROCESS PEAKS — BUILD PEAK TIME SERIES
% =========================================================

% Extract peak center frequencies and amplitudes over time
% Cluster peaks by frequency proximity to form "peak tracks"

all_peaks_t  = [];
all_peaks_cf = [];
all_peaks_amp = [];
all_peaks_bw  = [];

for wi = 1:n_wins
    pk = results.peaks{wi};
    if ~isempty(pk)
        for pi = 1:size(pk,1)
            all_peaks_t(end+1)   = win_centers(wi);
            all_peaks_cf(end+1)  = pk(pi, 1);
            all_peaks_amp(end+1) = pk(pi, 2);
            all_peaks_bw(end+1)  = pk(pi, 3);
        end
    end
end

%% =========================================================
%  SECTION 6: VISUALIZATION
% =========================================================

if plotUSE

    fprintf('=== Generating plots ===\n');

    fig = figure('Name','SPRiNT Results','Position',[50 50 1400 900],'Color','w');

    % ---- (a) Raw signal ----
    ax1 = subplot(4,3,[1 2 3]);
    plot(t, signal, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Simulated EEG-like Signal');
    xlim([0 T]);
    grid on; box off;

    % ---- (b) Spectrogram (log power) ----
    ax2 = subplot(4,3,[4 5 6]);
    imagesc(win_centers, freqs, log10(psd_matrix));
    axis xy;
    colormap(ax2, 'jet');
    colorbar;
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('Short-time Power Spectral Density (log_{10})');
    ylim(sprint.freq_range);
    clim_vals = prctile(log10(psd_matrix(:)), [5 95]);
    clim(clim_vals);

    % ---- (c) Aperiodic exponent over time ----
    ax3 = subplot(4,3,7);
    % Ground truth exponent (downsampled to window centers)
    % true_exp_at_wins = interp1(t, true_exp, win_centers);
    % plot(win_centers, true_exp_at_wins, 'k--', 'LineWidth', 2, 'DisplayName', 'Ground truth');
    hold on;
    plot(results.time, results.exponent, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', 'SPRiNT');
    legend('Location','best'); grid on; box off;
    xlabel('Time (s)'); ylabel('\chi (exponent)');
    title('Aperiodic Exponent');
    ylim([0.5 3.5]);

    % ---- (d) Aperiodic offset over time ----
    ax4 = subplot(4,3,8);
    plot(results.time, results.offset, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    % yline(true_offset, 'k--', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Offset (log power)');
    title('Aperiodic Offset');
    grid on; box off;

    % ---- (e) R-squared over time ----
    ax5 = subplot(4,3,9);
    plot(results.time, results.r_squared, 'g-', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('R²');
    title('Goodness of Fit (R²)');
    ylim([0 1]); grid on; box off;

    % ---- (f) Peak scatter: CF vs Time ----
    ax6 = subplot(4,3,[10 11 12]);
    scatter(all_peaks_t, all_peaks_cf, 20, all_peaks_amp, 'filled', 'MarkerFaceAlpha', 0.7);
    colormap(ax6, 'hot');
    cb = colorbar; cb.Label.String = 'Amplitude (log power)';
    xlabel('Time (s)'); ylabel('Peak CF (Hz)');
    title('Detected Oscillatory Peaks Over Time');
    ylim(sprint.freq_range);
    xlim([0 T]);
    % Overlay ground truth windows
    % hold on;
    % patch([15 45 45 15], [8 8 12 12], 'cyan', 'FaceAlpha', 0.15, 'EdgeColor','cyan', 'DisplayName','True alpha');
    % patch([30 60 60 30], [16 16 24 24], 'magenta', 'FaceAlpha', 0.15, 'EdgeColor','magenta', 'DisplayName','True beta');
    % legend('Detected peaks','True alpha window','True beta window','Location','northeast');
    % grid on; box off;

    sgtitle('SPRiNT: Time-Resolved Spectral Parameterization', 'FontSize', 14, 'FontWeight', 'bold');

    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x');

    fprintf('  Plot complete.\n\n');

    %% =========================================================
    %  SECTION 7: EXAMPLE — PLOT A FEW INDIVIDUAL WINDOW FITS
    % =========================================================

    example_wins = round(linspace(5, n_wins-5, 6));

    fig2 = figure('Name','Example Window Fits','Position',[100 100 1200 700],'Color','w');

    for ei = 1:length(example_wins)
        wi = example_wins(ei);
        psd      = psd_matrix(:, wi)';
        log_psd  = log10(psd);

        % Reconstruct fit
        % ap_params = [results.offset(wi), results.exponent(wi)];

        if strcmp(sprint.aperiodic_mode, 'knee')
            ap_params = [results.offset(wi), results.knee(wi), results.exponent(wi)];
        else
            ap_params = [results.offset(wi), results.exponent(wi)];
        end


        ap_fit    = aperiodic_model(log_freqs, ap_params, sprint.aperiodic_mode);
        pk        = results.peaks{wi};
        if ~isempty(pk)
            full_fit = ap_fit + sum(gaussian_peak_matrix(freqs, pk), 1);
        else
            full_fit = ap_fit;
        end

        subplot(2, 3, ei);
        plot(freqs, log_psd, 'k', 'LineWidth', 1.2, 'DisplayName', 'Measured');
        hold on;
        plot(freqs, ap_fit,  'b--', 'LineWidth', 1.5, 'DisplayName', 'Aperiodic');
        plot(freqs, full_fit,'r-', 'LineWidth', 1.5, 'DisplayName', 'Full model');
        xlabel('Frequency (Hz)'); ylabel('log_{10} Power');
        title(sprintf('t = %.1f s | \\chi=%.2f | R²=%.2f', ...
            win_centers(wi), results.exponent(wi), results.r_squared(wi)));
        legend('Location','northeast','FontSize',7);
        grid on; box off;
        xlim(sprint.freq_range);
    end
    sgtitle('FOOOF Fits at Example Time Windows', 'FontSize', 13, 'FontWeight', 'bold');


end

end












%% =========================================================
%  LOCAL FUNCTIONS
% =========================================================

function ap = aperiodic_model(log_freqs, params, mode)
% Evaluate aperiodic model in log-log space
% params = [offset, (knee,) exponent]
if strcmp(mode, 'knee')
    offset = params(1);
    knee = params(2);
    exp_ = params(3);
    ap = offset - log10(knee + 10.^(log_freqs * exp_));
else % fixed
    offset = params(1);
    exp_ = params(2);
    ap = offset - exp_ * log_freqs;
end
end

% function params = fit_aperiodic(log_freqs, log_psd, mode)
% % Fit aperiodic component using robust linear/nonlinear regression
% if strcmp(mode, 'knee')
%     % Nonlinear fit
%     f0 = [mean(log_psd) + 2, 0, 1.5];
%     opts = optimset('Display','off','TolFun',1e-6,'TolX',1e-6,'MaxIter',500);
%     try
%         params = fminsearch(@(p) sum((log_psd - aperiodic_model(log_freqs,p,'knee')).^2), ...
%             f0, opts);
%         params(3) = max(params(3), 0); % exponent >= 0
%     catch
%         params = fit_aperiodic(log_freqs, log_psd, 'fixed');
%         params = [params(1), 0, params(2)];
%     end
% else
%     % Linear regression in log-log space: log_psd = offset - exp*log_freqs
%     X = [ones(length(log_freqs),1), -log_freqs(:)];
%     b = X \ log_psd(:);
%     params = [b(1), max(b(2), 0)];
% end
% end

function params = fit_aperiodic(log_freqs, log_psd, mode)
if strcmp(mode, 'knee')
    f0 = [mean(log_psd) + 2, 10, 1.5];   % better initial knee guess
    lb = [0,    0,   0  ];
    ub = [20, 1000, 10  ];   % cap knee to prevent runaway
    opts = optimoptions('fmincon','Display','off','TolFun',1e-6,'TolX',1e-6,'MaxIter',500);
    try
        params = fmincon(@(p) sum((log_psd - aperiodic_model(log_freqs,p,'knee')).^2), ...
            f0, [], [], [], [], lb, ub, [], opts);
    catch
        params = fit_aperiodic(log_freqs, log_psd, 'fixed');
        params = [params(1), 0, params(2)];
    end
else
    X = [ones(length(log_freqs),1), -log_freqs(:)];
    b = X \ log_psd(:);
    params = [b(1), max(b(2), 0)];
end
end

function g = gaussian_peak(freqs, peak_p)
% Evaluate a single Gaussian peak in log-power vs linear-freq space
% peak_p = [cf, amplitude, sigma_hz]
if isempty(peak_p), g = zeros(size(freqs)); return; end
cf = peak_p(1); amp = peak_p(2); sigma = peak_p(3);
g  = amp * exp(-(freqs - cf).^2 / (2 * sigma^2));
end

function G = gaussian_peak_matrix(freqs, peaks_found)
% Evaluate multiple Gaussian peaks; returns n_peaks x n_freqs matrix
if isempty(peaks_found)
    G = zeros(1, length(freqs));
    return;
end
G = zeros(size(peaks_found,1), length(freqs));
for i = 1:size(peaks_found,1)
    G(i,:) = gaussian_peak(freqs, peaks_found(i,:));
end
end

function peak_p = fit_gaussian_peak(freqs, log_psd_flat, cf_init, amp_init, bw_init, sprint)
% Fit a single Gaussian peak; returns [cf, amplitude, sigma] or []
% Constrain cf within ±2*bw of initial estimate, and within freq_range
cf_lb  = max(cf_init - 5, sprint.freq_range(1));
cf_ub  = min(cf_init + 5, sprint.freq_range(2));
bw_lb  = sprint.peak_width_limits(1) / (2*sqrt(2*log(2))); % FWHM->sigma
bw_ub  = sprint.peak_width_limits(2) / (2*sqrt(2*log(2)));

% Simple bounded optimization
lb = [cf_lb, 0,    bw_lb];
% ub = [cf_ub, amp_init*3, bw_ub];
ub = [cf_ub, max(amp_init*3, 1.0), bw_ub]; 
p0 = [cf_init, amp_init, bw_init];
p0 = min(max(p0, lb), ub);

cost_fn = @(p) sum((log_psd_flat - gaussian_peak(freqs, p)).^2);
opts = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',300);

try
    % p_fit = fminsearchbnd(cost_fn, p0, lb, ub, opts);
    opts = optimoptions('fmincon','Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',300);
    p_fit = fmincon(cost_fn, p0, [], [], [], [], lb, ub, [], opts);
    % Validate
    if p_fit(2) < sprint.min_peak_height
        peak_p = [];
    elseif p_fit(3) < bw_lb || p_fit(3) > bw_ub
        peak_p = [];
    else
        peak_p = p_fit;
    end
catch
    peak_p = [];
end
end