

figure;

subplot(2,1,1);
fooofresidual = specParm_results.FIXED.PeakSpec - specParm_results.FIXED.fit;
plot(specParm_results.General.Freqs,fooofresidual,...
    'k')
xlim([2 55])
hold on
locs = specParm_results.FIXED.peaks.Summary.Hz;
peak_vals = specParm_results.FIXED.peaks.Summary.AmpRAW;
widths = specParm_results.FIXED.peaks.Summary.BW;

% Re-plot the peaks (markers)
plot(locs, peak_vals, 'ro', 'MarkerFaceColor','r');

% Approximate left/right bounds from widths (in Hz)
left_edges  = locs - widths/2;
right_edges = locs + widths/2;

% Draw horizontal width bars at the peak height
for i = 1:numel(locs)
    plot([left_edges(i) right_edges(i)], [peak_vals(i) peak_vals(i)], ...
        'r-', 'LineWidth', 2);
end



% -------------------------------------------------------------------------

fooofresidual = fooofRes.fooofed_spectrum - fooofRes.ap_fit;

subplot(2,1,2);
plot(fooofRes.freqs, fooofresidual, 'k'); 
xlim([3 55])
hold on;

locs = fooofRes.peak_params(:,1);
peak_vals = fooofRes.peak_params(:,2);
widths = fooofRes.peak_params(:,3);

% Re-plot the peaks (markers)
plot(locs, peak_vals, 'ro', 'MarkerFaceColor','r');

% Approximate left/right bounds from widths (in Hz)
left_edges  = locs - widths/2;
right_edges = locs + widths/2;

% Draw horizontal width bars at the peak height
for i = 1:numel(locs)
    plot([left_edges(i) right_edges(i)], [peak_vals(i) peak_vals(i)], ...
        'r-', 'LineWidth', 2);
end


%%

figure;
hold on
plot(fooofRes.power_spectrum)
plot(fooofRes.fooofed_spectrum)
plot(fooofRes.ap_fit)
legend('PSD','Fooof','Fit')
