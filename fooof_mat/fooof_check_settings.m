% fooof_check_settings() - Check a struct of settings for the FOOOF model.
%
% Usage:
%  >> settings = fooof_check_settings(settings)
%
% Inputs:
%   settings        = struct, can optionally include:
%       settings.peak_width_limts
%       settings.max_n_peaks
%       settings.min_peak_height
%       settings.peak_threshold
%       settings.aperiodic_mode
%       settings.verbose
%
% Outputs:
%   settings        = struct, with all settings defined:
%       settings.peak_width_limts
%       settings.max_n_peaks
%       settings.min_peak_height
%       settings.peak_threshold
%       settings.aperiodic_mode
%       settings.verbose
%
% Notes:
%   This is a helper function, probably not called directly by the user.
%   Any settings not specified are set to default values

function settings = fooof_check_settings(settings)

    % Set defaults for all settings
    defaults = struct(...
        'peak_width_limits', [2, 12], ...
        'max_n_peaks', 100, ...
        'min_peak_height', 0.0, ...
        'peak_threshold', 2.0, ...
        'aperiodic_mode', 'fixed', ...
        'verbose', true);

    %%%%%%%%% JAT added
    defaults.peak_width_limits = py.tuple(num2cell(double(defaults.peak_width_limits(:)')));  % (low, high)

    defaults.max_n_peaks = int64(defaults.max_n_peaks);         % int
    defaults.min_peak_height = double(defaults.min_peak_height);    % float
    defaults.peak_threshold = double(defaults.peak_threshold);     % float
    defaults.aperiodic_mode = char(defaults.aperiodic_mode);       % str
    defaults.verbose = logical(defaults.verbose);           % bool

    %%%%%%%%% JAT ADDED

    % Overwrite any non-existent or nan settings with defaults
    for field = fieldnames(defaults)'
        if ~isfield(settings, field) || all(isnan(settings.(field{1})))
            settings.(field{1}) = defaults.(field{1});
        end
    end

end