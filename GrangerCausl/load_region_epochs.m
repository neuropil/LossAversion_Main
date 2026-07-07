function [labels, epochs, fs] = load_region_epochs(filepath)
%% [labels, epochs, fs] = load_region_epochs(filepath)
% >>> ADJUST-ME <<<  Only function that touches your .mat internals.
%
% CLASE layout (from data inspection):
%   file  -> a 135 x 4 table, one ROW per trial
%   col 4 -> 'TrialTablesZS' : per-trial struct of event fields
%   event -> e.g. CHOICE / RESPONSEON / RESPONSEBUTTON / OUTCOME / ITI,
%            each a struct whose .Volts is [nChan x nSamp] bipolar voltages
%            (already z-scored; that's what the "ZS" column is).
%
% Emits ONE epoch per (trial, event), labelled by the event field name, in
% TRIAL ORDER. So position within an event == trial index, which is how the
% driver aligns the same trial across region files.
%
% Returns:
%   labels : nEpoch x 1 cellstr of event names
%   epochs : nEpoch x 1 cell of [nChan x nSamp]
%   fs     : [] unless your file stores it

S  = load(filepath);
fn = fieldnames(S);

% locate the main variable (a table in MATLAB; struct-array fallback)
main = [];
for i = 1:numel(fn)
    v = S.(fn{i});
    if safe_istable(v) || isstruct(v), main = v; break; end
end
if isempty(main)
    error('load_region_epochs:noVar', ...
        'No table/struct in %s. Edit load_region_epochs.m.', filepath);
end

% --- pull out the per-trial ZS column as a struct array or cell ---
zsName = 'TrialTablesZS';
if safe_istable(main)                            % <-- MATLAB table branch
    vn = main.Properties.VariableNames;
    if any(strcmp(vn, zsName))
        col = main.(zsName);
    else
        col = main.(vn{min(4, numel(vn))});      % fall back to column 4
    end
else
    col = main;                                  % already the column
end

if iscell(col), getTrial = @(k) col{k};
else,           getTrial = @(k) col(k); end
nTrial = numel(col);

% --- walk trials x event-fields, grabbing every field that has .Volts ---
labels = {};
epochs = {};
for k = 1:nTrial
    tr = getTrial(k);
    if iscell(tr), tr = tr{1}; end
    if ~isstruct(tr), continue; end
    ev = fieldnames(tr);
    for e = 1:numel(ev)
        node = tr.(ev{e});
        if isstruct(node) && isfield(node,'Volts') && ~isempty(node.Volts)
            labels{end+1,1} = ev{e};          %#ok<AGROW>
            epochs{end+1,1} = double(node.Volts); %#ok<AGROW>
        end
    end
end

if isempty(labels)
    error('load_region_epochs:noVolts', ...
        'Found no event.Volts in %s. Check field names.', filepath);
end
fs = [];
end

function tf = safe_istable(v)
% istable exists on MATLAB; returns false gracefully where it doesn't (Octave)
if exist('istable','builtin') || exist('istable','file')
    tf = istable(v);
else
    tf = false;
end
end
