% Plot Recovery inversion L-curves without automatic corner detection.
%
% This script only reads existing L-curve summary CSV files and redraws the
% L-curve figures. It does not run inversions and does not modify model files.
% The previous automatic corner marker is intentionally removed to keep the
% figures clean for thesis use. All six w101 groups are plotted together.

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();
models_root = paths.models;
out_root = paths.lcurve_figures;

if ~exist(out_root, 'dir')
    mkdir(out_root);
end

% Six available w101 groups in the current experiment archive.
groups = {
    'L_curve_models_500_10_1e-06_to_5e-04'
    'L_curve_models_5000_10_1e-06_to_5e-03'
    'L_curve_models_20000_10_1e-06_to_5e-04'
    'L_curve_models_40000_10_1e-06_to_5e-04'
    'L_curve_models_50000_10_1e-06_to_5e-04'
    'L_curve_models_70000_10_1e-06_to_5e-04'
};

dpi_png = 400;
w501_max_for_w101_5000 = 4.6e-4;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1050 760]);
ax = axes(fig);
hold(ax, 'on');

colors = lines(numel(groups));
legend_entries = {};

for i = 1:numel(groups)
    group_name = groups{i};
    group_dir = fullfile(models_root, group_name);
    summary_file = resolve_summary_csv(group_dir);

    if isempty(summary_file)
        warning('Missing summary CSV in: %s', group_dir);
        continue;
    end

    T = readtable(summary_file);
    T = T(strcmp(string(T.status), "success"), :);

    if isempty(T)
        warning('No successful runs in: %s', summary_file);
        continue;
    end

    T = sortrows(T, 'w501');
    w101 = T.w101(1);
    w103 = T.w103(1);

    if w101 == 5000
        T = T(T.w501 <= w501_max_for_w101_5000, :);
    end

    plot(ax, T.logPhiReg_unweighted, T.logPhiData_unweighted, 'o-', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', colors(i, :) * 0.65, ...
        'Color', colors(i, :));

    legend_entries{end + 1} = sprintf('w101=%g, w103=%g', w101, w103); %#ok<SAGROW>
end

if isempty(legend_entries)
    warning('No curves were plotted.');
else
    xlabel(ax, 'log_{10}(\phi_{reg,unweighted})', 'Interpreter', 'tex');
    ylabel(ax, 'log_{10}(\phi_{data,unweighted})', 'Interpreter', 'tex');
    title(ax, 'Unweighted L-curves for different w101 values', 'Interpreter', 'none');
    legend(ax, legend_entries, 'Location', 'best', 'Interpreter', 'none');
    grid(ax, 'on');
    box(ax, 'on');

    out_base = fullfile(out_root, 'Lcurve_unweighted_all_w101_groups');
    savefig(fig, [out_base '.fig']);
    print(fig, out_base, '-dpng', sprintf('-r%d', dpi_png));
end

close(fig);
fprintf('Done. Figure saved to:\n%s\n', out_root);

function summary_file = resolve_summary_csv(group_dir)
    candidates = dir(fullfile(group_dir, 'Lcurve_summary_w101_*_w103_*_w501scan_*.csv'));

    if ~isempty(candidates)
        if numel(candidates) > 1
            warning('Multiple weighted summary CSV files found in %s. Using %s.', ...
                group_dir, candidates(1).name);
        end
        summary_file = fullfile(group_dir, candidates(1).name);
        return;
    end

    legacy_file = fullfile(group_dir, 'Lcurve_summary.csv');
    if isfile(legacy_file)
        summary_file = legacy_file;
    else
        summary_file = '';
    end
end
