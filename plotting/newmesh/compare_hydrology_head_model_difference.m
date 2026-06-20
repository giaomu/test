% Compare HydrologyHead from two saved ISSM/SHAKTI model files.
%
% Difference definition:
%   head_difference = AllOutlet model - reference model
%
% Run from the project root:
%   run('plotting/newmesh/compare_hydrology_head_model_difference.m')

clear; clc; close all;

%% User settings

% 两个模型文件。脚本会分别读取 md.results.TransientSolution 里的 HydrologyHead。
alloutlet_model_file = 'D:\Auniversity\ISSM\ISSMmodel\My_ISSM\Antarctic_Wide_Subglacial_Hydrology_Modeling\SHAKTI_AIS\outputs\models_newmesh\RecoveryNewMesh_SHAKTI_0d_to_360d_1800s_np15_AllOutlet.mat';
reference_model_file = 'D:\Auniversity\ISSM\ISSMmodel\My_ISSM\Antarctic_Wide_Subglacial_Hydrology_Modeling\SHAKTI_AIS\outputs\models_newmesh\RecoveryNewMesh_SHAKTI_0d_to_30d_900s_np15.mat';

% 比较哪一个时间步。NaN 表示自动使用各自模型的最后一步。
alloutlet_step = NaN;
reference_step = NaN;

% 输出文件夹。
output_dir = 'D:\Auniversity\ISSM\ISSMmodel\My_ISSM\Antarctic_Wide_Subglacial_Hydrology_Modeling\SHAKTI_AIS\outputs\figures_newmesh\compare_head_difference_AllOutlet_minus_reference';

% 色标设置。true 表示色标以 0 为中心，正负变化更容易看。
use_symmetric_color_axis = true;

% 为了避免极端单点把色标拉得太宽，默认用差值绝对值的 98 百分位。
% 如果想看完整极值范围，改成 100。
color_abs_percentile = 98;

% 是否保存每个节点的两个水头和差值。
save_node_table = true;

% 是否同时画出两个模型水文设置里的 spchead。
save_spchead_maps = true;

%% Main

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);

if ~isfile(alloutlet_model_file)
    error('AllOutlet model does not exist: %s', alloutlet_model_file);
end
if ~isfile(reference_model_file)
    error('Reference model does not exist: %s', reference_model_file);
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

fprintf('Loading AllOutlet model:\n  %s\n', alloutlet_model_file);
md_alloutlet = load_md(alloutlet_model_file);

fprintf('Loading reference model:\n  %s\n', reference_model_file);
md_reference = load_md(reference_model_file);

check_same_mesh(md_alloutlet, md_reference);

[head_alloutlet, alloutlet_step, alloutlet_time_yr] = get_head_at_step(md_alloutlet, alloutlet_step, alloutlet_model_file);
[head_reference, reference_step, reference_time_yr] = get_head_at_step(md_reference, reference_step, reference_model_file);

head_difference = head_alloutlet - head_reference;
stats = summarize_values(head_difference);

fprintf('\nHydrologyHead difference = AllOutlet - reference\n');
fprintf('  AllOutlet step : %d, time = %.6g yr = %.4f days\n', ...
    alloutlet_step, alloutlet_time_yr, alloutlet_time_yr * 365);
fprintf('  Reference step : %d, time = %.6g yr = %.4f days\n', ...
    reference_step, reference_time_yr, reference_time_yr * 365);
fprintf('  vertices       : %d\n', numel(head_difference));
fprintf('  min            : %.6g m\n', stats.min);
fprintf('  p01            : %.6g m\n', stats.p01);
fprintf('  p05            : %.6g m\n', stats.p05);
fprintf('  mean           : %.6g m\n', stats.mean);
fprintf('  median         : %.6g m\n', stats.median);
fprintf('  p95            : %.6g m\n', stats.p95);
fprintf('  p99            : %.6g m\n', stats.p99);
fprintf('  max            : %.6g m\n', stats.max);
fprintf('  std            : %.6g m\n\n', stats.std);

fig = figure('Color', 'w', 'Name', 'HydrologyHead difference', ...
    'Position', [80 80 1100 850]);
ax = axes(fig);

patch(ax, 'Faces', md_alloutlet.mesh.elements, ...
    'Vertices', [md_alloutlet.mesh.x(:), md_alloutlet.mesh.y(:)], ...
    'FaceVertexCData', double(head_difference(:)), ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none');

axis(ax, 'equal');
axis(ax, 'tight');
box(ax, 'on');
xlabel(ax, 'X (m)');
ylabel(ax, 'Y (m)');
colormap(ax, blue_white_red_colormap(256));
cb = colorbar(ax);
ylabel(cb, '\Delta HydrologyHead (m)');
title(ax, sprintf('HydrologyHead difference: AllOutlet step %d - reference step %d [m]', ...
    alloutlet_step, reference_step), 'Interpreter', 'none');

if use_symmetric_color_axis
    finite_abs = abs(head_difference(isfinite(head_difference)));
    clim_abs = local_percentile(finite_abs, color_abs_percentile);
    if ~isfinite(clim_abs) || clim_abs <= 0
        clim_abs = max(finite_abs);
    end
    if isfinite(clim_abs) && clim_abs > 0
        clim(ax, [-clim_abs, clim_abs]);
    end
end

base_name = sprintf('HydrologyHead_difference_AllOutlet_step_%03d_minus_reference_step_%03d', ...
    alloutlet_step, reference_step);
fig_file = fullfile(output_dir, [base_name '.fig']);
png_file = fullfile(output_dir, [base_name '.png']);
stats_file = fullfile(output_dir, [base_name '_stats.csv']);
node_file = fullfile(output_dir, [base_name '_nodes.csv']);

savefig(fig, fig_file);
exportgraphics(fig, png_file, 'Resolution', 300);

stats_table = struct2table(stats);
stats_table.alloutlet_step = alloutlet_step;
stats_table.reference_step = reference_step;
stats_table.alloutlet_time_yr = alloutlet_time_yr;
stats_table.reference_time_yr = reference_time_yr;
writetable(stats_table, stats_file);

if save_node_table
    node_id = (1:numel(head_difference)).';
    x = md_alloutlet.mesh.x(:);
    y = md_alloutlet.mesh.y(:);
    head_diff = head_difference(:);
    node_table = table(node_id, x, y, head_alloutlet(:), head_reference(:), head_diff, ...
        'VariableNames', {'node_id', 'x', 'y', 'head_alloutlet', 'head_reference', 'head_difference'});
    writetable(node_table, node_file);
end

fprintf('Saved figure: %s\n', fig_file);
fprintf('Saved PNG   : %s\n', png_file);
fprintf('Saved stats : %s\n', stats_file);
if save_node_table
    fprintf('Saved nodes : %s\n', node_file);
end

if save_spchead_maps
    save_spchead_map_outputs(md_alloutlet, md_reference, output_dir);
end

%% Local functions

function md = load_md(model_file)
    S = load(model_file, 'md');
    if ~isfield(S, 'md')
        error('File does not contain variable md: %s', model_file);
    end
    md = S.md;
end

function [head, step, time_yr] = get_head_at_step(md, step, model_file)
    try
        results = md.results;
    catch
        error('No md.results found in %s.', model_file);
    end

    if ~isfield(results, 'TransientSolution') || isempty(results.TransientSolution)
        error('No md.results.TransientSolution found in %s.', model_file);
    end

    sol = results.TransientSolution;
    nsteps = numel(sol);
    if isnan(step)
        step = nsteps;
    end
    if step < 1 || step > nsteps || step ~= round(step)
        error('Step must be an integer between 1 and %d for %s.', nsteps, model_file);
    end

    if ~isfield(sol(step), 'HydrologyHead')
        error('HydrologyHead is not available at step %d in %s.', step, model_file);
    end

    head = sol(step).HydrologyHead(:);
    if numel(head) ~= md.mesh.numberofvertices
        error('HydrologyHead length %d does not match number of vertices %d in %s.', ...
            numel(head), md.mesh.numberofvertices, model_file);
    end

    if isfield(sol(step), 'time')
        time_yr = sol(step).time;
    else
        time_yr = NaN;
    end
end

function check_same_mesh(a, b)
    if a.mesh.numberofvertices ~= b.mesh.numberofvertices
        error('The two models have different numbers of vertices.');
    end
    if a.mesh.numberofelements ~= b.mesh.numberofelements
        error('The two models have different numbers of elements.');
    end
    if ~isequal(size(a.mesh.elements), size(b.mesh.elements))
        error('The two models have different element array sizes.');
    end
    if ~isequal(a.mesh.elements, b.mesh.elements)
        error('The two models have different element connectivity.');
    end

    xy_diff = max(abs([a.mesh.x(:) - b.mesh.x(:); a.mesh.y(:) - b.mesh.y(:)]));
    if xy_diff > 1e-6
        error('The two models have different mesh coordinates. Max difference = %.6g m.', xy_diff);
    end
end

function stats = summarize_values(values)
    v = values(isfinite(values));
    if isempty(v)
        error('No finite values found.');
    end

    stats.min = min(v);
    stats.p01 = local_percentile(v, 1);
    stats.p05 = local_percentile(v, 5);
    stats.mean = mean(v);
    stats.median = local_percentile(v, 50);
    stats.p95 = local_percentile(v, 95);
    stats.p99 = local_percentile(v, 99);
    stats.max = max(v);
    stats.std = std(v);
end

function save_spchead_map_outputs(md_alloutlet, md_reference, output_dir)
    spc_alloutlet = get_spchead(md_alloutlet, 'AllOutlet');
    spc_reference = get_spchead(md_reference, 'reference');

    finite_alloutlet = isfinite(spc_alloutlet);
    finite_reference = isfinite(spc_reference);

    fprintf('\nspchead maps\n');
    fprintf('  AllOutlet finite spchead nodes : %d\n', sum(finite_alloutlet));
    fprintf('  Reference finite spchead nodes : %d\n', sum(finite_reference));

    finite_values = [spc_alloutlet(finite_alloutlet); spc_reference(finite_reference)];
    if isempty(finite_values)
        warning('No finite spchead nodes found in either model.');
        return;
    end

    color_min = min(finite_values);
    color_max = max(finite_values);
    if color_min == color_max
        color_min = color_min - 1;
        color_max = color_max + 1;
    end

    base_name = 'spchead_maps_AllOutlet_and_reference';
    fig_file = fullfile(output_dir, [base_name '.fig']);
    png_file = fullfile(output_dir, [base_name '.png']);
    node_file = fullfile(output_dir, [base_name '_nodes.csv']);

    fig = figure('Color', 'w', 'Name', 'spchead maps', ...
        'Position', [80 80 1350 760]);
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    draw_background_mesh(ax1, md_alloutlet, md_alloutlet.geometry.thickness(:));
    hold(ax1, 'on');
    if any(finite_alloutlet)
        scatter(ax1, md_alloutlet.mesh.x(finite_alloutlet), md_alloutlet.mesh.y(finite_alloutlet), ...
            18, spc_alloutlet(finite_alloutlet), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.25);
    end
    colormap(ax1, turbo);
    clim(ax1, [color_min, color_max]);
    cb1 = colorbar(ax1);
    ylabel(cb1, 'spchead (m)');
    title(ax1, 'AllOutlet model spchead', 'Interpreter', 'none');

    ax2 = nexttile;
    draw_background_mesh(ax2, md_reference, md_reference.geometry.thickness(:));
    hold(ax2, 'on');
    if any(finite_reference)
        scatter(ax2, md_reference.mesh.x(finite_reference), md_reference.mesh.y(finite_reference), ...
            18, spc_reference(finite_reference), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.25);
    end
    colormap(ax2, turbo);
    clim(ax2, [color_min, color_max]);
    cb2 = colorbar(ax2);
    ylabel(cb2, 'spchead (m)');
    title(ax2, 'Reference model spchead', 'Interpreter', 'none');

    savefig(fig, fig_file);
    exportgraphics(fig, png_file, 'Resolution', 300);

    node_id = (1:numel(spc_alloutlet)).';
    x = md_alloutlet.mesh.x(:);
    y = md_alloutlet.mesh.y(:);
    node_table = table(node_id, x, y, spc_alloutlet(:), spc_reference(:), ...
        finite_alloutlet(:), finite_reference(:), ...
        'VariableNames', {'node_id', 'x', 'y', 'spchead_alloutlet', 'spchead_reference', ...
        'finite_alloutlet', 'finite_reference'});
    writetable(node_table, node_file);

    fprintf('Saved spchead figure: %s\n', fig_file);
    fprintf('Saved spchead PNG   : %s\n', png_file);
    fprintf('Saved spchead nodes : %s\n', node_file);
end

function spc = get_spchead(md, model_label)
    try
        spc = md.hydrology.spchead(:);
    catch
        error('Cannot read md.hydrology.spchead from %s model.', model_label);
    end

    if numel(spc) ~= md.mesh.numberofvertices
        error('%s spchead length %d does not match number of vertices %d.', ...
            model_label, numel(spc), md.mesh.numberofvertices);
    end
end

function draw_background_mesh(ax, md, bg)
    patch(ax, 'Faces', md.mesh.elements, ...
        'Vertices', [md.mesh.x(:), md.mesh.y(:)], ...
        'FaceVertexCData', double(bg(:)), ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'none');
    axis(ax, 'equal');
    axis(ax, 'tight');
    box(ax, 'on');
    xlabel(ax, 'X (m)');
    ylabel(ax, 'Y (m)');
    colormap(ax, turbo);
    cb = colorbar(ax);
    ylabel(cb, 'Ice thickness (m)');
end

function p = local_percentile(values, pct)
    v = sort(values(isfinite(values)));
    if isempty(v)
        p = NaN;
        return;
    end
    pct = max(0, min(100, pct));
    if numel(v) == 1
        p = v(1);
        return;
    end

    pos = 1 + (pct / 100) * (numel(v) - 1);
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        p = v(lo);
    else
        w = pos - lo;
        p = (1 - w) * v(lo) + w * v(hi);
    end
end

function cmap = blue_white_red_colormap(n)
    if nargin < 1
        n = 256;
    end
    n = max(3, round(n));
    x = linspace(0, 1, n).';

    blue = [0.10, 0.25, 0.85];
    white = [1.00, 1.00, 1.00];
    red = [0.85, 0.10, 0.05];

    cmap = zeros(n, 3);
    lower = x <= 0.5;
    upper = ~lower;

    t = x(lower) / 0.5;
    cmap(lower, :) = (1 - t) .* blue + t .* white;

    t = (x(upper) - 0.5) / 0.5;
    cmap(upper, :) = (1 - t) .* white + t .* red;
end
