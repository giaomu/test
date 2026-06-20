%CHECK_NEWMESH_SHAKTI_STEP_DIAGNOSTICS 诊断一个 SHAKTI 时间步。
%
% 从项目根目录运行:
%   run('plotting/newmesh/check_newmesh_shakti_step_diagnostics.m')
%
% 输出位置:
%   outputs/figures_newmesh/<model_tag>/step_###_diagnostics/
%     spatial_fig/      每个物理量的空间分布 .fig
%     spatial_png/      每个物理量的空间分布 PNG 预览
%     hist_fig/         每个物理量的直方图 .fig
%     hist_png/         每个物理量的直方图 PNG 预览
%     flow_fig/         可选的水头等值线和流向 proxy .fig
%     flow_png/         可选的水头等值线和流向 proxy PNG 预览
%     tables/           统计表和极值点表
clear; clc; close all;

%% 用户设置
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

% 1) 输入模型。留空 model_file_override 时，按下面两个默认文件自动查找。
model_file_override = '';
preferred_model_file = fullfile(paths.models_newmesh, 'RecoveryNewMesh_SHAKTI_360d_to_1080d_1800s_np15_noslide.mat');
fallback_model_file  = fullfile(paths.models_newmesh, 'RecoveryNewMesh_Simulation.mat');

% 2) 诊断时间步和物理量。
target_step = 364;
plot_fields = { ...
    'HydrologyHead', ...
    'EffectivePressure', ...
    'HydrologyGapHeight', ...
    'HydrologyBasalFlux', ...
    'GapHeightMinusBumpHeight', ...
    'Log10HydrologyGapHeight', ...
    'Log10HydrologyBasalFlux', ...
    'WaterPressureOverburdenRatio', ...
    'HeadMinusSurface' ...
};

% GapHeightMinusBumpHeight 使用的基岩粗糙度高度。
% NaN 表示优先读取模型里的 md.hydrology.bump_height。
% 如果模型里没有这个字段，就把这里改成常数，例如 0.1 表示 0.1 m。
manual_bump_height_m = 0.1;

% 3) 空间图色标。这里只影响 colorbar，不改变统计表、直方图和 CSV。
% color_percentiles = [0, 99] 表示每张图按当前场的 0-99 百分位显示。
% 想看完整极值范围时，设 use_percentile_color_limits = false。
use_percentile_color_limits = 1;
color_percentiles = [0, 99.5];

% 计算色标范围时，是否排除出口轮廓内的点/单元；图上仍然正常显示这些区域。
exclude_outlet_from_color_limits = 0;
outlet_exp_file = fullfile(paths.exp, 'Recovery_outlets.exp');

% 是否在空间图上叠加湖区轮廓。
overlay_lake_contours = true;
lake_exp_file = fullfile(paths.exp, 'new_lake.exp');
lake_line_color = [0, 0, 0];
lake_line_width = 1.4;

% 4) 可选的水头和流向诊断。流向 proxy 用 -grad(HydrologyHead) 近似。
save_head_contour_map = true;
save_head_flow_proxy_map = true;
save_head_zero_map = true;
head_zero_tolerance = 1e-9;
flow_proxy_background_field = 'HydrologyGapHeight';
flow_proxy_arrow_stride = 5;
flow_proxy_arrow_scale = 1.3;
head_contour_levels = 25;

% 5) 空间图、直方图和表格输出。
extreme_count = 0;
histogram_bins = 100;
outdir = '';
save_spatial_fig = true;
save_histogram_fig = true;
save_png_preview = true;
save_stats_csv = true;
save_extreme_csv = true;
save_all_values_csv = false;
figure_visibility_during_save = 'off';
open_output_folder_when_done = true;

%% 读取模型
model_file = choose_model_file(model_file_override, preferred_model_file, fallback_model_file);
fprintf('Loading model: %s\n', model_file);
Sload = load(model_file, 'md');
md = Sload.md;
clear Sload

assert(isfield(md.results, 'TransientSolution') && ~isempty(md.results.TransientSolution), ...
    'No TransientSolution found in %s.', model_file);

sol = md.results.TransientSolution;
nsteps = numel(sol);
if target_step < 1 || target_step > nsteps || target_step ~= round(target_step)
    error('target_step must be an integer between 1 and %d.', nsteps);
end

S = sol(target_step);
fprintf('Diagnosing step %d/%d, t = %.6g yr = %.4f days\n', ...
    target_step, nsteps, S.time, S.time * 365);
fprintf('Mesh: %d vertices, %d elements\n', md.mesh.numberofvertices, md.mesh.numberofelements);

if isempty(outdir)
    [~, model_tag] = fileparts(model_file);
    outdir = fullfile(paths.figures_newmesh, model_tag, sprintf('step_%03d_diagnostics', target_step));
end

spatial_fig_dir = fullfile(outdir, 'spatial_fig');
spatial_png_dir = fullfile(outdir, 'spatial_png');
hist_fig_dir = fullfile(outdir, 'hist_fig');
hist_png_dir = fullfile(outdir, 'hist_png');
flow_fig_dir = fullfile(outdir, 'flow_fig');
flow_png_dir = fullfile(outdir, 'flow_png');
tables_dir = fullfile(outdir, 'tables');
values_dir = fullfile(outdir, 'values');

ensure_folder(outdir);
if save_spatial_fig, ensure_folder(spatial_fig_dir); end
if save_histogram_fig, ensure_folder(hist_fig_dir); end
if save_png_preview
    ensure_folder(spatial_png_dir);
    ensure_folder(hist_png_dir);
end
if save_head_contour_map || save_head_flow_proxy_map || save_head_zero_map
    if save_spatial_fig, ensure_folder(flow_fig_dir); end
    if save_png_preview, ensure_folder(flow_png_dir); end
end
if save_stats_csv || save_extreme_csv, ensure_folder(tables_dir); end
if save_all_values_csv, ensure_folder(values_dir); end

plot_specs = build_plot_specs(md, S, plot_fields, manual_bump_height_m);
assert(~isempty(plot_specs), 'None of the requested plot_fields are available.');

outlet_nodes_for_color_limits = false(md.mesh.numberofvertices, 1);
if use_percentile_color_limits && exclude_outlet_from_color_limits
    outlet_nodes_for_color_limits = read_contour_nodes(md, outlet_exp_file, 'outlet color-limit mask');
end

lake_contours = struct('x', {}, 'y', {});
if overlay_lake_contours
    lake_contours = read_exp_contours(lake_exp_file, 'lake overlay');
end

stats_rows = cell(0, 22);
extreme_rows = cell(0, 18);

%% Optional Head-Flow Diagnostics
if save_head_contour_map || save_head_flow_proxy_map || save_head_zero_map
    head = get_result_field(S, 'HydrologyHead');
    if isempty(head) || numel(head) ~= md.mesh.numberofvertices
        warning('HydrologyHead is not available as a vertex field. Skipping head-flow diagnostics.');
    else
        head = head(:);

        if save_head_zero_map
            fig = figure('Name', sprintf('HydrologyHead zero nodes step %03d', target_step), ...
                'Color', 'w', 'Position', [80 80 1080 850], ...
                'Visible', figure_visibility_during_save);
            n_head_zero = draw_head_zero_points_map(md, head, head_zero_tolerance);
            apply_spatial_color_limits(md, head, use_percentile_color_limits, ...
                color_percentiles, outlet_nodes_for_color_limits);
            draw_exp_contours(lake_contours, lake_line_color, lake_line_width);
            title(sprintf('HydrologyHead = 0 nodes | step %d/%d | n = %d | tol = %.1e', ...
                target_step, nsteps, n_head_zero, head_zero_tolerance), 'Interpreter', 'none');
            save_figure_outputs(fig, flow_fig_dir, flow_png_dir, ...
                sprintf('HydrologyHeadZeroPoints_step_%03d', target_step), ...
                save_spatial_fig, save_png_preview);
            close(fig);
        end

        if save_head_contour_map
            fig = figure('Name', sprintf('HydrologyHead contours step %03d', target_step), ...
                'Color', 'w', 'Position', [80 80 1080 850], ...
                'Visible', figure_visibility_during_save);
            draw_head_contour_map(md, head, head_contour_levels);
            apply_spatial_color_limits(md, head, use_percentile_color_limits, ...
                color_percentiles, outlet_nodes_for_color_limits);
            draw_exp_contours(lake_contours, lake_line_color, lake_line_width);
            title(sprintf('HydrologyHead contours | step %d/%d | t = %.4f hours', ...
                target_step, nsteps, S.time * 8760), 'Interpreter', 'none');
            save_figure_outputs(fig, flow_fig_dir, flow_png_dir, ...
                sprintf('HydrologyHeadContours_step_%03d', target_step), ...
                save_spatial_fig, save_png_preview);
            close(fig);
        end

        if save_head_flow_proxy_map
            bg_specs = build_plot_specs(md, S, {flow_proxy_background_field}, manual_bump_height_m);
            if isempty(bg_specs)
                warning('flow_proxy_background_field "%s" is not available. Using HydrologyHead as background.', ...
                    flow_proxy_background_field);
                bg_data = head;
                bg_label = 'HydrologyHead';
            else
                [bg_data, bg_label] = get_plot_data(md, S, bg_specs(1), manual_bump_height_m);
            end

            fig = figure('Name', sprintf('Head-gradient flow proxy step %03d', target_step), ...
                'Color', 'w', 'Position', [100 100 1080 850], ...
                'Visible', figure_visibility_during_save);
            draw_head_flow_proxy_map(md, head, bg_data, bg_label, ...
                flow_proxy_arrow_stride, flow_proxy_arrow_scale);
            apply_spatial_color_limits(md, bg_data, use_percentile_color_limits, ...
                color_percentiles, outlet_nodes_for_color_limits);
            draw_exp_contours(lake_contours, lake_line_color, lake_line_width);
            title(sprintf('-grad(HydrologyHead) flow proxy | step %d/%d | t = %.4f hours', ...
                target_step, nsteps, S.time * 8760), 'Interpreter', 'none');
            save_figure_outputs(fig, flow_fig_dir, flow_png_dir, ...
                sprintf('HeadGradientFlowProxy_step_%03d', target_step), ...
                save_spatial_fig, save_png_preview);
            close(fig);
        end
    end
end

%% Field Diagnostics
for p = 1:numel(plot_specs)
    [data, label] = get_plot_data(md, S, plot_specs(p), manual_bump_height_m);
    safe_key = make_safe_name(plot_specs(p).key);
    [location_type, ids, xs, ys] = data_locations(md, data);
    stats = field_stats(data);

    stats_rows(end+1, :) = { ...
        plot_specs(p).key, label, target_step, S.time, S.time * 365, S.time * 8760, ...
        location_type, stats.n, stats.n_nan, stats.n_inf, ...
        stats.min, stats.p01, stats.p02, stats.p05, stats.mean, stats.median, ...
        stats.p95, stats.p98, stats.p99, stats.p999, stats.max, stats.std}; %#ok<AGROW>

    field_extreme_rows = build_extreme_rows( ...
        plot_specs(p).key, label, target_step, S.time, S.time * 365, ...
        location_type, ids, xs, ys, data, stats, extreme_count);
    extreme_rows = [extreme_rows; field_extreme_rows]; %#ok<AGROW>

    if save_all_values_csv
        values_table = table(ids(:), xs(:), ys(:), data(:), ...
            'VariableNames', {'id', 'x', 'y', 'value'});
        writetable(values_table, fullfile(values_dir, sprintf('%s_step_%03d_values.csv', safe_key, target_step)));
    end

    fig = figure('Name', sprintf('%s step %03d spatial', plot_specs(p).key, target_step), ...
        'Color', 'w', 'Position', [80 80 980 830], ...
        'Visible', figure_visibility_during_save);
    draw_model_field(md, data);
    apply_spatial_color_limits(md, data, use_percentile_color_limits, ...
        color_percentiles, outlet_nodes_for_color_limits);
    draw_exp_contours(lake_contours, lake_line_color, lake_line_width);
    mark_extreme_points(xs, ys, data, extreme_count);
    title(sprintf('%s | step %d/%d | t = %.4f hours', ...
        label, target_step, nsteps, S.time * 8760), 'Interpreter', 'none');
    save_figure_outputs(fig, spatial_fig_dir, spatial_png_dir, ...
        sprintf('%s_step_%03d', safe_key, target_step), ...
        save_spatial_fig, save_png_preview);
    close(fig);

    fig = figure('Name', sprintf('%s step %03d histogram', plot_specs(p).key, target_step), ...
        'Color', 'w', 'Position', [120 120 980 720], ...
        'Visible', figure_visibility_during_save);
    draw_histogram(data, histogram_bins, label, stats);
    title(sprintf('%s histogram | step %d/%d | t = %.4f hours', ...
        label, target_step, nsteps, S.time * 8760), 'Interpreter', 'none');
    save_figure_outputs(fig, hist_fig_dir, hist_png_dir, ...
        sprintf('%s_step_%03d_hist', safe_key, target_step), ...
        save_histogram_fig, save_png_preview);
    close(fig);

    fprintf('Saved diagnostics for %s\n', plot_specs(p).key);
end

if save_stats_csv
    stats_table = cell2table(stats_rows, 'VariableNames', { ...
        'field', 'label', 'step', 'time_yr', 'time_days', 'time_hours', ...
        'location', 'n', 'n_nan', 'n_inf', ...
        'min', 'p01', 'p02', 'p05', 'mean', 'median', ...
        'p95', 'p98', 'p99', 'p999', 'max', 'std'});
    writetable(stats_table, fullfile(tables_dir, sprintf('field_stats_step_%03d.csv', target_step)));
end

if save_extreme_csv
    extreme_table = cell2table(extreme_rows, 'VariableNames', { ...
        'field', 'label', 'step', 'time_yr', 'time_days', ...
        'location', 'extreme_type', 'rank', 'id', 'x', 'y', 'value', ...
        'field_min', 'field_median', 'field_mean', 'field_p98', 'field_p99', 'field_max'});
    writetable(extreme_table, fullfile(tables_dir, sprintf('extreme_points_step_%03d.csv', target_step)));
end

fprintf('\nDone. Output folder:\n  %s\n', outdir);
if open_output_folder_when_done
    open_output_folder(outdir);
end


function model_file = choose_model_file(model_file_override, preferred_model_file, fallback_model_file)
    if ~isempty(model_file_override)
        model_file = model_file_override;
    elseif exist(preferred_model_file, 'file') == 2
        model_file = preferred_model_file;
    elseif exist(fallback_model_file, 'file') == 2
        model_file = fallback_model_file;
    else
        error('Cannot find SHAKTI model. Checked:\n  %s\n  %s', preferred_model_file, fallback_model_file);
    end
end


function specs = build_plot_specs(md, S, plot_fields, manual_bump_height_m)
    specs = struct('key', {}, 'label', {}, 'derived', {});
    D = build_derived_fields(md, S, manual_bump_height_m);

    for i = 1:numel(plot_fields)
        key = plot_fields{i};
        if isfield(S, key)
            specs(end+1).key = key; %#ok<AGROW>
            specs(end).label = key;
            specs(end).derived = false;
        elseif isfield(D, key)
            specs(end+1).key = key; %#ok<AGROW>
            specs(end).label = derived_label(key);
            specs(end).derived = true;
        else
            warning('Requested field "%s" is not available and will be skipped.', key);
        end
    end
end


function [data, label] = get_plot_data(md, S, spec, manual_bump_height_m)
    if spec.derived
        D = build_derived_fields(md, S, manual_bump_height_m);
        data = D.(spec.key);
        label = spec.label;
    else
        raw = get_result_field(S, spec.key);
        [data, suffix] = normalize_plot_data(md, raw);
        label = [spec.label suffix];
    end
    data = data(:);
end


function [data, suffix] = normalize_plot_data(md, raw)
    suffix = '';
    data = raw;
    nv = md.mesh.numberofvertices;
    ne = md.mesh.numberofelements;

    if isempty(raw) || ~isnumeric(raw)
        data = [];
        return;
    end

    if numel(raw) == nv || numel(raw) == ne
        data = raw(:);
    elseif ismatrix(raw) && size(raw, 2) == 2 && (size(raw, 1) == nv || size(raw, 1) == ne)
        data = hypot(raw(:,1), raw(:,2));
        suffix = ' magnitude';
    elseif ismatrix(raw) && size(raw, 1) == 2 && (size(raw, 2) == nv || size(raw, 2) == ne)
        data = hypot(raw(1,:), raw(2,:))';
        suffix = ' magnitude';
    else
        data = raw(:);
    end
end


function D = build_derived_fields(md, S, manual_bump_height_m)
    D = struct();

    gap_height = get_result_field(S, 'HydrologyGapHeight');
    if ~isempty(gap_height)
        [gap_height, ~] = normalize_plot_data(md, gap_height);
        D.Log10HydrologyGapHeight = safe_log10_positive(gap_height);

        bump_height = get_bump_height_for_gap_height(md, gap_height, manual_bump_height_m);
        if ~isempty(bump_height)
            D.GapHeightMinusBumpHeight = gap_height(:) - bump_height(:);
        end
    end

    basal_flux = get_result_field(S, 'HydrologyBasalFlux');
    if ~isempty(basal_flux)
        [basal_flux, ~] = normalize_plot_data(md, basal_flux);
        D.Log10HydrologyBasalFlux = safe_log10_positive(basal_flux);
    end

    head = get_result_field(S, 'HydrologyHead');
    if isempty(head)
        return;
    end

    rho_i = md.materials.rho_ice;
    rho_w = md.materials.rho_freshwater;
    g = md.constants.g;
    head = head(:);
    base = md.geometry.base(:);
    surface = md.geometry.surface(:);
    thickness = md.geometry.thickness(:);

    water_pressure = rho_w * g * (head - base);
    overburden_pressure = rho_i * g * thickness;

    D.WaterPressureOverburdenRatio = safe_divide(water_pressure, overburden_pressure);
    D.HeadMinusSurface = head - surface;
end


function bump_height = get_bump_height_for_gap_height(md, gap_height, manual_bump_height_m)
    bump_height = [];

    if ~isnan(manual_bump_height_m)
        bump_height = manual_bump_height_m + zeros(numel(gap_height), 1);
        return;
    end

    if isfield(md, 'hydrology') && isfield(md.hydrology, 'bump_height')
        model_bump_height = md.hydrology.bump_height(:);
        if numel(model_bump_height) == numel(gap_height)
            bump_height = model_bump_height;
        end
    end
end


function label = derived_label(key)
    switch key
        case 'Log10HydrologyGapHeight'
            label = 'log10(HydrologyGapHeight) [m]';
        case 'GapHeightMinusBumpHeight'
            label = 'HydrologyGapHeight - bump height [m]';
        case 'Log10HydrologyBasalFlux'
            label = 'log10(HydrologyBasalFlux)';
        case 'WaterPressureOverburdenRatio'
            label = 'WaterPressure / overburden [-]';
        case 'HeadMinusSurface'
            label = 'HydrologyHead - surface [m]';
        otherwise
            label = key;
    end
end


function val = get_result_field(S, name)
    if isfield(S, name)
        val = S.(name);
    else
        val = [];
    end
end


function out = safe_divide(a, b)
    out = nan(size(a));
    ok = isfinite(a) & isfinite(b) & b ~= 0;
    out(ok) = a(ok) ./ b(ok);
end


function out = safe_log10_positive(a)
    out = nan(size(a));
    ok = isfinite(a) & a > 0;
    out(ok) = log10(a(ok));
end


function [location_type, ids, xs, ys] = data_locations(md, data)
    data = data(:);
    nv = md.mesh.numberofvertices;
    ne = md.mesh.numberofelements;

    if numel(data) == nv
        location_type = 'vertex';
        ids = (1:nv)';
        xs = md.mesh.x(:);
        ys = md.mesh.y(:);
    elseif numel(data) == ne
        location_type = 'element';
        ids = (1:ne)';
        xs = mean(md.mesh.x(md.mesh.elements), 2);
        ys = mean(md.mesh.y(md.mesh.elements), 2);
    else
        location_type = sprintf('array_%d', numel(data));
        ids = (1:numel(data))';
        xs = nan(size(ids));
        ys = nan(size(ids));
    end
end


function stats = field_stats(data)
    data = data(:);
    finite = data(isfinite(data));
    stats.n = numel(data);
    stats.n_nan = sum(isnan(data));
    stats.n_inf = sum(isinf(data));

    if isempty(finite)
        stats.min = NaN;
        stats.p01 = NaN;
        stats.p02 = NaN;
        stats.p05 = NaN;
        stats.mean = NaN;
        stats.median = NaN;
        stats.p95 = NaN;
        stats.p98 = NaN;
        stats.p99 = NaN;
        stats.p999 = NaN;
        stats.max = NaN;
        stats.std = NaN;
        return;
    end

    stats.min = min(finite);
    pct = prctile(finite, [1 2 5 95 98 99 99.9]);
    stats.p01 = pct(1);
    stats.p02 = pct(2);
    stats.p05 = pct(3);
    stats.mean = mean(finite);
    stats.median = median(finite);
    stats.p95 = pct(4);
    stats.p98 = pct(5);
    stats.p99 = pct(6);
    stats.p999 = pct(7);
    stats.max = max(finite);
    stats.std = std(finite);
end


function draw_head_contour_map(md, head, nlevels)
    draw_model_field(md, head);
    hold on;
    [X, Y, Z] = interpolate_vertex_field_to_grid(md, head, 360);
    contour(X, Y, Z, nlevels, 'k', 'LineWidth', 0.45);
    mark_spchead_points(md);
    xlabel('X (m)');
    ylabel('Y (m)');
end


function n_head_zero = draw_head_zero_points_map(md, head, tolerance)
    if nargin < 3 || isempty(tolerance)
        tolerance = 1e-9;
    end

    draw_model_field(md, head);
    hold on;

    head = head(:);
    head_zero = isfinite(head) & abs(head) <= tolerance;
    n_head_zero = sum(head_zero);

    if any(head_zero)
        plot(md.mesh.x(head_zero), md.mesh.y(head_zero), 'ro', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 4.5, 'LineWidth', 0.8);
    end

    text(0.02, 0.02, sprintf('HydrologyHead = 0 nodes: %d', n_head_zero), ...
        'Units', 'normalized', 'Interpreter', 'none', ...
        'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 5);
    xlabel('X (m)');
    ylabel('Y (m)');
end


function draw_head_flow_proxy_map(md, head, background_data, background_label, arrow_stride, arrow_scale)
    draw_model_field(md, background_data);
    hold on;

    [xc, yc, flow_x, flow_y] = element_head_flow_proxy(md, head);
    mag = hypot(flow_x, flow_y);
    ok = isfinite(xc) & isfinite(yc) & isfinite(mag) & mag > 0;
    idx = find(ok);
    idx = idx(1:max(1, round(arrow_stride)):end);

    if ~isempty(idx)
        arrow_length = typical_mesh_spacing(md) * arrow_scale;
        u = flow_x(idx) ./ mag(idx) * arrow_length;
        v = flow_y(idx) ./ mag(idx) * arrow_length;
        quiver(xc(idx), yc(idx), u, v, 0, 'k', 'LineWidth', 0.75, 'MaxHeadSize', 0.9);
    end

    mark_spchead_points(md);
    text(0.02, 0.02, sprintf('background: %s', background_label), ...
        'Units', 'normalized', 'Interpreter', 'none', ...
        'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 5);
    xlabel('X (m)');
    ylabel('Y (m)');
end


function [xc, yc, flow_x, flow_y] = element_head_flow_proxy(md, head)
    elems = md.mesh.elements;
    x = md.mesh.x(:);
    y = md.mesh.y(:);
    h = head(:);

    x1 = x(elems(:,1)); x2 = x(elems(:,2)); x3 = x(elems(:,3));
    y1 = y(elems(:,1)); y2 = y(elems(:,2)); y3 = y(elems(:,3));
    h1 = h(elems(:,1)); h2 = h(elems(:,2)); h3 = h(elems(:,3));

    twice_area = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);
    dh_dx = (h1 .* (y2 - y3) + h2 .* (y3 - y1) + h3 .* (y1 - y2)) ./ twice_area;
    dh_dy = (h1 .* (x3 - x2) + h2 .* (x1 - x3) + h3 .* (x2 - x1)) ./ twice_area;

    xc = mean(x(elems), 2);
    yc = mean(y(elems), 2);
    flow_x = -dh_dx;
    flow_y = -dh_dy;

    bad = ~isfinite(twice_area) | abs(twice_area) < eps;
    flow_x(bad) = NaN;
    flow_y(bad) = NaN;
end


function [X, Y, Z] = interpolate_vertex_field_to_grid(md, data, ngrid)
    x = md.mesh.x(:);
    y = md.mesh.y(:);
    data = data(:);

    xv = linspace(min(x), max(x), ngrid);
    yv = linspace(min(y), max(y), ngrid);
    [X, Y] = meshgrid(xv, yv);

    F = scatteredInterpolant(x, y, double(data), 'linear', 'none');
    Z = F(X, Y);

    TR = triangulation(md.mesh.elements, x, y);
    inside = ~isnan(pointLocation(TR, [X(:), Y(:)]));
    Z(~reshape(inside, size(Z))) = NaN;
end


function spacing = typical_mesh_spacing(md)
    elems = md.mesh.elements;
    x = md.mesh.x(:);
    y = md.mesh.y(:);

    e12 = hypot(x(elems(:,2)) - x(elems(:,1)), y(elems(:,2)) - y(elems(:,1)));
    e23 = hypot(x(elems(:,3)) - x(elems(:,2)), y(elems(:,3)) - y(elems(:,2)));
    e31 = hypot(x(elems(:,1)) - x(elems(:,3)), y(elems(:,1)) - y(elems(:,3)));
    lengths = [e12; e23; e31];
    lengths = lengths(isfinite(lengths) & lengths > 0);

    if isempty(lengths)
        spacing = 1;
    else
        spacing = median(lengths);
    end
end


function mark_spchead_points(md)
    try
        spchead = md.hydrology.spchead(:);
    catch
        return;
    end

    spc = isfinite(spchead);
    if any(spc)
        plot(md.mesh.x(spc), md.mesh.y(spc), 'wo', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 4.5, 'LineWidth', 0.8);
    end
end


function rows = build_extreme_rows(field, label, step, time_yr, time_days, location_type, ids, xs, ys, data, stats, n_extreme)
    rows = cell(0, 18);
    data = data(:);
    finite_idx = find(isfinite(data));
    if isempty(finite_idx) || n_extreme <= 0
        return;
    end

    [~, order_high] = sort(data(finite_idx), 'descend');
    [~, order_low] = sort(data(finite_idx), 'ascend');
    n = min(n_extreme, numel(finite_idx));

    for r = 1:n
        idx = finite_idx(order_high(r));
        rows(end+1, :) = extreme_row(field, label, step, time_yr, time_days, ...
            location_type, 'high', r, ids(idx), xs(idx), ys(idx), data(idx), stats); %#ok<AGROW>
    end

    for r = 1:n
        idx = finite_idx(order_low(r));
        rows(end+1, :) = extreme_row(field, label, step, time_yr, time_days, ...
            location_type, 'low', r, ids(idx), xs(idx), ys(idx), data(idx), stats); %#ok<AGROW>
    end
end


function row = extreme_row(field, label, step, time_yr, time_days, location_type, extreme_type, rank, id, x, y, value, stats)
    row = {field, label, step, time_yr, time_days, location_type, extreme_type, rank, ...
        id, x, y, value, stats.min, stats.median, stats.mean, stats.p98, stats.p99, stats.max};
end


function draw_model_field(md, data)
    data = data(:);
    if numel(data) == md.mesh.numberofvertices
        patch('Faces', md.mesh.elements, ...
              'Vertices', [md.mesh.x(:), md.mesh.y(:)], ...
              'FaceVertexCData', double(data), ...
              'FaceColor', 'interp', ...
              'EdgeColor', 'none');
    elseif numel(data) == md.mesh.numberofelements
        patch('Faces', md.mesh.elements, ...
              'Vertices', [md.mesh.x(:), md.mesh.y(:)], ...
              'FaceVertexCData', double(data), ...
              'FaceColor', 'flat', ...
              'EdgeColor', 'none');
    else
        text(0.5, 0.5, sprintf('Cannot plot data length %d', numel(data)), ...
            'HorizontalAlignment', 'center');
        axis off;
        return;
    end

    axis equal tight; box on;
    xlabel('X (m)');
    ylabel('Y (m)');
    colorbar;
    colormap(turbo);
end


function apply_spatial_color_limits(md, data, enabled, percentiles, excluded_nodes)
    if ~enabled
        return;
    end

    values = color_limit_values(md, data, excluded_nodes);
    values = values(isfinite(values));
    if isempty(values)
        return;
    end

    percentiles = percentiles(:)';
    if numel(percentiles) ~= 2
        warning('color_percentiles must contain exactly two numbers. Color limits were not changed.');
        return;
    end

    percentiles = sort(max(0, min(100, percentiles)));
    limits = prctile(values, percentiles);
    if numel(limits) == 2 && all(isfinite(limits)) && limits(2) > limits(1)
        clim(limits);
    end
end


function values = color_limit_values(md, data, excluded_nodes)
    data = data(:);
    nv = md.mesh.numberofvertices;
    ne = md.mesh.numberofelements;

    if nargin < 3 || isempty(excluded_nodes)
        excluded_nodes = false(nv, 1);
    end
    excluded_nodes = logical(excluded_nodes(:));

    if numel(data) == nv
        keep = true(nv, 1);
        keep(excluded_nodes) = false;
        values = data(keep);
    elseif numel(data) == ne
        excluded_elements = any(excluded_nodes(md.mesh.elements), 2);
        values = data(~excluded_elements);
    else
        values = data;
    end
end


function nodes = read_contour_nodes(md, exp_file, label)
    nodes = false(md.mesh.numberofvertices, 1);
    if exist(exp_file, 'file') ~= 2
        warning('Cannot find %s file: %s', label, exp_file);
        return;
    end

    try
        nodes = logical(ContourToNodes(md.mesh.x(:), md.mesh.y(:), exp_file, 1));
        nodes = nodes(:);
        if numel(nodes) ~= md.mesh.numberofvertices
            warning('ContourToNodes returned an unexpected node count for %s. Ignoring this mask.', exp_file);
            nodes = false(md.mesh.numberofvertices, 1);
        end
    catch ME
        warning('Failed to read %s from %s:\n%s', label, exp_file, ME.message);
        nodes = false(md.mesh.numberofvertices, 1);
    end
end


function contours = read_exp_contours(exp_file, label)
    contours = struct('x', {}, 'y', {});
    if exist(exp_file, 'file') ~= 2
        warning('Cannot find %s file: %s', label, exp_file);
        return;
    end

    fid = fopen(exp_file, 'r');
    if fid < 0
        warning('Cannot open %s file: %s', label, exp_file);
        return;
    end

    cleanup = onCleanup(@() fclose(fid));
    while ~feof(fid)
        line = next_exp_data_line(fid);
        if isempty(line)
            continue;
        end

        header_values = sscanf(line, '%f');
        if isempty(header_values)
            continue;
        end

        npoints = header_values(1);
        if npoints < 2 || npoints ~= round(npoints)
            continue;
        end

        xy = nan(npoints, 2);
        ok = true;
        for i = 1:npoints
            point_line = next_exp_data_line(fid);
            if isempty(point_line)
                ok = false;
                break;
            end

            point_values = sscanf(point_line, '%f');
            if numel(point_values) < 2
                ok = false;
                break;
            end

            xy(i, :) = point_values(1:2)';
        end

        if ok && all(isfinite(xy(:)))
            contours(end+1).x = xy(:, 1); %#ok<AGROW>
            contours(end).y = xy(:, 2);
        end
    end
end


function line = next_exp_data_line(fid)
    line = '';
    while ~feof(fid)
        candidate = strtrim(fgetl(fid));
        if isempty(candidate) || startsWith(candidate, '#')
            continue;
        end
        line = candidate;
        return;
    end
end


function draw_exp_contours(contours, line_color, line_width)
    if isempty(contours)
        return;
    end

    hold on;
    for i = 1:numel(contours)
        plot(contours(i).x, contours(i).y, '-', ...
            'Color', line_color, 'LineWidth', line_width);
    end
end


function mark_extreme_points(xs, ys, data, n_extreme)
    data = data(:);
    if n_extreme <= 0 || all(isnan(xs)) || all(isnan(ys))
        return;
    end

    finite_idx = find(isfinite(data));
    if isempty(finite_idx)
        return;
    end

    [~, order_high] = sort(data(finite_idx), 'descend');
    [~, order_low] = sort(data(finite_idx), 'ascend');
    n = min(n_extreme, numel(finite_idx));
    high_idx = finite_idx(order_high(1:n));
    low_idx = finite_idx(order_low(1:n));

    hold on;
    plot(xs(low_idx), ys(low_idx), 'wo', 'MarkerSize', 5, 'LineWidth', 0.9);
    plot(xs(high_idx), ys(high_idx), 'k.', 'MarkerSize', 11);
end


function draw_histogram(data, nbins, label, stats)
    finite = data(isfinite(data));
    if isempty(finite)
        text(0.5, 0.5, 'No finite values', 'HorizontalAlignment', 'center');
        axis off;
        return;
    end

    histogram(finite, nbins, 'EdgeColor', 'none');
    grid on;
    xlabel(label, 'Interpreter', 'none');
    ylabel('Count');
    hold on;
    xline(stats.median, 'k-', 'median', 'LabelOrientation', 'horizontal');
    xline(stats.mean, 'b-', 'mean', 'LabelOrientation', 'horizontal');
    xline(stats.p98, 'Color', [0.85 0.33 0.10], 'LineStyle', '--', 'Label', 'p98', ...
        'LabelOrientation', 'horizontal');
    xline(stats.p99, 'Color', [0.49 0.18 0.56], 'LineStyle', '--', 'Label', 'p99', ...
        'LabelOrientation', 'horizontal');
    xline(stats.max, 'r-', 'max', 'LabelOrientation', 'horizontal');

    txt = sprintf(['n = %d\nmin = %.6g\np02 = %.6g\nmedian = %.6g\nmean = %.6g\n' ...
                   'p98 = %.6g\np99 = %.6g\np99.9 = %.6g\nmax = %.6g\nstd = %.6g'], ...
        stats.n, stats.min, stats.p02, stats.median, stats.mean, ...
        stats.p98, stats.p99, stats.p999, stats.max, stats.std);
    text(0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'Interpreter', 'none', 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], ...
        'Margin', 6);
end


function save_figure_outputs(fig, fig_dir, png_dir, base_name, save_fig_file, save_png_file)
    if save_fig_file
        save_visible_fig(fig, fullfile(fig_dir, [base_name '.fig']));
    end
    if save_png_file
        saveas(fig, fullfile(png_dir, [base_name '.png']));
    end
end


function ensure_folder(folder)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end


function safe = make_safe_name(name)
    safe = regexprep(name, '[^A-Za-z0-9_]+', '_');
    safe = regexprep(safe, '_+', '_');
    safe = regexprep(safe, '^_|_$', '');
    if isempty(safe)
        safe = 'field';
    end
end


function save_visible_fig(fig, fig_file)
    original_visibility = get(fig, 'Visible');
    set(fig, 'Visible', 'on');
    savefig(fig, fig_file);
    set(fig, 'Visible', original_visibility);
end


function open_output_folder(folder)
    try
        if ispc
            winopen(folder);
        elseif ismac
            system(sprintf('open "%s"', folder));
        else
            system(sprintf('xdg-open "%s" >/dev/null 2>&1 &', folder));
        end
    catch ME
        fprintf('Could not open output folder automatically: %s\n', ME.message);
    end
end
