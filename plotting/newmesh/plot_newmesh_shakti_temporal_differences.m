%PLOT_NEWMESH_SHAKTI_TEMPORAL_DIFFERENCES 绘制 SHAKTI 物理量的时间差分空间图。
%
% 从项目根目录运行:
%   run('plotting/newmesh/plot_newmesh_shakti_temporal_differences.m')
%
% 差分定义:
%   delta_field = field(step k) - field(step k - difference_lag_steps)
%
% 输出:
%   outputs/figures_newmesh/<模型名>/temporal_differences/
%     diff_png/<field_name>/diff_step_###_minus_###_<field_name>.png
%     diff_fig/<field_name>/diff_step_###_minus_###_<field_name>.fig

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

% 2) 差分和输出时间步。
% difference_lag_steps = 1 表示 step k - step k-1。
% difference_lag_steps = 5 表示 step k - step k-5。
difference_lag_steps = 1;

% plot_steps 留空时，按 plot_every_n_steps 自动输出。
% 例如 plot_every_n_steps = 5 表示每 5 个保存步输出一次差分图。
% 如果只想画指定目标步，可以写 plot_steps = [20 40 80 154];
plot_steps = [];
plot_every_n_steps = 12;
include_final_step = true;

% 3) 要画差分的物理量。后两个是由 HydrologyHead 推导出来的诊断量。
plot_fields = { ...
    'HydrologyHead', ...
    'EffectivePressure', ...
    'HydrologyGapHeight', ...
    'HydrologyBasalFlux', ...
    'WaterPressureOverburdenRatio', ...
    'HeadMinusSurface' ...
};

% 4) 色标设置。差分图默认用以 0 为中心的对称色标。
% 'percentile' 用 abs(delta) 的百分位决定色标，避免极端点压扁主体结构。
% 'full'       用 abs(delta) 的最大值决定色标。
% 'fixed'      使用 fixed_symmetric_caxis。
color_limit_mode = 'percentile';
color_abs_percentile = 100;
fixed_symmetric_caxis = [-1, 1];

% 计算色标范围时，是否排除出口轮廓内的点/单元；图上仍然正常显示这些区域。
exclude_outlet_from_color_limits = true;
outlet_exp_file = fullfile(paths.exp, 'Recovery_outlets.exp');

% 是否叠加湖区轮廓和出口轮廓。
overlay_lake_contours = true;
lake_exp_file = fullfile(paths.exp, 'recovery_active_lakes.exp');
lake_line_color = [1.0, 0.85, 0.05];
lake_line_width = 1.4;

overlay_outlet_contours = true;
outlet_line_color = [0.95, 0.05, 0.05];
outlet_line_width = 1.1;

% 5) 输出设置。
outdir = '';
save_png = true;
save_fig = false;
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
fprintf('Transient solutions saved: %d\n', nsteps);
fprintf('Mesh: %d vertices, %d elements\n', md.mesh.numberofvertices, md.mesh.numberofelements);

assert(difference_lag_steps >= 1 && difference_lag_steps == round(difference_lag_steps), ...
    'difference_lag_steps must be a positive integer.');
assert(difference_lag_steps < nsteps, ...
    'difference_lag_steps=%d must be smaller than saved steps=%d.', difference_lag_steps, nsteps);

target_steps = resolve_target_steps(plot_steps, plot_every_n_steps, include_final_step, ...
    difference_lag_steps, nsteps);
plot_specs = build_plot_specs(md, sol, plot_fields);
assert(~isempty(plot_specs), 'None of the requested plot_fields are available.');

if isempty(outdir)
    [~, model_tag] = fileparts(model_file);
    outdir = fullfile(paths.figures_newmesh, model_tag, 'temporal_differences');
end

png_root = fullfile(outdir, 'diff_png');
fig_root = fullfile(outdir, 'diff_fig');
ensure_folder(outdir);
if save_png, ensure_folder(png_root); end
if save_fig, ensure_folder(fig_root); end

png_dirs = prepare_field_output_dirs(png_root, plot_specs, save_png);
fig_dirs = prepare_field_output_dirs(fig_root, plot_specs, save_fig);

outlet_nodes_for_color_limits = false(md.mesh.numberofvertices, 1);
if exclude_outlet_from_color_limits
    outlet_nodes_for_color_limits = read_contour_nodes(md, outlet_exp_file, 'outlet color-limit mask');
end

lake_contours = struct('x', {}, 'y', {});
if overlay_lake_contours
    lake_contours = read_exp_contours(lake_exp_file, 'lake overlay');
end

outlet_contours = struct('x', {}, 'y', {});
if overlay_outlet_contours
    outlet_contours = read_exp_contours(outlet_exp_file, 'outlet overlay');
end

%% 批量输出差分空间图
fprintf('\nStart plotting temporal differences...\n');
for ii = 1:numel(target_steps)
    k = target_steps(ii);
    k0 = k - difference_lag_steps;
    t_days = get_solution_time_days(sol(k));
    t0_days = get_solution_time_days(sol(k0));

    fprintf('[%d/%d] step %d - %d | %.4f d - %.4f d\n', ...
        ii, numel(target_steps), k, k0, t_days, t0_days);

    for p = 1:numel(plot_specs)
        [data_now, label] = get_plot_data(md, sol(k), plot_specs(p));
        [data_prev, ~] = get_plot_data(md, sol(k0), plot_specs(p));

        if isempty(data_now) || isempty(data_prev) || numel(data_now) ~= numel(data_prev)
            warning('Skipping %s at step %d because data are missing or size changed.', plot_specs(p).key, k);
            continue;
        end

        delta_data = data_now(:) - data_prev(:);
        safe_key = make_safe_name(plot_specs(p).key);

        fig = figure('Name', sprintf('%s diff step %03d minus %03d', plot_specs(p).key, k, k0), ...
            'Color', 'w', 'Position', [80 80 980 830], ...
            'Visible', figure_visibility_during_save);

        draw_model_field(md, delta_data);
        colormap(blue_white_red_colormap(256));
        apply_difference_color_limits(md, delta_data, color_limit_mode, color_abs_percentile, ...
            fixed_symmetric_caxis, outlet_nodes_for_color_limits);
        draw_exp_contours(lake_contours, lake_line_color, lake_line_width);
        draw_exp_contours(outlet_contours, outlet_line_color, outlet_line_width);
        title(sprintf('Delta %s | step %d - %d | %.2f d - %.2f d', ...
            label, k, k0, t_days, t0_days), 'Interpreter', 'none');

        base_name = sprintf('diff_step_%03d_minus_%03d_%s', k, k0, safe_key);
        if save_png
            saveas(fig, fullfile(png_dirs.(plot_specs(p).key), [base_name '.png']));
        end
        if save_fig
            save_visible_fig(fig, fullfile(fig_dirs.(plot_specs(p).key), [base_name '.fig']));
        end

        close(fig);
    end
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


function target_steps = resolve_target_steps(plot_steps, plot_every_n_steps, include_final_step, lag_steps, nsteps)
    first_valid = lag_steps + 1;

    if isempty(plot_steps)
        target_steps = first_valid:plot_every_n_steps:nsteps;
        if include_final_step && target_steps(end) ~= nsteps
            target_steps(end+1) = nsteps;
        end
    else
        target_steps = unique(plot_steps(:)');
    end

    if isempty(target_steps)
        error('No target steps selected.');
    end
    if any(target_steps < first_valid) || any(target_steps > nsteps) || any(target_steps ~= round(target_steps))
        error('plot_steps must contain integer indices between %d and %d.', first_valid, nsteps);
    end
end


function specs = build_plot_specs(md, sol, plot_fields)
    specs = struct('key', {}, 'label', {}, 'derived', {});
    D = build_derived_fields(md, sol(1));

    for i = 1:numel(plot_fields)
        key = plot_fields{i};
        if isfield(sol(1), key)
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


function [data, label] = get_plot_data(md, S, spec)
    if spec.derived
        D = build_derived_fields(md, S);
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


function D = build_derived_fields(md, S)
    D = struct();
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


function label = derived_label(key)
    switch key
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


function t_days = get_solution_time_days(S)
    if isfield(S, 'time')
        t_days = S.time * 365;
    elseif isfield(S, 'Time')
        t_days = S.Time * 365;
    else
        t_days = NaN;
    end
end


function dirs = prepare_field_output_dirs(root_dir, specs, enabled)
    dirs = struct();
    if ~enabled
        return;
    end

    for i = 1:numel(specs)
        key = specs(i).key;
        dirs.(key) = fullfile(root_dir, make_safe_name(key));
        ensure_folder(dirs.(key));
    end
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
end


function apply_difference_color_limits(md, data, mode, abs_percentile, fixed_caxis, excluded_nodes)
    switch lower(mode)
        case 'percentile'
            values = color_limit_values(md, data, excluded_nodes);
            values = abs(values(isfinite(values)));
            if isempty(values)
                return;
            end
            max_abs = prctile(values, max(0, min(100, abs_percentile)));
            if ~isfinite(max_abs) || max_abs <= 0
                max_abs = max(values);
            end
            if isfinite(max_abs) && max_abs > 0
                clim([-max_abs, max_abs]);
            end
        case 'full'
            values = color_limit_values(md, data, excluded_nodes);
            values = abs(values(isfinite(values)));
            if ~isempty(values)
                max_abs = max(values);
                if isfinite(max_abs) && max_abs > 0
                    clim([-max_abs, max_abs]);
                end
            end
        case 'fixed'
            if numel(fixed_caxis) == 2 && all(isfinite(fixed_caxis)) && fixed_caxis(2) > fixed_caxis(1)
                clim(fixed_caxis);
            end
        otherwise
            error('Unknown color_limit_mode "%s". Use percentile, full, or fixed.', mode);
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


function cmap = blue_white_red_colormap(n)
    if nargin < 1
        n = 256;
    end
    cmap = interp1([-1 0 1], ...
        [0.1 0.3 0.9; 1 1 1; 0.9 0.1 0.1], ...
        linspace(-1, 1, n));
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


function save_visible_fig(fig, fig_file)
    old_visibility = get(fig, 'Visible');
    set(fig, 'Visible', 'on');
    drawnow;
    savefig(fig, fig_file);
    set(fig, 'Visible', old_visibility);
end


function ensure_folder(folder)
    if exist(folder, 'dir') ~= 7
        mkdir(folder);
    end
end


function name = make_safe_name(name)
    name = regexprep(name, '[^\w\-]', '_');
end


function open_output_folder(folder)
    if ispc
        winopen(folder);
    elseif ismac
        system(sprintf('open "%s"', folder));
    else
        system(sprintf('xdg-open "%s" >/dev/null 2>&1 &', folder));
    end
end
