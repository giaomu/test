%CHECK_NEWMESH_SHAKTI_TIMESERIES 按指定时间步批量输出 SHAKTI 场的空间图。
%
% 从项目根目录运行:
%   run('plotting/newmesh/check_newmesh_shakti_timeseries.m')
%
% 输出:
%   outputs/figures_newmesh/<输入模型文件名>/selected_steps/
%     fields_png/<field_name>/step_###_<field_name>.png
%     fields_fig/<field_name>/step_###_<field_name>.fig

clear; clc; close all;

%% 用户设置
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

% 如果要手动指定模型文件，就填完整路径；留空则按下面两个文件自动查找。
model_file_override = '';
preferred_model_file = fullfile(paths.models_newmesh, 'RecoveryNewMesh_SHAKTI_360d_to_1080d_1800s_np15_noslide.mat');
fallback_model_file  = fullfile(paths.models_newmesh, 'RecoveryNewMesh_Simulation.mat');

% 要输出的保存时间步。留空时按 plot_every_n_steps 自动抽样。
% 例如: plot_steps = [1 24 48 72];
plot_steps = [];
plot_every_n_steps = 12;
include_final_step = true;

% 要画的物理量。前四个通常是 ISSM 直接输出，后两个由 HydrologyHead 推导。
plot_fields = { ...
    'HydrologyHead', ...
    'EffectivePressure', ...
    'HydrologyGapHeight', ...
    'HydrologyBasalFlux', ...
    'WaterPressureOverburdenRatio', ...
    'HeadMinusSurface' ...
};

% 输出目录。留空时自动使用 outputs/figures_newmesh/<模型名>/selected_steps。
outdir = '';
field_png_root = '';
field_fig_root = '';

% 批量出图默认保存 PNG；如果需要能在 MATLAB 里二次编辑，再打开 FIG。
save_individual_field_png = true;
save_individual_field_fig = false;
figure_visibility_during_save = 'off';
open_output_folder_when_done = true;

% 色标范围。默认 [0,100] 使用真实最小值到最大值。
% 如果主体结构被极值压扁，可以临时改成 [2,98] 或 [0,99.9]。
color_percentiles = [0, 100];
% 色标计算方式：
%   'per_step' 每个时间步单独自适应，避免后面异常值影响前面的图。
%   'global'   同一物理量所有输出时间步共用一个色标，方便跨时间比较。
color_scale_mode = 'per_step';

%% 读取模型
model_file = choose_model_file(model_file_override, preferred_model_file, fallback_model_file);
if isempty(outdir)
    [~, model_tag] = fileparts(model_file);
    outdir = fullfile(paths.figures_newmesh, model_tag, 'selected_steps');
end

field_png_root = fullfile(outdir, 'fields_png');
field_fig_root = fullfile(outdir, 'fields_fig');
ensure_folder(outdir);
if save_individual_field_png, ensure_folder(field_png_root); end
if save_individual_field_fig, ensure_folder(field_fig_root); end

fprintf('Loading model: %s\n', model_file);
Sload = load(model_file, 'md');
md = Sload.md;
clear Sload

assert(isfield(md.results, 'TransientSolution') && ~isempty(md.results.TransientSolution), ...
    'No TransientSolution found in %s.', model_file);

sol = md.results.TransientSolution;
nsteps = numel(sol);
fprintf('Transient solutions: %d\n', nsteps);
fprintf('Mesh: %d vertices, %d elements\n', md.mesh.numberofvertices, md.mesh.numberofelements);

plot_indices = resolve_plot_steps(plot_steps, plot_every_n_steps, include_final_step, nsteps);
plot_specs = build_plot_specs(md, sol, plot_fields);
assert(~isempty(plot_specs), 'None of the requested plot_fields are available.');

global_caxes = build_global_caxes(md, sol, plot_specs, plot_indices, color_percentiles);
field_png_dirs = prepare_field_output_dirs(field_png_root, plot_specs, save_individual_field_png);
field_fig_dirs = prepare_field_output_dirs(field_fig_root, plot_specs, save_individual_field_fig);

%% 批量输出空间图
for ii = 1:numel(plot_indices)
    k = plot_indices(ii);
    fprintf('Plotting step %d/%d, t = %.6g yr = %.4f days\n', ...
        k, nsteps, sol(k).time, sol(k).time * 365);

    for p = 1:numel(plot_specs)
        [data, label] = get_plot_data(md, sol(k), plot_specs(p));
        safe_key = make_safe_name(plot_specs(p).key);
        cax = choose_color_axis(data, global_caxes.(plot_specs(p).key), color_percentiles, color_scale_mode);

        fig = figure('Name', sprintf('%s step %03d', plot_specs(p).key, k), ...
            'Color', 'w', 'Position', [80 80 920 800], ...
            'Visible', figure_visibility_during_save);
        draw_model_field(md, data, cax);
        title(sprintf('%s | step %d/%d | t = %.4f hours', ...
            label, k, nsteps, sol(k).time * 8760), 'Interpreter', 'none');

        if save_individual_field_png
            png_file = fullfile(field_png_dirs.(plot_specs(p).key), ...
                sprintf('step_%03d_%s.png', k, safe_key));
            saveas(fig, png_file);
        end

        if save_individual_field_fig
            fig_file = fullfile(field_fig_dirs.(plot_specs(p).key), ...
                sprintf('step_%03d_%s.fig', k, safe_key));
            save_visible_fig(fig, fig_file);
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


function plot_indices = resolve_plot_steps(plot_steps, plot_every_n_steps, include_final_step, nsteps)
    if isempty(plot_steps)
        plot_indices = 1:plot_every_n_steps:nsteps;
        if include_final_step && plot_indices(end) ~= nsteps
            plot_indices(end+1) = nsteps;
        end
    else
        plot_indices = unique(plot_steps(:)');
    end

    if any(plot_indices < 1) || any(plot_indices > nsteps) || any(plot_indices ~= round(plot_indices))
        error('plot_steps must contain integer step indices between 1 and %d.', nsteps);
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


function cax = choose_color_axis(data, global_cax, pct, color_scale_mode)
    switch lower(color_scale_mode)
        case 'per_step'
            cax = compute_color_axis(data, pct);
        case 'global'
            cax = global_cax;
        otherwise
            error('Unknown color_scale_mode "%s". Use "per_step" or "global".', color_scale_mode);
    end
end


function caxes = build_global_caxes(md, sol, specs, plot_indices, pct)
    caxes = struct();
    for p = 1:numel(specs)
        vals = [];
        for k = plot_indices
            [data, ~] = get_plot_data(md, sol(k), specs(p));
            vals = [vals; data(isfinite(data))]; %#ok<AGROW>
        end

        caxes.(specs(p).key) = compute_color_axis(vals, pct);
    end
end


function cax = compute_color_axis(data, pct)
    vals = data(:);
    vals = vals(isfinite(vals));

    if isempty(vals) || isempty(pct)
        cax = [];
        return;
    end

    cax = prctile(vals, pct);
    if ~all(isfinite(cax)) || cax(1) == cax(2)
        cax = [min(vals), max(vals)];
    end
    if cax(1) == cax(2)
        cax = cax + [-1, 1] * max(abs(cax(1)), 1) * 1e-6;
    end
end


function field_dirs = prepare_field_output_dirs(root_dir, specs, should_create)
    field_dirs = struct();
    if ~should_create
        return;
    end
    for p = 1:numel(specs)
        safe_key = make_safe_name(specs(p).key);
        d = fullfile(root_dir, safe_key);
        ensure_folder(d);
        field_dirs.(specs(p).key) = d;
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


function draw_model_field(md, data, cax)
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
    if ~isempty(cax) && all(isfinite(cax)) && cax(1) < cax(2)
        caxis(cax);
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
