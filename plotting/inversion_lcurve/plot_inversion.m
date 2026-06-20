%% plot_inversion.m
%  批量输出反演结果�?(19 �?
%  全部使用 MATLAB 原生绘图（patch / scatter / histogram），不使�?plotmodel
%  速度相关图使�?log10 色标
%  输出文件夹根据模型文件名自动命名�?模型�?_inversionplot
%  -----------------------------------------------------------------------
%  新增功能�?%    - �?md.results.StressbalanceSolution 回写到顶�?md
%    - 利用 BedMachine mask �?shelf 节点摩擦强制置零
%    - 另存 synced 模型�?outputs/models/
%  -----------------------------------------------------------------------

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

%% ===== 第一步：用户参数 =================================================
model_file      = fullfile(paths.models, 'L_curve_models_70000_10_1e-06_to_5e-04', 'Lcurve_run14_w501_5p500e-05.mat');
%model_file      = fullfile(paths.models, 'Recovery_Inversion.mat');
bedmachine_file = paths.bedmachine_file;
dpi_png         = 400;
marker_size     = 4;
show_mesh_edges = false;
residual_floor  = 1.0;          % 相对残差分母下限 (m/yr)
vel_log_floor   = 1e-1;         % 速度�?log10 时的下限，避�?log(0)
fric_log_floor  = 1e-2;         % 摩擦系数�?log10 时的下限

% --- 是否启用回写并另�?synced 模型 ---
enable_sync_save = 0;

%% ===== 第二步：�?model_file 自动解析三个权重，并创建输出目录 ============
assert(exist(model_file, 'file') == 2, 'model_file 不存�? %s', model_file);

[model_dir, model_name, ~] = fileparts(model_file);
[~, parent_dir_name] = fileparts(model_dir);

% 从父目录解析第一个、第二个权重
% 例如:
%   parent_dir_name = 'L_curve_models_40000_10_1e-06_to_5e-03'
% 解析得到:
%   coeff101_name = '40000'
%   coeff103_name = '10'
[coeff101_name, coeff103_name] = local_parse_w101_w103_from_folder(parent_dir_name);

% 从模型文件名解析第三个权�?% 例如:
%   model_name = 'Lcurve_run20_w501_1p571e-04'
% 解析得到:
%   coeff501_name = '1p571e-04'
coeff501_name = local_parse_w501_from_filename(model_name);

fprintf('Parsed weights from model path:\n');
fprintf('   w101 = %s\n', coeff101_name);
fprintf('   w103 = %s\n', coeff103_name);
fprintf('   w501 = %s\n', coeff501_name);

% 自动生成输出目录
outroot = fullfile(paths.figures, sprintf('L_curve_models_%s_%s_%s_inversion', ...
                  coeff101_name, coeff103_name, coeff501_name));

figdir  = fullfile(outroot, 'fig');
pngdir  = fullfile(outroot, 'png');

% --- 误差图阈值（超过该阈值的节点叠加红点标记�?--
abs_err_threshold  = 1000;   % �?11: |V_model - V_obs| > 此�?(m/yr) 标红
rel_err_threshold  = 0.5;   % �?12: |rel_res|         > 此�?      标红

if ~exist(outroot,'dir'), mkdir(outroot); end
if ~exist(figdir, 'dir'), mkdir(figdir);  end
if ~exist(pngdir, 'dir'), mkdir(pngdir);  end

fprintf('Output root: %s\n', outroot);

%% ===== 第三步：读取模型并提�?mesh ======================================
S  = load(model_file);
md = local_extract_md(S);
clear S;

x        = md.mesh.x(:);
y        = md.mesh.y(:);
elements = md.mesh.elements;
if size(elements,2) > 3
    faces = elements(:,1:3);
else
    faces = elements;
end
nv = md.mesh.numberofvertices;
ne = md.mesh.numberofelements;
fprintf('Vertices : %d\n', nv);
fprintf('Elements : %d\n', ne);

%% ===== 第四步：读取 BedMachine mask 并插值到节点 ========================
[node_type, ~] = local_read_bedmachine_mask(bedmachine_file, x, y);

is_ocean    = node_type == 0;
is_land     = node_type == 1;
is_grounded = node_type == 2;
is_shelf    = node_type == 3;
is_vostok   = node_type == 4;
is_noice    = is_ocean | is_land;

fprintf('Grounded : %d\n', sum(is_grounded));
fprintf('Shelf    : %d\n', sum(is_shelf));
fprintf('No-ice   : %d\n', sum(is_noice));
fprintf('Vostok   : %d\n', sum(is_vostok));

%% ===== 第四.五步：反演结果回�?+ shelf 修正 + 另存 synced 模型 ==========
if enable_sync_save
    fprintf('\n--- Sync-back module ---\n');

    % 检查是否存在反演结�?    has_sol = local_has_member(md,'results') && local_has_member(md.results,'StressbalanceSolution');
    if ~has_sol
        error(['enable_sync_save=true, 但当前模型没�?md.results.StressbalanceSolution。\n' ...
               '无法生成 synced 模型。请检查输入文件�?]);
    end

    sol = md.results.StressbalanceSolution;

    % --- 回写 FrictionCoefficient ---
    if local_has_member(sol, 'FrictionCoefficient')
        md.friction.coefficient = double(sol.FrictionCoefficient(:));
        fprintf('   [OK] FrictionCoefficient -> md.friction.coefficient\n');
    else
        warning('FrictionCoefficient not found in results, skipped.');
    end

    % --- 回写 Vx ---
    vx_synced = false;
    if local_has_member(sol, 'Vx')
        md.initialization.vx = double(sol.Vx(:));
        vx_synced = true;
        fprintf('   [OK] Vx -> md.initialization.vx\n');
    else
        warning('Vx not found in results, skipped.');
    end

    % --- 回写 Vy ---
    vy_synced = false;
    if local_has_member(sol, 'Vy')
        md.initialization.vy = double(sol.Vy(:));
        vy_synced = true;
        fprintf('   [OK] Vy -> md.initialization.vy\n');
    else
        warning('Vy not found in results, skipped.');
    end

    % --- 回写 Vel ---
    if local_has_member(sol, 'Vel')
        md.initialization.vel = double(sol.Vel(:));
        fprintf('   [OK] Vel -> md.initialization.vel\n');
    elseif vx_synced && vy_synced
        md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);
        fprintf('   [OK] Vel computed from synced vx/vy -> md.initialization.vel\n');
    else
        warning('Vel not available and cannot be computed, skipped.');
    end

    % --- vz 清零 ---
    md.initialization.vz = zeros(nv, 1);
    fprintf('   [OK] vz set to zeros\n');

    % --- shelf 节点摩擦强制归零 (双保�? ---
    if ~isempty(md.friction.coefficient)
        n_before = sum(md.friction.coefficient(is_shelf) ~= 0);
        md.friction.coefficient(is_shelf) = 0;
        fprintf('   [OK] Shelf friction zeroed (%d nodes corrected)\n', n_before);
    end

    % --- 从文件名解析 w501 字符�?---
    coeff501_name = local_parse_w501_from_filename(model_name);
    fprintf('   Parsed w501 from filename: %s\n', coeff501_name);

    % --- 另存 synced 模型 ---
    synced_fname = sprintf('%s_%s_%s_synced.mat', coeff101_name, coeff103_name, coeff501_name);
synced_path  = fullfile(paths.models, synced_fname);
    fprintf('   Saving synced model: %s\n', synced_path);
    save(synced_path, 'md', '-v7.3');
    fprintf('   [OK] Synced model saved.\n');

    fprintf('--- Sync-back complete ---\n\n');
end

%% ===== 第五步：准备常用数据�?===========================================

% ---------- 观测速度 ----------
vel_obs = local_try_field(md, {'inversion','vel_obs'}, {'initialization','vel'});
vx_obs  = local_try_field(md, {'inversion','vx_obs'},  {'initialization','vx'});
vy_obs  = local_try_field(md, {'inversion','vy_obs'},  {'initialization','vy'});

% ---------- 模型速度 ----------
vel_model = []; vx_model = []; vy_model = [];
if local_has_member(md,'results') && local_has_member(md.results,'StressbalanceSolution')
    sol = md.results.StressbalanceSolution;
    if local_has_member(sol,'Vel'), vel_model = double(sol.Vel(:)); end
    if local_has_member(sol,'Vx'),  vx_model  = double(sol.Vx(:));  end
    if local_has_member(sol,'Vy'),  vy_model  = double(sol.Vy(:));  end
end

% ---------- 摩擦系数 ----------
fric_coef = [];
if local_has_member(md,'friction') && local_has_member(md.friction,'coefficient')
    fric_coef = double(md.friction.coefficient(:));
end
% 鲁棒回退：如果顶层为空，尝试�?results 读取
if isempty(fric_coef) && local_has_member(md,'results') && local_has_member(md.results,'StressbalanceSolution')
    sol_tmp = md.results.StressbalanceSolution;
    if local_has_member(sol_tmp, 'FrictionCoefficient')
        fric_coef = double(sol_tmp.FrictionCoefficient(:));
        warning('Fallback: using results.StressbalanceSolution.FrictionCoefficient');
    end
end

% ---------- Stressbalance BC mask ----------
bc_mask = false(nv,1);
if local_has_member(md,'stressbalance')
    sb = md.stressbalance;
    for fn = {'spcvx','spcvy','spcvz'}
        if local_has_member(sb, fn{1})
            v = double(sb.(fn{1}));
            if numel(v) >= nv
                bc_mask = bc_mask | ~isnan(v(1:nv));
            end
        end
    end
end

% ---------- 截取�?nv 长度 ----------
vel_obs   = local_clip(vel_obs,   nv);
vx_obs    = local_clip(vx_obs,    nv);
vy_obs    = local_clip(vy_obs,    nv);
vel_model = local_clip(vel_model, nv);
vx_model  = local_clip(vx_model,  nv);
vy_model  = local_clip(vy_model,  nv);
fric_coef = local_clip(fric_coef, nv);

% ---------- 统一准备速度误差量（供图 11 / 12 / 17 / 18 共用�?----------
abs_err = [];
rel_res = [];

if ~isempty(vel_model) && ~isempty(vel_obs)
    abs_err = abs(vel_model - vel_obs);

    rel_res = (vel_model - vel_obs) ./ max(abs(vel_obs), residual_floor);

    % 相对差图与相对差直方图都不看 shelf / no-ice
    rel_res(is_shelf) = NaN;
    rel_res(is_noice) = NaN;
end

%% ========================================================================
%  �? �? (19 �?
%  ========================================================================
edge_flag = show_mesh_edges;

% ---- 01 BedMachine 节点类型�?-------------------------------------------
fig = local_new_fig();
local_plot_categorical(fig, faces, x, y, node_type);
title('BedMachine Node Types');
local_save_both(fig, figdir, pngdir, '01_node_types', dpi_png);

% ---- 02 Stressbalance 边界条件�?----------------------------------------
fig = local_new_fig();
local_plot_background(faces, x, y, [0.85 0.85 0.85]);
hold on;
idx_bc = find(bc_mask);
scatter(x(idx_bc), y(idx_bc), marker_size, 'r', 'filled');
hold off;
axis equal tight; box on; xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Stressbalance BC Nodes (N = %d)', numel(idx_bc)));
local_save_both(fig, figdir, pngdir, '02_stressbalance_BC_nodes', dpi_png);

% ---- 03 bed elevation ---------------------------------------------------
local_safe_plot('md.geometry.bed', md, {'geometry','bed'}, ...
    faces, x, y, nv, edge_flag, 'Bed Elevation (m)', ...
    figdir, pngdir, '03_bed_elevation', dpi_png, false, 0);

% ---- 04 base elevation --------------------------------------------------
local_safe_plot('md.geometry.base', md, {'geometry','base'}, ...
    faces, x, y, nv, edge_flag, 'Base Elevation (m)', ...
    figdir, pngdir, '04_base_elevation', dpi_png, false, 0);

% ---- 05 surface elevation -----------------------------------------------
local_safe_plot('md.geometry.surface', md, {'geometry','surface'}, ...
    faces, x, y, nv, edge_flag, 'Surface Elevation (m)', ...
    figdir, pngdir, '05_surface_elevation', dpi_png, false, 0);

% ---- 06 thickness -------------------------------------------------------
local_safe_plot('md.geometry.thickness', md, {'geometry','thickness'}, ...
    faces, x, y, nv, edge_flag, 'Ice Thickness (m)', ...
    figdir, pngdir, '06_thickness', dpi_png, false, 0);

% ---- 07 观测速度�?(log) -------------------------------------------------
if ~isempty(vel_obs)
    fig = local_new_fig();
    local_plot_log_velocity(fig, faces, x, y, vel_obs, edge_flag, vel_log_floor);
    title('Observed Velocity Magnitude (m/yr)  [log_{10}]');
    local_save_both(fig, figdir, pngdir, '07_observed_velocity', dpi_png);
else
    warning('vel_obs not found, skipping 07.');
end

% ---- 08 反演后模拟速度�?(log) -------------------------------------------
if ~isempty(vel_model)
    fig = local_new_fig();
    local_plot_log_velocity(fig, faces, x, y, vel_model, edge_flag, vel_log_floor);
    title('Modelled Velocity Magnitude (m/yr)  [log_{10}]');
    local_save_both(fig, figdir, pngdir, '08_modelled_velocity', dpi_png);
else
    warning('vel_model not found, skipping 08.');
end

% ---- 09 反演摩擦系数图（线性色标）---------------------------------------
if ~isempty(fric_coef)
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, fric_coef, edge_flag);
    title('Friction Coefficient');
    local_save_both(fig, figdir, pngdir, '09_friction_coefficient', dpi_png);
else
    warning('friction.coefficient not found, skipping 09.');
end

% ---- 10 反演摩擦系数图（log10 色标�?------------------------------------
if ~isempty(fric_coef)
    fig = local_new_fig();
    local_plot_log_field(fig, faces, x, y, fric_coef, edge_flag, fric_log_floor);
    title('Friction Coefficient  [log_{10}]');
    local_save_both(fig, figdir, pngdir, '10_friction_coefficient_log', dpi_png);
else
    warning('friction.coefficient not found, skipping 10.');
end

% ---- 11 速度绝对误差�?(log) --------------------------------------------
if ~isempty(abs_err)
    fig = local_new_fig();
    local_plot_log_velocity(fig, faces, x, y, abs_err, edge_flag, vel_log_floor);

    % 阈值标�?    idx_over = find(abs_err > abs_err_threshold);
    n_over   = numel(idx_over);
    n_valid  = sum(~isnan(abs_err));
    if n_over > 0
        hold on;
        scatter(x(idx_over), y(idx_over), marker_size, 'r', 'filled');
        hold off;
    end

    title(sprintf(['Velocity Absolute Error |V_{model}-V_{obs}| (m/yr)  [log_{10}]\n' ...
                   'Red: |err| > %g m/yr  (N = %d / %d, %.2f%%)'], ...
                   abs_err_threshold, n_over, n_valid, 100*n_over/max(n_valid,1)));
    local_save_both(fig, figdir, pngdir, '11_velocity_absolute_error', dpi_png);

    fprintf('--- Figure 11 | abs_err > %g: %d / %d nodes (%.2f%%)\n', ...
            abs_err_threshold, n_over, n_valid, 100*n_over/max(n_valid,1));
else
    warning('Cannot compute absolute error, skipping 11.');
end

% ---- 11b 非冰架速度绝对误差图（线性色标）--------------------------------
if ~isempty(abs_err)
    abs_err_no_shelf = abs_err;
    abs_err_no_shelf(is_shelf) = NaN;

    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, abs_err_no_shelf, edge_flag);

    idx_over = find(abs_err_no_shelf > abs_err_threshold);
    n_over   = numel(idx_over);
    n_valid  = sum(~isnan(abs_err_no_shelf));
    if n_over > 0
        hold on;
        scatter(x(idx_over), y(idx_over), marker_size, 'r', 'filled');
        hold off;
    end

    title(sprintf(['Velocity Absolute Error without Shelf |V_{model}-V_{obs}| (m/yr)\n' ...
                   'Linear scale; red: |err| > %g m/yr  (N = %d / %d, %.2f%%)'], ...
                   abs_err_threshold, n_over, n_valid, 100*n_over/max(n_valid,1)));
    local_save_both(fig, figdir, pngdir, '11b_velocity_absolute_error_no_shelf_linear', dpi_png);

    fprintf('--- Figure 11b | no-shelf abs_err > %g: %d / %d nodes (%.2f%%)\n', ...
            abs_err_threshold, n_over, n_valid, 100*n_over/max(n_valid,1));
else
    warning('Cannot compute no-shelf absolute error, skipping 11b.');
end

% ---- 12 速度相对误差�?--------------------------------------------------
if ~isempty(rel_res)
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, rel_res, edge_flag);
    colormap(fig, local_diverging_cmap());
    v_clean = rel_res(~isnan(rel_res));
    if ~isempty(v_clean)
        mx = max(abs(v_clean));
        if mx > 0, caxis([-mx mx]); end
    end

    % 阈值标记（用绝对值比较，正负偏差都标�?    idx_over = find(abs(rel_res) > rel_err_threshold);
    n_over   = numel(idx_over);
    n_valid  = sum(~isnan(rel_res));
    if n_over > 0
        hold on;
        scatter(x(idx_over), y(idx_over), marker_size, 'r', 'filled');
        hold off;
    end

    title(sprintf(['Relative Residual: (V_{model}-V_{obs}) / max(|V_{obs}|, floor)\n' ...
                   'Red: |rel\\_res| > %g  (N = %d / %d, %.2f%%)'], ...
                   rel_err_threshold, n_over, n_valid, 100*n_over/max(n_valid,1)));
    local_save_both(fig, figdir, pngdir, '12_velocity_relative_error', dpi_png);

    fprintf('--- Figure 12 | |rel_res| > %g: %d / %d nodes (%.2f%%)\n', ...
            rel_err_threshold, n_over, n_valid, 100*n_over/max(n_valid,1));
else
    warning('Cannot compute relative error, skipping 12.');
end

% ---- 13 Vx 误差�?-------------------------------------------------------
if ~isempty(vx_model) && ~isempty(vx_obs)
    vx_err = vx_model - vx_obs;
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, vx_err, edge_flag);
    colormap(fig, local_diverging_cmap());
    v_clean = vx_err(~isnan(vx_err));
    if ~isempty(v_clean)
        mx = max(abs(v_clean));
        if mx > 0, caxis([-mx mx]); end
    end
    title('Vx Error: Vx_{model} - Vx_{obs} (m/yr)');
    local_save_both(fig, figdir, pngdir, '13_vx_error', dpi_png);
else
    warning('Vx data not found, skipping 13.');
end

% ---- 14 Vy 误差�?-------------------------------------------------------
if ~isempty(vy_model) && ~isempty(vy_obs)
    vy_err = vy_model - vy_obs;
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, vy_err, edge_flag);
    colormap(fig, local_diverging_cmap());
    v_clean = vy_err(~isnan(vy_err));
    if ~isempty(v_clean)
        mx = max(abs(v_clean));
        if mx > 0, caxis([-mx mx]); end
    end
    title('Vy Error: Vy_{model} - Vy_{obs} (m/yr)');
    local_save_both(fig, figdir, pngdir, '14_vy_error', dpi_png);
else
    warning('Vy data not found, skipping 14.');
end

% ---- 15 shelf / no-ice 摩擦检查图 ---------------------------------------
if ~isempty(fric_coef)
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, fric_coef, edge_flag);
    hold on;
    h1 = scatter(x(is_shelf),  y(is_shelf),  marker_size, 'b', 'filled');
    h2 = scatter(x(is_noice),  y(is_noice),  marker_size, 'r', 'filled');
    hold off;
    legend([h1 h2], {'shelf nodes','no-ice nodes'}, 'Location','best');
    title('Friction Coefficient Check on Shelf / No-Ice Nodes');
    local_save_both(fig, figdir, pngdir, '15_friction_check_shelf_noice', dpi_png);

    fc_s = fric_coef(is_shelf);
    fc_n = fric_coef(is_noice);
    fprintf('\n--- Friction on shelf nodes ---\n');
    fprintf('  min=%.4g  max=%.4g  mean=%.4g\n', min(fc_s), max(fc_s), mean(fc_s));
    fprintf('--- Friction on no-ice nodes ---\n');
    fprintf('  min=%.4g  max=%.4g  mean=%.4g\n', min(fc_n), max(fc_n), mean(fc_n));
else
    warning('friction.coefficient not found, skipping 15.');
end

% ---- 16 摩擦系数直方�?--------------------------------------------------
if ~isempty(fric_coef)
    fig = local_new_fig();
    fc_valid = fric_coef(~isnan(fric_coef));

    subplot(2,1,1);
    histogram(fc_valid, 100);
    xlabel('Friction Coefficient'); ylabel('Count');
    title('Friction Coefficient Histogram (all nodes)');
    grid on;

    subplot(2,1,2);
    fc_gr = fric_coef(is_grounded & ~isnan(fric_coef));
    histogram(fc_gr, 100);
    xlabel('Friction Coefficient'); ylabel('Count');
    title('Friction Coefficient Histogram (grounded nodes only)');
    grid on;

    local_save_both(fig, figdir, pngdir, '16_friction_histogram', dpi_png);
else
    warning('friction.coefficient not found, skipping 16.');
end

% ---- 17 速度绝对差直方图 ------------------------------------------------
if ~isempty(abs_err)
    fig = local_new_fig();

    abs_err_valid = abs_err(~isnan(abs_err));
    abs_err_gr    = abs_err(is_grounded & ~isnan(abs_err));

    subplot(2,1,1);
    histogram(abs_err_valid, 100);
    xlabel('|V_{model} - V_{obs}| (m/yr)');
    ylabel('Count');
    title('Velocity Absolute Error Histogram (all valid nodes)');
    grid on;

    subplot(2,1,2);
    histogram(abs_err_gr, 100);
    xlabel('|V_{model} - V_{obs}| (m/yr)');
    ylabel('Count');
    title('Velocity Absolute Error Histogram (grounded nodes only)');
    grid on;

    local_save_both(fig, figdir, pngdir, '17_velocity_absolute_error_histogram', dpi_png);

    fprintf('\n--- Velocity absolute error histogram stats ---\n');
    fprintf('  all valid : min=%.4g  max=%.4g  mean=%.4g  median=%.4g\n', ...
        min(abs_err_valid), max(abs_err_valid), mean(abs_err_valid), median(abs_err_valid));
    if ~isempty(abs_err_gr)
        fprintf('  grounded  : min=%.4g  max=%.4g  mean=%.4g  median=%.4g\n', ...
            min(abs_err_gr), max(abs_err_gr), mean(abs_err_gr), median(abs_err_gr));
    end
else
    warning('abs_err not found, skipping 17.');
end

% ---- 18 速度相对差直方图 ------------------------------------------------
if ~isempty(rel_res)
    fig = local_new_fig();

    rel_res_valid = rel_res(~isnan(rel_res));
    rel_res_gr    = rel_res(is_grounded & ~isnan(rel_res));

    subplot(2,1,1);
    histogram(rel_res_valid, 100);
    xlabel('(V_{model}-V_{obs}) / max(|V_{obs}|, floor)');
    ylabel('Count');
    title('Velocity Relative Residual Histogram (all valid ice nodes)');
    grid on;

    subplot(2,1,2);
    histogram(rel_res_gr, 100);
    xlabel('(V_{model}-V_{obs}) / max(|V_{obs}|, floor)');
    ylabel('Count');
    title('Velocity Relative Residual Histogram (grounded nodes only)');
    grid on;

    local_save_both(fig, figdir, pngdir, '18_velocity_relative_error_histogram', dpi_png);

    fprintf('\n--- Velocity relative residual histogram stats ---\n');
    fprintf('  valid ice : min=%.4g  max=%.4g  mean=%.4g  median=%.4g\n', ...
        min(rel_res_valid), max(rel_res_valid), mean(rel_res_valid), median(rel_res_valid));
    if ~isempty(rel_res_gr)
        fprintf('  grounded  : min=%.4g  max=%.4g  mean=%.4g  median=%.4g\n', ...
            min(rel_res_gr), max(rel_res_gr), mean(rel_res_gr), median(rel_res_gr));
    end
else
    warning('rel_res not found, skipping 18.');
end

%% ===== 完成 =============================================================
fprintf('\nDone. 18 figures generated.\n');
fprintf('FIG -> %s\n', figdir);
fprintf('PNG -> %s\n', pngdir);
close all;

%% ========================================================================
%  �?�?�?�?%  ========================================================================

function md = local_extract_md(S)
    if isfield(S,'md'), md = S.md; return; end
    fns = fieldnames(S);
    for k = 1:numel(fns)
        obj = S.(fns{k});
        if local_has_member(obj,'mesh'), md = obj; return; end
    end
    error('Cannot find md object in MAT file.');
end

function tf = local_has_member(obj, name)
    if isstruct(obj),     tf = isfield(obj, name);
    elseif isobject(obj), tf = isprop(obj, name);
    else,                 tf = false;
    end
end

function v = local_try_field(md, path1, path2)
    v = local_deep_get(md, path1);
    if isempty(v) && nargin >= 3
        v = local_deep_get(md, path2);
        if ~isempty(v)
            warning('Fallback: using %s', strjoin(path2,'.'));
        end
    end
    if ~isempty(v), v = double(v(:)); end
end

function v = local_deep_get(obj, path)
    v = [];
    cur = obj;
    for k = 1:numel(path)
        if local_has_member(cur, path{k})
            cur = cur.(path{k});
        else
            return;
        end
    end
    v = cur;
end

function v = local_clip(v, nv)
    if ~isempty(v) && numel(v) >= nv
        v = v(1:nv);
    end
end

function fig = local_new_fig()
    fig = figure('Visible','off','Color','w');
end

function [node_type, info] = local_read_bedmachine_mask(ncfile, xq, yq)
    x_bm = double(ncread(ncfile,'x'));
    y_bm = double(ncread(ncfile,'y'));
    mask_bm = double(ncread(ncfile,'mask'))';
    [x_bm, ix] = sort(x_bm);
    [y_bm, iy] = sort(y_bm);
    mask_bm = mask_bm(iy, ix);
    F = griddedInterpolant({y_bm, x_bm}, mask_bm, 'nearest','nearest');
    node_type = int32(F(yq, xq));
    info.x_bm = x_bm;  info.y_bm = y_bm;
end

function s = local_parse_w501_from_filename(model_name)
% 从文件名中解�?w501 字符�?% 例如 'Lcurve_run22_w501_7p501e-06' -> '7p501e-06'
    tokens = regexp(model_name, 'w501_([0-9p\-eE]+)', 'tokens');
    if ~isempty(tokens)
        s = tokens{1}{1};
    else
        % 回退：取最后一个下划线后的内容
        parts = strsplit(model_name, '_');
        s = parts{end};
        warning('Could not parse w501 from filename, using last segment: %s', s);
    end
end

function local_plot_continuous(fig, faces, x, y, data, show_edges)
    figure(fig);
    if show_edges, ec = [0.3 0.3 0.3]; else, ec = 'none'; end
    patch('Faces',faces,'Vertices',[x y],...
          'FaceVertexCData',data(:),...
          'FaceColor','interp','EdgeColor',ec);
    axis equal tight; box on;
    xlabel('x (m)'); ylabel('y (m)');
    colormap(fig, parula); colorbar;
end

function local_plot_log_velocity(fig, faces, x, y, data, show_edges, vfloor)
% log10 patch 绘图，colorbar 标注�?10^n（用于速度类正值场�?    figure(fig);
    if show_edges, ec = [0.3 0.3 0.3]; else, ec = 'none'; end
    log_data = log10(max(data(:), vfloor));
    patch('Faces',faces,'Vertices',[x y],...
          'FaceVertexCData',log_data,...
          'FaceColor','interp','EdgeColor',ec);
    axis equal tight; box on;
    xlabel('x (m)'); ylabel('y (m)');
    colormap(fig, parula);
    cb = colorbar;
    lo = floor(min(log_data(~isnan(log_data))));
    hi = ceil(max(log_data(~isnan(log_data))));
    if lo >= hi, hi = lo + 1; end
    tks = lo:hi;
    cb.Ticks = tks;
    cb.TickLabels = arrayfun(@(t) sprintf('10^{%d}',t), tks, 'uni',false);
    caxis([lo hi]);
end

function local_plot_log_field(fig, faces, x, y, data, show_edges, floor_val)
% 通用 log10 patch 绘图（用于摩擦系数等可能含零/负的场）
    figure(fig);
    if show_edges, ec = [0.3 0.3 0.3]; else, ec = 'none'; end
    log_data = log10(max(abs(data(:)), floor_val));
    patch('Faces',faces,'Vertices',[x y],...
          'FaceVertexCData',log_data,...
          'FaceColor','interp','EdgeColor',ec);
    axis equal tight; box on;
    xlabel('x (m)'); ylabel('y (m)');
    colormap(fig, parula);
    cb = colorbar;
    lo = floor(min(log_data(~isnan(log_data))));
    hi = ceil(max(log_data(~isnan(log_data))));
    if lo >= hi, hi = lo + 1; end
    tks = lo:hi;
    cb.Ticks = tks;
    cb.TickLabels = arrayfun(@(t) sprintf('10^{%d}',t), tks, 'uni',false);
    caxis([lo hi]);
end

function local_plot_background(faces, x, y, color)
    patch('Faces',faces,'Vertices',[x y],...
          'FaceColor',color,'EdgeColor','none');
end

function local_plot_categorical(fig, faces, x, y, node_type)
    figure(fig);
    face_type = mode(double(node_type(faces)),2);
    cmap = [0 0.4 0.8; 0.76 0.60 0.42; 0.6 0.6 0.6; 0 0.8 0.8; 0.6 0.2 0.8];
    patch('Faces',faces,'Vertices',[x y],...
          'FaceVertexCData',face_type(:),...
          'FaceColor','flat','EdgeColor','none');
    colormap(fig, cmap);
    caxis([-0.5 4.5]);
    cb = colorbar;
    cb.Ticks = 0:4;
    cb.TickLabels = {'ocean','land','grounded','shelf','vostok'};
    axis equal tight; box on;
    xlabel('x (m)'); ylabel('y (m)');
end

function local_safe_plot(label, md, path, faces, x, y, nv, edge_flag, ttl, ...
                         figdir, pngdir, basename, dpi_png, use_log, vfloor)
    data = local_deep_get(md, path);
    if isempty(data)
        warning('%s not found, skipping %s.', label, basename);
        return;
    end
    data = double(data(:));
    if numel(data) >= nv, data = data(1:nv); end
    fig = local_new_fig();
    if use_log
        local_plot_log_velocity(fig, faces, x, y, data, edge_flag, vfloor);
    else
        local_plot_continuous(fig, faces, x, y, data, edge_flag);
    end
    title(ttl);
    local_save_both(fig, figdir, pngdir, basename, dpi_png);
end

function cmap = local_diverging_cmap()
    n = 128;
    r = [linspace(0,1,n)'; linspace(1,1,n)'];
    g = [linspace(0,1,n)'; linspace(1,0,n)'];
    b = [linspace(1,1,n)'; linspace(1,0,n)'];
    cmap = [r g b];
end

function local_save_both(fig, figdir, pngdir, basename, dpi_png)
    savefig(fig, fullfile(figdir,[basename '.fig']));
    print(fig, fullfile(pngdir, basename), '-dpng', sprintf('-r%d',dpi_png));
end
function [coeff101_name, coeff103_name] = local_parse_w101_w103_from_folder(folder_name)
% 从父目录名中解析 w101 �?w103
%
% 例如:
%   folder_name = 'L_curve_models_40000_10_1e-06_to_5e-03'
% 返回:
%   coeff101_name = '40000'
%   coeff103_name = '10'

    tokens = regexp(folder_name, '^L_curve_models_([^_]+)_([^_]+)_.+$', 'tokens');

    if ~isempty(tokens)
        coeff101_name = tokens{1}{1};
        coeff103_name = tokens{1}{2};
    else
        error(['无法从父目录名中解析 w101/w103�?s\n' ...
               '期望格式类似：L_curve_models_40000_10_1e-06_to_5e-03'], folder_name);
    end
end
