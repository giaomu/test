 % =========================================================================
%  runme_Recovery.m
%  Recovery Glacier — Full Workflow (ISSM + SHAKTI)
%
%  Steps:
%    1. Mesh          — 自适应网格生成
%    2. Parameterize  — 物理参数化
%    3. Inversion     — 摩擦系数反演
%    4. Hydrology     — SHAKTI 水文模型配置
%    5. Simulation    — 瞬态水文模拟
% =========================================================================
clear; clc; close all;
paths = shaktiais_paths();
if ~exist(paths.models, 'dir'), mkdir(paths.models); end

% -------------------------------------------------------------------------
%  0. 全局配置
% -------------------------------------------------------------------------
steps = [1 2];
org   = organizer('repository', paths.models, 'prefix', 'Recovery_', 'steps', steps);
clear steps

% --- 路径 ---
domain_file    = fullfile(paths.exp, 'Recovery_all_Noice.exp');
velocity_file  = paths.velocity_file;
bedmachine_file = paths.bedmachine_file;
par_file       = fullfile(paths.par, 'Recovery_reset_bed.par');
outlet_exp     = fullfile(paths.exp, 'Recovery_outlets.exp');

% --- 初始化模型路径 ---
% 非空且文件存在：直接加载该 .mat 进入 step 2
% 留空 ''       ：退回到 loadmodel(org, 'Mesh')（需要 step 1 已经产出）
init_model_path = fullfile(paths.models, 'Recovery_Initialize.mat');

% --- 网格参数 ---
mesh_hinit = 3000;   % 初始均匀尺寸 (m)
mesh_hmin  = 500;   % 自适应最小尺寸 (m)
mesh_hmax  = 20000;  % 自适应最大尺寸 (m)
mesh_err   = 0.5;    % 插值误差容忍度 (m/yr)

% --- 反演参数 ---
inv_maxsteps      = 200;
inv_maxiter       = 400;
inv_coeff_abs     = 70000;
inv_coeff_log     = 10;
inv_coeff_reg     = 5.5e-5;
inv_friction_min  = 0.05;
inv_friction_max  = 3500;

% --- 水文参数 ---
hydro_head_fraction   = 0.8;   % 初始水头 = base + fraction * (rho_i/rho_w) * H
hydro_gap_height_init = 0.01; % 初始 gap (m)
hydro_gap_height_max  = 0.05;  % cap gap growth during coarse-grid spin-up (m)
hydro_bump_spacing    = 2.0;   % bump 间距 (m), SHAKTI 论文示例 lr=2 m
hydro_bump_height     = 0.1;   % bump 高度 (m), SHAKTI 论文示例 br=0.1 m；0 = 不考虑 opening-by-sliding
hydro_min_thickness   = 25;    % 最小冰厚安全阈值 (m)
hydro_storage         = 1e-4;  % englacial storage coefficient; damps coarse-grid head jumps
hydro_relaxation      = 0.05;   % under-relaxation for nonlinear Shakti head updates
hydro_open_spc_thickness = 500;  % 高基岩开放排水点的薄冰阈值 (m)
hydro_open_spc_base_min  = 300; % 高基岩开放排水点的基岩高程阈值 (m)
hydro_open_spc_rings     = 1;   % 保留参数；当前 spchead=base 只使用 H/base 阈值

% --- Simulation parameters ---
sim_dt_seconds    = 300;        % time step (s), 5 min
sim_n_steps       = 5;          % diagnostic step count
sim_final_days    = sim_n_steps * sim_dt_seconds / 86400;
sim_output_freq   = 1;          % save every step
if sim_n_steps == 1
    sim_result_name = sprintf('Recovery_SHAKTI_OneStep_%ds', sim_dt_seconds);
else
    sim_result_name = sprintf('Recovery_SHAKTI_%dSteps_%ds', sim_n_steps, sim_dt_seconds);
end
n_cores           = 15;         % 并行核数

% =========================================================================
%  STEP 1. Mesh — 自适应网格生成
% =========================================================================
if perform(org, 'Mesh') % {{{

    fprintf('\n===== STEP 1: Mesh =====\n');
    md = model();

    % 1.1 生成初始均匀网格
    md = triangle(md, domain_file, mesh_hinit);

    % 1.2 读取速度数据用于自适应
    x_raw = double(ncread(velocity_file, 'x'));
    y_raw = double(ncread(velocity_file, 'y'));
    [x_sorted, ix] = sort(x_raw);
    [y_sorted, iy] = sort(y_raw);

    vx_raw = double(ncread(velocity_file, 'VX'));
    vy_raw = double(ncread(velocity_file, 'VY'));

    % 清理 → 转置 → 按排序索引重排
    vx_raw(isnan(vx_raw) | vx_raw > 1e5) = 0;
    vy_raw(isnan(vy_raw) | vy_raw > 1e5) = 0;
    vx_grid = vx_raw';  vx_grid = vx_grid(iy, ix);
    vy_grid = vy_raw';  vy_grid = vy_grid(iy, ix);

    vel_grid = sqrt(vx_grid.^2 + vy_grid.^2);

    % 1.3 插值到网格并自适应加密
    vel_on_mesh = InterpFromGridToMesh(x_sorted, y_sorted, vel_grid, ...
                                       md.mesh.x, md.mesh.y, 0);
    vel_on_mesh(isnan(vel_on_mesh)) = 0;

    md = bamg(md, 'hmin', mesh_hmin, 'hmax', mesh_hmax, ...
              'field', vel_on_mesh, 'err', mesh_err, 'gradation', 1.5);

    % 1.4 南极投影 (EPSG:3031)
    md.mesh.epsg = 3031;
    [md.mesh.lat, md.mesh.long] = xy2ll(md.mesh.x, md.mesh.y, -1, 0, 71);

    fprintf('   Mesh: %d nodes, %d elements\n', ...
            md.mesh.numberofvertices, md.mesh.numberofelements);

    savemodel(org, md);

end % }}}

% =========================================================================
%  STEP 2. Parameterize — 物理参数化
% =========================================================================
if perform(org, 'Parameterize') % {{{

    fprintf('\n===== STEP 2: Parameterize =====\n');    
    % ---- 选择模型加载源 ----
    if ~isempty(init_model_path)
        assert(exist(init_model_path, 'file') == 2, ...
            sprintf('初始化模型不存在: %s', init_model_path));
        fprintf('   Loading initialized model: %s\n', init_model_path);
        S  = load(init_model_path);
        md = S.md;
        clear S
    else
        fprintf('   Loading organizer Mesh product\n');
        md = loadmodel(org, 'Mesh');
    end

    % 2.1 流动方程 & 掩码
    md = setflowequation(md, 'SSA', 'all');

    % 2.2 应用 .par 文件
    md = parameterize(md, par_file);
    md.miscellaneous.name = 'Recovery';

    fprintf('   Nodes: %d | BC nodes: %d\n', ...
            md.mesh.numberofvertices, sum(~isnan(md.stressbalance.spcvx)));

    savemodel(org, md);

end % }}}

% =========================================================================
%  STEP 3. Inversion — Basal friction inversion (M1QN3)
% =========================================================================
if perform(org, 'Inversion') % {{{

    fprintf('\n===== STEP 3: Inversion =====\n');
    md = loadmodel(org, 'Parameterize');
    nv = md.mesh.numberofvertices;
    ne = md.mesh.numberofelements;

    % ------------------------------------------------------------
    % 0. 基本检查
    % ------------------------------------------------------------
    assert(~isempty(md.initialization.vx) && ~isempty(md.initialization.vy), ...
        'md.initialization.vx / vy 为空，请先检查 par 文件');
    assert(any(abs(md.initialization.vx) > 0) || any(abs(md.initialization.vy) > 0), ...
        '初始速度全为零，请检查 par 文件');

    % ------------------------------------------------------------
    % 1. 读 BedMachine Antarctica V04.1 + Cliff suppression
    %    —— 与 par 节 1 / 1b 的节点分类完全对齐
    % ------------------------------------------------------------
    bm_file = paths.bedmachine_file;

    masks = classify_bedmachine_nodes(md, bm_file);
    node_type   = masks.node_type;
    is_ocean    = masks.is_ocean;
    is_land     = masks.is_land;
    is_grounded = masks.is_grounded;
    is_shelf    = masks.is_shelf;
    is_vostok   = masks.is_vostok;
    is_noice    = masks.is_noice;
    is_ice      = masks.is_ice;

    % 2.1 friction law
    if ~isa(md.friction, 'friction')
        warning('md.friction 不是 friction()；回填为 friction().');
        md.friction = friction();
    end
    if ~isfield(md.friction, 'coefficient') || isempty(md.friction.coefficient) || numel(md.friction.coefficient) ~= nv
        md.friction.coefficient = 20 * ones(nv,1);
    end
    if ~isfield(md.friction, 'p') || isempty(md.friction.p) || numel(md.friction.p) ~= ne
        md.friction.p = ones(ne,1);
    end
    if ~isfield(md.friction, 'q') || isempty(md.friction.q) || numel(md.friction.q) ~= ne
        md.friction.q = ones(ne,1);
    end
    if ~isfield(md.friction, 'coupling') || isempty(md.friction.coupling)
        md.friction.coupling = 0;
    end

    % 2.2 inversion 观测速
    need_obs = false;
    if ~isfield(md, 'inversion') || isempty(md.inversion)
        need_obs = true;
    else
        if ~isfield(md.inversion, 'vx_obs') || isempty(md.inversion.vx_obs) || numel(md.inversion.vx_obs) ~= nv
            need_obs = true;
        end
        if ~isfield(md.inversion, 'vy_obs') || isempty(md.inversion.vy_obs) || numel(md.inversion.vy_obs) ~= nv
            need_obs = true;
        end
    end
    if need_obs
        warning('md.inversion.vx_obs / vy_obs 缺失；回填为 md.initialization.vx / vy。');
        md.inversion.vx_obs = md.initialization.vx;
        md.inversion.vy_obs = md.initialization.vy;
    end
    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2);

    % 2.3 stressbalance 边界条件（缺失时用 CS 后的 is_ice 回填）
    need_spc = false;
    if isempty(md.stressbalance.spcvx) || numel(md.stressbalance.spcvx) ~= nv
        need_spc = true;
    end
    if isempty(md.stressbalance.spcvy) || numel(md.stressbalance.spcvy) ~= nv
        need_spc = true;
    end
    if isempty(md.stressbalance.spcvz) || numel(md.stressbalance.spcvz) ~= nv
        need_spc = true;
    end
    if need_spc
        warning('stressbalance SPC 缺失；回填为边界冰节点速度约束（基于 BedMachine + CS）。');
        md.stressbalance.spcvx = NaN * ones(nv,1);
        md.stressbalance.spcvy = NaN * ones(nv,1);
        md.stressbalance.spcvz = NaN * ones(nv,1);
        md.stressbalance.referential  = NaN * ones(nv,6);
        md.stressbalance.loadingforce = zeros(nv,3);
        pos_bc = find(logical(md.mesh.vertexonboundary(:)) & is_ice);
        md.stressbalance.spcvx(pos_bc) = md.initialization.vx(pos_bc);
        md.stressbalance.spcvy(pos_bc) = md.initialization.vy(pos_bc);
        md.stressbalance.spcvz(pos_bc) = 0;
    end

    % ------------------------------------------------------------
    % 3. 反演框架
    % ------------------------------------------------------------
    md.inversion = m1qn3inversion(md.inversion);
    md.inversion.iscontrol = 1;
    md.transient.amr_frequency = 0;
    md.verbose = verbose('solution', true, 'control', true, 'convergence', false);

    % ------------------------------------------------------------
    % 4. 摩擦处理（冰架锁死 / 非冰 + 过渡节点激活）
    %    —— is_shelf / is_noice 均为 CS 后版本，与 par 对齐
    % ------------------------------------------------------------

    % 4.1 真·冰架节点：摩擦强制为 0
    md.friction.coefficient(is_shelf) = 0;

    % 4.2 非冰 + CS 吞下来的过渡节点激活：
    %     par 节 5 把这批节点摩擦设为 0。这会在冰/非冰交界处产生巨大的
    %     |∇c| 伪梯度，被 501 正则误算为惩罚，从而压扁冰侧边界摩擦。
    %     这里把"当前摩擦 ≤ 0 或非有限"的非冰节点抬回 20，
    %     让 501 有合理起点、不在交界处产生人造惩罚。
    %     幂等：若已被抬起过，此处不会重复写入。
    reactivate_noice = is_noice & (~isfinite(md.friction.coefficient) | md.friction.coefficient <= 0);
    md.friction.coefficient(reactivate_noice) = max(inv_friction_min, 20);

    fprintf('   Friction: shelf forced to 0 = %d, non-ice reactivated to %g = %d\n', ...
        sum(is_shelf), max(inv_friction_min, 20), sum(reactivate_noice));

    % ------------------------------------------------------------
    % 5. 代价函数
    %    101 = 绝对速度误差
    %    103 = 对数速度误差
    %    501 = 控制参数梯度正则化
    % ------------------------------------------------------------
    md.inversion.cost_functions = [101, 103, 501];
    md.inversion.cost_functions_coefficients = zeros(nv, 3);
    md.inversion.cost_functions_coefficients(:,1) = inv_coeff_abs;
    md.inversion.cost_functions_coefficients(:,2) = inv_coeff_log;
    md.inversion.cost_functions_coefficients(:,3) = inv_coeff_reg;

    % 真·冰架节点：三列全关（不由速度失配驱动，也不参与梯度正则）
    md.inversion.cost_functions_coefficients(is_shelf, :) = 0;

    % 非冰 + 过渡节点：速度失配关，保留 501 平滑
    %md.inversion.cost_functions_coefficients(is_noice, 1:2) = 0;

    fprintf('   Cost masks: shelf(all off)=%d, no-ice(101/103 off)=%d\n', ...
        sum(is_shelf), sum(is_noice));

    % ------------------------------------------------------------
    % 6. 控制变量与上下界
    % ------------------------------------------------------------
    md.inversion.control_parameters      = {'FrictionCoefficient'};
    md.inversion.maxsteps                = inv_maxsteps;
    md.inversion.maxiter                 = inv_maxiter;
    md.inversion.control_scaling_factors = 1;
    minp = inv_friction_min * ones(nv,1);
    maxp = inv_friction_max * ones(nv,1);

    % 真·冰架节点：min=max=0 完全锁死
    minp(is_shelf) = 0;
    maxp(is_shelf) = 0;
    md.inversion.min_parameters = minp;
    md.inversion.max_parameters = maxp;

    % ------------------------------------------------------------
    % 7. 求解器精度
    % ------------------------------------------------------------
    md.stressbalance.restol = 0.01;
    md.stressbalance.reltol = 0.1;
    md.stressbalance.abstol = NaN;

    % ------------------------------------------------------------
    % 8. 求解
    % ------------------------------------------------------------
    md.cluster = generic('name', oshostname(), 'np', n_cores);
    md = solve(md, 'sb');

    % ------------------------------------------------------------
    % 9. 后处理：回写反演结果
    % ------------------------------------------------------------
    assert(isfield(md.results, 'StressbalanceSolution'), '反演失败：无 StressbalanceSolution');
    md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient;

    % 真·冰架节点再次强制为 0（双保险）
    md.friction.coefficient(is_shelf) = 0;

    md.initialization.vx  = md.results.StressbalanceSolution.Vx;
    md.initialization.vy  = md.results.StressbalanceSolution.Vy;
    md.initialization.vz  = zeros(nv,1);
    md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);

    fprintf('   Inversion complete.\n');
    savemodel(org, md);

end % }}}

% =========================================================================
%  STEP 4. Hydrology — SHAKTI 参数配置
% =========================================================================
if perform(org, 'Hydrology') % {{{

    fprintf('\n===== STEP 4: Hydrology Setup =====\n');

    % ------------------------------------------------------------
    %  4.0 加载模型
    % ------------------------------------------------------------
    lcurve_model_path = fullfile(paths.models, 'Recovery_Inversion.mat');
    assert(exist(lcurve_model_path, 'file') == 2, ...
        sprintf('候选模型不存在: %s', lcurve_model_path));

    fprintf('   Loading L-curve candidate: %s\n', lcurve_model_path);
    S  = load(lcurve_model_path);
    md = S.md;
    clear S

    nv = md.mesh.numberofvertices;
    ne = md.mesh.numberofelements;

    % ------------------------------------------------------------
    %  4.0.1 重新用 BedMachine mask 识别节点类型（与 Step 3 一致）
    % ------------------------------------------------------------
    bm_file = paths.bedmachine_file;

    masks = classify_bedmachine_nodes(md, bm_file);
    node_type   = masks.node_type;
    is_ocean    = masks.is_ocean;
    is_land     = masks.is_land;
    is_grounded = masks.is_grounded;
    is_shelf    = masks.is_shelf;
    is_vostok   = masks.is_vostok;
    is_noice    = masks.is_noice;
    is_ice      = masks.is_ice;

    assert(isfield(md.results, 'StressbalanceSolution'), ...
        '候选模型 md.results.StressbalanceSolution 缺失，无法回写反演结果');

    md.friction.coefficient       = md.results.StressbalanceSolution.FrictionCoefficient;
    md.friction.coefficient(is_shelf) = 0;                   % shelf nodes have no basal friction
    md.initialization.vx          = md.results.StressbalanceSolution.Vx;
    md.initialization.vy          = md.results.StressbalanceSolution.Vy;
    md.initialization.vz          = zeros(nv, 1);
    md.initialization.vel         = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);

    fprintf('   Inversion result written back to md.friction and md.initialization.\n');
    fprintf('   Mean friction (grounded): %.2f\n', ...
        mean(md.friction.coefficient(is_grounded)));

    md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);
    % ------------------------------------------------------------
    %  4.1 关闭报错阈值 & 清理反演开关
    % ------------------------------------------------------------
    md.settings.solver_residue_threshold = NaN;
    md.inversion.iscontrol = 0;

    if isprop(md.mask, 'ocean_levelset')
        pos_clear = md.mask.ocean_levelset < 0 | md.mask.ice_levelset > 0;
        md.friction.coefficient(pos_clear) = 0;
    end
    assert(mean(md.friction.coefficient) > 0, '摩擦系数全为零');

    % ------------------------------------------------------------
    %  4.2 安全修正：极薄冰
    % ------------------------------------------------------------
    pos_thin = md.geometry.thickness < hydro_min_thickness;
    if any(pos_thin)
        md.geometry.thickness(pos_thin) = hydro_min_thickness;
        md.geometry.surface(pos_thin)   = md.geometry.base(pos_thin) + hydro_min_thickness;
        fprintf('   Fixed %d nodes with H < %d m\n', sum(pos_thin), hydro_min_thickness);
    end

    % ------------------------------------------------------------
    %  4.3 初始化 SHAKTI
    % ------------------------------------------------------------
    md.hydrology            = hydrologyshakti();
    md.hydrology.relaxation = hydro_relaxation;
    md.hydrology.storage = hydro_storage * ones(nv, 1);
    md.hydrology.gap_height_min = 1e-3;
    md.hydrology.gap_height_max = hydro_gap_height_max;
    md.hydrology.melt_flag = 0;
    md.hydrology.requested_outputs = {'default'};

    assert(hydro_bump_spacing > 0, 'hydro_bump_spacing must be positive.');
    assert(hydro_bump_height >= 0, 'hydro_bump_height must be nonnegative.');
    if hydro_bump_height == 0
        fprintf('   Warning: hydro_bump_height=0; opening-by-sliding is disabled.\n');
    elseif hydro_bump_height > hydro_bump_spacing
        fprintf(['   Warning: hydro_bump_height > hydro_bump_spacing; ' ...
                 'check that bump height/spacing were not swapped.\n']);
    end

    rho_i = md.materials.rho_ice;
    rho_w = md.materials.rho_freshwater;

    % 初始水头
    md.hydrology.head = md.geometry.base + ...
        hydro_head_fraction * (rho_i / rho_w) * md.geometry.thickness;

    % 分布式排水几何
    md.hydrology.gap_height   = hydro_gap_height_init * ones(ne, 1);
    md.hydrology.bump_spacing = hydro_bump_spacing     * ones(ne, 1);
    md.hydrology.bump_height  = hydro_bump_height      * ones(ne, 1);
    md.hydrology.reynolds     = 1000 * ones(ne, 1);

    % 外部供水 = 0（冬季模式）
    md.hydrology.englacial_input = zeros(nv, 1);
    md.hydrology.moulin_input    = zeros(nv, 1);
    md.hydrology.neumannflux     = zeros(ne, 1);

    % ------------------------------------------------------------
    %  4.4 Dirichlet 边界条件 (spchead)
    % ------------------------------------------------------ ------
    spc_cfg.thin_h        = hydro_open_spc_thickness;
    spc_cfg.high_base     = hydro_open_spc_base_min;
    spc_cfg.open_rings    = hydro_open_spc_rings;
    md.hydrology.spchead = build_spchead(md, bm_file, masks, outlet_exp, spc_cfg);
    n_spchead = sum(~isnan(md.hydrology.spchead));
    fprintf('   Prescribed-head nodes: %d\n', n_spchead);
    validate_spchead_threshold(md, spc_cfg);
    plot_spchead_diagnostics(md, masks, fullfile(paths.figures, 'Recovery_step4_diagnostics'));
    assert(n_spchead > 0, 'md.hydrology.spchead 全为 NaN，没有设置水头边界条件。');
    assert(all(isfinite(md.hydrology.head)), 'md.hydrology.head contains NaN or Inf.');
    assert(all(isfinite(md.hydrology.gap_height)), 'md.hydrology.gap_height contains NaN or Inf.');
    assert(all(isfinite(md.hydrology.storage(:))) && all(md.hydrology.storage(:) >= 0), ...
        'md.hydrology.storage must be finite and nonnegative.');
    assert(all(isfinite(md.hydrology.bump_spacing(:))) && all(md.hydrology.bump_spacing(:) > 0), ...
        'md.hydrology.bump_spacing must be finite and positive.');
    assert(all(isfinite(md.hydrology.bump_height(:))) && all(md.hydrology.bump_height(:) >= 0), ...
        'md.hydrology.bump_height must be finite and nonnegative.');
    if any(md.hydrology.storage(:) == 0)
        fprintf('   Warning: hydrology.storage contains zero values; running zero-storage SHAKTI sensitivity.\n');
    end
    assert(all(isfinite(md.basalforcings.groundedice_melting_rate)), ...
        'md.basalforcings.groundedice_melting_rate contains NaN or Inf.');

    % ------------------------------------------------------------
    %  4.5 摩擦耦合（保持 coupling=0，与原脚本一致）
    % ------------------------------------------------------------
    C    = md.friction.coefficient;
    g    = md.constants.g;
    Neff = rho_i * g * md.geometry.thickness - ...
           rho_w * g * (md.hydrology.head - md.geometry.base);
    Neff(Neff < 0) = 0;
    md.friction.effective_pressure = Neff;
    md.friction.coupling           = 0;

    % ------------------------------------------------------------
    %  4.6 Transient 开关：仅水文
    % ------------------------------------------------------------
    md.transient                 = deactivateall(md.transient);
    md.transient.ishydrology     = 1;
    md.transient.isstressbalance = 0;

    savemodel(org, md);

end % }}}

% =========================================================================
%  STEP 5. Simulation — 瞬态水文模拟
% =========================================================================
if perform(org, 'Simulation') % {{{

    fprintf('\n===== STEP 5: Transient Simulation =====\n');
    md = loadmodel(org, 'Hydrology');

    % 5.2 时间步长
    md.timestepping.start_time = 0.0;
    md.timestepping.time_step  = sim_dt_seconds / md.constants.yts;
    md.timestepping.final_time = sim_final_days / 365;
    assert(sim_n_steps >= 1 && mod(sim_n_steps, 1) == 0, 'sim_n_steps must be a positive integer.');
    fprintf('   Time step: %.1f s, number of steps: %d, final_time: %.9g yr\n', ...
        sim_dt_seconds, sim_n_steps, md.timestepping.final_time);

    %md.hydrology.melt_flag = 0;


    % 5.3 输出 & 求解器
    md.settings.output_frequency = sim_output_freq;
    md.stressbalance.restol  = 0.05;
    md.stressbalance.reltol  = 0.05;
    md.stressbalance.abstol  = NaN;
    md.stressbalance.maxiter = 200;

    % 5.4 日志
    md.verbose.solution    = 1;
    md.verbose.module      = 0;
    md.verbose.convergence = 0;

    % 5.5 求解
    md.cluster = generic('name', oshostname(), 'np', n_cores);
    solve_failed = false;
    solve_error_identifier = '';
    solve_error_message = '';
    solve_error_report = '';
    solve_error_stack = struct([]);

    try
        md_solved = solve(md, 'Transient');
        md = md_solved;
    catch ME
        solve_failed = true;
        solve_error_identifier = ME.identifier;
        solve_error_message = ME.message;
        solve_error_report = getReport(ME, 'extended', 'hyperlinks', 'off');
        solve_error_stack = ME.stack;
        fprintf(2, '\nTransient solve failed; saving model for diagnostics anyway.\n');
        fprintf(2, 'Error: %s\n', solve_error_message);
    end

    has_transient_results = false;
    try
        has_transient_results = isfield(md.results, 'TransientSolution') && ...
            ~isempty(md.results.TransientSolution);
    catch
        has_transient_results = false;
    end
    if has_transient_results
        last_time = md.results.TransientSolution(end).time;
    else
        last_time = NaN;
    end

    % 先保存当前模型和错误信息，再判断是否提前停止。
    % 这样即使中途发散、没有完成全部时间步，也能保留可诊断的模型文件。
    md.miscellaneous.name = sim_result_name;
    diagnostic_file = fullfile(paths.models, [sim_result_name '.mat']);
    save(diagnostic_file, 'md', 'solve_failed', 'solve_error_identifier', ...
        'solve_error_message', 'solve_error_report', 'solve_error_stack', ...
        'has_transient_results', 'last_time', '-v7.3');
    fprintf('Diagnostic model saved to %s (last_time %.6g yr, solve_failed=%d)\n', ...
        diagnostic_file, last_time, solve_failed);

    if ~has_transient_results
        fprintf(2, 'Warning: no TransientSolution was returned; saved pre-solve model and error info.\n');
    elseif last_time < md.timestepping.final_time - md.timestepping.time_step
        fprintf(2, 'Warning: Transient solve stopped early at %.6g yr, before final_time %.6g yr.\n', ...
            last_time, md.timestepping.final_time);
    end

    % 5.6 保存
    try
        savemodel(org, md);
    catch ME_save
        fprintf(2, 'Warning: savemodel failed after diagnostic save: %s\n', ME_save.message);
    end
    fprintf('Diagnostic model saved to %s\n', diagnostic_file);
    if solve_failed
        fprintf('Done with failed solve; diagnostic model was saved.\n');
    else
        fprintf('Done.\n');
    end

end % }}}


function masks = classify_bedmachine_nodes(md, bedmachine_file)

    x_bm    = double(ncread(bedmachine_file, 'x'));
    y_bm    = double(ncread(bedmachine_file, 'y'));
    mask_bm = double(ncread(bedmachine_file, 'mask'));

    [x_bm, ix] = sort(x_bm);
    [y_bm, iy] = sort(y_bm);
    mask_bm = mask_bm';
    mask_bm = mask_bm(iy, ix);

    Fmask = griddedInterpolant({y_bm, x_bm}, mask_bm, 'nearest', 'nearest');
    node_type = round(Fmask(md.mesh.y, md.mesh.x));

    is_ocean_raw    = (node_type == 0);
    is_land_raw     = (node_type == 1);
    is_grounded_raw = (node_type == 2);
    is_shelf_raw    = (node_type == 3);
    is_vostok_raw   = (node_type == 4);

    is_ice_raw   = is_grounded_raw | is_shelf_raw | is_vostok_raw;
    is_noice_raw = is_ocean_raw | is_land_raw;

    fprintf('   Raw BedMachine: ocean=%d, land=%d, grounded=%d, shelf=%d, vostok=%d\n', ...
        sum(is_ocean_raw), sum(is_land_raw), sum(is_grounded_raw), sum(is_shelf_raw), sum(is_vostok_raw));

    if any(is_vostok_raw)
        warning('Lake Vostok nodes detected; keeping them as grounded-like ice nodes.');
    end

    is_ice = is_ice_raw;
    nonice_elem = any(is_noice_raw(md.mesh.elements), 2);
    is_ice(md.mesh.elements(nonice_elem,:)) = false;

    is_noice    = ~is_ice;
    is_ocean    = is_ocean_raw;
    is_land     = is_land_raw;
    is_grounded = is_grounded_raw & is_ice;
    is_shelf    = is_shelf_raw & is_ice;
    is_vostok   = is_vostok_raw & is_ice;

    n_demoted = sum(is_ice_raw & ~is_ice);
    fprintf('   After cliff suppression: ice=%d, non-ice=%d (demoted from ice = %d, shelf kept = %d)\n', ...
        sum(is_ice), sum(is_noice), n_demoted, sum(is_shelf));

    if sum(is_ice) == 0
        error('Cliff suppression removed all ice nodes. Check BedMachine mask and model mesh alignment.');
    end

    masks.node_type   = node_type;
    masks.is_ocean    = is_ocean;
    masks.is_land     = is_land;
    masks.is_grounded = is_grounded;
    masks.is_shelf    = is_shelf;
    masks.is_vostok   = is_vostok;
    masks.is_noice    = is_noice;
    masks.is_ice      = is_ice;
end


function spchead = build_spchead(md, bedmachine_file, masks, outlet_file, spc_cfg)

    nv = md.mesh.numberofvertices;
    spchead = NaN(nv,1);
    on_boundary = logical(md.mesh.vertexonboundary(:));

    if nargin < 5 || isempty(spc_cfg)
        spc_cfg.thin_h     = 50;
        spc_cfg.high_base  = 400;
        spc_cfg.open_rings = 1;
    end

    % 先读取显式 outlet；真正赋值放在 base 阈值之后，用 head=0 覆盖。
    outlet_nodes = false(nv,1);
    if nargin >= 4 && ~isempty(outlet_file) && exist(outlet_file, 'file')
        outlet_nodes = logical(ContourToNodes(md.mesh.x, md.mesh.y, outlet_file, 1));
        outlet_nodes = outlet_nodes(:);
        if nargin >= 3 && isfield(masks, 'is_ice')
            outlet_nodes = outlet_nodes & masks.is_ice(:);
        end
        fprintf('   Outlet contour nodes: %d\n', sum(outlet_nodes));
    end

    % -----------------------------
    % 1) 纯几何阈值：薄冰/浅基岩区域 -> base
    % -----------------------------
    threshold_base_nodes = ...
        (md.geometry.thickness(:) <= spc_cfg.thin_h) & ...
        (md.geometry.base(:) >= spc_cfg.high_base);
    threshold_base_new = threshold_base_nodes & isnan(spchead);
    spchead(threshold_base_nodes) = md.geometry.base(threshold_base_nodes);
    fprintf(['   Threshold spchead=base nodes: candidates=%d, newly prescribed=%d ' ...
             '(H<=%.1f m, base>=%.1f m)\n'], ...
        sum(threshold_base_nodes), sum(threshold_base_new), ...
        spc_cfg.thin_h, spc_cfg.high_base);

    % -----------------------------
    % 2) 海洋边界 -> 0
    % -----------------------------
    if isprop(md,'mask') && isprop(md.mask,'ocean_levelset') && ~isempty(md.mask.ocean_levelset)
        ocean_nodes = on_boundary & (md.mask.ocean_levelset(:) < 0);
        spchead(ocean_nodes) = 0;
    end

    % 显式 outlet 最后设为海平面水头 head=0。
    spchead(outlet_nodes) = 0;
    fprintf('   Outlet spchead=0 nodes: %d\n', sum(outlet_nodes));

    % -----------------------------
    % 3) 冰架节点 -> 0
    % 不要求在边界上
    % -----------------------------
    shelf_nodes = false(nv,1);

    % 优先用 BedMachine mask 识别冰架
    if nargin >= 3 && isfield(masks, 'is_shelf')
        shelf_nodes = masks.is_shelf(:);
    elseif nargin >= 2 && ~isempty(bedmachine_file) && exist(bedmachine_file,'file')
        bm_mask = interpBedmachineAntarctica(md.mesh.x, md.mesh.y, 'mask', 'nearest', bedmachine_file);
        shelf_nodes = (bm_mask == 3);
    else
        % 如果没有 BedMachine，可退化为：
        % 海洋侧且有冰的点，近似视为浮冰/冰架
        if isprop(md,'mask') && isprop(md.mask,'ocean_levelset') && ~isempty(md.mask.ocean_levelset) ...
           && isprop(md.mask,'ice_levelset')   && ~isempty(md.mask.ice_levelset)
            shelf_nodes = (md.mask.ocean_levelset(:) < 0) & (md.mask.ice_levelset(:) < 0);
        end
    end

    spchead(shelf_nodes) = 0;

    % -----------------------------
    % 4) 如果海洋边界 + 冰架 一个都没找到，再用无冰边界 -> base
    % -----------------------------
    if ~any(~isnan(spchead)) && isprop(md,'mask') && isprop(md.mask,'ice_levelset') && ~isempty(md.mask.ice_levelset)
        icefree_nodes = on_boundary & (md.mask.ice_levelset(:) > 0);
        spchead(icefree_nodes) = md.geometry.base(icefree_nodes);
    end

    % -----------------------------
    % 5) 如果还没有，再用薄冰边界 -> base
    % -----------------------------
    if ~any(~isnan(spchead))
        hlim = 30;
        thin_nodes = on_boundary & (md.geometry.thickness(:) <= hlim);
        spchead(thin_nodes) = md.geometry.base(thin_nodes);
    end
end


function out = expand_nodes_by_elements(md, seed_nodes, nrings)
    out = logical(seed_nodes(:));
    nrings = max(0, round(nrings));
    for k = 1:nrings
        touched_elements = any(out(md.mesh.elements), 2);
        out(md.mesh.elements(touched_elements,:)) = true;
    end
end


function plot_spchead_diagnostics(md, masks, step4_fig_dir)
    spc = md.hydrology.spchead(:);
    finite_spc = isfinite(spc);

    tol_zero = 1e-9;
    tol_base = 1e-6;
    spc_zero = finite_spc & abs(spc) <= tol_zero;
    spc_base = finite_spc & abs(spc - md.geometry.base(:)) <= tol_base & ~spc_zero;
    spc_other = finite_spc & ~spc_zero & ~spc_base;

    fprintf('   spchead=0 nodes: %d\n', sum(spc_zero));
    fprintf('   spchead=base nodes: %d\n', sum(spc_base));
    fprintf('   spchead=other finite nodes: %d\n', sum(spc_other));

    x = md.mesh.x(:);
    y = md.mesh.y(:);

    if nargin >= 2 && isfield(masks, 'is_shelf') && isfield(masks, 'is_grounded') && isfield(masks, 'is_noice')
        bg = zeros(md.mesh.numberofvertices, 1);
        bg(masks.is_noice(:)) = 1;
        bg(masks.is_grounded(:)) = 2;
        bg(masks.is_shelf(:)) = 3;
        bg_label = '节点类型: 1=noice, 2=grounded, 3=shelf';
    else
        bg = md.geometry.thickness(:);
        bg_label = '冰厚 H [m]';
    end

    fig = figure('Name', 'STEP 4 spchead diagnostics', ...
        'Color', 'w', 'Position', [80 80 1350 850]);
    tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    draw_spchead_panel(md, bg, bg_label, spc_zero, spc_base, spc_other);
    title('spchead 位置总览', 'Interpreter', 'none');

    nexttile;
    draw_spchead_panel(md, md.geometry.thickness(:), '冰厚 H [m]', spc_zero, spc_base, spc_other);
    title('spchead 叠加在冰厚上', 'Interpreter', 'none');

    if ~exist(step4_fig_dir, 'dir'), mkdir(step4_fig_dir); end
    saveas(fig, fullfile(step4_fig_dir, 'Recovery_step4_spchead_diagnostics.png'));
end


function validate_spchead_threshold(md, spc_cfg)
    spc = md.hydrology.spchead(:);
    threshold_nodes = ...
        (md.geometry.thickness(:) <= spc_cfg.thin_h) & ...
        (md.geometry.base(:) >= spc_cfg.high_base);

    threshold_finite = threshold_nodes & isfinite(spc);
    threshold_missing = threshold_nodes & ~isfinite(spc);
    threshold_zero = threshold_finite & abs(spc) <= 1e-9;
    threshold_base = threshold_finite & abs(spc - md.geometry.base(:)) <= 1e-6 & ~threshold_zero;
    threshold_other = threshold_finite & ~threshold_zero & ~threshold_base;

    fprintf(['   Threshold coverage check: candidates=%d, finite=%d, ' ...
             'base=%d, zero=%d, other=%d, missing=%d\n'], ...
        sum(threshold_nodes), sum(threshold_finite), sum(threshold_base), ...
        sum(threshold_zero), sum(threshold_other), sum(threshold_missing));

    if any(threshold_missing)
        ids = find(threshold_missing);
        show = ids(1:min(10, numel(ids)));
        fprintf('   First missing threshold spchead nodes:\n');
        for k = 1:numel(show)
            id = show(k);
            fprintf('      node %6d  H=%.2f  base=%.2f  x=%.1f  y=%.1f\n', ...
                id, md.geometry.thickness(id), md.geometry.base(id), md.mesh.x(id), md.mesh.y(id));
        end
    end

    assert(~any(threshold_missing), ...
        'Some H/base threshold nodes do not have finite spchead. Check build_spchead ordering.');
end


function draw_spchead_panel(md, bg, bg_label, spc_zero, spc_base, spc_other)
    patch('Faces', md.mesh.elements, ...
          'Vertices', [md.mesh.x(:), md.mesh.y(:)], ...
          'FaceVertexCData', double(bg(:)), ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    axis equal tight; box on; hold on;
    xlabel('X (m)');
    ylabel('Y (m)');
    colorbar;
    title(bg_label, 'Interpreter', 'none');

    plot(md.mesh.x(spc_zero), md.mesh.y(spc_zero), 'bo', ...
        'MarkerSize', 4, 'LineWidth', 1);
    plot(md.mesh.x(spc_base), md.mesh.y(spc_base), 'ro', ...
        'MarkerSize', 4, 'LineWidth', 1);
    if any(spc_other)
        plot(md.mesh.x(spc_other), md.mesh.y(spc_other), 'ko', ...
            'MarkerSize', 5, 'LineWidth', 1.2);
    end

    if any(spc_other)
        legend({'background', 'spchead=0', 'spchead=base', 'spchead=other'}, ...
            'Location', 'bestoutside');
    else
        legend({'background', 'spchead=0', 'spchead=base'}, ...
            'Location', 'bestoutside');
    end
end
