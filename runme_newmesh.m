 % =========================================================================
%  runme_newmesh.m
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
if ~exist(paths.models_newmesh, 'dir'), mkdir(paths.models_newmesh); end
if ~exist(paths.figures_newmesh, 'dir'), mkdir(paths.figures_newmesh); end

% -------------------------------------------------------------------------
%  0. 全局配置
% -------------------------------------------------------------------------
% 常用:
%   steps = [5];  Generate RecoveryNewMesh_Hydrology_ShelfGroundedRingNoSlidingCoupled.mat
%   steps = [6];  基于 Hydrology 模型进行瞬态水文模拟
%   steps = [7];  接地线一层接地冰 spchead 敏感性实验
steps = [7];
org   = organizer('repository', paths.models_newmesh, 'prefix', 'RecoveryNewMesh_', 'steps', steps);
clear steps

% --- 路径 ---
domain_file    = fullfile(paths.exp, 'Recovery_all.exp');
velocity_file  = paths.newmesh_velocity_file;
bedmachine_file = paths.bedmachine_file;
par_file       = fullfile(paths.par, 'Recovery_newmesh_bedmachinev4.par');
outlet_exp     = fullfile(paths.exp, 'Recovery_outlets.exp');

% --- 初始化模型路径 ---
% 非空且文件存在：直接加载该 .mat 进入 step 2
% 留空 ''       ：退回到 loadmodel(org, 'Mesh')（需要 step 1 已经产出）
init_model_path = '';

% --- 网格参数 ---
mesh_hinit = 2500;   % 初始均匀尺寸 (m)
mesh_hmin  = 1000;   % 自适应最小尺寸 (m)
mesh_hmax  = 4000;  % 自适应最大尺寸 (m)
mesh_err   = 0.5;    % 插值误差容忍度 (m/yr)

% --- 反演参数 ---
inv_maxsteps      = 200;
inv_maxiter       = 400;
inv_coeff_abs = 5000;
inv_coeff_log = 1.5e-3;
inv_coeff_reg = 1.5e-5;
inv_friction_min  = 0.05;
inv_friction_max  = 3500;

% --- 水文参数 ---
hydro_head_fraction   = 0.8;   % 初始水头 = base + fraction * (rho_i/rho_w) * H
hydro_gap_height_init = 0.01; % 初始 gap (m)
hydro_bump_spacing    = 2.0;   % bump 间距 (m), SHAKTI 论文示例 lr=2 m
hydro_bump_height     = 0.1;   % bump 高度 (m), SHAKTI 论文示例 br=0.1 m；0 = 不考虑 opening-by-sliding
hydro_min_thickness   = 25;    % 最小冰厚安全阈值 (m)
hydro_storage         = 0;  % englacial storage coefficient; damps coarse-grid head jumps
hydro_relaxation      = 1   % under-relaxation for nonlinear Shakti head updates
hydro_open_spc_thickness = 300;  % 高基岩开放排水点的薄冰阈值 (m)
hydro_open_spc_base_min  = -5000;  % 高基岩开放排水点的基岩高程阈值 (m)
hydro_open_spc_rings     = 1;    % unused; current spchead=base uses only H/base thresholds
% --- Step 6: Simulation settings ---
sim_start_time    = 0.0;
sim_dt_seconds    = 30 * 60;     % 30 min time step (s)
sim_n_steps       = 48 * 15;     % diagnostic step count
sim_final_days    = sim_n_steps * sim_dt_seconds / 86400;
sim_output_freq   = 1;           % save every step
sim_stress_restol = 0.05;
sim_stress_reltol = 0.05;
sim_stress_abstol = NaN;
sim_stress_maxiter = 200;
sim_verbose_solution = 1;
sim_verbose_module = 0;
sim_verbose_convergence = 0;
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
    [x_sorted, y_sorted, vx_grid, vy_grid] = read_velocity_grid(velocity_file);

    % 清理 → 转置 → 按排序索引重排
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
    md.miscellaneous.name = 'RecoveryNewMesh';

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
    lcurve_model_path = fullfile(paths.models_newmesh, 'RecoveryNewMesh_Inversion.mat');
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
    plot_spchead_diagnostics(md, masks);
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
    assert(all(isfinite(md.basalforcings.geothermalflux)) && all(md.basalforcings.geothermalflux >= 0), ...
        'md.basalforcings.geothermalflux contains NaN/Inf or negative values.');

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
%  STEP 5. Hydrology_ShelfGroundedRingNoSlidingCoupled
%          Reset spchead, set shelf + one grounded ring to head=0,
%          enable stressbalance-hydrology coupling, and keep bump_height=0.
% =========================================================================
% --- Step 5 settings ---
step5_model_name = 'Hydrology_ShelfGroundedRingNoSlidingCoupled';
step5_disable_sliding_opening = true;
step5_enable_two_way_coupling = 0;
step5_friction_coupling = 4;

if perform(org, step5_model_name) % {{{

    fprintf('\n===== STEP 5: Reset spchead; shelf + grounded ring set to head=0 =====\n');
    md = loadmodel(org, 'Hydrology');

    [md, spc_report] = reset_all_spchead(md);
    fprintf('   Original finite spchead nodes: %d\n', spc_report.n_finite_before);
    fprintf('   Cleared spchead=0 nodes inherited from Step 4: %d\n', spc_report.n_zero_before);
    fprintf('   Cleared spchead=base nodes inherited from Step 4: %d\n', spc_report.n_base_before);
    fprintf('   Cleared other finite spchead nodes inherited from Step 4: %d\n', spc_report.n_other_before);
    fprintf('   Remaining finite spchead nodes after reset: %d\n', spc_report.n_finite_after);

    masks = classify_bedmachine_nodes(md, bedmachine_file);
    [md, shelf_report] = set_shelf_plus_grounded_ring_spchead_zero(md, masks);
    fprintf('   Shelf nodes set to spchead=0: %d\n', shelf_report.n_shelf_nodes);
    fprintf('   One-ring grounded nodes set to spchead=0: %d\n', shelf_report.n_grounded_ring_nodes);
    fprintf('   Total target nodes set to spchead=0: %d\n', shelf_report.n_target_nodes);
    fprintf('   Newly prescribed target nodes: %d\n', shelf_report.n_new_zero);
    fprintf('   Target nodes already spchead=0: %d\n', shelf_report.n_already_zero);
    fprintf('   Total finite spchead nodes after target update: %d\n', shelf_report.n_finite_after);
    fprintf('   Total spchead=0 nodes after target update: %d\n', shelf_report.n_zero_after);

    assert(shelf_report.n_shelf_nodes > 0, ...
        'No shelf nodes were identified. Check BedMachine mask classification.');
    assert(shelf_report.n_zero_after > 0, ...
        'No spchead=0 nodes remain after setting target nodes.');
    assert(shelf_report.n_finite_after > 0, ...
        'md.hydrology.spchead is all NaN after Step 5 updates.');

    if step5_disable_sliding_opening
        md.hydrology.bump_height = zeros(size(md.hydrology.bump_height));
        fprintf('   Sliding opening disabled: hydrology.bump_height set to 0 on all elements.\n');
    end

    if step5_enable_two_way_coupling
        md.friction.coupling = step5_friction_coupling;
        md.transient.ishydrology = 1;
        md.transient.isstressbalance = 1;
        fprintf('   Two-way coupling enabled: transient hydrology + stressbalance, friction.coupling=%d.\n', ...
            md.friction.coupling);
    end

    plot_step5_spchead_outlet_diagnostics(md, masks, outlet_exp);
    savemodel(org, md);

end % }}}

% =========================================================================
%  STEP 6. Simulation - 瞬态水文模拟
% =========================================================================
% --- Step 6 settings ---
% Set to an organizer model name, for example:
%   'Hydrology'
%   'Hydrology_ShelfGroundedRingNoSlidingCoupled'
step6_input_model_name = 'Hydrology_ShelfGroundedRingNoSlidingCoupled';
% spchead=0 setup before solving:
%   'model'                 keep spchead from the input model
%   'outlet_exp'            clear spchead, then set outlet EXP nodes to 0
%   'model_plus_outlet_exp' keep input spchead, and also set outlet EXP nodes to 0
step6_spchead_zero_mode = 'outlet_exp';
step6_spchead_outlet_exp = outlet_exp;
% Output model name in Models/. Leave empty '' to auto-generate from steps/dt.
% You can include or omit the .mat suffix.
step6_output_model_name = 'RecoveryNewMesh_SHAKTI_720Steps_1800s_outlets_coupled';

if perform(org, 'Simulation') % {{{

    fprintf('\n===== STEP 6: Transient Simulation =====\n');
    sim_result_name = resolve_step6_output_model_name(step6_output_model_name, ...
        sim_n_steps, sim_dt_seconds);
    fprintf('   Output model: %s.mat\n', sim_result_name);
    fprintf('   Loading input model: RecoveryNewMesh_%s.mat\n', step6_input_model_name);
    md = loadmodel(org, step6_input_model_name);
    [md, step6_spc_report] = apply_step6_spchead_zero_mode(md, ...
        step6_spchead_zero_mode, step6_spchead_outlet_exp);
    fprintf(['   Step 6 spchead mode: %s | finite before=%d, finite after=%d, ' ...
             'outlet nodes set to 0=%d\n'], ...
        step6_spc_report.mode, step6_spc_report.n_finite_before, ...
        step6_spc_report.n_finite_after, step6_spc_report.n_outlet_zero);

    % 5.2 时间步长
    md.timestepping.start_time = sim_start_time;
    md.timestepping.time_step  = sim_dt_seconds / md.constants.yts;
    md.timestepping.final_time = sim_final_days / 365;
    assert(sim_n_steps >= 1 && mod(sim_n_steps, 1) == 0, 'sim_n_steps must be a positive integer.');
    fprintf('   Time step: %.1f s, number of steps: %d, final_time: %.9g yr\n', ...
        sim_dt_seconds, sim_n_steps, md.timestepping.final_time);

    %md.hydrology.melt_flag = 0;


    % 5.3 输出 & 求解器
    md.settings.output_frequency = sim_output_freq;
    md.stressbalance.restol  = sim_stress_restol;
    md.stressbalance.reltol  = sim_stress_reltol;
    md.stressbalance.abstol  = sim_stress_abstol;
    md.stressbalance.maxiter = sim_stress_maxiter;

    % 5.4 日志
    md.verbose.solution    = sim_verbose_solution;
    md.verbose.module      = sim_verbose_module;
    md.verbose.convergence = sim_verbose_convergence;

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
    diagnostic_file = fullfile(paths.models_newmesh, [sim_result_name '.mat']);
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


% =========================================================================
%  STEP 7. 接地线一层接地冰 spchead 敏感性实验
%          冰架节点固定 head=0；冰架相邻的接地冰 ring 使用
%          head = bed + alpha * rho_i/rho_w * H。
%          alpha = Pw/Pi，alpha 越大越接近浮托，排水越弱。
% =========================================================================
% --- Step 7 settings ---
step7_step_name = 'GroundingRingHeadSensitivity';
step7_input_model_name = 'Hydrology_ShelfGroundedRingNoSliding';
step7_alpha_list = [0.85 0.9 0.95 0.99];

% 是否额外加入一组“接地 ring 也固定 head=0”的基准实验。
% true  = 在 alpha 敏感性实验之外，先生成/运行 Head0，用来和原始强排水边界对照；
% false = 只运行 step7_alpha_list 中的冰覆压比例实验。
step7_include_head0_baseline = 0;

% 从冰架节点向接地冰内部扩展的层数。
% 1 表示只选“紧邻冰架的一层接地冰节点”作为过渡水头边界；
% 数值越大，被固定 spchead 的接地冰范围越宽，边界约束越强。
step7_grounded_ring_count = 1;

% 冰架/海洋侧节点的固定水头值。
% 0 表示海平面水头基准；这些节点作为海洋连通边界处理，
% 不使用 alpha 公式计算。
step7_shelf_head_value = 0;

% --- Step 7 两个小步开关 ---
% step7_save_setup_outputs:
%   true  = 保存每个 case（Head0 和 alpha）的 setup 模型、诊断图和 CSV；
%   false = 不保存 setup/诊断，只在内存中构造边界条件。
% 注意：即使设为 false，只要 step7_run_transient=true，程序仍会在内存中
% 生成当前 case 的 md_case，因为 transient 求解必须基于这个边界条件。
step7_save_setup_outputs = 1;

% step7_run_transient:
%   true  = 对每个 case 挨个运行 transient，并保存瞬态结果模型；
%   false = 不运行 transient。
step7_run_transient = 0;
step7_output_root = fullfile(paths.figures_newmesh, 'step7_grounded_ring_head_sensitivity');

% --- Step 7 瞬态模拟设置 ---
% 默认继承 Step 6 的设置；如果只想调整第七步，在这里改 step7_* 变量即可。
% step7_sim_dt_seconds: 每个水文时间步长度（秒）
% step7_sim_n_steps:    总步数；总时长 = step7_sim_n_steps * step7_sim_dt_seconds
% step7_sim_output_freq:每隔多少步保存一次结果，1 表示每步都保存
step7_sim_start_time = 0;
step7_sim_dt_seconds = 30 * 60;
step7_sim_n_steps = 48 * 60;
step7_sim_final_days = step7_sim_n_steps * step7_sim_dt_seconds / 86400;
step7_sim_output_freq = 12;

% SHAKTI 非线性水头求解松弛系数。
% 1 与之前跑通的 cluster 脚本一致；NaN 表示保留输入模型中的原值。
step7_hydro_relaxation = 1;

% --- Step 7 stressbalance 收敛设置 ---
% 当前第七步输入模型通常 isstressbalance=0；如果以后打开耦合，
% 这些阈值会传给 stressbalance 求解器。
step7_stress_restol = sim_stress_restol;
step7_stress_reltol = sim_stress_reltol;
step7_stress_abstol = sim_stress_abstol;
step7_stress_maxiter = sim_stress_maxiter;

% --- Step 7 日志和并行设置 ---
step7_verbose_solution = sim_verbose_solution;
step7_verbose_module = sim_verbose_module;
step7_verbose_convergence = sim_verbose_convergence;
step7_n_cores = n_cores;

if perform(org, step7_step_name) % {{{

    fprintf('\n===== STEP 7: Grounded-ring spchead sensitivity =====\n');
    assert(step7_save_setup_outputs || step7_run_transient, ...
        'Step 7 has nothing to do: enable step7_save_setup_outputs or step7_run_transient.');

    fprintf('   Loading input model: RecoveryNewMesh_%s.mat\n', step7_input_model_name);
    md_base = loadmodel(org, step7_input_model_name);

    if step7_save_setup_outputs && ~exist(step7_output_root, 'dir')
        mkdir(step7_output_root);
    end

    % 用 BedMachine 分类重新识别冰架和接地冰节点，保证 Step 7
    % 与 Step 5 使用同一套 shelf/grounded ring 定义。
    masks = classify_bedmachine_nodes(md_base, bedmachine_file);
    shelf_nodes = logical(masks.is_shelf(:));
    grounded_nodes = logical(masks.is_grounded(:));

    % 从冰架节点向外扩展指定层数，然后只保留其中的接地冰节点。
    % ring_count=1 表示“紧邻冰架的一层接地冰节点”。
    expanded_from_shelf = expand_nodes_by_elements(md_base, shelf_nodes, step7_grounded_ring_count);
    grounded_ring_nodes = expanded_from_shelf & grounded_nodes;

    assert(any(shelf_nodes), 'Step 7 did not identify any shelf nodes.');
    assert(any(grounded_ring_nodes), 'Step 7 did not identify any grounded ring nodes.');

    fprintf('   Shelf nodes fixed to head=%.3g: %d\n', step7_shelf_head_value, sum(shelf_nodes));
    fprintf('   Grounded ring count: %d, grounded ring nodes: %d\n', ...
        step7_grounded_ring_count, sum(grounded_ring_nodes));
    fprintf('   Step 7 save setup outputs: %d\n', step7_save_setup_outputs);
    fprintf('   Step 7 run transient: %d\n', step7_run_transient);
    if step7_run_transient
        fprintf('   Step 7 transient dt: %.1f s, steps: %d, final time: %.3f days\n', ...
            step7_sim_dt_seconds, step7_sim_n_steps, step7_sim_final_days);
        fprintf('   Step 7 hydrology relaxation: %.4g\n', step7_hydro_relaxation);
        fprintf('   Step 7 stressbalance tolerances: restol=%.4g, reltol=%.4g, abstol=%.4g, maxiter=%d\n', ...
            step7_stress_restol, step7_stress_reltol, step7_stress_abstol, step7_stress_maxiter);
    end

    rho_i = md_base.materials.rho_ice;
    rho_w = md_base.materials.rho_freshwater;

    step7_summary_rows = struct([]);

    n_step7_cases = double(logical(step7_include_head0_baseline)) + numel(step7_alpha_list);
    assert(n_step7_cases > 0, ...
        'Step 7 has no experiment cases: enable step7_include_head0_baseline or set step7_alpha_list.');
    step7_case_list = repmat(struct('case_name', '', 'alpha', NaN, ...
        'use_head0_baseline', false), 1, n_step7_cases);
    step7_case_index = 0;

    if step7_include_head0_baseline
        step7_case_index = step7_case_index + 1;
        step7_case_list(step7_case_index) = struct( ...
            'case_name', 'Head0', ...
            'alpha', NaN, ...
            'use_head0_baseline', true);
    end
    for ia = 1:numel(step7_alpha_list)
        alpha_i = step7_alpha_list(ia);
        assert(alpha_i > 0 && alpha_i <= 1, 'Each Step 7 alpha must satisfy 0 < alpha <= 1.');
        alpha_name_i = step7_alpha_name(alpha_i);
        step7_case_index = step7_case_index + 1;
        step7_case_list(step7_case_index) = struct( ...
            'case_name', ['Alpha' alpha_name_i], ...
            'alpha', alpha_i, ...
            'use_head0_baseline', false);
    end

    for icase = 1:numel(step7_case_list)
        case_info = step7_case_list(icase);
        alpha = case_info.alpha;

        % alpha 是 Pw/Pi，而不是 N/Pi。
        % alpha 越大，水压越接近冰覆压，接地线排水边界越弱。
        % N/Pi = 1 - alpha。
        if case_info.use_head0_baseline
            alpha_name = 'Head0';
            setup_model_name = 'Hydrology_Shelf0GroundedRingHead0_NoSliding';
        else
            alpha_name = step7_alpha_name(alpha);
            setup_model_name = sprintf('Hydrology_Shelf0GroundedRingAlpha%s_NoSliding', alpha_name);
        end
        setup_file = fullfile(paths.models_newmesh, ['RecoveryNewMesh_' setup_model_name '.mat']);
        case_outdir = fullfile(step7_output_root, case_info.case_name);
        if step7_save_setup_outputs && ~exist(case_outdir, 'dir')
            mkdir(case_outdir);
        end

        md_case = md_base;
        if ~isnan(step7_hydro_relaxation)
            assert(step7_hydro_relaxation >= 0, 'step7_hydro_relaxation must be nonnegative or NaN.');
            if has_model_member(md_case.hydrology, 'relaxation')
                md_case.hydrology.relaxation = step7_hydro_relaxation;
            else
                warning('md_case.hydrology.relaxation is missing; cannot set Step 7 relaxation.');
            end
        end

        H = md_case.geometry.thickness(:);
        bed = md_case.geometry.base(:);

        % 接地冰 ring 的固定水头：
        %   head = bed + alpha * rho_i/rho_w * H
        % 等价于 Pw/Pi=alpha，N/Pi=1-alpha。
        if case_info.use_head0_baseline
            % Head0 基准实验：直接把扩展的一层接地冰 ring 固定为 head=0，
            % 用来复现“接地线过渡节点排水很强”的原始边界条件。
            head_ring = zeros(size(bed));
        else
            head_ring = bed + alpha .* (rho_i ./ rho_w) .* H;
        end

        % 只固定两类节点：
        %   1) 冰架节点：head=0，代表海洋/冰架侧水头基准；
        %   2) 接地冰 ring：Head0 case 直接固定为 0；alpha case 使用上面按比例计算的过渡水头。
        % 其他接地冰节点保持 NaN，让 SHAKTI 自己求解。
        spc = NaN(md_case.mesh.numberofvertices, 1);
        spc(shelf_nodes) = step7_shelf_head_value;
        spc(grounded_ring_nodes) = head_ring(grounded_ring_nodes);
        md_case.hydrology.spchead = spc;

        % 让固定水头节点的初始 head 与 spchead 一致，避免初始旧水头和强约束不一致。
        finite_spc = isfinite(spc);
        md_case.hydrology.head(finite_spc) = spc(finite_spc);

        md_case.miscellaneous.name = ['RecoveryNewMesh_' setup_model_name];
        step7_report = summarize_step7_case(md_case, masks, shelf_nodes, grounded_ring_nodes, ...
            alpha, case_info.case_name, case_info.use_head0_baseline);

        if case_info.use_head0_baseline
            fprintf(['\n   Head0 baseline: shelf head=%.3g, grounded ring spchead=0, ' ...
                     'finite spchead=%d, grounded ring nodes=%d\n'], ...
                step7_shelf_head_value, step7_report.n_finite_spchead, ...
                step7_report.n_grounded_ring_nodes);
            fprintf('      Equivalent grounded ring Pw/Pi median %.3f, N/Pi median %.3f\n', ...
                step7_report.ring_Pw_over_Pi_median, step7_report.ring_N_over_Pi_median);
        else
            fprintf(['\n   Alpha %s: Pw/Pi=%.3f, N/Pi=%.3f, finite spchead=%d, ' ...
                     'grounded ring nodes=%d\n'], ...
                alpha_name, alpha, 1-alpha, step7_report.n_finite_spchead, ...
                step7_report.n_grounded_ring_nodes);
        end
        fprintf('      Grounded ring spchead median %.3f m, min %.3f m, max %.3f m\n', ...
            step7_report.ring_spchead_median, step7_report.ring_spchead_min, ...
            step7_report.ring_spchead_max);

        if step7_save_setup_outputs
            % 小步 7A：保存每个 case 的 setup 模型；变量名保持为 md，
            % 方便后续 load(model,'md')。同时保存诊断图和 CSV。
            md = md_case;
            save(setup_file, 'md', 'step7_report', '-v7.3');
            clear md
            fprintf('      Saved setup model: %s\n', setup_file);

            % 每个 case 都保存一套诊断图和 CSV，先检查边界条件是否合理，
            % 再决定是否打开 step7_run_transient 进行长时间模拟。
            save_step7_diagnostics(md_case, masks, shelf_nodes, grounded_ring_nodes, ...
                alpha, case_info.case_name, case_info.use_head0_baseline, case_outdir, step7_report);
        else
            fprintf('      Setup output saving disabled; using this case setup in memory only.\n');
        end

        if step7_save_setup_outputs
            if isempty(step7_summary_rows)
                step7_summary_rows = step7_report;
            else
                step7_summary_rows(end+1) = step7_report; %#ok<SAGROW>
            end
        end

        if step7_run_transient
            % 小步 7B：直接对当前 case 的 setup 模型跑 transient。
            % 该分支与 setup 输出保存相互独立；即使不保存 setup 文件，
            % 这里也会使用上面刚刚在内存中构造好的 md_case。
            result_model_name = sprintf('SHAKTI_%s_%dSteps_%ds_noslide', ...
                case_info.case_name, step7_sim_n_steps, step7_sim_dt_seconds);
            result_file = fullfile(paths.models_newmesh, ['RecoveryNewMesh_' result_model_name '.mat']);
            md_result = run_step7_transient_case(md_case, result_model_name, ...
                step7_sim_start_time, step7_sim_dt_seconds, step7_sim_n_steps, step7_sim_output_freq, ...
                step7_stress_restol, step7_stress_reltol, step7_stress_abstol, ...
                step7_stress_maxiter, step7_verbose_solution, step7_verbose_module, ...
                step7_verbose_convergence, step7_n_cores);
            md = md_result;
            save(result_file, 'md', '-v7.3');
            clear md
            fprintf('      Saved transient result model: %s\n', result_file);
        end
    end

    if step7_save_setup_outputs
        summary_file = fullfile(step7_output_root, 'step7_sensitivity_summary.csv');
        writetable(struct2table(step7_summary_rows), summary_file);
        fprintf('\n   Step 7 summary saved: %s\n', summary_file);
        fprintf('   Step 7 diagnostics folder: %s\n', step7_output_root);
    else
        fprintf('\n   Step 7 setup outputs were disabled; no setup models or diagnostics were saved.\n');
    end
    fprintf('Done.\n');
end % }}}




function [x_sorted, y_sorted, vx_grid, vy_grid] = read_velocity_grid(velocity_file)

    x_raw = double(ncread(velocity_file, 'x'));
    y_raw = double(ncread(velocity_file, 'y'));
    [x_sorted, ix] = sort(x_raw);
    [y_sorted, iy] = sort(y_raw);

    vx_raw = double(read_first_nc_variable(velocity_file, {'VX', 'vx', 'Vx'}));
    vy_raw = double(read_first_nc_variable(velocity_file, {'VY', 'vy', 'Vy'}));

    vx_raw(~isfinite(vx_raw) | abs(vx_raw) > 1e5) = 0;
    vy_raw(~isfinite(vy_raw) | abs(vy_raw) > 1e5) = 0;

    if size(vx_raw, 1) == numel(x_raw) && size(vx_raw, 2) == numel(y_raw)
        vx_grid = vx_raw';
        vy_grid = vy_raw';
    elseif size(vx_raw, 1) == numel(y_raw) && size(vx_raw, 2) == numel(x_raw)
        vx_grid = vx_raw;
        vy_grid = vy_raw;
    else
        error('Velocity grid size does not match x/y dimensions in %s.', velocity_file);
    end

    vx_grid = vx_grid(iy, ix);
    vy_grid = vy_grid(iy, ix);
end


function data = read_first_nc_variable(nc_file, names)
    for i = 1:numel(names)
        try
            data = ncread(nc_file, names{i});
            return;
        catch
        end
    end
    error('None of these variables were found in %s: %s', nc_file, strjoin(names, ', '));
end


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
    outlet_nodes = false(nv,1);

    if nargin < 5 || isempty(spc_cfg)
        spc_cfg.thin_h     = 50;
        spc_cfg.high_base  = 400;
        spc_cfg.open_rings = 1;
    end

    % 先读取显式 outlet；真正赋值放在海洋边界之后，用 head=0 覆盖。
    outlet_nodes = false(nv,1);
    if nargin >= 4 && ~isempty(outlet_file) && exist(outlet_file, 'file')
        outlet_nodes = logical(ContourToNodes(md.mesh.x, md.mesh.y, outlet_file, 1));
        outlet_nodes = outlet_nodes(:);
        if nargin >= 3 && isfield(masks, 'is_ice')
            outlet_nodes = outlet_nodes & masks.is_ice(:);
        end
        fprintf('   Outlet contour nodes: %d\n', sum(outlet_nodes));
    end

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
    % 1) 海洋边界 -> 0
    % -----------------------------
    if isprop(md,'mask') && isprop(md.mask,'ocean_levelset') && ~isempty(md.mask.ocean_levelset)
        ocean_nodes = on_boundary & (md.mask.ocean_levelset(:) < 0);
        spchead(ocean_nodes) = 0;
    end

    % 显式 outlet 最后设为海平面水头 head=0。
    spchead(outlet_nodes) = 0;
    fprintf('   Outlet spchead=0 nodes: %d\n', sum(outlet_nodes));

    % -----------------------------
    % 2) 冰架节点 -> 0
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
    fprintf('   Shelf spchead=0 nodes: %d\n', sum(shelf_nodes));

    if ~any(~isnan(spchead)) && isprop(md,'mask') && isprop(md.mask,'ice_levelset') && ~isempty(md.mask.ice_levelset)
        icefree_nodes = on_boundary & (md.mask.ice_levelset(:) > 0);
        spchead(icefree_nodes) = md.geometry.base(icefree_nodes);
    end

    if ~any(~isnan(spchead))
        hlim = 30;
        thin_nodes = on_boundary & (md.geometry.thickness(:) <= hlim);
        spchead(thin_nodes) = md.geometry.base(thin_nodes);
    end
end


function [md, report] = remove_atmospheric_spchead(md)
    spc = md.hydrology.spchead(:);
    base = md.geometry.base(:);

    tol_zero = 1e-9;
    tol_base = 1e-6;
    finite_spc = isfinite(spc);
    spc_zero = finite_spc & abs(spc) <= tol_zero;
    spc_base = finite_spc & abs(spc - base) <= tol_base & ~spc_zero;
    spc_other = finite_spc & ~spc_zero & ~spc_base;

    report.n_finite_before = sum(finite_spc);
    report.n_zero_kept = sum(spc_zero);
    report.n_base_removed = sum(spc_base);
    report.n_other_kept = sum(spc_other);

    spc(spc_base) = NaN;
    md.hydrology.spchead = spc;

    report.n_finite_after = sum(isfinite(md.hydrology.spchead(:)));
end


function [md, report] = reset_all_spchead(md)
    spc = md.hydrology.spchead(:);
    base = md.geometry.base(:);

    tol_zero = 1e-9;
    tol_base = 1e-6;
    finite_spc = isfinite(spc);
    spc_zero = finite_spc & abs(spc) <= tol_zero;
    spc_base = finite_spc & abs(spc - base) <= tol_base & ~spc_zero;
    spc_other = finite_spc & ~spc_zero & ~spc_base;

    report.n_finite_before = sum(finite_spc);
    report.n_zero_before = sum(spc_zero);
    report.n_base_before = sum(spc_base);
    report.n_other_before = sum(spc_other);

    md.hydrology.spchead = NaN(md.mesh.numberofvertices, 1);
    report.n_finite_after = sum(isfinite(md.hydrology.spchead(:)));
end


function validate_no_atmospheric_spchead(md)
    spc = md.hydrology.spchead(:);
    base = md.geometry.base(:);

    tol_zero = 1e-9;
    tol_base = 1e-6;
    finite_spc = isfinite(spc);
    spc_zero = finite_spc & abs(spc) <= tol_zero;
    spc_base_nonzero = finite_spc & abs(spc - base) <= tol_base & ~spc_zero;

    assert(~any(spc_base_nonzero), ...
        'Nonzero spchead=base nodes remain after removing atmospheric spchead.');
end


function [md, report] = set_shelf_plus_grounded_ring_spchead_zero(md, masks)
    spc = md.hydrology.spchead(:);
    assert(isfield(masks, 'is_shelf'), 'masks.is_shelf is required.');
    assert(isfield(masks, 'is_grounded'), 'masks.is_grounded is required.');
    shelf_nodes = logical(masks.is_shelf(:));
    grounded_nodes = logical(masks.is_grounded(:));
    assert(numel(shelf_nodes) == md.mesh.numberofvertices, ...
        'masks.is_shelf length does not match md.mesh.numberofvertices.');
    assert(numel(grounded_nodes) == md.mesh.numberofvertices, ...
        'masks.is_grounded length does not match md.mesh.numberofvertices.');

    expanded_from_shelf = expand_nodes_by_elements(md, shelf_nodes, 1);
    grounded_ring_nodes = expanded_from_shelf & grounded_nodes;
    target_nodes = shelf_nodes | grounded_ring_nodes;

    finite_before = isfinite(spc);
    zero_before = finite_before & abs(spc) <= 1e-9;
    new_zero = target_nodes & ~zero_before;
    already_zero = target_nodes & zero_before;
    overwritten_finite_nonzero = target_nodes & finite_before & ~zero_before;
    overwritten_nan = target_nodes & ~finite_before;

    spc(target_nodes) = 0;
    md.hydrology.spchead = spc;

    finite_after = isfinite(spc);
    zero_after = finite_after & abs(spc) <= 1e-9;

    report.n_shelf_nodes = sum(shelf_nodes);
    report.n_grounded_ring_nodes = sum(grounded_ring_nodes);
    report.n_target_nodes = sum(target_nodes);
    report.n_new_zero = sum(new_zero);
    report.n_already_zero = sum(already_zero);
    report.n_overwritten_finite_nonzero = sum(overwritten_finite_nonzero);
    report.n_overwritten_nan = sum(overwritten_nan);
    report.n_finite_after = sum(finite_after);
    report.n_zero_after = sum(zero_after);
end


function [md, report] = apply_step6_spchead_zero_mode(md, mode, outlet_file)
    mode = lower(strtrim(char(mode)));
    spc = md.hydrology.spchead(:);

    report.mode = mode;
    report.n_finite_before = sum(isfinite(spc));
    report.n_outlet_nodes = 0;
    report.n_outlet_zero = 0;

    switch mode
        case 'model'
            report.n_finite_after = report.n_finite_before;
            return;

        case {'outlet_exp', 'model_plus_outlet_exp'}
            assert(exist(outlet_file, 'file') == 2, ...
                sprintf('Step 6 outlet EXP file does not exist: %s', outlet_file));
            outlet_nodes = logical(ContourToNodes(md.mesh.x, md.mesh.y, outlet_file, 1));
            outlet_nodes = outlet_nodes(:);
            assert(numel(outlet_nodes) == md.mesh.numberofvertices, ...
                'Outlet EXP node mask length does not match md.mesh.numberofvertices.');
            assert(any(outlet_nodes), ...
                sprintf('Step 6 outlet EXP selected zero nodes: %s', outlet_file));

            if strcmp(mode, 'outlet_exp')
                spc = NaN(md.mesh.numberofvertices, 1);
            end
            spc(outlet_nodes) = 0;
            md.hydrology.spchead = spc;

            report.n_outlet_nodes = sum(outlet_nodes);
            report.n_outlet_zero = sum(outlet_nodes);
            report.n_finite_after = sum(isfinite(spc));

        otherwise
            error(['Unknown step6_spchead_zero_mode: %s. Use ''model'', ' ...
                   '''outlet_exp'', or ''model_plus_outlet_exp''.'], mode);
    end
end


function name = resolve_step6_output_model_name(user_name, sim_n_steps, sim_dt_seconds)
    name = strtrim(char(user_name));
    if isempty(name)
        if sim_n_steps == 1
            name = sprintf('RecoveryNewMesh_SHAKTI_OneStep_%ds', sim_dt_seconds);
        else
            name = sprintf('RecoveryNewMesh_SHAKTI_%dSteps_%ds', sim_n_steps, sim_dt_seconds);
        end
        return;
    end

    [folder, base, ext] = fileparts(name);
    assert(isempty(folder), ...
        'step6_output_model_name should be a file name only, not a path.');
    if strcmpi(ext, '.mat')
        name = base;
    else
        name = [base ext];
    end
    assert(~isempty(name), 'step6_output_model_name cannot be empty.');
end


function out = expand_nodes_by_elements(md, seed_nodes, nrings)
    out = logical(seed_nodes(:));
    nrings = max(0, round(nrings));
    for k = 1:nrings
        touched_elements = any(out(md.mesh.elements), 2);
        out(md.mesh.elements(touched_elements,:)) = true;
    end
end


function name = step7_alpha_name(alpha)
    % 把 alpha=0.95 转成文件名里的 095，方便排序和批量识别。
    name = sprintf('%03d', round(alpha * 100));
end


function report = summarize_step7_case(md, masks, shelf_nodes, grounded_ring_nodes, alpha, case_name, use_head0_baseline)
    % 汇总当前 case 的边界条件统计量。
    % 这里的 Pw/Pi 和 N/Pi 只用于诊断固定水头节点的物理含义，
    % 不会改变模型本身。
    spc = md.hydrology.spchead(:);
    H = md.geometry.thickness(:);
    bed = md.geometry.base(:);
    rho_i = md.materials.rho_ice;
    rho_w = md.materials.rho_freshwater;
    g = md.constants.g;

    Pi = rho_i .* g .* H;
    Pw = rho_w .* g .* (spc - bed);
    N = Pi - Pw;
    Pw_over_Pi = Pw ./ Pi;
    N_over_Pi = N ./ Pi;
    head_float = bed + (rho_i ./ rho_w) .* H;
    head_float_minus_spc = head_float - spc;

    ring = grounded_ring_nodes(:);
    finite_spc = isfinite(spc);
    zero_spc = finite_spc & abs(spc) <= 1e-9;

    report = struct();
    report.case_name = {case_name};
    report.use_head0_baseline = logical(use_head0_baseline);
    report.alpha = alpha;
    if use_head0_baseline
        report.alpha_name = {'Head0'};
        report.target_Pw_over_Pi = NaN;
        report.target_N_over_Pi = NaN;
    else
        report.alpha_name = {step7_alpha_name(alpha)};
        report.target_Pw_over_Pi = alpha;
        report.target_N_over_Pi = 1 - alpha;
    end
    report.n_shelf_nodes = sum(shelf_nodes(:));
    report.n_grounded_ring_nodes = sum(ring);
    report.n_finite_spchead = sum(finite_spc);
    report.n_zero_spchead = sum(zero_spc);
    report.n_grounded = sum(masks.is_grounded(:));
    report.n_shelf = sum(masks.is_shelf(:));
    if has_model_member(md.hydrology, 'relaxation')
        report.hydrology_relaxation = md.hydrology.relaxation;
    else
        report.hydrology_relaxation = NaN;
    end

    report.ring_spchead_min = finite_min(spc(ring));
    report.ring_spchead_median = finite_median(spc(ring));
    report.ring_spchead_max = finite_max(spc(ring));
    report.ring_head_float_minus_spchead_min = finite_min(head_float_minus_spc(ring));
    report.ring_head_float_minus_spchead_median = finite_median(head_float_minus_spc(ring));
    report.ring_head_float_minus_spchead_max = finite_max(head_float_minus_spc(ring));
    report.ring_Pw_over_Pi_min = finite_min(Pw_over_Pi(ring));
    report.ring_Pw_over_Pi_median = finite_median(Pw_over_Pi(ring));
    report.ring_Pw_over_Pi_max = finite_max(Pw_over_Pi(ring));
    report.ring_N_over_Pi_min = finite_min(N_over_Pi(ring));
    report.ring_N_over_Pi_median = finite_median(N_over_Pi(ring));
    report.ring_N_over_Pi_max = finite_max(N_over_Pi(ring));
end


function save_step7_diagnostics(md, masks, shelf_nodes, grounded_ring_nodes, alpha, case_name, use_head0_baseline, outdir, report)
    % 保存 Step 7 的诊断输出：
    %   - step7_spchead_nodes.csv: 每个固定水头节点的几何、水头、Pw/Pi、N/Pi
    %   - step7_summary.csv: 当前 case 的汇总统计
    %   - PNG 图: spchead 位置、水头值、Pw/Pi、N/Pi 和离浮托水头的差值
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    spc = md.hydrology.spchead(:);
    H = md.geometry.thickness(:);
    bed = md.geometry.base(:);
    rho_i = md.materials.rho_ice;
    rho_w = md.materials.rho_freshwater;
    g = md.constants.g;

    Pi = rho_i .* g .* H;
    Pw = rho_w .* g .* (spc - bed);
    N = Pi - Pw;
    Pw_over_Pi = Pw ./ Pi;
    N_over_Pi = N ./ Pi;
    head_float = bed + (rho_i ./ rho_w) .* H;
    head_float_minus_spc = head_float - spc;

    finite_spc = isfinite(spc);
    target_nodes = finite_spc;

    % 保存节点级 CSV，后续可直接检查每个固定水头节点的几何和压力状态。
    T = table();
    ids = find(target_nodes);
    T.node_id = ids(:);
    T.x = reshape(md.mesh.x(target_nodes), [], 1);
    T.y = reshape(md.mesh.y(target_nodes), [], 1);
    T.is_shelf = shelf_nodes(target_nodes);
    T.is_grounded_ring = grounded_ring_nodes(target_nodes);
    T.spchead_m = spc(target_nodes);
    T.bed_m = bed(target_nodes);
    T.thickness_m = H(target_nodes);
    T.head_float_m = head_float(target_nodes);
    T.head_float_minus_spchead_m = head_float_minus_spc(target_nodes);
    T.Pw_over_Pi = Pw_over_Pi(target_nodes);
    T.N_over_Pi = N_over_Pi(target_nodes);
    T.alpha = alpha * ones(height(T), 1);
    T.case_name = repmat({case_name}, height(T), 1);
    T.use_head0_baseline = repmat(logical(use_head0_baseline), height(T), 1);
    writetable(T, fullfile(outdir, 'step7_spchead_nodes.csv'));

    writetable(struct2table(report), fullfile(outdir, 'step7_summary.csv'));

    % 保存空间诊断图：位置、水头值、ring 上 Pw/Pi 和 N/Pi。
    x = md.mesh.x(:);
    y = md.mesh.y(:);
    bg = zeros(md.mesh.numberofvertices, 1);
    bg(masks.is_noice(:)) = 1;
    bg(masks.is_grounded(:)) = 2;
    bg(masks.is_shelf(:)) = 3;

    if use_head0_baseline
        ring_title = 'spchead 节点位置: shelf=0, ring head=0';
        ring_legend = 'grounded ring spchead=0';
    else
        ring_title = sprintf('spchead 节点位置: shelf=0, ring alpha=%.3f', alpha);
        ring_legend = 'grounded ring spchead=head\_alpha';
    end

    fig = figure('Visible', 'off', 'Name', ['STEP 7 ' case_name], ...
        'Color', 'w', 'Position', [80 80 1350 950]);
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    patch('Faces', md.mesh.elements, 'Vertices', [x, y], ...
        'FaceVertexCData', bg, 'FaceColor', 'interp', 'EdgeColor', 'none');
    axis equal tight; box on; hold on; colorbar;
    plot(x(shelf_nodes), y(shelf_nodes), 'b.', 'MarkerSize', 5);
    plot(x(grounded_ring_nodes), y(grounded_ring_nodes), 'k.', 'MarkerSize', 6);
    title(ring_title, 'Interpreter', 'none');
    legend({'node class', 'shelf spchead=0', ring_legend}, ...
        'Location', 'bestoutside');

    nexttile;
    scatter(x(finite_spc), y(finite_spc), 12, spc(finite_spc), 'filled');
    axis equal tight; box on; colorbar;
    title('有限 spchead 值 [m]', 'Interpreter', 'none');
    xlabel('X (m)'); ylabel('Y (m)');

    nexttile;
    scatter(x(grounded_ring_nodes), y(grounded_ring_nodes), 14, ...
        Pw_over_Pi(grounded_ring_nodes), 'filled');
    axis equal tight; box on; colorbar;
    title('接地冰 ring: Pw/Pi', 'Interpreter', 'none');
    xlabel('X (m)'); ylabel('Y (m)');

    nexttile;
    scatter(x(grounded_ring_nodes), y(grounded_ring_nodes), 14, ...
        N_over_Pi(grounded_ring_nodes), 'filled');
    axis equal tight; box on; colorbar;
    title('接地冰 ring: N/Pi', 'Interpreter', 'none');
    xlabel('X (m)'); ylabel('Y (m)');

    png_file = fullfile(outdir, sprintf('step7_%s_spchead_diagnostics.png', case_name));
    saveas(fig, png_file);
    close(fig);

    fig = figure('Visible', 'off', 'Name', ['STEP 7 head offset ' case_name], ...
        'Color', 'w', 'Position', [120 120 900 760]);
    scatter(x(grounded_ring_nodes), y(grounded_ring_nodes), 14, ...
        head_float_minus_spc(grounded_ring_nodes), 'filled');
    axis equal tight; box on; colorbar;
    title('接地冰 ring: 浮托水头 - 固定 spchead [m]', 'Interpreter', 'none');
    xlabel('X (m)'); ylabel('Y (m)');
    saveas(fig, fullfile(outdir, sprintf('step7_%s_head_float_minus_spchead.png', case_name)));
    close(fig);
end


function md_out = run_step7_transient_case(md_in, result_model_name, ...
    sim_start_time, sim_dt_seconds, sim_n_steps, sim_output_freq, ...
    sim_stress_restol, sim_stress_reltol, sim_stress_abstol, ...
    sim_stress_maxiter, sim_verbose_solution, sim_verbose_module, ...
    sim_verbose_convergence, n_cores)

    % 对一个 alpha setup 模型运行 transient。这个函数只在
    % step7_run_transient=true 时调用。
    md_in.miscellaneous.name = ['RecoveryNewMesh_' result_model_name];
    md_in.timestepping.start_time = sim_start_time;
    md_in.timestepping.time_step  = sim_dt_seconds / md_in.constants.yts;
    md_in.timestepping.final_time = (sim_n_steps * sim_dt_seconds / 86400) / 365;
    md_in.settings.output_frequency = sim_output_freq;
    md_in.stressbalance.restol = sim_stress_restol;
    md_in.stressbalance.reltol = sim_stress_reltol;
    md_in.stressbalance.abstol = sim_stress_abstol;
    md_in.stressbalance.maxiter = sim_stress_maxiter;
    md_in.verbose.solution = sim_verbose_solution;
    md_in.verbose.module = sim_verbose_module;
    md_in.verbose.convergence = sim_verbose_convergence;
    md_in.cluster = generic('name', oshostname(), 'np', n_cores);

    try
        md_out = solve(md_in, 'Transient');
    catch ME
        fprintf(2, 'Step 7 transient solve failed for %s: %s\n', result_model_name, ME.message);
        md_out = md_in;
    end

    try
        has_transient_results = isfield(md_out.results, 'TransientSolution') && ...
            ~isempty(md_out.results.TransientSolution);
    catch
        has_transient_results = false;
    end

    if has_transient_results
        last_time = md_out.results.TransientSolution(end).time;
        if last_time < md_out.timestepping.final_time - md_out.timestepping.time_step
            fprintf(2, ['Warning: Step 7 transient case %s stopped early at %.9g yr, ' ...
                'before final_time %.9g yr.\n'], ...
                result_model_name, last_time, md_out.timestepping.final_time);
        end
    else
        fprintf(2, 'Warning: Step 7 transient case %s returned no TransientSolution.\n', ...
            result_model_name);
    end
end


function tf = has_model_member(s, name)
    tf = (isstruct(s) && isfield(s, name)) || (isobject(s) && isprop(s, name));
end


function v = finite_min(values)
    values = values(isfinite(values));
    if isempty(values), v = NaN; else, v = min(values); end
end


function v = finite_max(values)
    values = values(isfinite(values));
    if isempty(values), v = NaN; else, v = max(values); end
end


function v = finite_median(values)
    values = values(isfinite(values));
    if isempty(values), v = NaN; else, v = median(values); end
end


function plot_spchead_diagnostics(md, masks)
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

end


function plot_step5_spchead_outlet_diagnostics(md, masks, outlet_file)
    spc = md.hydrology.spchead(:);
    finite_spc = isfinite(spc);

    tol_zero = 1e-9;
    tol_base = 1e-6;
    spc_zero = finite_spc & abs(spc) <= tol_zero;
    spc_base = finite_spc & abs(spc - md.geometry.base(:)) <= tol_base & ~spc_zero;
    spc_other = finite_spc & ~spc_zero & ~spc_base;

    outlet_nodes = false(md.mesh.numberofvertices, 1);
    if nargin >= 3 && ~isempty(outlet_file) && exist(outlet_file, 'file')
        outlet_nodes = logical(ContourToNodes(md.mesh.x, md.mesh.y, outlet_file, 1));
        outlet_nodes = outlet_nodes(:);
    end

    fprintf('   STEP 5 diagnostic spchead nodes: finite=%d, zero=%d, base=%d, other=%d\n', ...
        sum(finite_spc), sum(spc_zero), sum(spc_base), sum(spc_other));
    fprintf('   STEP 5 outlet EXP nodes: %d\n', sum(outlet_nodes));

    fig = figure('Name', 'STEP 5 spchead and outlet diagnostics', ...
        'Color', 'w', 'Position', [80 80 1100 850]);
    bg = zeros(md.mesh.numberofvertices, 1);
    bg(masks.is_noice(:)) = 1;
    bg(masks.is_grounded(:)) = 2;
    bg(masks.is_shelf(:)) = 3;
    bg_label = 'Node class: 1=no ice, 2=grounded, 3=shelf';
    draw_spchead_panel(md, bg, bg_label, spc_zero, spc_base, spc_other);
    hold on;

    if any(outlet_nodes)
        plot(md.mesh.x(outlet_nodes), md.mesh.y(outlet_nodes), 'm.', ...
            'MarkerSize', 7, 'DisplayName', 'outlet exp nodes');
    end
    draw_exp_overlay(outlet_file);
    title('STEP 5 finite spchead over node classes and outlet EXP', 'Interpreter', 'none');

    drawnow;
    fprintf('   STEP 5 spchead/outlet diagnostic figure displayed.\n');
end


function draw_exp_overlay(exp_file)
    if nargin < 1 || isempty(exp_file) || exist(exp_file, 'file') ~= 2
        return;
    end

    try
        exp_struct = expread(exp_file);
    catch ME
        fprintf(2, 'Warning: could not read outlet EXP for plotting: %s\n', ME.message);
        return;
    end

    if ~isstruct(exp_struct)
        return;
    end

    for k = 1:numel(exp_struct)
        if isfield(exp_struct(k), 'x') && isfield(exp_struct(k), 'y')
            x = exp_struct(k).x;
            y = exp_struct(k).y;
        elseif isfield(exp_struct(k), 'X') && isfield(exp_struct(k), 'Y')
            x = exp_struct(k).X;
            y = exp_struct(k).Y;
        else
            continue;
        end

        x = x(:);
        y = y(:);
        valid = isfinite(x) & isfinite(y);
        if nnz(valid) < 2
            continue;
        end
        patch(x(valid), y(valid), 'm', ...
            'FaceAlpha', 0.08, 'EdgeColor', 'm', 'LineWidth', 1.8, ...
            'DisplayName', 'outlet exp area');
    end
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
