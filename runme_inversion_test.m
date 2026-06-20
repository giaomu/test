% =========================================================================
%  runme_Lcurve_sweep_multiw101.m
%  Recovery Glacier — L-curve coarse sweep for regularization weight w501
%
%  新版功能：
%    - 外层扫描多个 w101
%    - 内层扫描 w501_list
%    - 每个 w101 自动生成独立输出目录
%    - 每个 w101 自动生成独立 summary / plot
%    - 同时生成总 summary
%
%  输出目录格式：
%    outputs/models/L_curve_models_<w101>_<w103>_<w501min>_to_<w501max>
%
% =========================================================================
clear; clc; close all;
paths = shaktiais_paths();
if ~exist(paths.models, 'dir'), mkdir(paths.models); end

issm_root = paths.issm_root;
issm_bin  = fullfile(issm_root, 'bin');
if exist(issm_bin, 'dir') == 7
    addpath(issm_bin);
else
    error('ISSM bin folder not found: %s', issm_bin);
end

% -------------------------------------------------------------------------
%  0. 全局配置
% -------------------------------------------------------------------------

% --- 基准模型加载方式 ---
% 方式 A: 用 organizer
org = organizer('repository', paths.models, 'prefix', 'Recovery_', 'steps', [2]);
base_model_loader = @() loadmodel(org, 'Parameterize');

% 方式 B: 直接读 mat 文件 (如需切换，注释掉上面，启用下面)
% base_model_loader = @() load_and_extract(fullfile(paths.models, 'Recovery_Parameterize.mat'));

% --- BedMachine 路径 ---
bm_file = paths.bedmachine_file;

% --- 反演固定参数 ---
w101_list = [20000 40000 50000 70000];   % <<< 在这里写多个 w101
w103 = 10;
 
inv_maxsteps      = 200;
inv_maxiter       = 400;
inv_friction_min  = 0.05;   

inv_friction_max  = 3500;
n_cores           = 15;

% --- warm start：用已有单次反演结果的摩擦系数作为初始摩擦 ---
use_warm_start_friction = 0;
warm_start_model_file   = fullfile(paths.models, 'Recovery_Inversion.mat');

% --- w501 扫描列表 (对数粗扫) ---

w501_list = [ ...
    1.000000e-06
    1.800000e-06
    3.200000e-06
    5.600000e-06
    1.000000e-05
    1.500000e-05
    2.000000e-05
    2.500000e-05
    3.000000e-05
    3.500000e-05
    4.000000e-05
    4.500000e-05
    5.000000e-05
    5.500000e-05
    6.000000e-05
    7.000000e-05
    8.500000e-05
    1.000000e-04
    1.300000e-04
    1.800000e-04
    2.800000e-04
    5.000000e-04 ];

%{
w501_list = [ ...
    1.000000e-06
    1.304947e-06
    1.702886e-06
    2.222176e-06
    2.899821e-06
    3.784113e-06
    4.938066e-06
    6.443913e-06
    8.408964e-06
    1.097325e-05
    1.431951e-05
    1.868620e-05
    2.438449e-05
    3.182047e-05
    4.152402e-05
    5.418664e-05
    7.071068e-05
    9.227367e-05
    1.204122e-04
    1.571316e-04
    2.050483e-04
    2.675772e-04
    3.491740e-04
    4.556535e-04
    5.946036e-04
    7.759260e-04
    1.012542e-03
    1.321314e-03
    1.724244e-03
    2.250047e-03 
    2.936192e-03
    3.831574e-03
    5.000000e-03 ];
%}
%{
w501_list = [ ...
    5.000000e-03
    6.524734e-03
    8.514431e-03
    1.111088e-02
    1.449911e-02
    1.892056e-02
    2.469033e-02
    3.221957e-02
    4.204482e-02
    5.486626e-02
    7.159754e-02
    9.343099e-02
    1.219225e-01
    1.591023e-01
    2.076201e-01
    2.709332e-01
    3.535534e-01
    4.613684e-01
    6.020612e-01
    7.856578e-01 ];
%}
n_w101 = numel(w101_list);
n_w501 = numel(w501_list);
n_runs_total = n_w101 * n_w501;

% --- 目标经验比例 10:5:2 ---
target_ratio = [10, 5, 2] / 17;

% -------------------------------------------------------------------------
%  1. 节点分类 (BedMachine mask, 只算一次)
% -------------------------------------------------------------------------
fprintf('Loading base model for node classification...\n');
md_tmp = base_model_loader();
[is_grounded, is_shelf, is_noice, is_vostok] = ...
    local_extract_node_types_from_bedmachine(bm_file, md_tmp.mesh.x, md_tmp.mesh.y, md_tmp.mesh.elements);
fprintf('   [After cliff suppression] grounded=%d  shelf=%d  no-ice=%d  vostok=%d\n', ...
    sum(is_grounded), sum(is_shelf), sum(is_noice), sum(is_vostok));
clear md_tmp

warm_start_friction = [];
if use_warm_start_friction
    warm_start_friction = local_load_warm_start_friction(warm_start_model_file);
    fprintf('Warm-start friction loaded from: %s\n', warm_start_model_file);
    fprintf('   min=%.4g  max=%.4g  mean=%.4g\n', ...
        min(warm_start_friction), max(warm_start_friction), mean(warm_start_friction));
end

% -------------------------------------------------------------------------
%  2. 初始化总汇总表
% -------------------------------------------------------------------------
col_names = { ...
    'run_id', 'w101', 'w103', 'w501', 'status', 'error_message', ...
    'runtime_sec', 'n_iter', 'InversionStopFlag', 'grad_norm_final', ...
    'J101', 'J103', 'J501', 'Jtotal', 'Jdata', ...
    'phi_data_unweighted', 'phi_reg_unweighted', 'logPhiData_unweighted', 'logPhiReg_unweighted', ...
    'ratio101', 'ratio103', 'ratio501', ...
    'target_ratio101', 'target_ratio103', 'target_ratio501', 'ratio_distance', ...
    'logPhiData', 'logPhiReg', ...
    'model_file', 'log_file'};

summary_all = local_empty_table(n_runs_total, col_names, target_ratio);

% -------------------------------------------------------------------------
%  3. 主循环：外层扫 w101，内层扫 w501
% -------------------------------------------------------------------------
fprintf('\n========================================\n');
fprintf('  Multi-w101 L-curve sweep: %d w101 values x %d w501 values = %d runs\n', ...
        n_w101, n_w501, n_runs_total);
fprintf('========================================\n\n');

global_run_id = 0;

for i101 = 1:n_w101

    w101 = w101_list(i101);

    out_dir = fullfile(paths.models, sprintf('L_curve_models_%s_%s_%s_to_%s', ...
        local_num_to_str(w101), ...
        local_num_to_str(w103), ...
        local_num_to_expstr(w501_list(1)), ...
        local_num_to_expstr(w501_list(end))));

    log_dir = fullfile(out_dir, 'logs');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    if ~exist(log_dir, 'dir'), mkdir(log_dir); end

    fprintf('\n========================================\n');
    fprintf('  Sweeping w101 = %g\n', w101);
    fprintf('  Output dir: %s\n', out_dir);
    fprintf('========================================\n\n');

    summary_one = local_empty_table(n_w501, col_names, target_ratio);

    for k = 1:n_w501

        global_run_id = global_run_id + 1;
        w501 = w501_list(k);

        safe_name = local_make_safe_filename(w501);
        model_fname = sprintf('Lcurve_run%02d_w501_%s.mat', k, safe_name);
        log_fname   = sprintf('Lcurve_run%02d_w501_%s.txt', k, safe_name);
        model_path  = fullfile(out_dir, model_fname);
        log_path    = fullfile(log_dir, log_fname);

        fprintf('--- Global %03d/%03d | w101 = %.4g | w501 = %.4e ---\n', ...
                global_run_id, n_runs_total, w101, w501);

        row = local_empty_row(global_run_id, w101, w103, w501, ...
                              model_fname, log_fname, target_ratio, col_names);

        try
            % ==============================================================
            % 3.1 从基准模型重新加载
            % ==============================================================
            md = base_model_loader();
            nv = md.mesh.numberofvertices;
            ne = md.mesh.numberofelements;

            if use_warm_start_friction
                if numel(warm_start_friction) ~= nv
                    error('warm-start friction length (%d) does not match current nv (%d).', ...
                          numel(warm_start_friction), nv);
                end
                md.friction.coefficient = warm_start_friction(:);
            end

            % ==============================================================
            % 3.2 反演框架设置
            % ==============================================================
            md.miscellaneous.name = sprintf('Lcurve_w101_%s_run%02d_w501_%s', ...
                local_num_to_str(w101), k, safe_name);

            md.inversion = m1qn3inversion(md.inversion);
            md.inversion.iscontrol = 1;
            md.transient.amr_frequency = 0;
            md.verbose = verbose('solution', true, 'control', true, 'convergence', true);

            % ==============================================================
            % 3.3 代价函数与权重
            % ==============================================================
            md.inversion.cost_functions = [101, 103, 501];
            md = local_apply_cost_masks(md, is_shelf, is_noice, w101, w103, w501, ...
                                        inv_friction_min, inv_friction_max);
            %{
            fc = md.friction.coefficient;
            fprintf('Before solve: min = %.4g, max = %.4g\n', min(fc), max(fc));
            fprintf('Before solve: n_zero = %d\n', sum(fc==0));
            fprintf('Before solve: n_20   = %d\n', sum(abs(fc-20)<1e-12));
            % ===== 绘制进入 solve 前的摩擦系数初值场 =====
            figure('Name','Friction coefficient before solve', ...
                   'Color','w', 'Position',[100 100 850 700]);
            plotmodel(md, ...
                'data', fc, ...
                'title', sprintf(['Friction coefficient before solve\n' ...
                                  'min=%.4g, max=%.4g, n_{zero}=%d, n_{20}=%d'], ...
                                  min(fc), max(fc), sum(fc==0), sum(abs(fc-20)<1e-12)), ...
                'edgecolor', 'none', ...
                'colorbar', 'on', ...
                'colormap', 'parula');
            drawnow;
            saveas(gcf, fullfile(out_dir, 'friction_before_solve.fig'));

            % ===== 绘制 101 和 103 都不参与的节点 =====
            coeff = md.inversion.cost_functions_coefficients;
            
            flag_no101_103 = zeros(md.mesh.numberofvertices,1);
            idx_no101_103 = find(coeff(:,1)==0 & coeff(:,2)==0);
            flag_no101_103(idx_no101_103) = 1;
            
            fprintf('Nodes with both 101 and 103 OFF: %d\n', numel(idx_no101_103));
            
            figure('Name','Nodes excluded from 101 and 103', ...
                   'Color','w', 'Position',[120 120 850 700]);
            
            plotmodel(md, ...
                'data', flag_no101_103, ...
                'title', sprintf('Nodes with 101 and 103 excluded (N = %d)', numel(idx_no101_103)), ...
                'edgecolor', 'none', ...
                'colorbar', 'on', ...
                'colormap', 'parula');
            
            caxis([0 1]);
            drawnow;
            
            saveas(gcf, fullfile(out_dir, sprintf('no101_no103_nodes_w101_%g_w501_%g.fig', w101, w501)));
            %}

            % ==============================================================
            % 3.4 控制变量
            % ==============================================================
            md.inversion.control_parameters      = {'FrictionCoefficient'};
            md.inversion.maxsteps                = inv_maxsteps;
            md.inversion.maxiter                 = inv_maxiter;
            md.inversion.control_scaling_factors = 1;

            % ==============================================================
            % 3.5 求解器精度
            % ==============================================================
            md.stressbalance.restol = 0.01;
            md.stressbalance.reltol = 0.1;
            md.stressbalance.abstol = NaN;

            % ==============================================================
            % 3.6 求解 (evalc 捕获日志, tic/toc 计时)
            % ==============================================================
            md.cluster = generic('name', oshostname(), 'np', n_cores);

            tic_start = tic;
            cmdout = evalc('md = solve(md,''sb'');');
            elapsed = toc(tic_start);

            % 保存日志
            fid = fopen(log_path, 'w');
            if fid > 0
                fwrite(fid, cmdout);
                fclose(fid);
            end

            % ==============================================================
            % 3.7 提取结果
            % ==============================================================
            assert(isfield(md.results, 'StressbalanceSolution'), 'No StressbalanceSolution');

            Jhist = md.results.StressbalanceSolution.J;
            n_iter = size(Jhist, 1);
            J101_val   = Jhist(end, 1);
            J103_val   = Jhist(end, 2);
            J501_val   = Jhist(end, 3);
            Jtotal_val = Jhist(end, 4);
            Jdata_val  = J101_val + J103_val;

            % 未加权量
            phi_data_uw = J101_val / w101 + J103_val / w103;
            phi_reg_uw  = J501_val / w501;

            % 梯度范数 (从日志解析)
            [~, ~, grad_norm, ~, ~, ~] = local_parse_final_m1qn3_line(cmdout);

            % StopFlag
            stop_flag = NaN;
            if isfield(md.results.StressbalanceSolution, 'InversionStopFlag')
                stop_flag = md.results.StressbalanceSolution.InversionStopFlag;
            end

            % 归一化比例
            Jsum = J101_val + J103_val + J501_val;
            if Jsum > 0
                r101 = J101_val / Jsum;
                r103 = J103_val / Jsum;
                r501 = J501_val / Jsum;
            else
                r101 = NaN; r103 = NaN; r501 = NaN;
            end

            rdist = sqrt((r101 - target_ratio(1))^2 + ...
                         (r103 - target_ratio(2))^2 + ...
                         (r501 - target_ratio(3))^2);

            % ---- 反演后回写（方便直接衔接 step 4）----
            md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient;
            md.friction.coefficient(is_shelf) = 0;
            md.initialization.vx  = md.results.StressbalanceSolution.Vx;
            md.initialization.vy  = md.results.StressbalanceSolution.Vy;
            md.initialization.vz  = zeros(nv,1);
            md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);

            % 保存模型
            save(model_path, 'md', '-v7.3');

            % 填表
            row.status                = {'success'};
            row.runtime_sec           = elapsed;
            row.n_iter                = n_iter;
            row.InversionStopFlag     = stop_flag;
            row.grad_norm_final       = grad_norm;
            row.J101                  = J101_val;
            row.J103                  = J103_val;
            row.J501                  = J501_val;
            row.Jtotal                = Jtotal_val;
            row.Jdata                 = Jdata_val;
            row.phi_data_unweighted   = phi_data_uw;
            row.phi_reg_unweighted    = phi_reg_uw;
            row.logPhiData_unweighted = log10(max(phi_data_uw, eps));
            row.logPhiReg_unweighted  = log10(max(phi_reg_uw, eps));
            row.ratio101              = r101;
            row.ratio103              = r103;
            row.ratio501              = r501;
            row.ratio_distance        = rdist;
            row.logPhiData            = log10(max(Jdata_val, eps));
            row.logPhiReg             = log10(max(J501_val, eps));

            fprintf('   OK  %.1fs | %d iters | J=[%.2e %.2e %.2e] tot=%.2e | rdist=%.4f\n', ...
                    elapsed, n_iter, J101_val, J103_val, J501_val, Jtotal_val, rdist);

        catch ME
            row.status        = {'failed'};
            row.error_message = {ME.message};
            fprintf('   FAILED: %s\n', ME.message);
        end

        % 当前 w101 的汇总
        summary_one(k,:) = row;

        % 总汇总
        summary_all(global_run_id,:) = row;

        % 即时保存
        save(fullfile(out_dir, 'Lcurve_summary.mat'), 'summary_one');
        writetable(summary_one, fullfile(out_dir, 'Lcurve_summary.csv'));

        save(fullfile(paths.models, 'Lcurve_summary_all.mat'), 'summary_all');
        writetable(summary_all, fullfile(paths.models, 'Lcurve_summary_all.csv'));

    end

    % ---------------------------------------------------------------------
    %  4. 当前 w101 的后处理
    % ---------------------------------------------------------------------
    local_postprocess_one_summary(summary_one, out_dir, w101, w103);

end

fprintf('\n========================================\n');
fprintf('  All sweeps finished.\n');
fprintf('  Overall summary saved to:\n');
fprintf('     %s\n', fullfile(paths.models, 'Lcurve_summary_all.mat'));
fprintf('     %s\n', fullfile(paths.models, 'Lcurve_summary_all.csv'));
fprintf('========================================\n');


% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [is_grounded, is_shelf, is_noice, is_vostok] = ...
        local_extract_node_types_from_bedmachine(bm_file, mesh_x, mesh_y, elements)
% 读 BedMachine Antarctica mask，做与 par / step 3 一致的 Cliff suppression。
% 返回的四个分类都是 CS 后的版本：
%   is_grounded : CS 后仍为 grounded 的节点
%   is_shelf    : CS 后仍为 shelf 的节点（未被降级）
%   is_noice    : 原始非冰 ∪ CS 吞下来的过渡节点
%   is_vostok   : CS 后仍为 Lake Vostok 的节点

    x_bm    = double(ncread(bm_file, 'x'));
    y_bm    = double(ncread(bm_file, 'y'));
    mask_bm = double(ncread(bm_file, 'mask'));

    [x_bm, ix] = sort(x_bm);
    [y_bm, iy] = sort(y_bm);
    mask_bm = mask_bm';
    mask_bm = mask_bm(iy, ix);

    Fmask     = griddedInterpolant({y_bm, x_bm}, mask_bm, 'nearest', 'nearest');
    node_type = round(Fmask(mesh_y, mesh_x));

    % ---- 原始分类 ----
    is_ocean        = (node_type == 0);
    is_land         = (node_type == 1);
    is_grounded_raw = (node_type == 2);
    is_shelf_raw    = (node_type == 3);
    is_vostok_raw   = (node_type == 4);

    is_ice_raw   = is_grounded_raw | is_shelf_raw | is_vostok_raw;
    is_noice_raw = is_ocean | is_land;

    % ---- Cliff suppression：含任何原始非冰节点的单元，整单元降级为非冰 ----
    is_ice = is_ice_raw;
    nonice_elem = any(is_noice_raw(elements), 2);
    is_ice(elements(nonice_elem,:)) = false;

    % ---- CS 后分类输出 ----
    is_grounded = is_grounded_raw & is_ice;
    is_shelf    = is_shelf_raw    & is_ice;
    is_vostok   = is_vostok_raw   & is_ice;
    is_noice    = ~is_ice;
end


function md = local_apply_cost_masks(md, is_shelf, is_noice, w101, w103, w501, ...
                                      friction_min, friction_max)
% 等价于 runme step 3 的节 4 + 节 5 + 节 6
% is_shelf / is_noice 必须已经是 cliff-suppression 之后的分类

    nv = md.mesh.numberofvertices;

    % ---- 节 5：代价函数权重 ----
    md.inversion.cost_functions_coefficients = zeros(nv, 3);
    md.inversion.cost_functions_coefficients(:,1) = w101;
    md.inversion.cost_functions_coefficients(:,2) = w103;
    md.inversion.cost_functions_coefficients(:,3) = w501;

    % 真·冰架：三列全关
    md.inversion.cost_functions_coefficients(is_shelf, :) = 0;

    % 非冰 + 过渡节点：101/103 关，501 保留（平滑）
    %md.inversion.cost_functions_coefficients(is_noice, 1:2) = 0;

    % ---- 节 4：摩擦处理 ----
    md.friction.coefficient(is_shelf) = 0;

    reactivate = is_noice & (~isfinite(md.friction.coefficient) | md.friction.coefficient <= 0);
    md.friction.coefficient(reactivate) = max(friction_min, 20);

    % ---- 节 6：控制变量上下界 ----
    minp = friction_min * ones(nv, 1);
    maxp = friction_max * ones(nv, 1);
    minp(is_shelf) = 0;
    maxp(is_shelf) = 0;
    md.inversion.min_parameters = minp;
    md.inversion.max_parameters = maxp;
end


function friction = local_load_warm_start_friction(model_file)
% Load an existing inversion result and return a nodal friction coefficient.
% Prefer the solved FrictionCoefficient stored in StressbalanceSolution;
% fall back to md.friction.coefficient if the result field is unavailable.

    assert(exist(model_file, 'file') == 2, ...
        'warm-start model file does not exist: %s', model_file);

    S = load(model_file);
    md = local_extract_md_from_struct(S);
    clear S

    friction = [];
    try
        if ~isempty(md.results.StressbalanceSolution.FrictionCoefficient)
            friction = double(md.results.StressbalanceSolution.FrictionCoefficient(:));
        end
    catch
    end

    if isempty(friction)
        try
            if ~isempty(md.friction.coefficient)
                friction = double(md.friction.coefficient(:));
            end
        catch
        end
    end

    if isempty(friction)
        error('No warm-start friction coefficient found in: %s', model_file);
    end

    if any(~isfinite(friction))
        error('Warm-start friction contains NaN/Inf: %s', model_file);
    end
end


function md = local_extract_md_from_struct(S)
    if isfield(S, 'md')
        md = S.md;
        return;
    end

    fns = fieldnames(S);
    for i = 1:numel(fns)
        obj = S.(fns{i});
        if isstruct(obj) && isfield(obj, 'mesh')
            md = obj;
            return;
        end
        if isobject(obj) && isprop(obj, 'mesh')
            md = obj;
            return;
        end
    end

    error('Cannot find md object in warm-start model file.');
end


function [iter_final, fx_final, grad_norm, J101_log, J103_log, J501_log] = ...
        local_parse_final_m1qn3_line(cmdout)
    iter_final = NaN; fx_final = NaN; grad_norm = NaN;
    J101_log = NaN; J103_log = NaN; J501_log = NaN;
    try
        lines = strsplit(cmdout, {'\n', '\r'});
        fx_lines = lines(contains(lines, 'f(x)='));
        if isempty(fx_lines), return; end
        last_line = strtrim(fx_lines{end});
        tokens = regexp(last_line, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
        vals = cellfun(@str2double, tokens);
        if numel(vals) >= 6
            iter_final = vals(1); fx_final = vals(2); grad_norm = vals(3);
            J101_log = vals(4); J103_log = vals(5); J501_log = vals(6);
        elseif numel(vals) >= 3
            iter_final = vals(1); fx_final = vals(2); grad_norm = vals(3);
        end
    catch
    end
end


function s = local_make_safe_filename(w501)
    raw = sprintf('%.3e', w501);
    s   = strrep(raw, '.', 'p');
    s   = strrep(s, '+', '');
end


function T = local_empty_table(n, col_names, target_ratio)
    nc = numel(col_names);
    C = cell(1, nc);
    for i = 1:nc
        nm = col_names{i};
        if ismember(nm, {'status','error_message','model_file','log_file'})
            C{i} = repmat({''}, n, 1);
        elseif strcmp(nm, 'target_ratio101')
            C{i} = repmat(target_ratio(1), n, 1);
        elseif strcmp(nm, 'target_ratio103')
            C{i} = repmat(target_ratio(2), n, 1);
        elseif strcmp(nm, 'target_ratio501')
            C{i} = repmat(target_ratio(3), n, 1);
        else
            C{i} = zeros(n, 1);
        end
    end
    T = table(C{:}, 'VariableNames', col_names);
end


function row = local_empty_row(run_id, w101, w103, w501, model_fname, log_fname, target_ratio, col_names)
    nc = numel(col_names);
    C = cell(1, nc);
    for i = 1:nc
        nm = col_names{i};
        switch nm
            case 'run_id',            C{i} = run_id;
            case 'w101',              C{i} = w101;
            case 'w103',              C{i} = w103;
            case 'w501',              C{i} = w501;
            case 'status',            C{i} = {'pending'};
            case 'error_message',     C{i} = {''};
            case 'model_file',        C{i} = {model_fname};
            case 'log_file',          C{i} = {log_fname};
            case 'target_ratio101',   C{i} = target_ratio(1);
            case 'target_ratio103',   C{i} = target_ratio(2);
            case 'target_ratio501',   C{i} = target_ratio(3);
            otherwise,                C{i} = NaN;
        end
    end
    row = table(C{:}, 'VariableNames', col_names);
end


function idx = local_find_lcurve_corner(logPhiData, logPhiReg)
    n = numel(logPhiData);
    if n < 3
        idx = NaN;
        return;
    end
    kappa = zeros(n, 1);
    for i = 2:(n-1)
        x1 = logPhiReg(i-1);  y1 = logPhiData(i-1);
        x2 = logPhiReg(i);    y2 = logPhiData(i);
        x3 = logPhiReg(i+1);  y3 = logPhiData(i+1);

        area2 = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
        d12 = sqrt((x2-x1)^2 + (y2-y1)^2);
        d23 = sqrt((x3-x2)^2 + (y3-y2)^2);
        d13 = sqrt((x3-x1)^2 + (y3-y1)^2);

        denom = d12 * d23 * d13;
        if denom > 0
            kappa(i) = area2 / denom;
        end
    end
    [~, idx] = max(kappa);
    if kappa(idx) == 0
        idx = NaN;
    end
end


function local_postprocess_one_summary(summary, out_dir, w101, w103)

    fprintf('\n========================================\n');
    fprintf('  Post-processing for w101 = %g, w103 = %g\n', w101, w103);
    fprintf('========================================\n\n');

    ok = strcmp(summary.status, 'success');
    if ~any(ok)
        fprintf('All runs failed for w101 = %g. No candidates.\n', w101);
        return
    end

    S = summary(ok, :);
    [S_sorted, sort_idx] = sortrows(S, 'w501');

    % --- 候选 A: 最接近目标比例 ---
    [~, idxA] = min(S.ratio_distance);
    fprintf('Candidate A (closest ratio):\n');
    fprintf('   run_id = %d | w501 = %.4e | ratio_distance = %.4f\n', ...
            S.run_id(idxA), S.w501(idxA), S.ratio_distance(idxA));
    fprintf('   J101=%.2e  J103=%.2e  J501=%.2e\n', ...
            S.J101(idxA), S.J103(idxA), S.J501(idxA));

    % --- 候选 B: L-curve 拐点 (加权版) ---
    idxB_local = local_find_lcurve_corner(S_sorted.logPhiData, S_sorted.logPhiReg);
    if ~isnan(idxB_local)
        idxB = sort_idx(idxB_local);
        fprintf('Candidate B (weighted L-curve corner):\n');
        fprintf('   run_id = %d | w501 = %.4e\n', S.run_id(idxB), S.w501(idxB));
        fprintf('   logPhiData=%.3f  logPhiReg=%.3f\n', ...
                S.logPhiData(idxB), S.logPhiReg(idxB));
    else
        fprintf('Candidate B: could not determine weighted L-curve corner.\n');
    end

    % --- 候选 C: L-curve 拐点 (未加权版) ---
    idxC_local = local_find_lcurve_corner(S_sorted.logPhiData_unweighted, S_sorted.logPhiReg_unweighted);
    if ~isnan(idxC_local)
        idxC = sort_idx(idxC_local);
        fprintf('Candidate C (unweighted L-curve corner):\n');
        fprintf('   run_id = %d | w501 = %.4e\n', S.run_id(idxC), S.w501(idxC));
        fprintf('   logPhiData_uw=%.3f  logPhiReg_uw=%.3f\n', ...
                S.logPhiData_unweighted(idxC), S.logPhiReg_unweighted(idxC));
    else
        fprintf('Candidate C: could not determine unweighted L-curve corner.\n');
    end

    % --- 图 1: 加权版 L-curve ---
    fig1 = figure('Visible', 'off', 'Name', 'L-curve (weighted)', 'Color', 'w', ...
                  'Position', [100 100 800 600]);
    plot(S_sorted.logPhiReg, S_sorted.logPhiData, 'o-', ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.4 0.8]);
    hold on;
    for j = 1:height(S_sorted)
        text(S_sorted.logPhiReg(j), S_sorted.logPhiData(j), ...
             sprintf('  %.1e', S_sorted.w501(j)), 'FontSize', 7);
    end
    if ~isnan(idxB_local)
        plot(S_sorted.logPhiReg(idxB_local), S_sorted.logPhiData(idxB_local), ...
             'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r');
    end
    xlabel('log_{10}(J_{reg})');
    ylabel('log_{10}(J_{data})');
    title(sprintf('L-curve (weighted): w101=%g, w103=%g', w101, w103));
    grid on; box on;
    saveas(fig1, fullfile(out_dir, 'Lcurve_plot.png'));
    close(fig1);

    % --- 图 2: 未加权版 L-curve ---
    fig2 = figure('Visible', 'off', 'Name', 'L-curve (unweighted)', 'Color', 'w', ...
                  'Position', [950 100 800 600]);
    plot(S_sorted.logPhiReg_unweighted, S_sorted.logPhiData_unweighted, 's-', ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.3 0.2]);
    hold on;
    for j = 1:height(S_sorted)
        text(S_sorted.logPhiReg_unweighted(j), S_sorted.logPhiData_unweighted(j), ...
             sprintf('  %.1e', S_sorted.w501(j)), 'FontSize', 7);
    end
    if ~isnan(idxC_local)
        plot(S_sorted.logPhiReg_unweighted(idxC_local), S_sorted.logPhiData_unweighted(idxC_local), ...
             'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r');
    end
    xlabel('log_{10}(\phi_{reg,unweighted})');
    ylabel('log_{10}(\phi_{data,unweighted})');
    title(sprintf('L-curve (unweighted): w101=%g, w103=%g', w101, w103));
    grid on; box on;
    saveas(fig2, fullfile(out_dir, 'Lcurve_plot_unweighted.png'));
    close(fig2);

    fprintf('\nDone for w101 = %g. Summary saved to:\n   %s\n   %s\n', ...
            w101, fullfile(out_dir, 'Lcurve_summary.mat'), fullfile(out_dir, 'Lcurve_summary.csv'));
end


function s = local_num_to_str(x)
% 适合 w101 / w103 这种一般整数权重
    if abs(x - round(x)) < 1e-12
        s = sprintf('%d', round(x));
    else
        s = strrep(sprintf('%.15g', x), '.', 'p');
    end
end


function s = local_num_to_expstr(x)
% 适合 w501 范围字符串
% 例如 1e-06 -> '1e-06', 5e-03 -> '5e-03'
    s = sprintf('%.0e', x);
end
