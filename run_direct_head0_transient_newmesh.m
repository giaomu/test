% =========================================================================
% run_direct_head0_transient_newmesh.m
%
% 直接对已保存的 Head0 setup 模型运行 SHAKTI 瞬态模拟。
% 用途：
%   1) 不重建 shelf / grounded ring；
%   2) 不重新写 spchead；
%   3) 不跑 Step 7 的 Head0 + alpha case 循环；
%   4) 只验证 RecoveryNewMesh_Hydrology_ShelfGroundedRingNoSliding.mat
%      这个旧 Head0 模型本身是否还能稳定跑瞬态。
%
% 当前瞬态设置与 runme_newmesh.m 的 Step 7 保持一致。
% =========================================================================

clear; clc; close all;

paths = shaktiais_paths();
if ~exist(paths.models_newmesh, 'dir'), mkdir(paths.models_newmesh); end

% -------------------------------------------------------------------------
% 1. 输入输出模型
% -------------------------------------------------------------------------
input_model_name = 'Hydrology_ShelfGroundedRingNoSliding';
input_model_file = fullfile(paths.models_newmesh, ['RecoveryNewMesh_' input_model_name '.mat']);
result_model_name = 'SHAKTI_DirectHead0_480Steps_1800s_noslide';
result_file = fullfile(paths.models_newmesh, ['RecoveryNewMesh_' result_model_name '.mat']);

fprintf('\n===== Direct Head0 transient check =====\n');
fprintf('Input model:  RecoveryNewMesh_%s.mat\n', input_model_name);
fprintf('Output model: RecoveryNewMesh_%s.mat\n', result_model_name);

% -------------------------------------------------------------------------
% 2. 读取旧 Head0 setup 模型
% -------------------------------------------------------------------------
assert(exist(input_model_file, 'file') == 2, ...
    'Input model file does not exist: %s', input_model_file);
S = load(input_model_file, 'md');
assert(isfield(S, 'md'), 'Input file does not contain variable md: %s', input_model_file);
md = S.md;
clear S

if has_transient_solutions(md)
    md = apply_last_shakti_state(md, md.results.TransientSolution(end));
    md.results = struct();
end

% -------------------------------------------------------------------------
% 3. 运行前诊断：确认这确实是关闭滑动开腔、单向水文的 Head0 setup
% -------------------------------------------------------------------------
local_validate_hydrology_model(md, input_model_file);

spc = md.hydrology.spchead(:);
finite_spc = isfinite(spc);
zero_spc = finite_spc & abs(spc) <= 1e-9;

if local_has_member(md.hydrology, 'bump_height') && ~isempty(md.hydrology.bump_height)
    bump_height = md.hydrology.bump_height(:);
    bump_min = min(bump_height);
    bump_max = max(bump_height);
    bump_all_zero = all(abs(bump_height) <= 1e-12);
else
    bump_min = NaN;
    bump_max = NaN;
    bump_all_zero = false;
end

fprintf('\n--- Pre-run diagnostics ---\n');
fprintf('bump_height: min = %.6g, max = %.6g, all_zero = %d\n', ...
    bump_min, bump_max, bump_all_zero);
fprintf('transient.ishydrology     = %d\n', md.transient.ishydrology);
fprintf('transient.isstressbalance = %d\n', md.transient.isstressbalance);
fprintf('friction.coupling         = %d\n', md.friction.coupling);
fprintf('finite spchead nodes      = %d\n', sum(finite_spc));
fprintf('spchead = 0 nodes         = %d\n', sum(zero_spc));
fprintf('spchead finite min/max    = %.6g / %.6g m\n', ...
    finite_min(spc(finite_spc)), finite_max(spc(finite_spc)));

% -------------------------------------------------------------------------
% 4. 与 Step 7 对应的瞬态模拟设置
% -------------------------------------------------------------------------
step7_sim_start_time = 0;
step7_sim_dt_seconds = 30 * 60;
step7_sim_n_steps = 48 * 10;
step7_sim_final_days = step7_sim_n_steps * step7_sim_dt_seconds / 86400;
step7_sim_output_freq = 12;

step7_stress_restol = 0.05;
step7_stress_reltol = 0.05;
step7_stress_abstol = NaN;
step7_stress_maxiter = 200;

step7_verbose_solution = 1;
step7_verbose_module = 0;
step7_verbose_convergence = 0;
step7_n_cores = 15;

% 与之前成功的 cluster 脚本保持一致：显式使用 SHAKTI 非线性松弛系数 1。
% 如果设为 NaN，则保留输入模型中的 md.hydrology.relaxation。
step7_hydro_relaxation = 1;
treat_early_stop_as_failure = true;

fprintf('\n--- Transient settings copied from Step 7 ---\n');
fprintf('dt = %.1f s, n_steps = %d, final_days = %.6g, output_frequency = %d\n', ...
    step7_sim_dt_seconds, step7_sim_n_steps, step7_sim_final_days, step7_sim_output_freq);
fprintf('stressbalance restol/reltol/abstol/maxiter = %.4g / %.4g / %.4g / %d\n', ...
    step7_stress_restol, step7_stress_reltol, step7_stress_abstol, step7_stress_maxiter);
fprintf('hydrology relaxation = %.6g\n', step7_hydro_relaxation);
fprintf('np = %d\n', step7_n_cores);

% -------------------------------------------------------------------------
% 5. 设置瞬态参数
% -------------------------------------------------------------------------
md.miscellaneous.name = ['RecoveryNewMesh_' result_model_name];
md.timestepping.start_time = step7_sim_start_time;
md.timestepping.time_step = step7_sim_dt_seconds / md.constants.yts;
md.timestepping.final_time = step7_sim_final_days / 365;
md.settings.output_frequency = step7_sim_output_freq;

if ~isnan(step7_hydro_relaxation)
    assert(step7_hydro_relaxation >= 0, 'step7_hydro_relaxation must be nonnegative or NaN.');
    if local_has_member(md.hydrology, 'relaxation')
        md.hydrology.relaxation = step7_hydro_relaxation;
    else
        warning('md.hydrology.relaxation is missing; cannot set relaxation.');
    end
end

md.stressbalance.restol = step7_stress_restol;
md.stressbalance.reltol = step7_stress_reltol;
md.stressbalance.abstol = step7_stress_abstol;
md.stressbalance.maxiter = step7_stress_maxiter;

md.verbose.solution = step7_verbose_solution;
md.verbose.module = step7_verbose_module;
md.verbose.convergence = step7_verbose_convergence;
md.cluster = generic('name', oshostname(), 'np', step7_n_cores);

% -------------------------------------------------------------------------
% 6. 求解。失败时也保存诊断模型和错误信息。
% -------------------------------------------------------------------------
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
    fprintf(2, '\nDirect Head0 transient solve failed; saving diagnostics anyway.\n');
    fprintf(2, 'Error: %s\n', solve_error_message);
end

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

if ~has_transient_results
    solve_failed = true;
    if isempty(solve_error_message)
        solve_error_message = 'No TransientSolution was returned.';
    end
elseif treat_early_stop_as_failure && ...
        last_time < md.timestepping.final_time - md.timestepping.time_step
    solve_failed = true;
    if isempty(solve_error_message)
        solve_error_message = sprintf( ...
            'Transient solve stopped early at %.9g yr, before final_time %.9g yr.', ...
            last_time, md.timestepping.final_time);
    end
end

save(result_file, 'md', ...
    'input_model_name', 'input_model_file', 'result_model_name', ...
    'step7_sim_start_time', 'step7_sim_dt_seconds', 'step7_sim_n_steps', ...
    'step7_sim_final_days', 'step7_sim_output_freq', ...
    'step7_stress_restol', 'step7_stress_reltol', 'step7_stress_abstol', ...
    'step7_stress_maxiter', 'step7_verbose_solution', 'step7_verbose_module', ...
    'step7_verbose_convergence', 'step7_n_cores', 'step7_hydro_relaxation', ...
    'treat_early_stop_as_failure', ...
    'solve_failed', 'solve_error_identifier', 'solve_error_message', ...
    'solve_error_report', 'solve_error_stack', ...
    'has_transient_results', 'last_time', '-v7.3');

fprintf('\nSaved diagnostic result: %s\n', result_file);
fprintf('solve_failed = %d, has_transient_results = %d, last_time = %.9g yr\n', ...
    solve_failed, has_transient_results, last_time);

if ~has_transient_results
    fprintf(2, 'Warning: no TransientSolution was returned.\n');
elseif last_time < md.timestepping.final_time - md.timestepping.time_step
    fprintf(2, 'Warning: solve stopped early at %.9g yr, before final_time %.9g yr.\n', ...
        last_time, md.timestepping.final_time);
end

if solve_failed
    fprintf('Done with failed solve; diagnostic model and error info were saved.\n');
    error('Direct Head0 transient failed or stopped early. Diagnostic model was saved: %s', result_file);
else
    fprintf('Done.\n');
end


function tf = local_has_member(s, name)
    tf = (isstruct(s) && isfield(s, name)) || (isobject(s) && isprop(s, name));
end


function local_validate_hydrology_model(md, input_file)
    required = {'mesh', 'geometry', 'hydrology', 'basalforcings', ...
        'materials', 'constants', 'transient'};
    for i = 1:numel(required)
        if ~local_has_member(md, required{i})
            error('Input model is missing md.%s: %s', required{i}, input_file);
        end
    end

    nv = md.mesh.numberofvertices;
    ne = md.mesh.numberofelements;
    check_len(md.geometry.thickness, nv, 'md.geometry.thickness');
    check_len(md.geometry.base, nv, 'md.geometry.base');
    check_len(md.hydrology.head, nv, 'md.hydrology.head');
    check_len(md.hydrology.spchead, nv, 'md.hydrology.spchead');
    check_len(md.hydrology.gap_height, ne, 'md.hydrology.gap_height');

    assert(all(isfinite(md.hydrology.head(:))), 'md.hydrology.head contains NaN/Inf.');
    assert(all(isfinite(md.hydrology.gap_height(:))), 'md.hydrology.gap_height contains NaN/Inf.');
    assert(all(isfinite(md.hydrology.storage(:))) && all(md.hydrology.storage(:) >= 0), ...
        'md.hydrology.storage must be finite and nonnegative.');
    assert(all(isfinite(md.basalforcings.geothermalflux(:))) && ...
        all(md.basalforcings.geothermalflux(:) >= 0), ...
        'md.basalforcings.geothermalflux contains NaN/Inf or negative values.');

    if ~local_has_member(md.hydrology, 'melt_flag')
        error('md.hydrology.melt_flag is missing.');
    end
    if md.hydrology.melt_flag ~= 0
        warning('md.hydrology.melt_flag is %g. For friction/geothermal melt, expected 0.', ...
            md.hydrology.melt_flag);
    end
    if md.transient.ishydrology ~= 1
        error('md.transient.ishydrology must be 1.');
    end

    fprintf('Hydrology model validated: %d vertices, %d elements.\n', nv, ne);
    fprintf('gap_height: min %.4g, max %.4g m\n', ...
        min(md.hydrology.gap_height(:)), max(md.hydrology.gap_height(:)));
    fprintf('geothermal flux: min %.4g, max %.4g W/m^2\n', ...
        min(md.basalforcings.geothermalflux(:)), max(md.basalforcings.geothermalflux(:)));
end


function check_len(x, n, name)
    if numel(x) ~= n
        error('%s length is %d, expected %d.', name, numel(x), n);
    end
end


function tf = has_transient_solutions(md)
    try
        tf = isfield(md.results, 'TransientSolution') && ...
            ~isempty(md.results.TransientSolution);
    catch
        tf = false;
    end
end


function md = apply_last_shakti_state(md, S)
    if isfield(S, 'HydrologyHead')
        md.hydrology.head = S.HydrologyHead;
    else
        warning('Cannot update md.hydrology.head: HydrologyHead missing in last solution.');
    end

    if isfield(S, 'HydrologyGapHeight')
        md.hydrology.gap_height = S.HydrologyGapHeight;
    else
        warning('Cannot update md.hydrology.gap_height: HydrologyGapHeight missing in last solution.');
    end
end


function v = finite_min(values)
    values = values(isfinite(values));
    if isempty(values), v = NaN; else, v = min(values); end
end


function v = finite_max(values)
    values = values(isfinite(values));
    if isempty(values), v = NaN; else, v = max(values); end
end
