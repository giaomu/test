% =========================================================================
% runme_outletshead_sensitivity_cluster.m
%
% Cluster runner for one outlet-head sensitivity SHAKTI transient simulation.
%
% Put this script in:
%   /share/home/u10109/data/mumu/ISSM-Linux-MATLAB_1773134435/SHAKTI/outletshead_sensitivity_test
%
% Required input models:
%   Models/RecoveryNewMesh_Hydrology_Shelf0GroundedRingAlpha085_NoSliding.mat
%   Models/RecoveryNewMesh_Hydrology_Shelf0GroundedRingAlpha090_NoSliding.mat
%   Models/RecoveryNewMesh_Hydrology_Shelf0GroundedRingAlpha095_NoSliding.mat
%   Models/RecoveryNewMesh_Hydrology_Shelf0GroundedRingAlpha099_NoSliding.mat
%
% Main output:
%   Results/RecoveryNewMesh_SHAKTI_<case>_<start>d_to_<end>d_<dt>s_<tag>/
%   Results/.../segments/RecoveryNewMesh_SHAKTI_<case>_<seg_start>d_to_<seg_end>d_<dt>s_<tag>.mat
% =========================================================================
clear; clc; close all;

%% ========================================================================
%% User settings
%% ========================================================================

% --- Server paths ---
workdir     = '/share/home/u10109/data/mumu/ISSM-Linux-MATLAB_1773134435/SHAKTI/outletshead_sensitivity_test';
issm_root   = '/share/home/u10109/data/mumu/issm_local';
issm_src    = fullfile(issm_root, 'src', 'ISSM');
issm_prefix = fullfile(issm_root, 'prefix', 'issm');
issm_env_sh = fullfile(issm_root, 'env', 'issm_env.sh');
matlab_bin  = '/share/apps/MATLAB/R2023b/bin/matlab';

% --- Input and output files ---
% Change only case_model_name to switch among the four Step 7 models.
case_model_name = 'RecoveryNewMesh_Hydrology_Shelf0GroundedRingAlpha095_NoSliding.mat';
input_models_dir = fullfile(workdir, 'Models');
input_model_file = fullfile(input_models_dir, case_model_name);
runs_root        = fullfile(workdir, 'Results');
case_tag         = local_model_tag_from_file_name(case_model_name);

% --- Simulation settings ---
sim_dt_seconds  = 1800;       % Hydrology time step: 30 min.
sim_n_steps     = 2 * 24 *360 * 5;    % Step 7 setting: 2880 steps = 60 days.
sim_output_freq = 2 * 24 * 10;         % Save every 12 steps, i.e. every 6 hours.
checkpoint_save_days = 360;    % MAT checkpoint interval in model days.
sim_start_day = NaN;          % NaN = infer from the input model's last saved time.
hydro_relaxation = 1;         % SHAKTI nonlinear under-relaxation; NaN = keep input value.
run_tag = 'np15';             % Output suffix; empty = auto script or Slurm job tag.

sim_duration_days = sim_n_steps * sim_dt_seconds / 86400;
work_files_root = fullfile(workdir, 'workfiles');
execution_root = fullfile(work_files_root, 'execution');

% --- Cluster settings ---
requested_cores = 15;
require_slurm   = false;      % Set true to refuse non-Slurm runs.

% --- Runtime maintenance ---
repair_issm_links = true;     % Keep true if the ISSM install needs symlinks.

%% ========================================================================
%% Bootstrap ISSM environment if needed
%% ========================================================================

if ~strcmp(getenv('RUNME_OUTLETSHEAD_SENS_BOOTSTRAPPED'), '1')
    assert_exists(workdir,     'dir',  'workdir does not exist: %s');
    assert_exists(issm_env_sh, 'file', 'ISSM environment script does not exist: %s');
    assert_exists(matlab_bin,  'file', 'MATLAB executable does not exist: %s');

    if can_use_current_matlab(issm_root, issm_src, issm_prefix)
        setenv('RUNME_OUTLETSHEAD_SENS_BOOTSTRAPPED', '1');
    else
        this_script = get_this_script_fullpath();
        assert_exists(this_script, 'file', 'Cannot locate current script: %s');

        fprintf('Bootstrap: restarting MATLAB inside the ISSM environment.\n');
        fprintf('Script: %s\n', this_script);
        fprintf('Environment: %s\n', issm_env_sh);
        drawnow;

        wrapper_sh = fullfile(workdir, ...
            ['bootstrap_outletshead_sensitivity_' datestr(now, 'yyyymmdd_HHMMSSFFF') '.sh']);
        fid = fopen(wrapper_sh, 'w');
        if fid < 0
            error('Cannot create bootstrap wrapper: %s', wrapper_sh);
        end

        matlab_cmd = ['run(''' strrep(this_script, '''', '''''') ''')'];
        fprintf(fid, '#!/bin/bash\n');
        fprintf(fid, 'set -euo pipefail\n');
        fprintf(fid, 'source %s\n', shell_quote(issm_env_sh));
        fprintf(fid, 'export RUNME_OUTLETSHEAD_SENS_BOOTSTRAPPED=1\n');
        fprintf(fid, 'cd %s\n', shell_quote(workdir));
        fprintf(fid, 'if command -v stdbuf >/dev/null 2>&1; then\n');
        fprintf(fid, '  exec stdbuf -oL -eL %s -batch "%s"\n', shell_quote(matlab_bin), matlab_cmd);
        fprintf(fid, 'else\n');
        fprintf(fid, '  exec %s -batch "%s"\n', shell_quote(matlab_bin), matlab_cmd);
        fprintf(fid, 'fi\n');
        fclose(fid);

        system(sprintf('chmod +x %s', shell_quote(wrapper_sh)));
        st = system(sprintf('bash %s', shell_quote(wrapper_sh)), '-echo');
        if exist(wrapper_sh, 'file') == 2
            delete(wrapper_sh);
        end
        if st ~= 0
            error('Bootstrap subprocess failed. Check the output above.');
        end
        return;
    end
end

%% ========================================================================
%% ISSM runtime setup
%% ========================================================================

assert_exists(workdir,          'dir',  'workdir does not exist: %s');
assert_exists(issm_src,         'dir',  'ISSM source dir does not exist: %s');
assert_exists(issm_prefix,      'dir',  'ISSM prefix dir does not exist: %s');
assert_exists(input_model_file, 'file', 'Input model file does not exist: %s');

cd(workdir);
ensure_dir(work_files_root);
ensure_dir(execution_root);

issm_bin = fullfile(issm_prefix, 'bin');
issm_lib = fullfile(issm_prefix, 'lib');
issm_etc = fullfile(issm_src, 'etc');

gcc_libdir   = '/share/apps/gcc/11.4.0/lib64';
openmpi_lib  = '/share/apps/openmpi4.1.6/lib';
openblas_lib = '/share/apps/OpenBLAS/lib';

triangle_lib = fullfile(issm_src, 'externalpackages', 'triangle', 'install', 'lib');
petsc_lib    = fullfile(issm_src, 'externalpackages', 'petsc', 'install', 'lib');
m1qn3_lib    = fullfile(issm_src, 'externalpackages', 'm1qn3', 'install');
local_libs   = fullfile(issm_root, 'local_libs');

wrapper_matlab_dir  = fullfile(issm_src, 'src', 'wrappers', 'matlab');
wrapper_matlab_libs = fullfile(wrapper_matlab_dir, '.libs');

setenv('ISSM_ROOT',   issm_root);
setenv('ISSM_DIR',    issm_src);
setenv('ISSM_SRC',    issm_src);
setenv('ISSM_PREFIX', issm_prefix);

setenv('OMP_NUM_THREADS',      '1');
setenv('OPENBLAS_NUM_THREADS', '1');
setenv('MKL_NUM_THREADS',      '1');
setenv('OMPI_MCA_btl_openib_warn_no_device_params_found', '0');

ensure_dir(issm_lib);

if repair_issm_links
    link_shared_libs(triangle_lib, {'libtriangle.so', 'libtriangle.so.*'}, issm_lib);
    link_shared_libs(petsc_lib, ...
        {'libpetsc.so','libpetsc.so.*', ...
         'libscalapack.so','libscalapack.so.*', ...
         'libparmetis.so','libparmetis.so.*', ...
         'libmetis.so','libmetis.so.*'}, issm_lib);
    link_shared_libs(m1qn3_lib, ...
        {'libm1qn3.so','libm1qn3.so.*','libddot.so','libddot.so.*'}, issm_lib);
end

ld_candidates = {gcc_libdir, openmpi_lib, openblas_lib, ...
                 issm_lib, petsc_lib, triangle_lib, m1qn3_lib, local_libs};
ld_candidates = ld_candidates(cellfun(@(p) exist(p, 'dir') == 7, ld_candidates));
old_ld = getenv('LD_LIBRARY_PATH');
if ~isempty(old_ld)
    ld_candidates = [ld_candidates, strsplit(old_ld, ':')]; %#ok<AGROW>
end
ld_candidates = ld_candidates(~cellfun(@isempty, ld_candidates));
[~, ia] = unique(ld_candidates, 'stable');
setenv('LD_LIBRARY_PATH', strjoin(ld_candidates(sort(ia)), ':'));

assert_exists(issm_bin, 'dir', 'ISSM bin dir does not exist: %s');

if repair_issm_links
    ensure_required_exe(fullfile(issm_bin, 'issm.exe'),     fullfile(issm_src, 'issm.exe'));
    ensure_optional_exe(fullfile(issm_bin, 'issm_slc.exe'), fullfile(issm_src, 'issm_slc.exe'));
    ensure_optional_exe(fullfile(issm_bin, 'kriging.exe'),  fullfile(issm_src, 'kriging.exe'));

    env_sh = fullfile(issm_etc, 'environment.sh');
    assert_exists(env_sh, 'file', 'environment.sh does not exist: %s');
    prefix_etc = fullfile(issm_prefix, 'etc');
    ensure_dir(prefix_etc);
    ensure_symlink(env_sh, fullfile(prefix_etc, 'environment.sh'));
end

dirs_to_add = {issm_bin, issm_lib};
if exist(fullfile(issm_prefix, 'share'), 'dir') == 7
    dirs_to_add{end + 1} = fullfile(issm_prefix, 'share');
end
if exist(wrapper_matlab_dir, 'dir') == 7
    dirs_to_add{end + 1} = wrapper_matlab_dir;
end
if exist(wrapper_matlab_libs, 'dir') == 7
    dirs_to_add{end + 1} = wrapper_matlab_libs;
end

for i = 1:numel(dirs_to_add)
    if exist(dirs_to_add{i}, 'dir') == 7
        addpath(genpath(dirs_to_add{i}), '-begin');
    end
end
rehash;

if isempty(which('issmversion')) || isempty(which('generic'))
    error('ISSM MATLAB path loading failed.');
end

mexname = ['IssmConfig_matlab.' mexext];
mex_candidates = {};
w = which('IssmConfig_matlab');
if ~isempty(w), mex_candidates{end + 1} = w; end
fixed_paths = { ...
    fullfile(issm_prefix, 'lib', mexname), fullfile(issm_prefix, 'bin', mexname), ...
    fullfile(wrapper_matlab_libs, mexname), fullfile(wrapper_matlab_dir, mexname)};
for i = 1:numel(fixed_paths)
    if exist(fixed_paths{i}, 'file') == 2
        mex_candidates{end + 1} = fixed_paths{i};
    end
end
[~, ia] = unique(mex_candidates, 'stable');
mex_candidates = mex_candidates(sort(ia));
if isempty(mex_candidates)
    error('Cannot find IssmConfig_matlab MEX.');
end
addpath(fileparts(mex_candidates{1}), '-begin');
rehash;
IssmConfig('_HAVE_MPI_');

%% ========================================================================
%% Slurm / core count and logging
%% ========================================================================

slurm_job_id = getenv('SLURM_JOB_ID');
n_cores = local_read_slurm_cores();

if isempty(slurm_job_id)
    if require_slurm
        error('This script must run inside a Slurm job.');
    end
    warning('SLURM_JOB_ID is empty. Using requested_cores=%d.', requested_cores);
    n_cores = requested_cores;
elseif isnan(n_cores) || n_cores < 1
    warning('Could not read Slurm task count. Using requested_cores=%d.', requested_cores);
    n_cores = requested_cores;
end

if n_cores ~= requested_cores
    warning('Requested %d cores, but this run will use n_cores=%d.', requested_cores, n_cores);
end

%% ========================================================================
%% Load model, configure transient run, and solve
%% ========================================================================

S = load(input_model_file, 'md');
if ~isfield(S, 'md')
    error('Input file does not contain variable md: %s', input_model_file);
end
md = S.md;
clear S

inferred_start_day = infer_start_day_from_model(md);
if isnan(sim_start_day)
    sim_start_day = inferred_start_day;
end
sim_end_day = sim_start_day + sim_duration_days;
sim_period_tag = sprintf('%sd_to_%sd_%ds', ...
    local_num_token(sim_start_day), local_num_token(sim_end_day), sim_dt_seconds);
sim_result_base = sprintf('%s_SHAKTI_%s', case_tag, sim_period_tag);

if isempty(run_tag)
    run_tag = local_auto_run_tag();
end
run_tag = local_safe_tag(run_tag);
if isempty(run_tag)
    sim_result_name = sim_result_base;
else
    sim_result_name = sprintf('%s_%s', sim_result_base, run_tag);
end
run_dir = fullfile(runs_root, sim_result_name);
models_dir = run_dir;
segment_dir = fullfile(run_dir, 'segments');
logs_dir = fullfile(run_dir, 'logs');
checkpoint_dir = fullfile(run_dir, 'checkpoints');
ensure_dir(models_dir);
ensure_dir(segment_dir);
ensure_dir(logs_dir);
ensure_dir(checkpoint_dir);

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if isempty(slurm_job_id)
    job_tag = 'nojob';
else
    job_tag = slurm_job_id;
end
master_logfile = fullfile(logs_dir, ...
    sprintf('%s_job%s_%s.txt', sim_result_name, job_tag, timestamp));
diary(master_logfile);

logmsg('============================================================');
logmsg('Outlet-head sensitivity SHAKTI transient simulation started');
logmsg('workdir     : %s', workdir);
logmsg('case model  : %s', case_model_name);
logmsg('model tag   : %s', case_tag);
logmsg('input model : %s', input_model_file);
logmsg('models_dir  : %s', models_dir);
logmsg('segments    : %s', segment_dir);
logmsg('logs_dir    : %s', logs_dir);
logmsg('checkpoints : %s', checkpoint_dir);
logmsg('work files  : %s', work_files_root);
logmsg('execution   : %s', execution_root);
logmsg('Slurm job   : %s', job_tag);
logmsg('cores       : %d', n_cores);
logmsg('SLURM_JOB_NODELIST  : %s', getenv('SLURM_JOB_NODELIST'));
logmsg('SLURM_NNODES        : %s', getenv('SLURM_NNODES'));
logmsg('SLURM_NTASKS        : %s', getenv('SLURM_NTASKS'));
logmsg('SLURM_CPUS_ON_NODE  : %s', getenv('SLURM_CPUS_ON_NODE'));
logmsg('SLURM_CPUS_PER_TASK : %s', getenv('SLURM_CPUS_PER_TASK'));
logmsg('auto inferred start day: %g', inferred_start_day);
logmsg('start day   : %g', sim_start_day);
logmsg('end day     : %g', sim_end_day);
logmsg('dt seconds  : %g', sim_dt_seconds);
logmsg('steps       : %d', sim_n_steps);
logmsg('output freq : %d', sim_output_freq);
logmsg('checkpoint days: %g', checkpoint_save_days);
logmsg('hydrology relaxation: %g', hydro_relaxation);
logmsg('result name : %s', sim_result_name);
logmsg('============================================================');

if has_transient_solutions(md)
    md = apply_last_shakti_state(md, md.results.TransientSolution(end));
    md.results = struct();
end

local_validate_hydrology_model(md, input_model_file);

if ~isnan(hydro_relaxation)
    assert(hydro_relaxation >= 0, 'hydro_relaxation must be nonnegative or NaN.');
    md.hydrology.relaxation = hydro_relaxation;
end

md.miscellaneous.name = sim_result_name;
md.timestepping.time_step = sim_dt_seconds / md.constants.yts;
md.settings.output_frequency = sim_output_freq;

md.stressbalance.restol  = 0.05;
md.stressbalance.reltol  = 0.05;
md.stressbalance.abstol  = NaN;
md.stressbalance.maxiter = 200;

md.verbose.solution    = 1;
md.verbose.module      = 0;
md.verbose.convergence = 0;

solve_failed = false;
solve_error_identifier = '';
solve_error_message = '';
solve_error_report = '';
solve_error_stack = struct([]);
elapsed = 0;
all_transient = struct([]);

checkpoint_steps = max(1, round(checkpoint_save_days * 86400 / sim_dt_seconds));
checkpoint_steps = min(checkpoint_steps, sim_n_steps);
if checkpoint_save_days <= 0
    checkpoint_steps = sim_n_steps;
end
logmsg('checkpoint steps: %d', checkpoint_steps);

step_done = 0;
days_per_year = md.constants.yts / (24 * 3600);
current_time = sim_start_day / days_per_year;
segment_id = 0;
dt_yr = sim_dt_seconds / md.constants.yts;
dt_seconds = sim_dt_seconds;

while step_done < sim_n_steps
    segment_id = segment_id + 1;
    segment_steps = min(checkpoint_steps, sim_n_steps - step_done);
    segment_start_step = step_done + 1;
    segment_end_step = step_done + segment_steps;
    global_start_step = segment_start_step;
    global_end_step = segment_end_step;

    segment_start_day = sim_start_day + (segment_start_step - 1) * sim_dt_seconds / 86400;
    segment_end_day = sim_start_day + segment_end_step * sim_dt_seconds / 86400;
    segment_period_tag = sprintf('%sd_to_%sd_%ds', ...
        local_num_token(segment_start_day), local_num_token(segment_end_day), sim_dt_seconds);
    segment_base_name = sprintf('%s_SHAKTI_%s', case_tag, segment_period_tag);
    if isempty(run_tag)
        segment_name = segment_base_name;
    else
        segment_name = sprintf('%s_%s', segment_base_name, run_tag);
    end
    segment_work_dir = fullfile(work_files_root, segment_name);

    md.miscellaneous.name = segment_name;
    md.timestepping.start_time = current_time;
    md.timestepping.final_time = current_time + segment_steps * dt_yr;

    md.cluster = generic('name', oshostname(), 'np', n_cores);
    if local_has_member(md.cluster, 'codepath'),      md.cluster.codepath      = issm_src;       end
    if local_has_member(md.cluster, 'executionpath'), md.cluster.executionpath = execution_root; end
    if local_has_member(md.cluster, 'srcpath'),       md.cluster.srcpath       = issm_src;       end

    ensure_dir(segment_work_dir);
    cleanup_stale_artifacts(segment_work_dir, execution_root, md.miscellaneous.name);

    logmsg('------------------------------------------------------------');
    logmsg('Segment %03d: steps %d-%d, day %g -> %g, start %.10g yr, final %.10g yr', ...
        segment_id, segment_start_step, segment_end_step, ...
        segment_start_day, segment_end_day, ...
        md.timestepping.start_time, md.timestepping.final_time);
    logmsg('Segment name: %s', segment_name);
    logmsg('Segment work files: %s', segment_work_dir);

    previous_dir = pwd;
    cd(segment_work_dir);
    cwd_cleanup = onCleanup(@() cd(previous_dir));
    try
        tic_start = tic;
        md = solve(md, 'Transient');
        elapsed = elapsed + toc(tic_start);
    catch ME
        solve_failed = true;
        solve_error_identifier = ME.identifier;
        solve_error_message = ME.message;
        solve_error_report = getReport(ME, 'extended', 'hyperlinks', 'off');
        solve_error_stack = ME.stack;
        fprintf(2, '\nTransient solve failed in segment %d; saving checkpoint anyway.\n', segment_id);
        fprintf(2, 'Error: %s\n', solve_error_message);
    end
    delete(cwd_cleanup);
    clear cwd_cleanup
    cd(previous_dir);

    segment_has_results = false;
    try
        segment_has_results = isfield(md.results, 'TransientSolution') && ...
            ~isempty(md.results.TransientSolution);
    catch
        segment_has_results = false;
    end

    if segment_has_results
        segment_solutions = md.results.TransientSolution(:)';
        segment_n_saved = numel(segment_solutions);

        md.miscellaneous.name = segment_name;
        md.results.TransientSolution = segment_solutions;
        segment_model_file = fullfile(segment_dir, [segment_name '.mat']);
        save(segment_model_file, 'md', 'segment_id', ...
            'segment_start_day', 'segment_end_day', ...
            'global_start_step', 'global_end_step', 'dt_seconds', ...
            'segment_n_saved', '-v7.3');
        logmsg('Segment model saved: %s', segment_model_file);

        all_transient = append_transient_solutions(all_transient, segment_solutions);
        last_solution = segment_solutions(end);
        current_time = last_solution.time;
        step_done = segment_end_step;
        md = apply_last_shakti_state(md, last_solution);
        md.results.TransientSolution = all_transient;
    elseif ~solve_failed
        solve_failed = true;
        solve_error_message = 'Segment finished without TransientSolution.';
    end

    md.miscellaneous.name = sim_result_name;
    has_transient_results = ~isempty(all_transient);
    if has_transient_results
        last_time = all_transient(end).time;
        n_saved = numel(all_transient);
    else
        last_time = NaN;
        n_saved = 0;
    end

    checkpoint_file = fullfile(checkpoint_dir, ...
        sprintf('%s_checkpoint_seg%03d_%sd_to_%sd_step%05d.mat', ...
        sim_result_name, segment_id, ...
        local_num_token(segment_start_day), local_num_token(segment_end_day), step_done));
    latest_checkpoint_file = fullfile(checkpoint_dir, [sim_result_name '_checkpoint_latest.mat']);
    save(checkpoint_file, 'md', 'solve_failed', 'solve_error_identifier', ...
        'solve_error_message', 'solve_error_report', 'solve_error_stack', ...
        'has_transient_results', 'last_time', 'n_saved', 'elapsed', ...
        'segment_id', 'step_done', 'segment_start_day', 'segment_end_day', ...
        'global_start_step', 'global_end_step', 'dt_seconds', '-v7.3');
    save(latest_checkpoint_file, 'md', 'solve_failed', 'solve_error_identifier', ...
        'solve_error_message', 'solve_error_report', 'solve_error_stack', ...
        'has_transient_results', 'last_time', 'n_saved', 'elapsed', ...
        'segment_id', 'step_done', 'segment_start_day', 'segment_end_day', ...
        'global_start_step', 'global_end_step', 'dt_seconds', '-v7.3');
    logmsg('Checkpoint saved: %s', checkpoint_file);

    if solve_failed
        break;
    end

end

md.miscellaneous.name = sim_result_name;
if ~isempty(all_transient)
    md.results.TransientSolution = all_transient;
end

has_transient_results = ~isempty(all_transient);
if has_transient_results
    last_time = all_transient(end).time;
    n_saved = numel(all_transient);
else
    last_time = NaN;
    n_saved = 0;
end

diagnostic_file = fullfile(models_dir, [sim_result_name '.mat']);
save(diagnostic_file, 'md', 'solve_failed', 'solve_error_identifier', ...
    'solve_error_message', 'solve_error_report', 'solve_error_stack', ...
    'has_transient_results', 'last_time', 'n_saved', 'elapsed', ...
    'sim_start_day', 'sim_end_day', 'sim_dt_seconds', 'sim_n_steps', '-v7.3');

logmsg('============================================================');
logmsg('Solve failed     : %d', solve_failed);
logmsg('Elapsed seconds  : %.3f', elapsed);
logmsg('Saved time steps : %d', n_saved);
logmsg('Last time yr     : %.10g', last_time);
logmsg('Output model     : %s', diagnostic_file);
if ~has_transient_results
    logmsg('WARNING: no TransientSolution was returned.');
elseif last_time < md.timestepping.final_time - md.timestepping.time_step
    logmsg('WARNING: solve stopped early before final_time %.10g yr.', md.timestepping.final_time);
end
logmsg('============================================================');

diary('off');

if solve_failed
    error('Transient solve failed. Diagnostic model was saved: %s', diagnostic_file);
end


%% ========================================================================
%% Helper functions
%% ========================================================================

function ok = can_use_current_matlab(issm_root, issm_src, issm_prefix)
    ok = false;
    try
        wdir  = fullfile(issm_src, 'src', 'wrappers', 'matlab');
        wlibs = fullfile(wdir, '.libs');
        dirs  = {fullfile(issm_prefix, 'bin'), fullfile(issm_prefix, 'lib'), wdir, wlibs};
        for i = 1:numel(dirs)
            if exist(dirs{i}, 'dir') == 7
                addpath(genpath(dirs{i}), '-begin');
            end
        end
        rehash;
        setenv('ISSM_ROOT', issm_root);
        setenv('ISSM_DIR', issm_src);
        setenv('ISSM_SRC', issm_src);
        setenv('ISSM_PREFIX', issm_prefix);
        if ~isempty(which('IssmConfig'))
            IssmConfig('_HAVE_MPI_');
            ok = true;
        elseif ~isempty(which('IssmConfig_matlab'))
            IssmConfig('_HAVE_MPI_');
            ok = true;
        end
    catch
    end
end


function p = get_this_script_fullpath()
    p = mfilename('fullpath');
    if isempty(p)
        s = dbstack('-completenames');
        if ~isempty(s)
            p = s(1).file;
        end
    end
    if isempty(p)
        p = which(mfilename);
    end
    if isempty(p)
        error('Cannot locate current script path.');
    end
    if exist(p, 'file') ~= 2 && exist([p '.m'], 'file') == 2
        p = [p '.m'];
    end
end


function q = shell_quote(s)
    q = ['''' strrep(s, '''', '''"''"''') ''''];
end


function assert_exists(pathstr, kind, errmsg)
    switch lower(kind)
        case 'file'
            ok = exist(pathstr, 'file') == 2;
        case 'dir'
            ok = exist(pathstr, 'dir') == 7;
        otherwise
            error('assert_exists: unknown kind %s', kind);
    end
    if ~ok
        error(errmsg, pathstr);
    end
end


function ensure_dir(pathstr)
    if exist(pathstr, 'dir') ~= 7
        mkdir(pathstr);
    end
end


function link_shared_libs(srcdir, patterns, destdir)
    if exist(srcdir, 'dir') ~= 7
        return;
    end
    for ip = 1:numel(patterns)
        files = dir(fullfile(srcdir, patterns{ip}));
        for k = 1:numel(files)
            ensure_symlink(fullfile(files(k).folder, files(k).name), ...
                           fullfile(destdir, files(k).name));
        end
    end
end


function ensure_required_exe(src, dst)
    if exist(src, 'file') ~= 2
        error('Required executable does not exist: %s', src);
    end
    ensure_symlink(src, dst);
end


function ensure_optional_exe(src, dst)
    if exist(src, 'file') == 2
        ensure_symlink(src, dst);
    end
end


function ensure_symlink(src, dst)
    if exist(src, 'file') ~= 2
        return;
    end
    if exist(dst, 'file') == 2 || exist(dst, 'dir') == 7
        [ok, out] = system(sprintf('readlink -f %s', shell_quote(dst)));
        if ok == 0 && strcmp(strtrim(out), src)
            return;
        end
        system(sprintf('rm -f %s', shell_quote(dst)));
    end
    [st, out] = system(sprintf('ln -s %s %s', shell_quote(src), shell_quote(dst)));
    if st ~= 0
        error('Failed to create symlink:\nln -s %s %s\n%s', src, dst, out);
    end
end


function n = local_read_slurm_cores()
    n = str2double(getenv('SLURM_NTASKS'));
    if isnan(n) || n < 1
        n = str2double(getenv('SLURM_NPROCS'));
    end
    if isnan(n) || n < 1
        n = str2double(getenv('SLURM_CPUS_ON_NODE'));
    end
end


function tf = local_has_member(s, name)
    tf = (isstruct(s) && isfield(s, name)) || (isobject(s) && isprop(s, name));
end


function tag = local_auto_run_tag()
    job_id = getenv('SLURM_JOB_ID');
    [~, script_name] = fileparts(mfilename('fullpath'));
    if isempty(script_name)
        script_name = 'run';
    end
    if isempty(job_id)
        tag = script_name;
    else
        tag = sprintf('%s_job%s', script_name, job_id);
    end
end


function tag = local_safe_tag(tag)
    tag = strtrim(char(tag));
    tag = regexprep(tag, '[^A-Za-z0-9_\\-]+', '_');
    tag = regexprep(tag, '_+', '_');
    tag = regexprep(tag, '^_|_$', '');
end


function tag = local_model_tag_from_file_name(model_name)
    [~, base] = fileparts(model_name);
    tag = local_safe_tag(base);
    if isempty(tag)
        tag = 'Case';
    end
end


function start_day = infer_start_day_from_model(md)
    start_day = 0;
    try
        days_per_year = md.constants.yts / (24 * 3600);
    catch
        days_per_year = 365;
    end

    if has_transient_solutions(md)
        t = md.results.TransientSolution(end).time;
        if isfinite(t)
            start_day = t * days_per_year;
            return;
        end
    end

    try
        t0 = md.timestepping.start_time;
        if isfinite(t0)
            start_day = t0 * days_per_year;
        end
    catch
        start_day = 0;
    end
end


function tf = has_transient_solutions(md)
    tf = false;
    try
        tf = isfield(md.results, 'TransientSolution') && ...
            ~isempty(md.results.TransientSolution);
    catch
        tf = false;
    end
end


function s = local_num_token(x)
    if abs(x - round(x)) < 1e-10
        s = sprintf('%d', round(x));
    else
        s = sprintf('%.6g', x);
        s = strrep(s, '.', 'p');
        s = strrep(s, '-', 'm');
    end
end


function cleanup_stale_artifacts(workdir, execution_root, run_name)
    exts = {'.bin', '.queue', '.toolkits', '.outbin', '.outlog', '.errlog'};
    for i = 1:numel(exts)
        p = fullfile(workdir, [run_name exts{i}]);
        if exist(p, 'file') == 2
            delete(p);
        end
    end
    if exist(execution_root, 'dir') == 7
        d = dir(fullfile(execution_root, [run_name '*']));
        for i = 1:numel(d)
            t = fullfile(d(i).folder, d(i).name);
            if d(i).isdir && ~ismember(d(i).name, {'.', '..'})
                system(sprintf('rm -rf %s', shell_quote(t)));
            elseif ~d(i).isdir
                delete(t);
            end
        end
    end
end


function local_validate_hydrology_model(md, input_file)
    required = { ...
        'mesh', 'geometry', 'hydrology', 'basalforcings', ...
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
    spc = md.hydrology.spchead(:);
    finite_spc = isfinite(spc);
    fprintf('transient.ishydrology     = %d\n', md.transient.ishydrology);
    if local_has_member(md.transient, 'isstressbalance')
        fprintf('transient.isstressbalance = %d\n', md.transient.isstressbalance);
    end
    if local_has_member(md, 'friction') && local_has_member(md.friction, 'coupling')
        fprintf('friction.coupling         = %g\n', md.friction.coupling);
    end
    if local_has_member(md.hydrology, 'bump_height')
        bh = md.hydrology.bump_height(:);
        fprintf('bump_height: min %.4g, max %.4g, all_zero %d\n', ...
            min(bh), max(bh), all(abs(bh) <= 1e-12));
    end
    fprintf('spchead nodes: %d\n', sum(finite_spc));
    fprintf('spchead zero nodes: %d\n', sum(finite_spc & abs(spc) <= 1e-9));
    fprintf('spchead negative nodes: %d\n', sum(finite_spc & spc < 0));
    if any(finite_spc)
        fprintf('spchead finite min %.4g, max %.4g m\n', ...
            min(spc(finite_spc)), max(spc(finite_spc)));
    end
    fprintf('gap_height: min %.4g, max %.4g m\n', ...
        min(md.hydrology.gap_height(:)), max(md.hydrology.gap_height(:)));
    fprintf('geothermal flux: min %.4g, max %.4g W/m^2\n', ...
        min(md.basalforcings.geothermalflux(:)), max(md.basalforcings.geothermalflux(:)));
end


function all_solutions = append_transient_solutions(all_solutions, new_solutions)
    if isempty(new_solutions)
        return;
    end
    new_solutions = new_solutions(:)';
    if isempty(all_solutions)
        all_solutions = new_solutions;
    else
        all_solutions = [all_solutions(:)' new_solutions]; %#ok<AGROW>
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


function check_len(x, n, name)
    if numel(x) ~= n
        error('%s length is %d, expected %d.', name, numel(x), n);
    end
end


function logmsg(varargin)
    fprintf('[%s] ', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(varargin{:});
    fprintf('\n');
    drawnow;
end
