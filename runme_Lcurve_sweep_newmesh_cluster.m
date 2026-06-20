% =========================================================================
% runme_Lcurve_sweep_newmesh_cluster.m
%
% New-mesh Recovery Glacier L-curve sweep on the cluster.
%
% Output layout follows the local runme_inversion_test.m convention:
%   outputs/models/L_curve_models_<w101>_<w103>_<w501min>_to_<w501max>/
%   outputs/models/Lcurve_summary_all.mat
%   outputs/models/Lcurve_summary_all.csv
% =========================================================================
clear; clc; close all;

%% ========================================================================
%% User settings
%% ========================================================================

% --- Server paths ---
workdir     = '/share/home/u10109/data/mumu/ISSM-Linux-MATLAB_1773134435/inversion_test';
issm_root   = '/share/home/u10109/data/mumu/issm_local';
issm_src    = fullfile(issm_root, 'src', 'ISSM');
issm_prefix = fullfile(issm_root, 'prefix', 'issm');
issm_env_sh = fullfile(issm_root, 'env', 'issm_env.sh');
matlab_bin  = '/share/apps/MATLAB/R2023b/bin/matlab';

% --- Input files ---
input_model_file = fullfile(workdir, 'outputs', 'models_newmesh', ...
    'RecoveryNewMesh_Parameterize.mat');
bm_file = fullfile(workdir, 'data', ...
    'NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc');

% --- Output roots ---
models_dir     = fullfile(workdir, 'outputs', 'models');
execution_root = fullfile(workdir, 'execution');

% --- L-curve weights ---
w101_list = [500 2000 5000 10000 20000 40000];
w103 = 10;
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

% --- Inversion settings ---
inv_maxsteps     = 200;
inv_maxiter      = 400;
inv_friction_min = 0.05;
inv_friction_max = 3500;
target_ratio     = [10, 5, 2] / 17;

% Match runme_newmesh.m: shelf is excluded, but no-ice data terms are not
% zeroed unless this switch is set true.
zero_noice_data_terms = false;

% --- Cluster settings ---
requested_cores = 32;
require_slurm   = false;  % Set true if this script must refuse non-Slurm runs.

% --- Safety / maintenance settings ---
repair_issm_links = true; % Keep true if the ISSM install needs runtime symlinks.
test_mode = false;        % If true, run a small 1 x 2 sweep for path testing.

if test_mode
    w101_list = w101_list(1);
    w501_list = w501_list(1:min(2, numel(w501_list)));
end

%% ========================================================================
%% Bootstrap ISSM environment if needed
%% ========================================================================

if ~strcmp(getenv('RUNME_LCURVE_NEWMESH_BOOTSTRAPPED'), '1')
    assert_exists(workdir,     'dir',  'workdir does not exist: %s');
    assert_exists(issm_env_sh, 'file', 'ISSM environment script does not exist: %s');
    assert_exists(matlab_bin,  'file', 'MATLAB executable does not exist: %s');

    if can_use_current_matlab(issm_root, issm_src, issm_prefix)
        setenv('RUNME_LCURVE_NEWMESH_BOOTSTRAPPED', '1');
    else
        this_script = get_this_script_fullpath();
        assert_exists(this_script, 'file', 'Cannot locate current script: %s');

        fprintf('Bootstrap: restarting MATLAB inside the ISSM environment.\n');
        fprintf('Script: %s\n', this_script);
        fprintf('Environment: %s\n', issm_env_sh);
        drawnow;

        wrapper_sh = fullfile(workdir, ...
            ['bootstrap_lcurve_newmesh_' datestr(now, 'yyyymmdd_HHMMSSFFF') '.sh']);
        fid = fopen(wrapper_sh, 'w');
        if fid < 0
            error('Cannot create bootstrap wrapper: %s', wrapper_sh);
        end

        matlab_cmd = ['run(''' strrep(this_script, '''', '''''') ''')'];
        fprintf(fid, '#!/bin/bash\n');
        fprintf(fid, 'set -euo pipefail\n');
        fprintf(fid, 'source %s\n', shell_quote(issm_env_sh));
        fprintf(fid, 'export RUNME_LCURVE_NEWMESH_BOOTSTRAPPED=1\n');
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
%% Path checks and ISSM runtime setup
%% ========================================================================

assert_exists(workdir,          'dir',  'workdir does not exist: %s');
assert_exists(issm_src,         'dir',  'ISSM source dir does not exist: %s');
assert_exists(issm_prefix,      'dir',  'ISSM prefix dir does not exist: %s');
assert_exists(input_model_file, 'file', 'Input model file does not exist: %s');
assert_exists(bm_file,          'file', 'BedMachine file does not exist: %s');

cd(workdir);
if ~exist(models_dir, 'dir'), mkdir(models_dir); end
if ~exist(execution_root, 'dir'), mkdir(execution_root); end

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

if ~exist(issm_lib, 'dir'), mkdir(issm_lib); end

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
    if ~exist(prefix_etc, 'dir'), mkdir(prefix_etc); end
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

fprintf('ISSM_DIR       = %s\n', getenv('ISSM_DIR'));
fprintf('ISSM_PREFIX    = %s\n', getenv('ISSM_PREFIX'));
fprintf('issmversion at = %s\n', which('issmversion'));
fprintf('generic at     = %s\n\n', which('generic'));

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
fprintf('IssmConfig MEX = %s\n', mex_candidates{1});
try
    have_mpi = IssmConfig('_HAVE_MPI_');
    fprintf('IssmConfig loaded, _HAVE_MPI_ = %g\n\n', have_mpi);
catch ME
    error('IssmConfig MEX load failed: %s', ME.message);
end

%% ========================================================================
%% Slurm / core count
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
%% Master log
%% ========================================================================

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if isempty(slurm_job_id)
    job_tag = 'nojob';
else
    job_tag = slurm_job_id;
end
master_logfile = fullfile(models_dir, ...
    sprintf('Lcurve_newmesh_multiw101_job%s_%s.txt', job_tag, timestamp));
diary(master_logfile);

logmsg('============================================================');
logmsg('New-mesh multi-w101 L-curve sweep started');
logmsg('workdir     : %s', workdir);
logmsg('input model : %s', input_model_file);
logmsg('BedMachine  : %s', bm_file);
logmsg('models_dir  : %s', models_dir);
logmsg('Slurm job   : %s', job_tag);
logmsg('cores       : %d', n_cores);
logmsg('w101 values : %s', mat2str(w101_list));
logmsg('w103        : %g', w103);
logmsg('w501 count  : %d', numel(w501_list));
logmsg('test_mode   : %d', test_mode);
logmsg('============================================================');

%% ========================================================================
%% Load base model and classify nodes
%% ========================================================================

logmsg('Loading base model...');
S0 = load(input_model_file, 'md');
if ~isfield(S0, 'md')
    error('Input file does not contain variable md: %s', input_model_file);
end
md_base = S0.md;
clear S0
local_validate_base_model(md_base, input_model_file);

logmsg('Classifying nodes from BedMachine mask...');
[is_grounded, is_shelf, is_noice, is_vostok] = ...
    local_extract_node_types_from_bedmachine( ...
        bm_file, md_base.mesh.x, md_base.mesh.y, md_base.mesh.elements);
logmsg('After cliff suppression: grounded=%d shelf=%d no-ice=%d vostok=%d', ...
    sum(is_grounded), sum(is_shelf), sum(is_noice), sum(is_vostok));

%% ========================================================================
%% Summary tables
%% ========================================================================

n_w101 = numel(w101_list);
n_w501 = numel(w501_list);
n_runs_total = n_w101 * n_w501;

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

logmsg('Total runs: %d w101 values x %d w501 values = %d runs', ...
    n_w101, n_w501, n_runs_total);

%% ========================================================================
%% Main sweep
%% ========================================================================

global_start = tic;
global_run_id = 0;

for i101 = 1:n_w101

    w101 = w101_list(i101);
    out_dir = fullfile(models_dir, sprintf('L_curve_models_%s_%s_%s_to_%s', ...
        local_num_to_str(w101), ...
        local_num_to_str(w103), ...
        local_num_to_expstr(w501_list(1)), ...
        local_num_to_expstr(w501_list(end))));

    log_dir = fullfile(out_dir, 'logs');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    if ~exist(log_dir, 'dir'), mkdir(log_dir); end

    logmsg('============================================================');
    logmsg('Sweeping w101 = %g', w101);
    logmsg('Output dir: %s', out_dir);
    logmsg('============================================================');

    summary_one = local_empty_table(n_w501, col_names, target_ratio);

    for k = 1:n_w501

        global_run_id = global_run_id + 1;
        w501 = w501_list(k);

        safe_name = local_make_safe_filename(w501);
        model_fname = sprintf('Lcurve_run%02d_w501_%s.mat', k, safe_name);
        log_fname   = sprintf('Lcurve_run%02d_w501_%s.txt', k, safe_name);
        model_path  = fullfile(out_dir, model_fname);
        log_path    = fullfile(log_dir, log_fname);

        logmsg('--- Global %03d/%03d | w101 = %.4g | w501 = %.4e ---', ...
            global_run_id, n_runs_total, w101, w501);

        row = local_empty_row(global_run_id, w101, w103, w501, ...
            model_fname, log_fname, target_ratio, col_names);

        try
            md = md_base;
            nv = md.mesh.numberofvertices;

            md.miscellaneous.name = sprintf('Lcurve_w101_%s_run%02d_w501_%s', ...
                local_num_to_str(w101), k, safe_name);

            md.inversion = m1qn3inversion(md.inversion);
            md.inversion.iscontrol = 1;
            md.transient.amr_frequency = 0;
            md.verbose = verbose('solution', true, 'control', true, 'convergence', true);

            md.inversion.cost_functions = [101, 103, 501];
            md = local_apply_cost_masks(md, is_shelf, is_noice, w101, w103, w501, ...
                inv_friction_min, inv_friction_max, zero_noice_data_terms);

            md.inversion.control_parameters      = {'FrictionCoefficient'};
            md.inversion.maxsteps                = inv_maxsteps;
            md.inversion.maxiter                 = inv_maxiter;
            md.inversion.control_scaling_factors = 1;

            md.stressbalance.restol = 0.01;
            md.stressbalance.reltol = 0.1;
            md.stressbalance.abstol = NaN;

            md.cluster = generic('name', oshostname(), 'np', n_cores);
            if local_has_member(md.cluster, 'codepath'),      md.cluster.codepath      = issm_src;       end
            if local_has_member(md.cluster, 'executionpath'), md.cluster.executionpath = execution_root; end
            if local_has_member(md.cluster, 'srcpath'),       md.cluster.srcpath       = issm_src;       end

            cleanup_stale_artifacts(workdir, execution_root, md.miscellaneous.name);

            tic_start = tic;
            cmdout = evalc('md = solve(md, ''sb'');');
            elapsed = toc(tic_start);

            fid = fopen(log_path, 'w');
            if fid > 0
                fwrite(fid, cmdout);
                fclose(fid);
            end

            assert(local_has_member(md, 'results') && local_has_member(md.results, 'StressbalanceSolution'), ...
                'No StressbalanceSolution');

            sol = md.results.StressbalanceSolution;
            Jhist = sol.J;
            n_iter = size(Jhist, 1);
            J101_val   = Jhist(end, 1);
            J103_val   = Jhist(end, 2);
            J501_val   = Jhist(end, 3);
            Jtotal_val = Jhist(end, 4);
            Jdata_val  = J101_val + J103_val;

            phi_data_uw = J101_val / w101 + J103_val / w103;
            phi_reg_uw  = J501_val / w501;

            [~, ~, grad_norm, ~, ~, ~] = local_parse_final_m1qn3_line(cmdout);

            stop_flag = NaN;
            if local_has_member(sol, 'InversionStopFlag')
                stop_flag = sol.InversionStopFlag;
            end

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

            md.friction.coefficient = sol.FrictionCoefficient;
            md.friction.coefficient(is_shelf) = 0;
            md.initialization.vx  = sol.Vx;
            md.initialization.vy  = sol.Vy;
            md.initialization.vz  = zeros(nv, 1);
            md.initialization.vel = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);

            save(model_path, 'md', '-v7.3');

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

            logmsg('OK %.1fs | %d iters | J=[%.2e %.2e %.2e] total=%.2e | rdist=%.4f', ...
                elapsed, n_iter, J101_val, J103_val, J501_val, Jtotal_val, rdist);

        catch ME
            row.status        = {'failed'};
            row.error_message = {ME.message};
            logmsg('FAILED: %s', ME.message);
        end

        summary_one(k, :) = row;
        summary_all(global_run_id, :) = row;

        save(fullfile(out_dir, 'Lcurve_summary.mat'), 'summary_one');
        writetable(summary_one, fullfile(out_dir, 'Lcurve_summary.csv'));

        save(fullfile(models_dir, 'Lcurve_summary_all.mat'), 'summary_all');
        writetable(summary_all, fullfile(models_dir, 'Lcurve_summary_all.csv'));
    end

    local_postprocess_one_summary(summary_one, out_dir, w101, w103);
end

elapsed_all = toc(global_start);
logmsg('============================================================');
logmsg('Done. Overall summary saved to:');
logmsg('  %s', fullfile(models_dir, 'Lcurve_summary_all.mat'));
logmsg('  %s', fullfile(models_dir, 'Lcurve_summary_all.csv'));
logmsg('Master log: %s', master_logfile);
logmsg('Total elapsed: %.2f s (%.2f h)', elapsed_all, elapsed_all / 3600);
logmsg('============================================================');
diary('off');

%% ========================================================================
%% Helper functions: bootstrap and environment
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
        end
    catch
        ok = false;
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
            error('Unknown assert_exists kind: %s', kind);
    end
    if ~ok
        error(errmsg, pathstr);
    end
end


function link_shared_libs(srcdir, patterns, destdir)
    if ~exist(srcdir, 'dir')
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
        error('Failed to create symlink: ln -s %s %s\n%s', src, dst, out);
    end
end


function n_cores = local_read_slurm_cores()
    candidates = {'SLURM_NTASKS', 'SLURM_NPROCS', 'SLURM_CPUS_ON_NODE'};
    n_cores = NaN;
    for i = 1:numel(candidates)
        v = str2double(getenv(candidates{i}));
        if isfinite(v) && v >= 1
            n_cores = v;
            return;
        end
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


function logmsg(varargin)
    fprintf('[%s] ', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(varargin{:});
    fprintf('\n');
    drawnow;
end

%% ========================================================================
%% Helper functions: model, inversion, and summaries
%% ========================================================================

function local_validate_base_model(md, model_file)
    nv = md.mesh.numberofvertices;
    checks = {
        'md.mesh.x',               md.mesh.x,                 nv
        'md.mesh.y',               md.mesh.y,                 nv
        'md.geometry.thickness',   md.geometry.thickness,     nv
        'md.inversion.vx_obs',     md.inversion.vx_obs,       nv
        'md.inversion.vy_obs',     md.inversion.vy_obs,       nv
        'md.inversion.vel_obs',    md.inversion.vel_obs,      nv
        'md.initialization.vx',    md.initialization.vx,      nv
        'md.initialization.vy',    md.initialization.vy,      nv
        'md.initialization.vel',   md.initialization.vel,     nv
        'md.friction.coefficient', md.friction.coefficient,   nv
        'md.mask.ice_levelset',    md.mask.ice_levelset,      nv
        'md.mask.ocean_levelset',  md.mask.ocean_levelset,    nv
        };

    for i = 1:size(checks, 1)
        name = checks{i, 1};
        vals = checks{i, 2};
        expected = checks{i, 3};
        if numel(vals) ~= expected
            error('%s length mismatch in %s: expected %d, got %d', ...
                name, model_file, expected, numel(vals));
        end
        if any(~isfinite(vals(:)))
            error('%s contains NaN/Inf in %s', name, model_file);
        end
    end
end


function [is_grounded, is_shelf, is_noice, is_vostok] = ...
        local_extract_node_types_from_bedmachine(bm_file, mesh_x, mesh_y, elements)
    x_bm    = double(ncread(bm_file, 'x'));
    y_bm    = double(ncread(bm_file, 'y'));
    mask_bm = double(ncread(bm_file, 'mask'));

    [x_bm, ix] = sort(x_bm);
    [y_bm, iy] = sort(y_bm);
    mask_bm = mask_bm';
    mask_bm = mask_bm(iy, ix);

    Fmask = griddedInterpolant({y_bm, x_bm}, mask_bm, 'nearest', 'nearest');
    node_type = round(Fmask(mesh_y(:), mesh_x(:)));

    is_ocean        = (node_type == 0);
    is_land         = (node_type == 1);
    is_grounded_raw = (node_type == 2);
    is_shelf_raw    = (node_type == 3);
    is_vostok_raw   = (node_type == 4);

    is_ice_raw   = is_grounded_raw | is_shelf_raw | is_vostok_raw;
    is_noice_raw = is_ocean | is_land;

    if size(elements, 2) > 3
        elements = elements(:, 1:3);
    end
    is_ice = is_ice_raw;
    nonice_elem = any(is_noice_raw(elements), 2);
    is_ice(elements(nonice_elem, :)) = false;

    is_grounded = is_grounded_raw & is_ice;
    is_shelf    = is_shelf_raw    & is_ice;
    is_vostok   = is_vostok_raw   & is_ice;
    is_noice    = ~is_ice;
end


function md = local_apply_cost_masks(md, is_shelf, is_noice, w101, w103, w501, ...
        friction_min, friction_max, zero_noice_data_terms)
    nv = md.mesh.numberofvertices;

    md.inversion.cost_functions_coefficients = zeros(nv, 3);
    md.inversion.cost_functions_coefficients(:, 1) = w101;
    md.inversion.cost_functions_coefficients(:, 2) = w103;
    md.inversion.cost_functions_coefficients(:, 3) = w501;

    md.inversion.cost_functions_coefficients(is_shelf, :) = 0;
    if zero_noice_data_terms
        md.inversion.cost_functions_coefficients(is_noice, 1:2) = 0;
    end

    md.friction.coefficient(is_shelf) = 0;
    reactivate = is_noice & (~isfinite(md.friction.coefficient) | md.friction.coefficient <= 0);
    md.friction.coefficient(reactivate) = max(friction_min, 20);

    minp = friction_min * ones(nv, 1);
    maxp = friction_max * ones(nv, 1);
    minp(is_shelf) = 0;
    maxp(is_shelf) = 0;
    md.inversion.min_parameters = minp;
    md.inversion.max_parameters = maxp;
end


function tf = local_has_member(s, name)
    tf = (isstruct(s) && isfield(s, name)) || (isobject(s) && isprop(s, name));
end


function [iter_final, fx_final, grad_norm, J101_log, J103_log, J501_log] = ...
        local_parse_final_m1qn3_line(cmdout)
    iter_final = NaN; fx_final = NaN; grad_norm = NaN;
    J101_log = NaN; J103_log = NaN; J501_log = NaN;
    try
        lines = strsplit(cmdout, {'\n', '\r'});
        fx_lines = lines(contains(lines, 'f(x)='));
        if isempty(fx_lines)
            return;
        end
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
    s = strrep(raw, '.', 'p');
    s = strrep(s, '+', '');
end


function T = local_empty_table(n, col_names, target_ratio)
    nc = numel(col_names);
    C = cell(1, nc);
    for i = 1:nc
        nm = col_names{i};
        if ismember(nm, {'status', 'error_message', 'model_file', 'log_file'})
            C{i} = repmat({''}, n, 1);
        elseif strcmp(nm, 'target_ratio101')
            C{i} = repmat(target_ratio(1), n, 1);
        elseif strcmp(nm, 'target_ratio103')
            C{i} = repmat(target_ratio(2), n, 1);
        elseif strcmp(nm, 'target_ratio501')
            C{i} = repmat(target_ratio(3), n, 1);
        else
            C{i} = NaN(n, 1);
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
            case 'run_id'
                C{i} = run_id;
            case 'w101'
                C{i} = w101;
            case 'w103'
                C{i} = w103;
            case 'w501'
                C{i} = w501;
            case 'status'
                C{i} = {'pending'};
            case 'error_message'
                C{i} = {''};
            case 'model_file'
                C{i} = {model_fname};
            case 'log_file'
                C{i} = {log_fname};
            case 'target_ratio101'
                C{i} = target_ratio(1);
            case 'target_ratio103'
                C{i} = target_ratio(2);
            case 'target_ratio501'
                C{i} = target_ratio(3);
            otherwise
                C{i} = NaN;
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
    for i = 2:(n - 1)
        x1 = logPhiReg(i - 1);  y1 = logPhiData(i - 1);
        x2 = logPhiReg(i);      y2 = logPhiData(i);
        x3 = logPhiReg(i + 1);  y3 = logPhiData(i + 1);

        area2 = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
        d12 = hypot(x2 - x1, y2 - y1);
        d23 = hypot(x3 - x2, y3 - y2);
        d13 = hypot(x3 - x1, y3 - y1);

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
    logmsg('Post-processing w101=%g, w103=%g', w101, w103);

    ok = strcmp(summary.status, 'success');
    if ~any(ok)
        logmsg('All runs failed for w101=%g. No candidates.', w101);
        return;
    end

    S = summary(ok, :);
    [S_sorted, sort_idx] = sortrows(S, 'w501');

    [~, idxA] = min(S.ratio_distance);
    logmsg('Candidate A closest ratio: run_id=%d w501=%.4e ratio_distance=%.4f', ...
        S.run_id(idxA), S.w501(idxA), S.ratio_distance(idxA));
    logmsg('  J101=%.2e J103=%.2e J501=%.2e', ...
        S.J101(idxA), S.J103(idxA), S.J501(idxA));

    idxB_local = local_find_lcurve_corner(S_sorted.logPhiData, S_sorted.logPhiReg);
    if ~isnan(idxB_local)
        idxB = sort_idx(idxB_local);
        logmsg('Candidate B weighted corner: run_id=%d w501=%.4e logPhiData=%.3f logPhiReg=%.3f', ...
            S.run_id(idxB), S.w501(idxB), S.logPhiData(idxB), S.logPhiReg(idxB));
    else
        logmsg('Candidate B weighted corner: not available.');
    end

    idxC_local = local_find_lcurve_corner(S_sorted.logPhiData_unweighted, S_sorted.logPhiReg_unweighted);
    if ~isnan(idxC_local)
        idxC = sort_idx(idxC_local);
        logmsg('Candidate C unweighted corner: run_id=%d w501=%.4e logPhiData_uw=%.3f logPhiReg_uw=%.3f', ...
            S.run_id(idxC), S.w501(idxC), ...
            S.logPhiData_unweighted(idxC), S.logPhiReg_unweighted(idxC));
    else
        logmsg('Candidate C unweighted corner: not available.');
    end

    fig1 = figure('Visible', 'off', 'Name', 'L-curve weighted', ...
        'Color', 'w', 'Position', [100 100 800 600]);
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

    fig2 = figure('Visible', 'off', 'Name', 'L-curve unweighted', ...
        'Color', 'w', 'Position', [950 100 800 600]);
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

    logmsg('Post-processing done for w101=%g. Summary: %s', ...
        w101, fullfile(out_dir, 'Lcurve_summary.csv'));
end


function s = local_num_to_str(x)
    if abs(x - round(x)) < 1e-12
        s = sprintf('%d', round(x));
    else
        s = strrep(sprintf('%.15g', x), '.', 'p');
    end
end


function s = local_num_to_expstr(x)
    s = sprintf('%.0e', x);
end
