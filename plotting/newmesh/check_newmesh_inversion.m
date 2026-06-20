%CHECK_NEWMESH_INVERSION Validate and plot Step 3 outputs for runme_newmesh.m.
%
% Run from the project root after Step 3:
%   run('plotting/newmesh/check_newmesh_inversion.m')

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

model_file = fullfile(paths.models_newmesh, 'RecoveryNewMesh_Inversion.mat');
parameterize_file = fullfile(paths.models_newmesh, 'RecoveryNewMesh_Parameterize.mat');
outdir = fullfile(paths.figures_newmesh, 'RecoveryNewMesh_inversion_check');
figdir = fullfile(outdir, 'fig');
pngdir = fullfile(outdir, 'png');

if ~exist(outdir, 'dir'), mkdir(outdir); end
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(pngdir, 'dir'), mkdir(pngdir); end

assert(exist(model_file, 'file') == 2, 'Model file not found: %s', model_file);

report_file = fullfile(outdir, 'newmesh_inversion_check_report.txt');
if exist(report_file, 'file'), delete(report_file); end
diary(report_file);

fprintf('===== RecoveryNewMesh Inversion Check =====\n');
fprintf('Model file: %s\n', model_file);
fprintf('Output dir: %s\n\n', outdir);

S = load(model_file, 'md');
md = S.md;
clear S;

has_parameterize = exist(parameterize_file, 'file') == 2;
if has_parameterize
    P = load(parameterize_file, 'md');
    md0 = P.md;
    clear P;
    fprintf('Parameterize reference: %s\n\n', parameterize_file);
else
    md0 = [];
    fprintf('Parameterize reference not found; before/after plots will be skipped.\n\n');
end

x = md.mesh.x(:);
y = md.mesh.y(:);
elements = md.mesh.elements;
if size(elements, 2) > 3
    faces = elements(:, 1:3);
else
    faces = elements;
end
nv = md.mesh.numberofvertices;
ne = md.mesh.numberofelements;

fprintf('Mesh: vertices=%d, elements=%d\n\n', nv, ne);

sol = [];
has_solution = local_has_member(md, 'results') && local_has_member(md.results, 'StressbalanceSolution') && ...
    ~isempty(md.results.StressbalanceSolution);
if has_solution
    sol = md.results.StressbalanceSolution;
end

vx_obs = local_field_or_empty(md, {'inversion', 'vx_obs'});
vy_obs = local_field_or_empty(md, {'inversion', 'vy_obs'});
vel_obs = local_field_or_empty(md, {'inversion', 'vel_obs'});
if isempty(vel_obs) && ~isempty(vx_obs) && ~isempty(vy_obs)
    vel_obs = sqrt(vx_obs(:).^2 + vy_obs(:).^2);
end

vx_model = local_solution_or_initialization(md, sol, 'Vx', 'vx');
vy_model = local_solution_or_initialization(md, sol, 'Vy', 'vy');
vel_model = local_solution_or_initialization(md, sol, 'Vel', 'vel');
if isempty(vel_model) && ~isempty(vx_model) && ~isempty(vy_model)
    vel_model = sqrt(vx_model(:).^2 + vy_model(:).^2);
end

friction = local_field_or_empty(md, {'friction', 'coefficient'});
friction_result = local_solution_field_or_empty(sol, 'FrictionCoefficient');
if isempty(friction) && ~isempty(friction_result)
    friction = friction_result;
end

assert(numel(vx_obs) == nv, 'inversion.vx_obs missing or wrong length.');
assert(numel(vy_obs) == nv, 'inversion.vy_obs missing or wrong length.');
assert(numel(vel_obs) == nv, 'inversion.vel_obs missing or wrong length.');
assert(numel(vx_model) == nv, 'Model Vx after inversion missing or wrong length.');
assert(numel(vy_model) == nv, 'Model Vy after inversion missing or wrong length.');
assert(numel(vel_model) == nv, 'Model speed after inversion missing or wrong length.');
assert(numel(friction) == nv, 'Friction coefficient after inversion missing or wrong length.');

if has_parameterize
    vx_before = local_field_or_empty(md0, {'initialization', 'vx'});
    vy_before = local_field_or_empty(md0, {'initialization', 'vy'});
    vel_before = local_field_or_empty(md0, {'initialization', 'vel'});
    friction_before = local_field_or_empty(md0, {'friction', 'coefficient'});
else
    vx_before = [];
    vy_before = [];
    vel_before = [];
    friction_before = [];
end

speed_residual = vel_model(:) - vel_obs(:);
abs_speed_residual = abs(speed_residual);
component_residual = hypot(vx_model(:) - vx_obs(:), vy_model(:) - vy_obs(:));
relative_speed_residual = speed_residual ./ max(abs(vel_obs(:)), 1.0);
speed_change_from_parameterize = local_same_length_difference(vel_model, vel_before, nv);
friction_change_from_parameterize = local_same_length_difference(friction, friction_before, nv);
friction_ratio_from_parameterize = local_safe_ratio(friction, friction_before, nv);

ice = md.mask.ice_levelset(:) < 0;
ocean = md.mask.ocean_levelset(:) < 0;
grounded = ice & ~ocean;
shelf = ice & ocean;
nonice = ~ice;
bc_mask = local_boundary_velocity_mask(md, nv);

checks = {};
checks = add_check(checks, 'inversion.vx_obs', vx_obs, nv);
checks = add_check(checks, 'inversion.vy_obs', vy_obs, nv);
checks = add_check(checks, 'inversion.vel_obs', vel_obs, nv);
checks = add_check(checks, 'model.vx_after_inversion', vx_model, nv);
checks = add_check(checks, 'model.vy_after_inversion', vy_model, nv);
checks = add_check(checks, 'model.vel_after_inversion', vel_model, nv);
checks = add_check(checks, 'friction.coefficient', friction, nv);
checks = add_check(checks, 'results.FrictionCoefficient', friction_result, nv);
checks = add_check(checks, 'speed_residual', speed_residual, nv);
checks = add_check(checks, 'relative_speed_residual', relative_speed_residual, nv);
checks = add_check(checks, 'component_residual', component_residual, nv);
checks = add_check(checks, 'mask.ice_levelset', md.mask.ice_levelset, nv);
checks = add_check(checks, 'mask.ocean_levelset', md.mask.ocean_levelset, nv);
if has_parameterize
    checks = add_check(checks, 'parameterize.initialization.vel', vel_before, nv);
    checks = add_check(checks, 'parameterize.friction.coefficient', friction_before, nv);
end

field_table = cell2table(checks, 'VariableNames', ...
    {'field', 'expected_length', 'actual_length', 'finite_count', 'nan_count', 'inf_count', 'min_value', 'max_value', 'mean_value', 'status'});
writetable(field_table, fullfile(outdir, 'field_summary.csv'));
disp(field_table);

metrics = table;
metrics.metric = string({
    'has_stressbalance_solution'
    'n_vertices'
    'n_elements'
    'n_grounded_nodes'
    'n_shelf_nodes'
    'n_nonice_nodes'
    'n_boundary_velocity_nodes'
    'speed_residual_mean_m_per_yr'
    'speed_residual_median_m_per_yr'
    'speed_residual_rms_m_per_yr'
    'speed_residual_max_abs_m_per_yr'
    'component_residual_rms_m_per_yr'
    'relative_speed_residual_mean_abs'
    'relative_speed_residual_max_abs'
    'friction_min'
    'friction_max'
    'friction_mean'
    'n_zero_friction_nodes'
    'n_negative_friction_nodes'
    'max_abs_result_minus_top_friction'
});
metrics.value = [
    double(has_solution)
    nv
    ne
    sum(grounded)
    sum(shelf)
    sum(nonice)
    sum(bc_mask)
    local_mean_finite(speed_residual)
    local_median_finite(speed_residual)
    sqrt(local_mean_finite(speed_residual.^2))
    local_max_finite(abs_speed_residual)
    sqrt(local_mean_finite(component_residual.^2))
    local_mean_finite(abs(relative_speed_residual))
    local_max_finite(abs(relative_speed_residual))
    local_min_finite(friction)
    local_max_finite(friction)
    local_mean_finite(friction)
    sum(friction == 0)
    sum(friction < 0)
    local_max_abs_difference(friction_result, friction, nv)
];
writetable(metrics, fullfile(outdir, 'inversion_metrics.csv'));

fprintf('\nInversion metrics:\n');
disp(metrics);

local_write_result_fields(md, fullfile(outdir, 'stressbalance_solution_fields.txt'));

make_mesh_edge_plot(x, y, faces, fullfile(figdir, 'mesh_edges.fig'), fullfile(pngdir, 'mesh_edges.png'));
make_node_plot(x, y, faces, vel_obs, 'velocity_obs_m_per_yr', 'Observed speed (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, vel_model, 'velocity_model_after_inversion_m_per_yr', 'Model speed after inversion (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, speed_residual, 'speed_residual_model_minus_obs_m_per_yr', 'Speed residual: model - observed (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, abs_speed_residual, 'speed_residual_abs_m_per_yr', 'Absolute speed residual (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, relative_speed_residual, 'speed_residual_relative', 'Relative speed residual', figdir, pngdir);
make_node_plot(x, y, faces, component_residual, 'velocity_component_residual_m_per_yr', 'Vector velocity residual magnitude (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, vx_model(:) - vx_obs(:), 'vx_residual_model_minus_obs_m_per_yr', 'Vx residual: model - observed (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, vy_model(:) - vy_obs(:), 'vy_residual_model_minus_obs_m_per_yr', 'Vy residual: model - observed (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, friction, 'friction_coefficient_after_inversion', 'Friction coefficient after inversion', figdir, pngdir);
make_node_plot(x, y, faces, log10(max(friction, 1e-2)), 'log10_friction_coefficient_after_inversion', 'log10 friction coefficient after inversion', figdir, pngdir);
make_node_plot(x, y, faces, double(grounded), 'grounded_mask', 'Grounded ice mask (1=grounded)', figdir, pngdir);
make_node_plot(x, y, faces, double(shelf), 'shelf_mask', 'Floating ice mask (1=shelf)', figdir, pngdir);
make_node_plot(x, y, faces, double(nonice), 'nonice_mask', 'Non-ice mask (1=nonice)', figdir, pngdir);
make_node_plot(x, y, faces, double(bc_mask), 'stressbalance_boundary_velocity_mask', 'Velocity boundary condition mask', figdir, pngdir);

if has_parameterize
    make_node_plot(x, y, faces, vel_before, 'velocity_before_inversion_m_per_yr', 'Speed before inversion (m/yr)', figdir, pngdir);
    make_node_plot(x, y, faces, speed_change_from_parameterize, 'speed_change_after_minus_before_m_per_yr', 'Speed change: after - before inversion (m/yr)', figdir, pngdir);
    make_node_plot(x, y, faces, friction_before, 'friction_coefficient_before_inversion', 'Friction coefficient before inversion', figdir, pngdir);
    make_node_plot(x, y, faces, friction_change_from_parameterize, 'friction_change_after_minus_before', 'Friction change: after - before inversion', figdir, pngdir);
    make_node_plot(x, y, faces, friction_ratio_from_parameterize, 'friction_ratio_after_div_before', 'Friction ratio: after / before inversion', figdir, pngdir);
end

plot_cost_coefficients(md, x, y, faces, figdir, pngdir);
plot_parameter_bounds(md, x, y, faces, figdir, pngdir);
make_histogram_plot(abs_speed_residual, 'hist_abs_speed_residual_m_per_yr', 'Absolute speed residual histogram', '|model - observed| (m/yr)', figdir, pngdir);
make_histogram_plot(friction, 'hist_friction_coefficient', 'Friction coefficient histogram', 'Friction coefficient', figdir, pngdir);
make_vector_plot(x, y, vx_obs, vy_obs, 'velocity_vectors_observed', 'Observed velocity vectors', figdir, pngdir);
make_vector_plot(x, y, vx_model, vy_model, 'velocity_vectors_after_inversion', 'Model velocity vectors after inversion', figdir, pngdir);

fprintf('\nSaved report and figures to: %s\n', outdir);
diary off;


function checks = add_check(checks, name, values, expected_len, allow_nan)
    if nargin < 5
        allow_nan = false;
    end
    values = values(:);
    finite_values = values(isfinite(values));
    actual_len = numel(values);
    n_nan = sum(isnan(values));
    n_inf = sum(isinf(values));
    if isempty(finite_values)
        minv = NaN;
        maxv = NaN;
        meanv = NaN;
    else
        minv = min(finite_values);
        maxv = max(finite_values);
        meanv = mean(finite_values);
    end
    len_ok = actual_len == expected_len;
    finite_ok = (n_nan == 0 && n_inf == 0) || (allow_nan && n_inf == 0);
    if len_ok && finite_ok
        status = "ok";
    else
        status = "check";
    end
    checks(end + 1, :) = {string(name), expected_len, actual_len, numel(finite_values), n_nan, n_inf, minv, maxv, meanv, status};
end


function data = local_field_or_empty(s, names)
    data = [];
    cur = s;
    for i = 1:numel(names)
        if isstruct(cur) && isfield(cur, names{i})
            cur = cur.(names{i});
        elseif isobject(cur) && isprop(cur, names{i})
            cur = cur.(names{i});
        else
            return;
        end
    end
    if isnumeric(cur)
        data = double(cur(:));
    end
end


function tf = local_has_member(s, name)
    tf = (isstruct(s) && isfield(s, name)) || (isobject(s) && isprop(s, name));
end


function data = local_solution_field_or_empty(sol, name)
    data = [];
    if isempty(sol)
        return;
    end
    if isstruct(sol) && isfield(sol, name)
        data = double(sol.(name)(:));
    elseif isobject(sol) && isprop(sol, name)
        data = double(sol.(name)(:));
    end
end


function data = local_solution_or_initialization(md, sol, sol_name, init_name)
    data = local_solution_field_or_empty(sol, sol_name);
    if isempty(data)
        data = local_field_or_empty(md, {'initialization', init_name});
    end
end


function diffv = local_same_length_difference(after, before, nv)
    if numel(after) == nv && numel(before) == nv
        diffv = after(:) - before(:);
    else
        diffv = NaN(nv, 1);
    end
end


function ratio = local_safe_ratio(after, before, nv)
    ratio = NaN(nv, 1);
    if numel(after) ~= nv || numel(before) ~= nv
        return;
    end
    ok = isfinite(after(:)) & isfinite(before(:)) & abs(before(:)) > 0;
    ratio(ok) = after(ok) ./ before(ok);
end


function val = local_max_abs_difference(a, b, nv)
    if numel(a) == nv && numel(b) == nv
        val = max(abs(a(:) - b(:)));
    else
        val = NaN;
    end
end


function val = local_mean_finite(values)
    values = values(:);
    values = values(isfinite(values));
    if isempty(values)
        val = NaN;
    else
        val = mean(values);
    end
end


function val = local_median_finite(values)
    values = values(:);
    values = values(isfinite(values));
    if isempty(values)
        val = NaN;
    else
        val = median(values);
    end
end


function val = local_min_finite(values)
    values = values(:);
    values = values(isfinite(values));
    if isempty(values)
        val = NaN;
    else
        val = min(values);
    end
end


function val = local_max_finite(values)
    values = values(:);
    values = values(isfinite(values));
    if isempty(values)
        val = NaN;
    else
        val = max(values);
    end
end


function bc_mask = local_boundary_velocity_mask(md, nv)
    bc_mask = false(nv, 1);
    if ~local_has_member(md, 'stressbalance')
        return;
    end
    names = {'spcvx', 'spcvy', 'spcvz'};
    for i = 1:numel(names)
        vals = local_field_or_empty(md, {'stressbalance', names{i}});
        if numel(vals) == nv
            bc_mask = bc_mask | ~isnan(vals(:));
        end
    end
end


function local_write_result_fields(md, out_file)
    fid = fopen(out_file, 'w');
    if fid < 0
        warning('Could not write result field list: %s', out_file);
        return;
    end
    cleaner = onCleanup(@() fclose(fid));
    if ~local_has_member(md, 'results') || ~local_has_member(md.results, 'StressbalanceSolution')
        fprintf(fid, 'No md.results.StressbalanceSolution found.\n');
        return;
    end
    sol = md.results.StressbalanceSolution;
    if isstruct(sol)
        f = fieldnames(sol);
    else
        f = properties(sol);
    end
    for i = 1:numel(f)
        try
            v = sol.(f{i});
            fprintf(fid, '%s: %s %s\n', f{i}, class(v), mat2str(size(v)));
        catch
            fprintf(fid, '%s\n', f{i});
        end
    end
end


function make_node_plot(x, y, faces, data, basename, title_text, figdir, pngdir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 980 780]);
    patch('Faces', faces, ...
          'Vertices', [x, y], ...
          'FaceVertexCData', data(:), ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    axis equal tight;
    box on;
    colormap(gca, parula);
    colorbar;
    title(title_text, 'Interpreter', 'none');
    xlabel('X (m)');
    ylabel('Y (m)');
    set(fig, 'Visible', 'on');
    savefig(fig, fullfile(figdir, [basename '.fig']));
    print(fig, fullfile(pngdir, basename), '-dpng', '-r300');
    close(fig);
end


function make_mesh_edge_plot(x, y, faces, fig_file, png_file)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 980 780]);
    patch('Faces', faces, ...
          'Vertices', [x, y], ...
          'FaceColor', 'none', ...
          'EdgeColor', [0.15 0.15 0.15], ...
          'LineWidth', 0.15);
    axis equal tight;
    box on;
    title('Mesh edges', 'Interpreter', 'none');
    xlabel('X (m)');
    ylabel('Y (m)');
    set(fig, 'Visible', 'on');
    savefig(fig, fig_file);
    print(fig, png_file, '-dpng', '-r300');
    close(fig);
end


function make_histogram_plot(values, basename, title_text, xlabel_text, figdir, pngdir)
    values = values(:);
    values = values(isfinite(values));
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 620]);
    if isempty(values)
        text(0.5, 0.5, 'No finite values', 'HorizontalAlignment', 'center');
        axis off;
    else
        histogram(values, 60);
        box on;
        grid on;
        xlabel(xlabel_text, 'Interpreter', 'none');
        ylabel('Count');
    end
    title(title_text, 'Interpreter', 'none');
    set(fig, 'Visible', 'on');
    savefig(fig, fullfile(figdir, [basename '.fig']));
    print(fig, fullfile(pngdir, basename), '-dpng', '-r300');
    close(fig);
end


function make_vector_plot(x, y, vx, vy, basename, title_text, figdir, pngdir)
    vx = vx(:);
    vy = vy(:);
    ok = isfinite(x) & isfinite(y) & isfinite(vx) & isfinite(vy);
    ids = find(ok);
    max_arrows = 3000;
    if numel(ids) > max_arrows
        ids = ids(round(linspace(1, numel(ids), max_arrows)));
    end
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 980 780]);
    quiver(x(ids), y(ids), vx(ids), vy(ids), 1.5, 'k');
    axis equal tight;
    box on;
    title(title_text, 'Interpreter', 'none');
    xlabel('X (m)');
    ylabel('Y (m)');
    set(fig, 'Visible', 'on');
    savefig(fig, fullfile(figdir, [basename '.fig']));
    print(fig, fullfile(pngdir, basename), '-dpng', '-r300');
    close(fig);
end


function plot_cost_coefficients(md, x, y, faces, figdir, pngdir)
    coeff = local_field_or_empty(md, {'inversion', 'cost_functions_coefficients'});
    if isempty(coeff)
        return;
    end
    nv = numel(x);
    ncol = numel(coeff) / nv;
    if abs(ncol - round(ncol)) > 0
        return;
    end
    coeff = reshape(coeff, nv, round(ncol));
    names = {'cost_coeff_101_abs_velocity', 'cost_coeff_103_log_velocity', 'cost_coeff_501_regularization'};
    titles = {'Cost coefficient 101', 'Cost coefficient 103', 'Cost coefficient 501'};
    for i = 1:min(size(coeff, 2), numel(names))
        make_node_plot(x, y, faces, coeff(:, i), names{i}, titles{i}, figdir, pngdir);
    end
end


function plot_parameter_bounds(md, x, y, faces, figdir, pngdir)
    minp = local_field_or_empty(md, {'inversion', 'min_parameters'});
    maxp = local_field_or_empty(md, {'inversion', 'max_parameters'});
    nv = numel(x);
    if numel(minp) == nv
        make_node_plot(x, y, faces, minp, 'inversion_min_friction_parameter', 'Inversion minimum friction parameter', figdir, pngdir);
    end
    if numel(maxp) == nv
        make_node_plot(x, y, faces, maxp, 'inversion_max_friction_parameter', 'Inversion maximum friction parameter', figdir, pngdir);
    end
end
