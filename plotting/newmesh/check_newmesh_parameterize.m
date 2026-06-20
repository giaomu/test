%CHECK_NEWMESH_PARAMETERIZE Validate Step 1/2 outputs for runme_newmesh.m.
%
% Run from the project root:
%   run('plotting/newmesh/check_newmesh_parameterize.m')

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

model_file = fullfile(paths.models_newmesh, 'RecoveryNewMesh_Parameterize.mat');
outdir = fullfile(paths.figures_newmesh, 'RecoveryNewMesh_parameterize_check');
figdir = fullfile(outdir, 'fig');
pngdir = fullfile(outdir, 'png');

if ~exist(outdir, 'dir'), mkdir(outdir); end
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(pngdir, 'dir'), mkdir(pngdir); end

assert(exist(model_file, 'file') == 2, 'Model file not found: %s', model_file);

report_file = fullfile(outdir, 'newmesh_parameterize_check_report.txt');
if exist(report_file, 'file'), delete(report_file); end
diary(report_file);

fprintf('===== RecoveryNewMesh Parameterize Check =====\n');
fprintf('Model file: %s\n', model_file);
fprintf('Output dir: %s\n\n', outdir);

S = load(model_file, 'md');
md = S.md;
clear S;

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

checks = {};
checks = add_check(checks, 'geometry.bed', md.geometry.bed, nv);
checks = add_check(checks, 'geometry.base', md.geometry.base, nv);
checks = add_check(checks, 'geometry.surface', md.geometry.surface, nv);
checks = add_check(checks, 'geometry.thickness', md.geometry.thickness, nv);
checks = add_check(checks, 'initialization.vx', md.initialization.vx, nv);
checks = add_check(checks, 'initialization.vy', md.initialization.vy, nv);
checks = add_check(checks, 'initialization.vel', md.initialization.vel, nv);
checks = add_check(checks, 'inversion.vx_obs', md.inversion.vx_obs, nv);
checks = add_check(checks, 'inversion.vy_obs', md.inversion.vy_obs, nv);
checks = add_check(checks, 'inversion.vel_obs', md.inversion.vel_obs, nv);
checks = add_check(checks, 'mask.ice_levelset', md.mask.ice_levelset, nv);
checks = add_check(checks, 'mask.ocean_levelset', md.mask.ocean_levelset, nv);
checks = add_check(checks, 'friction.coefficient', md.friction.coefficient, nv);
checks = add_check(checks, 'basalforcings.geothermalflux', md.basalforcings.geothermalflux, nv);
checks = add_check(checks, 'basalforcings.groundedice_melting_rate', md.basalforcings.groundedice_melting_rate, nv);
checks = add_check(checks, 'materials.rheology_B', md.materials.rheology_B, nv);
checks = add_check(checks, 'materials.rheology_n', md.materials.rheology_n, ne);
checks = add_check(checks, 'stressbalance.spcvx', md.stressbalance.spcvx, nv, true);
checks = add_check(checks, 'stressbalance.spcvy', md.stressbalance.spcvy, nv, true);
checks = add_check(checks, 'stressbalance.spcvz', md.stressbalance.spcvz, nv, true);

check_table = cell2table(checks, 'VariableNames', ...
    {'field', 'expected_length', 'actual_length', 'finite_count', 'nan_count', 'inf_count', 'min_value', 'max_value', 'mean_value', 'status'});
writetable(check_table, fullfile(outdir, 'field_summary.csv'));
disp(check_table);

bed = md.geometry.bed(:);
base = md.geometry.base(:);
surface = md.geometry.surface(:);
thickness = md.geometry.thickness(:);
vx = md.initialization.vx(:);
vy = md.initialization.vy(:);
vel = md.initialization.vel(:);
vx_obs = md.inversion.vx_obs(:);
vy_obs = md.inversion.vy_obs(:);
vel_obs = md.inversion.vel_obs(:);
friction = md.friction.coefficient(:);
geothermalflux = md.basalforcings.geothermalflux(:);
grounded_melt = md.basalforcings.groundedice_melting_rate(:);
ice = md.mask.ice_levelset(:) < 0;
ocean = md.mask.ocean_levelset(:) < 0;
grounded = ice & ~ocean;
shelf = ice & ocean;
nonice = ~ice;

derived = table;
derived.metric = string({
    'max_abs_base_minus_bed'
    'max_abs_surface_minus_base_minus_thickness'
    'min_thickness'
    'max_abs_init_vel_minus_sqrt_vx_vy'
    'max_abs_obs_vel_minus_sqrt_vx_vy'
    'max_abs_obs_minus_init_vx'
    'max_abs_obs_minus_init_vy'
    'n_ice_nodes'
    'n_grounded_nodes'
    'n_shelf_nodes'
    'n_nonice_nodes'
    'n_boundary_velocity_nodes'
    'n_zero_friction_nodes'
    'n_negative_friction_nodes'
    'min_geothermalflux_W_m2'
    'max_geothermalflux_W_m2'
    'mean_geothermalflux_W_m2'
    'max_abs_grounded_melt_rate'
});
derived.value = [
    max(abs(base - bed))
    max(abs(surface - base - thickness))
    min(thickness)
    max(abs(vel - sqrt(vx.^2 + vy.^2)))
    max(abs(vel_obs - sqrt(vx_obs.^2 + vy_obs.^2)))
    max(abs(vx_obs - vx))
    max(abs(vy_obs - vy))
    sum(ice)
    sum(grounded)
    sum(shelf)
    sum(nonice)
    sum(~isnan(md.stressbalance.spcvx(:)) | ~isnan(md.stressbalance.spcvy(:)))
    sum(friction == 0)
    sum(friction < 0)
    min(geothermalflux)
    max(geothermalflux)
    mean(geothermalflux)
    max(abs(grounded_melt))
];
writetable(derived, fullfile(outdir, 'derived_consistency.csv'));

fprintf('\nDerived consistency metrics:\n');
disp(derived);

make_mesh_edge_plot(x, y, faces, fullfile(figdir, 'mesh_edges.fig'), fullfile(pngdir, 'mesh_edges.png'));
element_area = local_triangle_area(x, y, faces);
node_area = local_element_to_node_mean(faces, element_area, nv);
element_max_edge = local_triangle_max_edge(x, y, faces);
node_max_edge = local_element_to_node_max(faces, element_max_edge, nv);
make_node_plot(x, y, faces, node_area, 'node_mean_element_area_m2', 'Mean connected element area (m^2)', figdir, pngdir);
make_node_plot(x, y, faces, node_max_edge, 'node_max_edge_length_m', 'Maximum connected edge length (m)', figdir, pngdir);

make_node_plot(x, y, faces, thickness, 'thickness_m', 'Ice thickness (m)', figdir, pngdir);
make_node_plot(x, y, faces, bed, 'bed_m', 'Bed elevation (m)', figdir, pngdir);
make_node_plot(x, y, faces, surface, 'surface_m', 'Surface elevation (m)', figdir, pngdir);
make_node_plot(x, y, faces, vel_obs, 'velocity_obs_m_per_yr', 'Observed speed (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, vx_obs, 'vx_obs_m_per_yr', 'Observed Vx (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, vy_obs, 'vy_obs_m_per_yr', 'Observed Vy (m/yr)', figdir, pngdir);
make_node_plot(x, y, faces, friction, 'friction_coefficient', 'Friction coefficient', figdir, pngdir);
make_node_plot(x, y, faces, geothermalflux, 'geothermalflux_W_m2', 'Geothermal flux (W/m^2)', figdir, pngdir);
make_node_plot(x, y, faces, double(ice), 'ice_mask', 'Ice mask (1=ice)', figdir, pngdir);
make_node_plot(x, y, faces, double(ocean), 'ocean_mask', 'Ocean/floating mask (1=ocean/floating)', figdir, pngdir);

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


function area = local_triangle_area(x, y, faces)
    x1 = x(faces(:,1)); y1 = y(faces(:,1));
    x2 = x(faces(:,2)); y2 = y(faces(:,2));
    x3 = x(faces(:,3)); y3 = y(faces(:,3));
    area = 0.5 * abs((x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1));
end


function max_edge = local_triangle_max_edge(x, y, faces)
    e12 = hypot(x(faces(:,1)) - x(faces(:,2)), y(faces(:,1)) - y(faces(:,2)));
    e23 = hypot(x(faces(:,2)) - x(faces(:,3)), y(faces(:,2)) - y(faces(:,3)));
    e31 = hypot(x(faces(:,3)) - x(faces(:,1)), y(faces(:,3)) - y(faces(:,1)));
    max_edge = max([e12, e23, e31], [], 2);
end


function node_mean = local_element_to_node_mean(faces, elem_values, nv)
    node_sum = accumarray(faces(:), repmat(elem_values(:), 3, 1), [nv, 1], @sum, 0);
    node_count = accumarray(faces(:), 1, [nv, 1], @sum, 0);
    node_mean = node_sum ./ max(node_count, 1);
end


function node_max = local_element_to_node_max(faces, elem_values, nv)
    node_max = accumarray(faces(:), repmat(elem_values(:), 3, 1), [nv, 1], @max, NaN);
end
