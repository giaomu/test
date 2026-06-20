%% plot_parameterize_fields.m
%  Batch plots for basic physical fields in Recovery_Parameterize.mat.
%  Uses native MATLAB plotting only (patch / scatter / histogram).
%  Saves each figure as both .fig and high-resolution .png.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

%% ===== User settings =====================================================
model_file      = fullfile(paths.models, 'Recovery_Parameterize.mat');
bedmachine_file = paths.bedmachine_file;
issm_root       = paths.issm_root;

dpi_png         = 400;
marker_size     = 4;
show_mesh_edges = false;
vel_log_floor   = 1e-1;
fric_log_floor  = 1e-2;
min_thickness   = 25;

issm_bin = fullfile(issm_root, 'bin');
if exist(issm_bin, 'dir') == 7
    addpath(issm_bin);
else
    warning('ISSM bin folder not found: %s', issm_bin);
end

%% ===== Output folders ====================================================
assert(exist(model_file, 'file') == 2, 'model_file does not exist: %s', model_file);

[~, model_name, ~] = fileparts(model_file);
outroot = fullfile(paths.figures, sprintf('%s_physics', model_name));
figdir  = fullfile(outroot, 'fig');
pngdir  = fullfile(outroot, 'png');

if ~exist(outroot,'dir'), mkdir(outroot); end
if ~exist(figdir, 'dir'), mkdir(figdir);  end
if ~exist(pngdir, 'dir'), mkdir(pngdir);  end

fprintf('Model file : %s\n', model_file);
fprintf('Output root: %s\n', outroot);

%% ===== Load model ========================================================
S  = load(model_file);
md = local_extract_md(S);
clear S;

x        = double(md.mesh.x(:));
y        = double(md.mesh.y(:));
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

edge_flag = show_mesh_edges;

%% ===== BedMachine mask on model nodes ===================================
assert(exist(bedmachine_file, 'file') == 2, ...
    'bedmachine_file does not exist: %s', bedmachine_file);

[node_type, ~] = local_read_bedmachine_mask(bedmachine_file, x, y);

is_ocean    = node_type == 0;
is_land     = node_type == 1;
is_grounded = node_type == 2;
is_shelf    = node_type == 3;
is_vostok   = node_type == 4;
is_noice    = is_ocean | is_land;
is_ice_bm   = is_grounded | is_shelf | is_vostok;

fprintf('BedMachine node counts:\n');
fprintf('  ocean    : %d\n', sum(is_ocean));
fprintf('  land     : %d\n', sum(is_land));
fprintf('  grounded : %d\n', sum(is_grounded));
fprintf('  shelf    : %d\n', sum(is_shelf));
fprintf('  vostok   : %d\n', sum(is_vostok));

%% ===== Frequently used model fields =====================================
ice_levelset   = local_clip(local_try_field(md, {'mask','ice_levelset'}), nv);
ocean_levelset = local_clip(local_try_field(md, {'mask','ocean_levelset'}), nv);

bed       = local_clip(local_try_field(md, {'geometry','bed'}), nv);
base      = local_clip(local_try_field(md, {'geometry','base'}), nv);
surface   = local_clip(local_try_field(md, {'geometry','surface'}), nv);
thickness = local_clip(local_try_field(md, {'geometry','thickness'}), nv);

vx_init  = local_clip(local_try_field(md, {'initialization','vx'}), nv);
vy_init  = local_clip(local_try_field(md, {'initialization','vy'}), nv);
vel_init = local_clip(local_try_field(md, {'initialization','vel'}), nv);
if isempty(vel_init) && ~isempty(vx_init) && ~isempty(vy_init)
    vel_init = sqrt(vx_init.^2 + vy_init.^2);
end

vx_obs  = local_clip(local_try_field(md, {'inversion','vx_obs'}), nv);
vy_obs  = local_clip(local_try_field(md, {'inversion','vy_obs'}), nv);
vel_obs = local_clip(local_try_field(md, {'inversion','vel_obs'}), nv);
if isempty(vel_obs) && ~isempty(vx_obs) && ~isempty(vy_obs)
    vel_obs = sqrt(vx_obs.^2 + vy_obs.^2);
end

fric_coef = local_clip(local_try_field(md, {'friction','coefficient'}), nv);
rheo_B    = local_clip(local_try_field(md, {'materials','rheology_B'}), nv);
rheo_n    = local_try_field(md, {'materials','rheology_n'});
if ~isempty(rheo_n), rheo_n = double(rheo_n(:)); end

grounded_melt = local_clip(local_try_field(md, {'basalforcings','groundedice_melting_rate'}), nv);
floating_melt = local_clip(local_try_field(md, {'basalforcings','floatingice_melting_rate'}), nv);
watercolumn   = local_clip(local_try_field(md, {'initialization','watercolumn'}), nv);
waterfraction = local_clip(local_try_field(md, {'initialization','waterfraction'}), nv);

cost_coeff = local_try_field(md, {'inversion','cost_functions_coefficients'});

bc_mask = false(nv,1);
if local_has_member(md, 'stressbalance')
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

%% ===== Plotting ==========================================================

% ---- 01 BedMachine node types -------------------------------------------
fig = local_new_fig();
local_plot_categorical(fig, faces, x, y, node_type);
title('BedMachine Node Types');
local_save_both(fig, figdir, pngdir, '01_bedmachine_node_types', dpi_png);

% ---- 02 md.mask.ice_levelset --------------------------------------------
local_safe_node_plot('md.mask.ice_levelset', ice_levelset, ...
    faces, x, y, nv, edge_flag, 'md.mask.ice\_levelset', ...
    figdir, pngdir, '02_md_mask_ice_levelset', dpi_png, false, 0, false);

% ---- 03 md.mask.ocean_levelset ------------------------------------------
local_safe_node_plot('md.mask.ocean_levelset', ocean_levelset, ...
    faces, x, y, nv, edge_flag, 'md.mask.ocean\_levelset', ...
    figdir, pngdir, '03_md_mask_ocean_levelset', dpi_png, false, 0, false);

% ---- 04 BedMachine ice mask vs md.mask.ice_levelset ----------------------
if ~isempty(ice_levelset)
    md_ice = ice_levelset < 0;
    mismatch = md_ice ~= is_ice_bm;
    fig = local_new_fig();
    local_plot_background(faces, x, y, [0.85 0.85 0.85]);
    hold on;
    h1 = scatter(x(is_ice_bm), y(is_ice_bm), marker_size, [0.20 0.55 0.90], 'filled');
    h2 = scatter(x(mismatch), y(mismatch), marker_size*2, 'r', 'filled');
    hold off;
    axis equal tight; box on; xlabel('x (m)'); ylabel('y (m)');
    legend([h1 h2], {'BedMachine ice nodes','mismatch nodes'}, 'Location','best');
    title(sprintf('BedMachine Ice Mask vs md.mask.ice\\_levelset (Mismatch N = %d)', sum(mismatch)));
    local_save_both(fig, figdir, pngdir, '04_mask_compare_bedmachine_vs_md', dpi_png);
else
    warning('md.mask.ice_levelset not found, skipping 04.');
end

% ---- 05 bed elevation ----------------------------------------------------
local_safe_node_plot('md.geometry.bed', bed, ...
    faces, x, y, nv, edge_flag, 'Bed Elevation (m)', ...
    figdir, pngdir, '05_bed_elevation', dpi_png, false, 0, false);

% ---- 06 base elevation ---------------------------------------------------
local_safe_node_plot('md.geometry.base', base, ...
    faces, x, y, nv, edge_flag, 'Base Elevation (m)', ...
    figdir, pngdir, '06_base_elevation', dpi_png, false, 0, false);

% ---- 07 surface elevation ------------------------------------------------
local_safe_node_plot('md.geometry.surface', surface, ...
    faces, x, y, nv, edge_flag, 'Surface Elevation (m)', ...
    figdir, pngdir, '07_surface_elevation', dpi_png, false, 0, false);

% ---- 08 thickness --------------------------------------------------------
local_safe_node_plot('md.geometry.thickness', thickness, ...
    faces, x, y, nv, edge_flag, 'Ice Thickness (m)', ...
    figdir, pngdir, '08_thickness', dpi_png, false, 0, false);

% ---- 09 surface - bed - thickness check ---------------------------------
if ~isempty(surface) && ~isempty(bed) && ~isempty(thickness)
    geom_res = surface - bed - thickness;
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, geom_res, edge_flag);
    colormap(fig, local_diverging_cmap());
    local_balanced_caxis(geom_res);
    title('Geometry Check: surface - bed - thickness (m)');
    local_save_both(fig, figdir, pngdir, '09_surface_bed_thickness_check', dpi_png);
else
    warning('Cannot compute surface - bed - thickness check, skipping 09.');
end

% ---- 10 base - bed check -------------------------------------------------
if ~isempty(base) && ~isempty(bed)
    base_res = base - bed;
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, base_res, edge_flag);
    colormap(fig, local_diverging_cmap());
    local_balanced_caxis(base_res);
    title('Geometry Check: base - bed (m)');
    local_save_both(fig, figdir, pngdir, '10_base_bed_check', dpi_png);
else
    warning('Cannot compute base - bed check, skipping 10.');
end

% ---- 11 thin ice nodes ---------------------------------------------------
if ~isempty(thickness)
    thin_nodes = thickness < min_thickness;
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, thickness, edge_flag);
    hold on;
    h = scatter(x(thin_nodes), y(thin_nodes), marker_size*2, 'r', 'filled');
    hold off;
    if any(thin_nodes)
        legend(h, sprintf('thickness < %g m', min_thickness), 'Location','best');
    end
    title(sprintf('Thin Ice Check: thickness < %g m (N = %d)', min_thickness, sum(thin_nodes)));
    local_save_both(fig, figdir, pngdir, '11_thin_ice_nodes_lt_25m', dpi_png);
else
    warning('thickness not found, skipping 11.');
end

% ---- 12 initialization vx ------------------------------------------------
local_safe_node_plot('md.initialization.vx', vx_init, ...
    faces, x, y, nv, edge_flag, 'Initial Vx (m/yr)', ...
    figdir, pngdir, '12_initial_vx', dpi_png, false, 0, true);

% ---- 13 initialization vy ------------------------------------------------
local_safe_node_plot('md.initialization.vy', vy_init, ...
    faces, x, y, nv, edge_flag, 'Initial Vy (m/yr)', ...
    figdir, pngdir, '13_initial_vy', dpi_png, false, 0, true);

% ---- 14 initialization velocity magnitude, linear scale ------------------
local_safe_node_plot('initial velocity magnitude', vel_init, ...
    faces, x, y, nv, edge_flag, 'Initial Velocity Magnitude (m/yr)', ...
    figdir, pngdir, '14_initial_velocity_magnitude_linear', dpi_png, false, 0, false);

% ---- 15 initialization velocity magnitude, log scale ---------------------
local_safe_node_plot('initial velocity magnitude', vel_init, ...
    faces, x, y, nv, edge_flag, 'Initial Velocity Magnitude (m/yr)  [log_{10}]', ...
    figdir, pngdir, '15_initial_velocity_magnitude_log', dpi_png, true, vel_log_floor, false);

% ---- 16 observed velocity magnitude, log scale ---------------------------
local_safe_node_plot('md.inversion.vel_obs', vel_obs, ...
    faces, x, y, nv, edge_flag, 'Observed Velocity Magnitude (m/yr)  [log_{10}]', ...
    figdir, pngdir, '16_observed_velocity_magnitude_log', dpi_png, true, vel_log_floor, false);

% ---- 17 observed - initial velocity check --------------------------------
if ~isempty(vel_obs) && ~isempty(vel_init)
    vel_diff = vel_obs - vel_init;
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, vel_diff, edge_flag);
    colormap(fig, local_diverging_cmap());
    local_balanced_caxis(vel_diff);
    title('Velocity Check: vel\_obs - initial vel (m/yr)');
    local_save_both(fig, figdir, pngdir, '17_velocity_obs_minus_initial', dpi_png);
else
    warning('Cannot compute vel_obs - initial vel, skipping 17.');
end

% ---- 18 friction coefficient, linear -------------------------------------
local_safe_node_plot('md.friction.coefficient', fric_coef, ...
    faces, x, y, nv, edge_flag, 'Friction Coefficient', ...
    figdir, pngdir, '18_friction_coefficient_linear', dpi_png, false, 0, false);

% ---- 19 friction coefficient, log ----------------------------------------
local_safe_node_plot('md.friction.coefficient', fric_coef, ...
    faces, x, y, nv, edge_flag, 'Friction Coefficient  [log_{10}]', ...
    figdir, pngdir, '19_friction_coefficient_log', dpi_png, true, fric_log_floor, false);

% ---- 20 friction check on shelf / no-ice nodes ---------------------------
if ~isempty(fric_coef)
    fig = local_new_fig();
    local_plot_continuous(fig, faces, x, y, fric_coef, edge_flag);
    hold on;
    h1 = scatter(x(is_shelf), y(is_shelf), marker_size, 'b', 'filled');
    h2 = scatter(x(is_noice), y(is_noice), marker_size, 'r', 'filled');
    hold off;
    legend([h1 h2], {'shelf nodes','no-ice nodes'}, 'Location','best');
    title('Friction Coefficient Check on Shelf / No-Ice Nodes');
    local_save_both(fig, figdir, pngdir, '20_friction_check_shelf_noice', dpi_png);
else
    warning('friction.coefficient not found, skipping 20.');
end

% ---- 21 rheology B -------------------------------------------------------
local_safe_node_plot('md.materials.rheology_B', rheo_B, ...
    faces, x, y, nv, edge_flag, 'Rheology B', ...
    figdir, pngdir, '21_rheology_B', dpi_png, false, 0, false);

% ---- 22 rheology n -------------------------------------------------------
if ~isempty(rheo_n)
    fig = local_new_fig();
    if numel(rheo_n) == ne
        rheo_n_node = local_element_to_vertex(faces, rheo_n, nv);
        local_plot_continuous(fig, faces, x, y, rheo_n_node, edge_flag);
    else
        local_plot_continuous(fig, faces, x, y, local_clip(rheo_n, nv), edge_flag);
    end
    title('Rheology n');
    local_save_both(fig, figdir, pngdir, '22_rheology_n', dpi_png);
else
    warning('materials.rheology_n not found, skipping 22.');
end

% ---- 23 stressbalance boundary condition nodes ---------------------------
fig = local_new_fig();
local_plot_background(faces, x, y, [0.85 0.85 0.85]);
hold on;
idx_bc = find(bc_mask);
scatter(x(idx_bc), y(idx_bc), marker_size, 'r', 'filled');
hold off;
axis equal tight; box on; xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Stressbalance BC Nodes (N = %d)', numel(idx_bc)));
local_save_both(fig, figdir, pngdir, '23_stressbalance_BC_nodes', dpi_png);

% ---- 24 grounded ice melting rate ----------------------------------------
local_safe_node_plot('md.basalforcings.groundedice_melting_rate', grounded_melt, ...
    faces, x, y, nv, edge_flag, 'Grounded Ice Melting Rate', ...
    figdir, pngdir, '24_groundedice_melting_rate', dpi_png, false, 0, false);

% ---- 25 floating ice melting rate ----------------------------------------
local_safe_node_plot('md.basalforcings.floatingice_melting_rate', floating_melt, ...
    faces, x, y, nv, edge_flag, 'Floating Ice Melting Rate', ...
    figdir, pngdir, '25_floatingice_melting_rate', dpi_png, false, 0, false);

% ---- 26 initial water column ---------------------------------------------
local_safe_node_plot('md.initialization.watercolumn', watercolumn, ...
    faces, x, y, nv, edge_flag, 'Initial Water Column', ...
    figdir, pngdir, '26_initial_watercolumn', dpi_png, false, 0, false);

% ---- 27 initial water fraction -------------------------------------------
local_safe_node_plot('md.initialization.waterfraction', waterfraction, ...
    faces, x, y, nv, edge_flag, 'Initial Water Fraction', ...
    figdir, pngdir, '27_initial_waterfraction', dpi_png, false, 0, false);

% ---- 28 cost function coefficients column 1 ------------------------------
if ~isempty(cost_coeff) && size(cost_coeff,1) >= nv && size(cost_coeff,2) >= 1
    local_safe_node_plot('md.inversion.cost_functions_coefficients(:,1)', cost_coeff(1:nv,1), ...
        faces, x, y, nv, edge_flag, 'Cost Function Coefficients Column 1', ...
        figdir, pngdir, '28_cost_coefficients_col1', dpi_png, false, 0, false);
else
    warning('cost_functions_coefficients column 1 not available, skipping 28.');
end

% ---- 29 cost function coefficients column 2 ------------------------------
if ~isempty(cost_coeff) && size(cost_coeff,1) >= nv && size(cost_coeff,2) >= 2
    local_safe_node_plot('md.inversion.cost_functions_coefficients(:,2)', cost_coeff(1:nv,2), ...
        faces, x, y, nv, edge_flag, 'Cost Function Coefficients Column 2', ...
        figdir, pngdir, '29_cost_coefficients_col2', dpi_png, false, 0, false);
else
    warning('cost_functions_coefficients column 2 not available, skipping 29.');
end

% ---- 30 summary histograms -----------------------------------------------
fig = local_new_fig();
local_plot_histogram_summary(thickness, vel_init, fric_coef, rheo_B, is_grounded);
local_save_both(fig, figdir, pngdir, '30_histograms_summary', dpi_png);

%% ===== Done ==============================================================
fprintf('\nDone. Parameterize field figures generated.\n');
fprintf('FIG -> %s\n', figdir);
fprintf('PNG -> %s\n', pngdir);
close all;

%% ========================================================================
%  Local functions
%  ========================================================================

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
    if isstruct(obj)
        tf = isfield(obj, name);
    elseif isobject(obj)
        tf = isprop(obj, name);
    else
        tf = false;
    end
end

function v = local_try_field(md, path)
    v = local_deep_get(md, path);
    if ~isempty(v), v = double(v); end
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
    if isempty(v), return; end
    v = double(v(:));
    if numel(v) >= nv
        v = v(1:nv);
    else
        warning('Field has only %d entries, expected %d. Skipping field.', numel(v), nv);
        v = [];
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
    node_type = node_type(:);

    info.x_bm = x_bm;
    info.y_bm = y_bm;
end

function local_safe_node_plot(label, data, faces, x, y, nv, edge_flag, ttl, ...
                              figdir, pngdir, basename, dpi_png, use_log, floor_val, diverging)
    if isempty(data)
        warning('%s not found, skipping %s.', label, basename);
        return;
    end
    data = double(data(:));
    if numel(data) < nv
        warning('%s has fewer than nv entries, skipping %s.', label, basename);
        return;
    end
    data = data(1:nv);

    fig = local_new_fig();
    if use_log
        local_plot_log_field(fig, faces, x, y, data, edge_flag, floor_val);
    else
        local_plot_continuous(fig, faces, x, y, data, edge_flag);
        if diverging
            colormap(fig, local_diverging_cmap());
            local_balanced_caxis(data);
        end
    end
    title(ttl);
    local_save_both(fig, figdir, pngdir, basename, dpi_png);
end

function local_plot_continuous(fig, faces, x, y, data, show_edges)
    figure(fig);
    if show_edges, ec = [0.3 0.3 0.3]; else, ec = 'none'; end
    patch('Faces',faces,'Vertices',[x y], ...
          'FaceVertexCData',data(:), ...
          'FaceColor','interp','EdgeColor',ec);
    axis equal tight; box on;
    xlabel('x (m)'); ylabel('y (m)');
    colormap(fig, parula); colorbar;
end

function local_plot_log_field(fig, faces, x, y, data, show_edges, floor_val)
    figure(fig);
    if show_edges, ec = [0.3 0.3 0.3]; else, ec = 'none'; end
    log_data = log10(max(abs(data(:)), floor_val));
    patch('Faces',faces,'Vertices',[x y], ...
          'FaceVertexCData',log_data, ...
          'FaceColor','interp','EdgeColor',ec);
    axis equal tight; box on;
    xlabel('x (m)'); ylabel('y (m)');
    colormap(fig, parula);
    cb = colorbar;

    clean = log_data(isfinite(log_data));
    if isempty(clean), return; end
    lo = floor(min(clean));
    hi = ceil(max(clean));
    if lo >= hi, hi = lo + 1; end
    tks = lo:hi;
    cb.Ticks = tks;
    cb.TickLabels = arrayfun(@(t) sprintf('10^{%d}',t), tks, 'UniformOutput',false);
    clim([lo hi]);
end

function local_plot_background(faces, x, y, color)
    patch('Faces',faces,'Vertices',[x y], ...
          'FaceColor',color,'EdgeColor','none');
end

function local_plot_categorical(fig, faces, x, y, node_type)
    figure(fig);
    face_type = mode(double(node_type(faces)),2);
    cmap = [0 0.4 0.8; ...
            0.76 0.60 0.42; ...
            0.6 0.6 0.6; ...
            0 0.8 0.8; ...
            0.6 0.2 0.8];
    patch('Faces',faces,'Vertices',[x y], ...
          'FaceVertexCData',face_type(:), ...
          'FaceColor','flat','EdgeColor','none');
    colormap(fig, cmap);
    clim([-0.5 4.5]);
    cb = colorbar;
    cb.Ticks = 0:4;
    cb.TickLabels = {'ocean','land','grounded','shelf','vostok'};
    axis equal tight; box on;
    xlabel('x (m)'); ylabel('y (m)');
end

function local_balanced_caxis(data)
    clean = data(isfinite(data));
    if isempty(clean), return; end
    mx = max(abs(clean));
    if mx > 0
        clim([-mx mx]);
    end
end

function cmap = local_diverging_cmap()
    n = 128;
    r = [linspace(0,1,n)'; linspace(1,1,n)'];
    g = [linspace(0,1,n)'; linspace(1,0,n)'];
    b = [linspace(1,1,n)'; linspace(1,0,n)'];
    cmap = [r g b];
end

function local_plot_histogram_summary(thickness, vel_init, fric_coef, rheo_B, is_grounded)
    subplot(2,2,1);
    local_hist_one(thickness, 'Thickness (m)', 'Thickness Histogram');

    subplot(2,2,2);
    local_hist_one(vel_init, 'Initial velocity (m/yr)', 'Initial Velocity Histogram');

    subplot(2,2,3);
    if ~isempty(fric_coef)
        fc = fric_coef(is_grounded & isfinite(fric_coef));
    else
        fc = [];
    end
    local_hist_one(fc, 'Friction coefficient', 'Friction Histogram (grounded nodes)');

    subplot(2,2,4);
    local_hist_one(rheo_B, 'Rheology B', 'Rheology B Histogram');
end

function local_hist_one(data, xlabel_text, ttl)
    if isempty(data)
        title([ttl ' (missing)']);
        axis off;
        return;
    end
    clean = data(isfinite(data));
    if isempty(clean)
        title([ttl ' (no finite values)']);
        axis off;
        return;
    end
    histogram(clean, 100);
    xlabel(xlabel_text);
    ylabel('Count');
    title(ttl);
    grid on;
end

function node_data = local_element_to_vertex(faces, elem_data, nv)
    elem_data = double(elem_data(:));
    node_sum = zeros(nv,1);
    node_cnt = zeros(nv,1);

    for j = 1:size(faces,2)
        idx = faces(:,j);
        node_sum = node_sum + accumarray(idx, elem_data, [nv 1], @sum, 0);
        node_cnt = node_cnt + accumarray(idx, 1, [nv 1], @sum, 0);
    end

    node_data = NaN(nv,1);
    ok = node_cnt > 0;
    node_data(ok) = node_sum(ok) ./ node_cnt(ok);
end

function local_save_both(fig, figdir, pngdir, basename, dpi_png)
    savefig(fig, fullfile(figdir,[basename '.fig']));
    print(fig, fullfile(pngdir, basename), '-dpng', sprintf('-r%d',dpi_png));
end
