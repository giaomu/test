% Fig. 4.4 spatial comparison for Recovery inversion results.
% Reads MAT models, plots fields, and writes PNG files only.

clear; clc; close all;

%% User parameters
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

model_dir = paths.models;
output_dir = fullfile(paths.figures, 'Fig44_inversion_spatial_compare_fixed');

horizontal_w101 = [500 5000 20000 50000];
horizontal_w501_targets = [1e-4 1e-5];
horizontal_w501_1e5_match_targets = [1e-5 1.1e-5 1e-5 1e-5];
vertical_w101_500_w501 = [1e-6 1e-5 1e-4 5e-4];
vertical_w101_5000_w501 = [1e-6 1e-4 5e-4];
vertical_w101_40000_w501 = [1e-6 1e-5 1e-4 5e-4];

output_field_figures = true;
output_abs_error = false;
show_mesh_edges = false;
dpi_png = 600;
auto_color_prctile = [0 100];
diff_color_prctile = 99;
horizontal_1e5_diff_clim = [-100 100];
vertical_500_diff_clim = [-200 200];
velocity_log_floor = 1e-3;
exact_relative_tolerance = 1e-8;

use_manual_file_map = false;
manual_file_map = {
%   500,  1e-6,  'D:\...\Lcurve_run01_w501_1p000e-06.mat'
%   500,  1e-5,  'D:\...\Lcurve_run05_w501_1p000e-05.mat'
%   500,  1e-4,  'D:\...\Lcurve_run18_w501_1p000e-04.mat'
%   500,  5e-4,  'D:\...\Lcurve_run22_w501_5p000e-04.mat'
%   5000, 1e-6,  'D:\...\Lcurve_run01_w501_1p000e-06.mat'
%   5000, 1e-4,  'D:\...\Lcurve_run18_w501_9p227e-05.mat'
%   5000, 5e-4,  'D:\...\Lcurve_run24_w501_4p557e-04.mat'
};

%% Model matching
fprintf('\nFig44 inversion spatial comparison\n');
fprintf('Model root : %s\n', model_dir);
fprintf('Output root: %s\n\n', output_dir);

ensure_dir(output_dir);
out_observed = fullfile(output_dir, 'observed');
out_horizontal_1e4 = fullfile(output_dir, 'horizontal_w501_1e-4');
out_horizontal_1e5 = fullfile(output_dir, 'horizontal_w501_1e-5');
out_vertical_500 = fullfile(output_dir, 'vertical_w101_500');
out_vertical_5000 = fullfile(output_dir, 'vertical_w101_5000');
out_vertical_40000 = fullfile(output_dir, 'vertical_w101_40000');
ensure_dir(out_observed);
ensure_dir(out_horizontal_1e4);
ensure_dir(out_horizontal_1e5);
ensure_dir(out_vertical_500);
ensure_dir(out_vertical_5000);
ensure_dir(out_vertical_40000);

targets = build_unique_targets(horizontal_w101, horizontal_w501_targets, horizontal_w501_1e5_match_targets, ...
    vertical_w101_500_w501, vertical_w101_5000_w501, vertical_w101_40000_w501);
if use_manual_file_map
    catalog = catalog_from_manual_map(manual_file_map);
else
    catalog = scan_model_catalog(model_dir);
end
fprintf('Cataloged weighted model files: %d\n', numel(catalog));
matches = match_targets(targets, catalog, exact_relative_tolerance);

%% Read data
records = struct('ok', {}, 'w101', {}, 'w501_target', {}, 'w501_actual', {}, ...
    'match_is_exact', {}, 'rel_diff', {}, 'file', {}, 'x', {}, 'y', {}, ...
    'faces', {}, 'vel_model', {}, 'vel_obs', {}, 'diff_vel', {}, ...
    'abs_error', {}, 'friction', {});

for i = 1:numel(matches)
    m = matches(i);
    if ~m.ok
        fprintf('[MISS] w101=%g, target w501=%s: no model file found.\n', m.w101, sci(m.w501_target));
        records(end + 1) = empty_record(m); %#ok<SAGROW>
        continue;
    end
    fprintf('[LOAD] w101=%g, target w501=%s, actual w501=%s\n       %s\n', ...
        m.w101, sci(m.w501_target), sci(m.w501_actual), m.file);
    try
        records(end + 1) = read_model_record(m); %#ok<SAGROW>
    catch ME
        warning('Failed to read model for w101=%g, w501=%s: %s', m.w101, sci(m.w501_target), ME.message);
        records(end + 1) = empty_record(m); %#ok<SAGROW>
    end
end

valid_records = records([records.ok]);
if isempty(valid_records)
    error('No selected model could be loaded. Check model_dir or manual_file_map.');
end

%% Plot settings
fprintf('\nColor limits are adaptive for each subplot.\n');
fprintf('  Velocity fields are displayed as log10(max(velocity, %g)).\n', velocity_log_floor);
fprintf('  Velocity differences remain linear and use the %.1fth absolute percentile for color limits.\n', diff_color_prctile);
fprintf('  Exceptions: horizontal w501~1e-5 diff uses [%g, %g]; vertical w101=500 diff uses [%g, %g].\n', ...
    horizontal_1e5_diff_clim(1), horizontal_1e5_diff_clim(2), vertical_500_diff_clim(1), vertical_500_diff_clim(2));
fprintf('  Friction coefficients remain linear.\n\n');

%% Plot figures
saved_files = {};
obs_rec = first_record_with_obs(valid_records);
if ~isempty(obs_rec)
    saved_files = [saved_files, plot_observed_velocity(obs_rec, out_observed, dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
else
    warning('No observed velocity was found in loaded models; observed velocity figure skipped.');
end

horizontal_1e4_recs = select_records(records, horizontal_w101, repmat(1e-4, size(horizontal_w101)));
horizontal_1e5_recs = select_records(records, horizontal_w101, horizontal_w501_1e5_match_targets);
horizontal_1e5_recs = override_w501_title(horizontal_1e5_recs, 1e-5);
vertical_500_recs = select_records(records, repmat(500, size(vertical_w101_500_w501)), vertical_w101_500_w501);
vertical_5000_recs = select_records(records, repmat(5000, size(vertical_w101_5000_w501)), vertical_w101_5000_w501);
vertical_40000_recs = select_records(records, repmat(40000, size(vertical_w101_40000_w501)), vertical_w101_40000_w501);

saved_files = [saved_files, plot_group_all(horizontal_1e4_recs, ...
    'Spatial comparison for different w101 (w501 ~ 1e-4)', ...
    fullfile(out_horizontal_1e4, 'fig44_horizontal_w501_1e-4_all.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
saved_files = [saved_files, plot_group_all(horizontal_1e5_recs, ...
    'Spatial comparison for different w101 (w501 ~ 1e-5)', ...
    fullfile(out_horizontal_1e5, 'fig44_horizontal_w501_1e-5_all.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor, horizontal_1e5_diff_clim)]; %#ok<AGROW>
saved_files = [saved_files, plot_group_all(vertical_500_recs, ...
    'Spatial comparison for different w501 (w101 = 500)', ...
    fullfile(out_vertical_500, 'fig44_vertical_w101_500_all.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor, vertical_500_diff_clim)]; %#ok<AGROW>
saved_files = [saved_files, plot_group_all(vertical_5000_recs, ...
    'Spatial comparison for different w501 (w101 = 5000)', ...
    fullfile(out_vertical_5000, 'fig44_vertical_w101_5000_all.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
saved_files = [saved_files, plot_group_all(vertical_40000_recs, ...
    'Spatial comparison for different w501 (w101 = 40000)', ...
    fullfile(out_vertical_40000, 'fig44_vertical_w101_40000_all.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>

if output_field_figures
    saved_files = [saved_files, plot_group_field(horizontal_1e4_recs, 'vel_model', ...
        'Horizontal comparison: modeled velocity log10 (w501 ~ 1e-4)', ...
        fullfile(out_horizontal_1e4, 'fig44_horizontal_w501_1e-4_velocity.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(horizontal_1e4_recs, 'diff_vel', ...
        'Horizontal comparison: velocity difference (w501 ~ 1e-4)', ...
        fullfile(out_horizontal_1e4, 'fig44_horizontal_w501_1e-4_diff.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(horizontal_1e4_recs, 'friction', ...
        'Horizontal comparison: basal friction coefficient (w501 ~ 1e-4)', ...
        fullfile(out_horizontal_1e4, 'fig44_horizontal_w501_1e-4_friction.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(horizontal_1e5_recs, 'vel_model', ...
        'Horizontal comparison: modeled velocity log10 (w501 ~ 1e-5)', ...
        fullfile(out_horizontal_1e5, 'fig44_horizontal_w501_1e-5_velocity.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(horizontal_1e5_recs, 'diff_vel', ...
        'Horizontal comparison: velocity difference (w501 ~ 1e-5)', ...
        fullfile(out_horizontal_1e5, 'fig44_horizontal_w501_1e-5_diff.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor, horizontal_1e5_diff_clim)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(horizontal_1e5_recs, 'friction', ...
        'Horizontal comparison: basal friction coefficient (w501 ~ 1e-5)', ...
        fullfile(out_horizontal_1e5, 'fig44_horizontal_w501_1e-5_friction.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_500_recs, 'vel_model', ...
        'Vertical comparison: modeled velocity log10 (w101 = 500)', ...
        fullfile(out_vertical_500, 'fig44_vertical_w101_500_velocity.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_500_recs, 'diff_vel', ...
        'Vertical comparison: velocity difference (w101 = 500)', ...
        fullfile(out_vertical_500, 'fig44_vertical_w101_500_diff.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor, vertical_500_diff_clim)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_500_recs, 'friction', ...
        'Vertical comparison: basal friction coefficient (w101 = 500)', ...
        fullfile(out_vertical_500, 'fig44_vertical_w101_500_friction.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_5000_recs, 'vel_model', ...
        'Vertical comparison: modeled velocity log10 (w101 = 5000)', ...
        fullfile(out_vertical_5000, 'fig44_vertical_w101_5000_velocity.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_5000_recs, 'diff_vel', ...
        'Vertical comparison: velocity difference (w101 = 5000)', ...
        fullfile(out_vertical_5000, 'fig44_vertical_w101_5000_diff.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_5000_recs, 'friction', ...
        'Vertical comparison: basal friction coefficient (w101 = 5000)', ...
        fullfile(out_vertical_5000, 'fig44_vertical_w101_5000_friction.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_40000_recs, 'vel_model', ...
        'Vertical comparison: modeled velocity log10 (w101 = 40000)', ...
        fullfile(out_vertical_40000, 'fig44_vertical_w101_40000_velocity.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_40000_recs, 'diff_vel', ...
        'Vertical comparison: velocity difference (w101 = 40000)', ...
        fullfile(out_vertical_40000, 'fig44_vertical_w101_40000_diff.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_40000_recs, 'friction', ...
        'Vertical comparison: basal friction coefficient (w101 = 40000)', ...
        fullfile(out_vertical_40000, 'fig44_vertical_w101_40000_friction.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
end

if output_abs_error
    saved_files = [saved_files, plot_group_field(horizontal_1e4_recs, 'abs_error', ...
        'Horizontal comparison: absolute velocity error (w501 ~ 1e-4)', ...
        fullfile(out_horizontal_1e4, 'fig44_horizontal_w501_1e-4_abs_error.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(horizontal_1e5_recs, 'abs_error', ...
        'Horizontal comparison: absolute velocity error (w501 ~ 1e-5)', ...
        fullfile(out_horizontal_1e5, 'fig44_horizontal_w501_1e-5_abs_error.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_500_recs, 'abs_error', ...
        'Vertical comparison: absolute velocity error (w101 = 500)', ...
        fullfile(out_vertical_500, 'fig44_vertical_w101_500_abs_error.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_5000_recs, 'abs_error', ...
        'Vertical comparison: absolute velocity error (w101 = 5000)', ...
        fullfile(out_vertical_5000, 'fig44_vertical_w101_5000_abs_error.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
    saved_files = [saved_files, plot_group_field(vertical_40000_recs, 'abs_error', ...
        'Vertical comparison: absolute velocity error (w101 = 40000)', ...
        fullfile(out_vertical_40000, 'fig44_vertical_w101_40000_abs_error.png'), dpi_png, show_mesh_edges, auto_color_prctile, diff_color_prctile, velocity_log_floor)]; %#ok<AGROW>
end

readme_path = fullfile(output_dir, 'README_fig44_outputs.md');
write_output_readme(readme_path, matches, saved_files, diff_color_prctile, horizontal_1e5_diff_clim, vertical_500_diff_clim);

fprintf('\nAll output PNG files:\n');
for i = 1:numel(saved_files), fprintf('  %s\n', saved_files{i}); end
fprintf('Output README:\n  %s\n', readme_path);
fprintf('\nDone.\n');

%% Helper functions
function targets = build_unique_targets(horizontal_w101, horizontal_w501_targets, horizontal_w501_1e5_match_targets, w501_500, w501_5000, w501_40000)
    [h_w101, h_w501] = ndgrid(horizontal_w101(:), horizontal_w501_targets(:));
    w101 = [h_w101(:); horizontal_w101(:); repmat(500, numel(w501_500), 1); ...
        repmat(5000, numel(w501_5000), 1); repmat(40000, numel(w501_40000), 1)];
    w501 = [h_w501(:); horizontal_w501_1e5_match_targets(:); w501_500(:); w501_5000(:); w501_40000(:)];
    targets = struct('w101', {}, 'w501_target', {});
    keys = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    for i = 1:numel(w101)
        key = sprintf('%g_%0.16g', w101(i), w501(i));
        if ~isKey(keys, key)
            keys(key) = true;
            targets(end + 1).w101 = w101(i); %#ok<AGROW>
            targets(end).w501_target = w501(i);
        end
    end
end

function catalog = scan_model_catalog(model_dir)
    catalog = struct('file', {}, 'name', {}, 'folder', {}, 'w101', {}, 'w501', {});
    if ~isfolder(model_dir), warning('model_dir does not exist: %s', model_dir); return; end
    files = dir(fullfile(model_dir, '**', '*.mat'));
    for i = 1:numel(files)
        name = files(i).name;
        if contains(lower(name), 'summary'), continue; end
        f = fullfile(files(i).folder, name);
        text = [files(i).folder filesep name];
        [w101, ~] = parse_weight(text, 'w101');
        if isnan(w101), w101 = parse_w101_from_folder(text); end
        [w501, ~] = parse_weight(text, 'w501');
        if isnan(w101) || isnan(w501), continue; end
        catalog(end + 1).file = f; %#ok<AGROW>
        catalog(end).name = name;
        catalog(end).folder = files(i).folder;
        catalog(end).w101 = w101;
        catalog(end).w501 = w501;
    end
end

function catalog = catalog_from_manual_map(manual_file_map)
    catalog = struct('file', {}, 'name', {}, 'folder', {}, 'w101', {}, 'w501', {});
    for i = 1:size(manual_file_map, 1)
        if size(manual_file_map, 2) < 3 || isempty(manual_file_map{i, 3}), continue; end
        f = manual_file_map{i, 3};
        [folder, name, ext] = fileparts(f);
        [w501, ~] = parse_weight([folder filesep name ext], 'w501');
        if isnan(w501), w501 = manual_file_map{i, 2}; end
        catalog(end + 1).file = f; %#ok<AGROW>
        catalog(end).name = [name ext];
        catalog(end).folder = folder;
        catalog(end).w101 = manual_file_map{i, 1};
        catalog(end).w501 = w501;
    end
end

function matches = match_targets(targets, catalog, exact_tol)
    matches = struct('ok', {}, 'w101', {}, 'w501_target', {}, 'w501_actual', {}, ...
        'match_is_exact', {}, 'rel_diff', {}, 'file', {});
    for i = 1:numel(targets)
        target = targets(i);
        matches(i).ok = false;
        matches(i).w101 = target.w101;
        matches(i).w501_target = target.w501_target;
        matches(i).w501_actual = NaN;
        matches(i).match_is_exact = false;
        matches(i).rel_diff = NaN;
        matches(i).file = '';
        idx101 = find(abs([catalog.w101] - target.w101) < max(1e-12, abs(target.w101) * 1e-12));
        if isempty(idx101), fprintf('[MISS] No files found for w101=%g.\n', target.w101); continue; end
        w501_vals = [catalog(idx101).w501];
        rel = abs(w501_vals - target.w501_target) ./ max(abs(target.w501_target), eps);
        [best_rel, j] = min(rel);
        best = catalog(idx101(j));
        matches(i).ok = true;
        matches(i).w501_actual = best.w501;
        matches(i).rel_diff = best_rel;
        matches(i).file = best.file;
        matches(i).match_is_exact = best_rel <= exact_tol;
        if matches(i).match_is_exact
            fprintf('[MATCH] w101=%g, w501=%s -> %s\n', target.w101, sci(best.w501), best.name);
        else
            fprintf('[NEAREST] w101=%g, target w501=%s, actual w501=%s, relative difference=%.3g\n          %s\n', ...
                target.w101, sci(target.w501_target), sci(best.w501), best_rel, best.file);
        end
    end
end

function rec = read_model_record(m)
    S = load(m.file);
    md = extract_md(S);
    clear S;
    x = double(md.mesh.x(:));
    y = double(md.mesh.y(:));
    faces = double(md.mesh.elements);
    if size(faces, 2) > 3, faces = faces(:, 1:3); end
    nv = numel(x);
    vel_model = get_model_velocity(md, nv);
    vel_obs = get_observed_velocity(md, nv);
    friction = get_friction(md, nv);
    non_shelf_mask = get_non_shelf_velocity_mask(md, nv);
    vel_model(~non_shelf_mask) = NaN;
    vel_obs(~non_shelf_mask) = NaN;
    diff_vel = NaN(size(vel_model));
    abs_error = NaN(size(vel_model));
    if ~isempty(vel_model) && ~isempty(vel_obs) && numel(vel_model) == numel(vel_obs)
        diff_vel = vel_model - vel_obs;
        abs_error = abs(diff_vel);
    else
        warning('Velocity difference skipped for %s because model/observed velocity is missing or mismatched.', m.file);
    end
    rec.ok = true;
    rec.w101 = m.w101;
    rec.w501_target = m.w501_target;
    rec.w501_actual = m.w501_actual;
    rec.match_is_exact = m.match_is_exact;
    rec.rel_diff = m.rel_diff;
    rec.file = m.file;
    rec.x = x;
    rec.y = y;
    rec.faces = faces;
    rec.vel_model = clean_vector(vel_model);
    rec.vel_obs = clean_vector(vel_obs);
    rec.diff_vel = clean_vector(diff_vel);
    rec.abs_error = clean_vector(abs_error);
    rec.friction = clean_vector(friction);
end

function rec = empty_record(m)
    rec.ok = false;
    rec.w101 = m.w101;
    rec.w501_target = m.w501_target;
    rec.w501_actual = m.w501_actual;
    rec.match_is_exact = false;
    rec.rel_diff = m.rel_diff;
    rec.file = m.file;
    rec.x = [];
    rec.y = [];
    rec.faces = [];
    rec.vel_model = [];
    rec.vel_obs = [];
    rec.diff_vel = [];
    rec.abs_error = [];
    rec.friction = [];
end

function md = extract_md(S)
    if isfield(S, 'md'), md = S.md; return; end
    if isfield(S, 'md_recovery'), md = S.md_recovery; return; end
    fns = fieldnames(S);
    for k = 1:numel(fns)
        obj = S.(fns{k});
        if isobject(obj) && isa(obj, 'model'), md = obj; return; end
    end
    for k = 1:numel(fns)
        obj = S.(fns{k});
        if has_member(obj, 'mesh') && has_member(obj, 'friction'), md = obj; return; end
    end
    error('Cannot find md or md_recovery model variable in MAT file.');
end

function vel = get_model_velocity(md, nv)
    vel = [];
    sol = get_stressbalance_solution(md);
    if ~isempty(sol)
        vel = deep_get(sol, {'Vel'});
        if isempty(vel)
            vx = deep_get(sol, {'Vx'});
            vy = deep_get(sol, {'Vy'});
            if ~isempty(vx) && ~isempty(vy), vel = sqrt(double(vx(:)).^2 + double(vy(:)).^2); end
        end
    end
    vel = clip_or_pad(vel, nv);
end

function vel_obs = get_observed_velocity(md, nv)
    vel_obs = deep_get(md, {'inversion', 'vel_obs'});
    if isempty(vel_obs)
        vx_obs = deep_get(md, {'inversion', 'vx_obs'});
        vy_obs = deep_get(md, {'inversion', 'vy_obs'});
        if ~isempty(vx_obs) && ~isempty(vy_obs), vel_obs = sqrt(double(vx_obs(:)).^2 + double(vy_obs(:)).^2); end
    end
    if isempty(vel_obs)
        vel_obs = deep_get(md, {'initialization', 'vel'});
        if ~isempty(vel_obs), warning('Fallback: using md.initialization.vel as observed velocity.'); end
    end
    vel_obs = clip_or_pad(vel_obs, nv);
end

function friction = get_friction(md, nv)
    friction = deep_get(md, {'friction', 'coefficient'});
    if isempty(friction)
        sol = get_stressbalance_solution(md);
        friction = deep_get(sol, {'FrictionCoefficient'});
        if ~isempty(friction), warning('Fallback: using results.StressbalanceSolution.FrictionCoefficient as friction.'); end
    end
    friction = clip_or_pad(friction, nv);
end

function keep = get_non_shelf_velocity_mask(md, nv)
    keep = true(nv, 1);
    ocean_levelset = deep_get(md, {'mask', 'ocean_levelset'});
    if ~isempty(ocean_levelset)
        ocean_levelset = clip_or_pad(ocean_levelset, nv);
        shelf = ocean_levelset < 0;
        ice_levelset = deep_get(md, {'mask', 'ice_levelset'});
        if ~isempty(ice_levelset)
            ice_levelset = clip_or_pad(ice_levelset, nv);
            shelf = shelf & (ice_levelset <= 0);
        end
        keep(shelf) = false;
        fprintf('       velocity mask: removed %d shelf nodes using md.mask.ocean_levelset.\n', sum(shelf));
        return;
    end
    floating = deep_get(md, {'mask', 'vertexonfloatingice'});
    if ~isempty(floating)
        floating = clip_or_pad(floating, nv);
        shelf = floating ~= 0;
        keep(shelf) = false;
        fprintf('       velocity mask: removed %d shelf nodes using md.mask.vertexonfloatingice.\n', sum(shelf));
        return;
    end
    warning('No shelf mask found. Velocity fields are plotted without shelf-node removal.');
end

function sol = get_stressbalance_solution(md)
    sol = [];
    if ~has_member(md, 'results'), return; end
    results = md.results;
    if ~has_member(results, 'StressbalanceSolution'), return; end
    sol = results.StressbalanceSolution;
    if iscell(sol)
        sol = sol{end};
    elseif numel(sol) > 1
        sol = sol(end);
    end
end

function tf = has_member(obj, name)
    if isstruct(obj)
        tf = isfield(obj, name);
    elseif isobject(obj)
        tf = isprop(obj, name);
    else
        tf = false;
    end
end

function v = deep_get(obj, path)
    v = [];
    if isempty(obj), return; end
    cur = obj;
    for k = 1:numel(path)
        if has_member(cur, path{k})
            cur = cur.(path{k});
        else
            return;
        end
    end
    v = cur;
end

function v = clip_or_pad(v, n)
    if isempty(v), v = NaN(n, 1); return; end
    v = double(v(:));
    if numel(v) >= n
        v = v(1:n);
    else
        v = [v; NaN(n - numel(v), 1)];
    end
end

function v = clean_vector(v)
    v = double(v(:));
    v(~isfinite(v)) = NaN;
end

function [value, raw] = parse_weight(text, weight_name)
    value = NaN;
    raw = '';
    expr = [weight_name '[_\-= ]*([0-9]+(?:[pP\.][0-9]+)?(?:[eE][\+\-_]?\d+)?|[0-9]*\.[0-9]+(?:[eE][\+\-_]?\d+)?)'];
    tok = regexp(text, expr, 'tokens', 'once');
    if isempty(tok), return; end
    raw = tok{1};
    value = weight_text_to_double(raw);
end

function w101 = parse_w101_from_folder(text)
    w101 = NaN;
    tok = regexp(text, 'L_curve_models_([0-9]+(?:[pP\.][0-9]+)?(?:[eE][\+\-_]?\d+)?)_', 'tokens', 'once');
    if ~isempty(tok), w101 = weight_text_to_double(tok{1}); end
end

function val = weight_text_to_double(s)
    s = char(s);
    s = strrep(s, 'P', 'p');
    s = strrep(s, 'p', '.');
    s = regexprep(s, '[, ]', '');
    s = regexprep(s, '([eE])_', '$1-');
    val = str2double(s);
end

function q = local_prctile(vals, pct)
    vals = sort(vals(:));
    vals = vals(isfinite(vals));
    if isempty(vals), q = NaN; return; end
    if exist('prctile', 'file') == 2, q = prctile(vals, pct); return; end
    pos = 1 + (numel(vals) - 1) * pct / 100;
    lo = floor(pos); hi = ceil(pos);
    if lo == hi, q = vals(lo); else, q = vals(lo) + (vals(hi) - vals(lo)) * (pos - lo); end
end

function rec = first_record_with_obs(records)
    rec = [];
    for i = 1:numel(records)
        if records(i).ok && ~isempty(records(i).vel_obs) && any(isfinite(records(i).vel_obs))
            rec = records(i);
            return;
        end
    end
end

function selected = select_records(records, w101_list, w501_targets)
    selected = repmat(empty_record(struct('w101', NaN, 'w501_target', NaN, ...
        'w501_actual', NaN, 'rel_diff', NaN, 'file', '')), 1, numel(w101_list));
    for i = 1:numel(w101_list)
        idx = [];
        for j = 1:numel(records)
            same101 = records(j).w101 == w101_list(i);
            same501 = abs(records(j).w501_target - w501_targets(i)) <= max(eps, abs(w501_targets(i)) * 1e-12);
            if same101 && same501, idx = j; break; end
        end
        if ~isempty(idx)
            selected(i) = records(idx);
        else
            m.w101 = w101_list(i); m.w501_target = w501_targets(i);
            m.w501_actual = NaN; m.rel_diff = NaN; m.file = '';
            selected(i) = empty_record(m);
        end
    end
end

function recs = override_w501_title(recs, display_w501)
    for i = 1:numel(recs)
        if recs(i).ok
            recs(i).w501_target = display_w501;
            recs(i).w501_actual = display_w501;
            recs(i).match_is_exact = true;
        end
    end
end

function saved = plot_observed_velocity(rec, out_dir, dpi, show_edges, auto_pct, diff_pct, velocity_log_floor)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 780]);
    ax = axes(fig);
    patch_field(ax, rec, 'vel_obs', field_cmap('vel_obs'), show_edges, auto_pct, diff_pct, velocity_log_floor);
    title(ax, 'Observed velocity log10', 'Interpreter', 'none');
    add_manual_colorbar(ax, field_cmap('vel_obs'), colorbar_label('vel_obs'));
    saved = save_png(fig, fullfile(out_dir, 'fig44_observed_velocity.png'), dpi);
end

function saved = plot_group_all(recs, main_title, out_file, dpi, show_edges, auto_pct, diff_pct, velocity_log_floor, fixed_diff_clim)
    if nargin < 9, fixed_diff_clim = []; end
    nrow = numel(recs);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 80 1420 max(760, 310 * nrow)]);
    tl = tiledlayout(fig, nrow, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    fields = {'vel_model', 'diff_vel', 'friction'};
    labels = {'Modeled velocity log10', 'Velocity difference', 'Basal friction coefficient'};
    for r = 1:nrow
        for c = 1:3
            ax = nexttile(tl);
            if recs(r).ok
                patch_field(ax, recs(r), fields{c}, field_cmap(fields{c}), show_edges, auto_pct, diff_pct, velocity_log_floor, fixed_diff_clim);
                title(ax, sprintf('%s\n%s', labels{c}, record_title(recs(r))), 'Interpreter', 'none', 'FontSize', 9);
                add_manual_colorbar(ax, field_cmap(fields{c}), colorbar_label(fields{c}));
            else
                plot_missing(ax, recs(r), labels{c});
            end
        end
    end
    apply_sgtitle(tl, main_title);
    saved = save_png(fig, out_file, dpi);
end

function saved = plot_group_field(recs, field_name, main_title, out_file, dpi, show_edges, auto_pct, diff_pct, velocity_log_floor, fixed_diff_clim)
    if nargin < 10, fixed_diff_clim = []; end
    n = numel(recs);
    if n == 3
        nrow = 1; ncol = 3; pos = [100 120 1420 520];
    elseif n == 4
        nrow = 2; ncol = 2; pos = [100 80 1180 980];
    else
        nrow = ceil(sqrt(n)); ncol = ceil(n / nrow); pos = [100 80 1180 900];
    end
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', pos);
    tl = tiledlayout(fig, nrow, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:n
        ax = nexttile(tl);
        if recs(i).ok
            patch_field(ax, recs(i), field_name, field_cmap(field_name), show_edges, auto_pct, diff_pct, velocity_log_floor, fixed_diff_clim);
            title(ax, record_title(recs(i)), 'Interpreter', 'none', 'FontSize', 10);
            add_manual_colorbar(ax, field_cmap(field_name), colorbar_label(field_name));
        else
            plot_missing(ax, recs(i), field_label(field_name));
        end
    end
    apply_sgtitle(tl, main_title);
    saved = save_png(fig, out_file, dpi);
end

function patch_field(ax, rec, field_name, cmap, show_edges, auto_pct, diff_pct, velocity_log_floor, fixed_diff_clim)
    if nargin < 9, fixed_diff_clim = []; end
    data = display_field_data(rec.(field_name), field_name, velocity_log_floor);
    if show_edges, edge_color = [0.5 0.5 0.5]; else, edge_color = 'none'; end
    if numel(data) == numel(rec.x)
        patch(ax, 'Faces', rec.faces, 'Vertices', [rec.x rec.y], ...
            'FaceVertexCData', data(:), 'FaceColor', 'interp', 'EdgeColor', edge_color);
    elseif numel(data) == size(rec.faces, 1)
        patch(ax, 'Faces', rec.faces, 'Vertices', [rec.x rec.y], ...
            'FaceVertexCData', data(:), 'FaceColor', 'flat', 'EdgeColor', edge_color);
    else
        warning('Field %s has incompatible size in %s.', field_name, rec.file);
        cla(ax); text(ax, 0.5, 0.5, 'field size mismatch', 'HorizontalAlignment', 'center');
    end
    axis(ax, 'equal'); axis(ax, 'tight'); box(ax, 'on');
    xlabel(ax, 'x (m)'); ylabel(ax, 'y (m)');
    colormap(ax, cmap);
    clim = adaptive_color_limits(data, field_name, auto_pct, diff_pct, fixed_diff_clim);
    if all(isfinite(clim)) && clim(1) < clim(2), caxis(ax, clim); end
    set(ax, 'Color', 'w');
end

function cax = add_manual_colorbar(ax, cmap, label_text)
    cax = [];
    clim = get(ax, 'CLim');
    if ~all(isfinite(clim)) || clim(1) >= clim(2), return; end
    fig = ancestor(ax, 'figure');
    old_units = get(ax, 'Units');
    set(ax, 'Units', 'normalized');
    pos = get(ax, 'Position');
    set(ax, 'Units', old_units);

    cbar_w = min(0.012, pos(3) * 0.05);
    cbar_h = pos(4) * 0.62;
    cbar_x = pos(1) + pos(3) - cbar_w - 0.006;
    cbar_y = pos(2) + pos(4) * 0.19;
    cax = axes('Parent', fig, 'Units', 'normalized', ...
        'Position', [cbar_x cbar_y cbar_w cbar_h]);
    vals = linspace(clim(1), clim(2), size(cmap, 1))';
    imagesc(cax, 1, vals, vals);
    set(cax, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right', ...
        'FontSize', 7, 'Box', 'on');
    colormap(cax, cmap);
    cax.YLabel.String = label_text;
    cax.YLabel.Interpreter = 'none';
    cax.YLabel.FontSize = 8;
    uistack(cax, 'top');
end

function data = display_field_data(data, field_name, velocity_log_floor)
    data = clean_vector(data);
    if ismember(field_name, {'vel_model', 'vel_obs'})
        data(data < velocity_log_floor) = velocity_log_floor;
        data = log10(data);
    end
end

function lim = adaptive_color_limits(data, field_name, pct, diff_pct, fixed_diff_clim)
    vals = data(isfinite(data));
    if isempty(vals), lim = [NaN NaN]; return; end
    if strcmp(field_name, 'diff_vel')
        if numel(fixed_diff_clim) == 2 && all(isfinite(fixed_diff_clim)) && fixed_diff_clim(1) < fixed_diff_clim(2)
            lim = fixed_diff_clim;
            return;
        end
        half = local_prctile(abs(vals), diff_pct);
        if ~isfinite(half) || half <= 0, half = max(abs(vals)); end
        if ~isfinite(half) || half <= 0, half = 1; end
        lim = [-half half];
        return;
    end
    lim = [local_prctile(vals, pct(1)), local_prctile(vals, pct(2))];
    lim = expand_degenerate_limits(lim, vals);
end

function lim = expand_degenerate_limits(lim, vals)
    if all(isfinite(lim)) && lim(1) < lim(2), return; end
    mn = min(vals); mx = max(vals);
    if mn == mx
        pad = max(abs(mx) * 0.05, 1);
        lim = [mn - pad, mx + pad];
    else
        lim = [mn mx];
    end
end

function plot_missing(ax, rec, label)
    cla(ax); axis(ax, 'off');
    text(ax, 0.5, 0.55, 'missing model', 'HorizontalAlignment', 'center', 'FontSize', 11);
    text(ax, 0.5, 0.43, sprintf('w101=%g, w501=%s', rec.w101, sci(rec.w501_target)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Interpreter', 'none');
    title(ax, label, 'Interpreter', 'none');
end

function ttl = record_title(rec)
    target = sci(rec.w501_target); actual = sci(rec.w501_actual);
    if rec.match_is_exact
        ttl = sprintf('w101=%g, w501=%s', rec.w101, actual);
    else
        ttl = sprintf('w101=%g, w501~%s (actual %s)', rec.w101, target, actual);
    end
end

function label = field_label(field_name)
    switch field_name
        case 'vel_model', label = 'Modeled velocity log10';
        case 'diff_vel', label = 'Velocity difference';
        case 'abs_error', label = 'Absolute velocity error';
        case 'friction', label = 'Basal friction coefficient';
        case 'vel_obs', label = 'Observed velocity log10';
        otherwise, label = field_name;
    end
end

function label = colorbar_label(field_name)
    switch field_name
        case {'vel_model', 'vel_obs'}, label = 'log10(m/yr)';
        case {'diff_vel', 'abs_error'}, label = 'm/yr';
        otherwise, label = 'coefficient';
    end
end

function cmap = field_cmap(field_name)
    switch field_name
        case 'diff_vel', cmap = diverging_cmap(256);
        otherwise, cmap = parula(256);
    end
end

function cmap = diverging_cmap(n)
    if nargin < 1, n = 256; end
    n1 = floor(n / 2); n2 = n - n1;
    blue = [49 99 174] / 255; white = [1 1 1]; red = [178 24 43] / 255;
    cmap = [interp_color(blue, white, n1); interp_color(white, red, n2)];
end

function c = interp_color(c1, c2, n)
    t = linspace(0, 1, n)';
    c = c1 .* (1 - t) + c2 .* t;
end

function saved = save_png(fig, out_file, dpi)
    ensure_dir(fileparts(out_file));
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    print(fig, out_file, '-dpng', sprintf('-r%d', dpi));
    close(fig);
    saved = {out_file};
    fprintf('[SAVE] %s\n', out_file);
end

function apply_sgtitle(tl, txt)
    try
        title(tl, txt, 'Interpreter', 'none', 'FontSize', 13, 'FontWeight', 'bold');
    catch
        sgtitle(txt, 'Interpreter', 'none', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

function ensure_dir(d)
    if ~exist(d, 'dir'), mkdir(d); end
end

function s = sci(x)
    if isempty(x) || ~isfinite(x)
        s = 'NaN';
    else
        s = sprintf('%.3g', x);
        if abs(x) < 0.01 || abs(x) >= 1e4, s = sprintf('%.3e', x); end
    end
end

function write_output_readme(readme_path, matches, saved_files, diff_color_prctile, horizontal_1e5_diff_clim, vertical_500_diff_clim)
    fid = fopen(readme_path, 'w');
    if fid < 0, warning('Cannot write README: %s', readme_path); return; end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# Fig44 inversion spatial comparison outputs\n\n');
    fprintf(fid, 'Generated by `plotting/fig44_inversion_spatial_compare.m`.\n\n');
    fprintf(fid, '## Folder structure\n\n');
    fprintf(fid, '- `observed/`: observed velocity reference figure.\n');
    fprintf(fid, '- `horizontal_w501_1e-4/`: horizontal comparison at target `w501 = 1e-4`.\n');
    fprintf(fid, '- `horizontal_w501_1e-5/`: horizontal comparison at target `w501 = 1e-5`.\n');
    fprintf(fid, '- `vertical_w101_500/`: vertical comparison for `w101 = 500`.\n');
    fprintf(fid, '- `vertical_w101_5000/`: vertical comparison for `w101 = 5000`.\n');
    fprintf(fid, '- `vertical_w101_40000/`: vertical comparison for `w101 = 40000`.\n');
    fprintf(fid, '\nVelocity magnitudes are plotted as `log10(m/yr)`. Velocity differences and friction coefficients are linear.\n');
    fprintf(fid, 'Each subplot uses its own adaptive colorbar range. Velocity-difference plots use the %.1fth absolute percentile by default.\n', diff_color_prctile);
    fprintf(fid, 'Exceptions: `horizontal_w501_1e-5` velocity differences use `[%g, %g]`; `vertical_w101_500` velocity differences use `[%g, %g]`.\n\n', ...
        horizontal_1e5_diff_clim(1), horizontal_1e5_diff_clim(2), vertical_500_diff_clim(1), vertical_500_diff_clim(2));
    fprintf(fid, '## Matched model files\n\n');
    for i = 1:numel(matches)
        if matches(i).ok
            fprintf(fid, '- target `w101=%g, w501=%s`; actual `w501=%s`; rel.diff `%.3g`; file `%s`\n', ...
                matches(i).w101, sci(matches(i).w501_target), sci(matches(i).w501_actual), matches(i).rel_diff, matches(i).file);
        else
            fprintf(fid, '- missing target `w101=%g, w501=%s`\n', matches(i).w101, sci(matches(i).w501_target));
        end
    end
    fprintf(fid, '\n## Saved PNG files\n\n');
    for i = 1:numel(saved_files), fprintf(fid, '- `%s`\n', saved_files{i}); end
    delete(cleanup);
    fprintf('[SAVE] %s\n', readme_path);
end
