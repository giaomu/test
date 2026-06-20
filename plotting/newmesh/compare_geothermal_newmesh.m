%COMPARE_GEOTHERMAL_NEWMESH Plot old-source and new-mesh geothermal flux.
%
% Run from the project root:
%   run('plotting/newmesh/compare_geothermal_newmesh.m')

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

old_model_file = fullfile(project_root, '..', 'AIS_GlaDS_Initialize.mat');
new_model_file = fullfile(paths.models_newmesh, 'RecoveryNewMesh_Parameterize.mat');
outdir = fullfile(paths.figures_newmesh, 'RecoveryNewMesh_geothermal_compare');
figdir = fullfile(outdir, 'fig');
pngdir = fullfile(outdir, 'png');

if ~exist(old_model_file, 'file')
    error('Old geothermal source model not found: %s', old_model_file);
end
if ~exist(new_model_file, 'file')
    error('New parameterized model not found: %s', new_model_file);
end
if ~exist(outdir, 'dir'), mkdir(outdir); end
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(pngdir, 'dir'), mkdir(pngdir); end

fprintf('Loading old model: %s\n', old_model_file);
S = load(old_model_file, 'md');
md_old = S.md;
clear S

fprintf('Loading new model: %s\n', new_model_file);
S = load(new_model_file, 'md');
md_new = S.md;
clear S

x_old = md_old.mesh.x(:);
y_old = md_old.mesh.y(:);
faces_old = md_old.mesh.elements;
if size(faces_old, 2) > 3, faces_old = faces_old(:, 1:3); end
ghf_old = normalize_geothermal_units(md_old.basalforcings.geothermalflux(:), 'old');

x_new = md_new.mesh.x(:);
y_new = md_new.mesh.y(:);
faces_new = md_new.mesh.elements;
if size(faces_new, 2) > 3, faces_new = faces_new(:, 1:3); end
ghf_new = normalize_geothermal_units(md_new.basalforcings.geothermalflux(:), 'new');

if numel(ghf_old) ~= numel(x_old)
    error('Old geothermalflux length does not match old mesh vertices.');
end
if numel(ghf_new) ~= numel(x_new)
    error('New geothermalflux length does not match new mesh vertices.');
end

pad = 0.08 * max(range(x_new), range(y_new));
if pad <= 0 || ~isfinite(pad), pad = 20000; end
xlim_focus = [min(x_new) - pad, max(x_new) + pad];
ylim_focus = [min(y_new) - pad, max(y_new) + pad];

old_keep = elements_in_box(faces_old, x_old, y_old, xlim_focus, ylim_focus);
if ~any(old_keep)
    error('No old mesh elements overlap the new mesh bounding box. Check coordinate systems.');
end

old_vals = ghf_old(faces_old(old_keep, :));
all_vals = [old_vals(:); ghf_new(:)];
all_vals = all_vals(isfinite(all_vals));
cax = prctile(all_vals, [2 98]);
if ~all(isfinite(cax)) || cax(1) == cax(2)
    cax = [min(all_vals), max(all_vals)];
end

report_file = fullfile(outdir, 'geothermal_compare_report.txt');
fid = fopen(report_file, 'w');
cleanup = onCleanup(@() fclose_if_open(fid));
write_report(fid, 'Old source model', ghf_old);
write_report(fid, 'Old source in focus box', old_vals(:));
write_report(fid, 'New interpolated model', ghf_new);
fprintf('Report saved: %s\n', report_file);

fig = figure('Name', 'Geothermal flux compare', 'Color', 'w', ...
    'Position', [80 80 1500 720], 'Visible', 'on');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
draw_field(faces_old(old_keep, :), x_old, y_old, ghf_old, cax);
hold on;
plot_new_mesh_outline(faces_new, x_new, y_new);
title('Old source geothermal flux focused on new mesh', 'Interpreter', 'none');
xlim(xlim_focus); ylim(ylim_focus);

nexttile;
draw_field(faces_new, x_new, y_new, ghf_new, cax);
title('New mesh interpolated geothermal flux', 'Interpreter', 'none');
xlim(xlim_focus); ylim(ylim_focus);

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Geothermal flux (W/m^2)';

savefig(fig, fullfile(figdir, 'geothermal_old_vs_new_focus.fig'));
saveas(fig, fullfile(pngdir, 'geothermal_old_vs_new_focus.png'));

fig2 = figure('Name', 'New mesh over old geothermal flux', 'Color', 'w', ...
    'Position', [120 120 1000 850], 'Visible', 'on');
draw_field(faces_old(old_keep, :), x_old, y_old, ghf_old, cax);
hold on;
patch('Faces', faces_new, 'Vertices', [x_new, y_new], ...
    'FaceColor', 'none', 'EdgeColor', [0 0 0], 'LineWidth', 0.15);
title('New mesh over old-source geothermal flux', 'Interpreter', 'none');
xlim(xlim_focus); ylim(ylim_focus);
cb2 = colorbar;
cb2.Label.String = 'Geothermal flux (W/m^2)';
savefig(fig2, fullfile(figdir, 'new_mesh_over_old_geothermal.fig'));
saveas(fig2, fullfile(pngdir, 'new_mesh_over_old_geothermal.png'));

fprintf('Saved figures to:\n');
fprintf('  %s\n', figdir);
fprintf('  %s\n', pngdir);


function ghf = normalize_geothermal_units(ghf, label)
    ghf = ghf(:);
    valid = isfinite(ghf) & ghf >= 0;
    if ~any(valid)
        error('%s geothermal flux has no finite nonnegative values.', label);
    end
    m = mean(ghf(valid), 'omitnan');
    if m > 1
        fprintf('%s geothermal flux appears to be mW/m^2; converting to W/m^2.\n', label);
        ghf = ghf / 1000;
    end
end


function keep = elements_in_box(faces, x, y, xlim_box, ylim_box)
    cx = mean(x(faces), 2);
    cy = mean(y(faces), 2);
    keep = cx >= xlim_box(1) & cx <= xlim_box(2) & ...
           cy >= ylim_box(1) & cy <= ylim_box(2);
end


function draw_field(faces, x, y, data, cax)
    patch('Faces', faces, 'Vertices', [x(:), y(:)], ...
          'FaceVertexCData', data(:), ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    axis equal tight; box on;
    xlabel('X (m)');
    ylabel('Y (m)');
    colormap(turbo);
    caxis(cax);
end


function plot_new_mesh_outline(faces, x_new, y_new)
    edges = [faces(:, [1 2]); faces(:, [2 3]); faces(:, [3 1])];
    edges_sorted = sort(edges, 2);
    [unique_edges, ~, ic] = unique(edges_sorted, 'rows');
    edge_counts = accumarray(ic, 1);
    boundary_edges = unique_edges(edge_counts == 1, :);

    xx = [x_new(boundary_edges(:,1))'; x_new(boundary_edges(:,2))'; nan(1, size(boundary_edges,1))];
    yy = [y_new(boundary_edges(:,1))'; y_new(boundary_edges(:,2))'; nan(1, size(boundary_edges,1))];
    plot(xx(:), yy(:), 'k-', 'LineWidth', 1.2);
end


function write_report(fid, label, values)
    values = values(:);
    values = values(isfinite(values));
    fprintf(fid, '%s\n', label);
    fprintf(fid, '  n    = %d\n', numel(values));
    fprintf(fid, '  min  = %.10g W/m^2\n', min(values));
    fprintf(fid, '  p02  = %.10g W/m^2\n', prctile(values, 2));
    fprintf(fid, '  mean = %.10g W/m^2\n', mean(values));
    fprintf(fid, '  p98  = %.10g W/m^2\n', prctile(values, 98));
    fprintf(fid, '  max  = %.10g W/m^2\n\n', max(values));
end


function fclose_if_open(fid)
    if fid > 0
        fclose(fid);
    end
end
