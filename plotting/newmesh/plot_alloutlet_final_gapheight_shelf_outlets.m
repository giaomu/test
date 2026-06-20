% Plot final-step HydrologyGapHeight from the AllOutlet run, then overlay:
%   1) all ice-shelf nodes from BedMachine mask
%   2) Recovery_outlets.exp contour
%
% Run from the project root:
%   run('plotting/newmesh/plot_alloutlet_final_gapheight_shelf_outlets.m')

clear; clc; close all;

%% User settings

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

model_file = fullfile(paths.models_newmesh, 'RecoveryNewMesh_SHAKTI_0d_to_360d_1800s_np15_AllOutlet.mat');
outlet_exp_file = fullfile(paths.exp, 'Recovery_outlets.exp');
bedmachine_file = paths.bedmachine_file;
output_dir = fullfile(paths.figures_newmesh, 'RecoveryNewMesh_AllOutlet_final_gapheight_shelf_outlets');

% NaN 表示自动使用最后一个 TransientSolution。
target_step = NaN;

figure_visible = 'on';
shelf_marker_size = 8;

%% Load model and fields

assert(exist(model_file, 'file') == 2, 'Model file not found: %s', model_file);
assert(exist(outlet_exp_file, 'file') == 2, 'Outlet exp file not found: %s', outlet_exp_file);
assert(exist(bedmachine_file, 'file') == 2, 'BedMachine file not found: %s', bedmachine_file);

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('Loading model: %s\n', model_file);
S = load(model_file, 'md');
md = S.md;
clear S

[gap_height, target_step, time_yr, nsteps] = get_gapheight_at_step(md, target_step, model_file);
fprintf('Using step %d/%d, time = %.6g yr = %.4f days\n', target_step, nsteps, time_yr, time_yr * 365);

fprintf('Reading BedMachine shelf mask...\n');
bm_mask = interpBedmachineAntarctica(md.mesh.x, md.mesh.y, 'mask', 'nearest', bedmachine_file);
shelf_nodes = (bm_mask(:) == 3);
fprintf('Shelf nodes: %d\n', sum(shelf_nodes));

outlets = read_exp_contours_local(outlet_exp_file);
assert(~isempty(outlets), 'No contours found in %s.', outlet_exp_file);

%% Plot

fig = figure('Name', 'Final gap height with shelf nodes and outlet contours', ...
    'Color', 'w', 'Position', [80 80 1200 900], ...
    'Visible', figure_visible);
ax = axes(fig);

draw_mesh_field(ax, md, gap_height);
colormap(ax, turbo);
cb = colorbar(ax);
ylabel(cb, 'Hydrology gap height (m)');
title(ax, sprintf('HydrologyGapHeight final step %d/%d | shelf nodes + Recovery\\_outlets.exp', ...
    target_step, nsteps), 'Interpreter', 'tex');

hold(ax, 'on');
if any(shelf_nodes)
    scatter(ax, md.mesh.x(shelf_nodes), md.mesh.y(shelf_nodes), shelf_marker_size, ...
        'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25);
end

for k = 1:numel(outlets)
    plot(ax, outlets(k).x, outlets(k).y, 'r-', 'LineWidth', 2.2);
    plot(ax, outlets(k).x, outlets(k).y, 'ko', 'MarkerSize', 3, 'LineWidth', 0.8);
end

legend(ax, {'HydrologyGapHeight', 'ice-shelf nodes', 'Recovery\_outlets.exp', 'EXP vertices'}, ...
    'Location', 'bestoutside', 'Interpreter', 'tex');

base_name = sprintf('HydrologyGapHeight_step_%03d_with_shelf_nodes_and_Recovery_outlets', target_step);
fig_file = fullfile(output_dir, [base_name '.fig']);
png_file = fullfile(output_dir, [base_name '.png']);
shelf_table_file = fullfile(output_dir, [base_name '_shelf_nodes.csv']);

savefig(fig, fig_file);
exportgraphics(fig, png_file, 'Resolution', 300);

node_id = find(shelf_nodes);
shelf_table = table(node_id(:), md.mesh.x(shelf_nodes), md.mesh.y(shelf_nodes), ...
    'VariableNames', {'node_id', 'x', 'y'});
writetable(shelf_table, shelf_table_file);

fprintf('Saved figure: %s\n', fig_file);
fprintf('Saved PNG   : %s\n', png_file);
fprintf('Saved shelf nodes table: %s\n', shelf_table_file);

%% Local functions

function [gap_height, step, time_yr, nsteps] = get_gapheight_at_step(md, step, model_file)
    try
        results = md.results;
    catch
        error('No md.results found in %s.', model_file);
    end

    if ~isfield(results, 'TransientSolution') || isempty(results.TransientSolution)
        error('No md.results.TransientSolution found in %s.', model_file);
    end

    sol = results.TransientSolution;
    nsteps = numel(sol);
    if isnan(step)
        step = nsteps;
    end
    if step < 1 || step > nsteps || step ~= round(step)
        error('target_step must be an integer between 1 and %d.', nsteps);
    end

    if isfield(sol(step), 'HydrologyGapHeight')
        gap_height = sol(step).HydrologyGapHeight(:);
    elseif isfield(sol(step), 'Hydrologygap_height')
        gap_height = sol(step).Hydrologygap_height(:);
    elseif isfield(sol(step), 'GapHeight')
        gap_height = sol(step).GapHeight(:);
    else
        error('No HydrologyGapHeight field found at step %d in %s.', step, model_file);
    end

    if isfield(sol(step), 'time')
        time_yr = sol(step).time;
    else
        time_yr = NaN;
    end
end

function draw_mesh_field(ax, md, data)
    data = data(:);
    if numel(data) == md.mesh.numberofvertices
        patch(ax, 'Faces', md.mesh.elements, ...
            'Vertices', [md.mesh.x(:), md.mesh.y(:)], ...
            'FaceVertexCData', double(data), ...
            'FaceColor', 'interp', ...
            'EdgeColor', 'none');
    elseif numel(data) == md.mesh.numberofelements
        patch(ax, 'Faces', md.mesh.elements, ...
            'Vertices', [md.mesh.x(:), md.mesh.y(:)], ...
            'FaceVertexCData', double(data), ...
            'FaceColor', 'flat', ...
            'EdgeColor', 'none');
    else
        error('Field length %d does not match vertices %d or elements %d.', ...
            numel(data), md.mesh.numberofvertices, md.mesh.numberofelements);
    end

    axis(ax, 'equal');
    axis(ax, 'tight');
    box(ax, 'on');
    xlabel(ax, 'X (m)');
    ylabel(ax, 'Y (m)');
end

function contours = read_exp_contours_local(exp_file)
    fid = fopen(exp_file, 'r');
    assert(fid > 0, 'Cannot open exp file: %s', exp_file);
    cleaner = onCleanup(@() fclose(fid));

    contours = struct('x', {}, 'y', {}, 'value', {});
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end

        line = strtrim(line);
        if isempty(line) || startsWith(line, '#')
            continue;
        end

        header = sscanf(line, '%f');
        if numel(header) < 1
            continue;
        end

        npoints = round(header(1));
        value = NaN;
        if numel(header) >= 2
            value = header(2);
        end

        coords = NaN(npoints, 2);
        nread = 0;
        while nread < npoints
            coord_line = fgetl(fid);
            if ~ischar(coord_line)
                break;
            end
            coord_line = strtrim(coord_line);
            if isempty(coord_line) || startsWith(coord_line, '#')
                continue;
            end
            xy = sscanf(coord_line, '%f');
            if numel(xy) >= 2
                nread = nread + 1;
                coords(nread, :) = xy(1:2).';
            end
        end

        if nread == npoints
            contours(end+1).x = coords(:, 1); %#ok<AGROW>
            contours(end).y = coords(:, 2);
            contours(end).value = value;
        else
            warning('Incomplete contour in %s. Expected %d points, read %d.', exp_file, npoints, nread);
            break;
        end
    end
end
