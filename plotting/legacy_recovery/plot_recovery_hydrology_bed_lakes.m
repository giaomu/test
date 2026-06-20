% Plot bed elevation from the configured Recovery_Hydrology.mat with lake outlines.
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
[~, this_folder] = fileparts(script_dir);
project_root = script_dir;
if strcmpi(this_folder, 'plotting')
    project_root = fileparts(fileparts(script_dir));
end
addpath(project_root);
paths = shaktiais_paths();

model_file = fullfile(paths.models, 'Recovery_Hydrology.mat');
lake_file  = fullfile(paths.exp, 'recovery_active_lakes.exp');
if exist(lake_file, 'file') ~= 2
    lake_file = fullfile(paths.root, '..', 'recovery_active_lakes.exp');
end
out_dir    = fullfile(paths.figures, 'Recovery_Hydrology_bed_lakes');
out_png    = fullfile(out_dir, 'Recovery_Hydrology_bed_lakes.png');
out_fig    = fullfile(out_dir, 'Recovery_Hydrology_bed_lakes.fig');

assert(exist(model_file, 'file') == 2, 'Model file not found: %s', model_file);
assert(exist(lake_file, 'file') == 2, 'Lake contour file not found: %s', lake_file);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

if exist('loadmodel', 'file') == 2
    md = loadmodel(model_file);
else
    S = load(model_file, '-mat');
    names = fieldnames(S);
    assert(numel(names) == 1, 'Expected exactly one model variable in %s', model_file);
    md = S.(names{1});
end

x = md.mesh.x(:);
y = md.mesh.y(:);
elements = md.mesh.elements(:, 1:3);
bed = md.geometry.base(:);

assert(numel(bed) == md.mesh.numberofvertices, ...
    'md.geometry.base length does not match mesh vertices.');

if exist('expread', 'file') == 2
    lake_contours = expread(lake_file);
else
    lake_contours = read_exp_contours_local(lake_file);
end

fig = figure('Color', 'w', 'Position', [100 100 1050 820]);
patch('Faces', elements, ...
      'Vertices', [x, y], ...
      'FaceVertexCData', bed, ...
      'FaceColor', 'interp', ...
      'EdgeColor', 'none');
hold on;

axis equal tight;
box on;
colormap(gca, parula);
c = colorbar;
c.Label.String = 'Model bed elevation: md.geometry.base (m)';
xlabel('X (m)');
ylabel('Y (m)');
%title('Recovery Hydrology model bed elevation with lake outlines', ...
%    'Interpreter', 'none');

for i = 1:numel(lake_contours)
    xi = lake_contours(i).x(:);
    yi = lake_contours(i).y(:);
    if isempty(xi) || isempty(yi), continue; end

    plot(xi, yi, 'k-', 'LineWidth', 2.0);
    plot(xi, yi, 'w--', 'LineWidth', 0.9);
end

savefig(fig, out_fig);
try
    exportgraphics(fig, out_png, 'Resolution', 300);
catch
    saveas(fig, out_png);
end

fprintf('Saved figure: %s\n', out_png);
fprintf('Saved MATLAB figure: %s\n', out_fig);


function contours = read_exp_contours_local(exp_file)
    fid = fopen(exp_file, 'r');
    assert(fid > 0, 'Could not open exp file: %s', exp_file);
    cleanup = onCleanup(@() fclose(fid));

    contours = struct('name', {}, 'x', {}, 'y', {});
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end

        if startsWith(strtrim(line), '## Name:')
            item.name = strtrim(extractAfter(strtrim(line), '## Name:'));

            fgetl(fid); % icon line
            fgetl(fid); % points-count header
            count_line = strsplit(strtrim(fgetl(fid)));
            npts = str2double(count_line{1});
            fgetl(fid); % coordinate header

            xy = nan(npts, 2);
            for k = 1:npts
                coord_line = strsplit(strtrim(fgetl(fid)));
                xy(k, 1) = str2double(coord_line{1});
                xy(k, 2) = str2double(coord_line{2});
            end

            item.x = xy(:, 1);
            item.y = xy(:, 2);
            contours(end + 1) = item; %#ok<AGROW>
        end
    end
end
