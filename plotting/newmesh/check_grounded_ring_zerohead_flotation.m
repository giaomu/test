%CHECK_GROUNDED_RING_ZEROHEAD_FLOTATION
% Diagnose what flotation/effective-pressure state is implied by setting
% the Step-5 one-ring grounded outlet nodes to spchead = head = 0.
%
% Run from the project root:
%   run('plotting/newmesh/check_grounded_ring_zerohead_flotation.m')

clear; clc; close all;

%% User settings
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

model_file_override = '';
default_model_file = fullfile(paths.models_newmesh, ...
    'RecoveryNewMesh_Hydrology_ShelfGroundedRingNoSlidingCoupled.mat');

% Node selection mode:
%   'step5_ring'          Reconstruct Step-5 shelf-expanded one grounded ring.
%   'actual_grounded_spc0' Use grounded nodes that are actually spchead=0.
selection_mode = 'step5_ring';

outdir = fullfile(paths.figures_newmesh, 'zerohead_grounded_ring_flotation');
save_png = true;
save_csv = true;

%% Load model
if ~isempty(model_file_override)
    model_file = model_file_override;
else
    model_file = default_model_file;
end
assert(exist(model_file, 'file') == 2, 'Model file not found: %s', model_file);

fprintf('Loading model: %s\n', model_file);
S = load(model_file, 'md');
md = S.md;
clear S

if ~exist(outdir, 'dir'), mkdir(outdir); end

%% Reconstruct masks
nv = md.mesh.numberofvertices;
elements = md.mesh.elements;

is_ice = true(nv, 1);
if has_member(md, 'mask') && has_member(md.mask, 'ice_levelset') && ~isempty(md.mask.ice_levelset)
    is_ice = md.mask.ice_levelset(:) < 0;
end

assert(has_member(md, 'mask') && has_member(md.mask, 'ocean_levelset') && ~isempty(md.mask.ocean_levelset), ...
    'md.mask.ocean_levelset is required.');

is_grounded = (md.mask.ocean_levelset(:) > 0) & is_ice;
is_shelf = (md.mask.ocean_levelset(:) < 0) & is_ice;

assert(numel(is_grounded) == nv && numel(is_shelf) == nv, ...
    'Mask lengths do not match md.mesh.numberofvertices.');

expanded_from_shelf = expand_nodes_by_elements_local(md, is_shelf, 1);
step5_grounded_ring = expanded_from_shelf & is_grounded;

spc = NaN(nv, 1);
if has_member(md, 'hydrology') && has_member(md.hydrology, 'spchead') && ~isempty(md.hydrology.spchead)
    spc = md.hydrology.spchead(:);
end
actual_grounded_spc0 = is_grounded & isfinite(spc) & abs(spc) <= 1e-9;

switch lower(selection_mode)
    case 'step5_ring'
        target = step5_grounded_ring;
    case 'actual_grounded_spc0'
        target = actual_grounded_spc0;
    otherwise
        error('Unknown selection_mode: %s', selection_mode);
end

assert(any(target), 'No target grounded outlet nodes selected.');

%% Compute flotation/effective-pressure state implied by head = 0
rhoi = md.materials.rho_ice;
rhow = md.materials.rho_freshwater;
g = md.constants.g;

H = md.geometry.thickness(:);
bed = md.geometry.base(:);

Pi = rhoi .* g .* H;
Pw_zero = rhow .* g .* (0 - bed);
N_zero = Pi - Pw_zero;

% These two ratios are the main diagnostics:
%   Pw_over_Pi_zero = 1 means water pressure is at flotation/overburden.
%   N_over_Pi_zero  = 0 means flotation; larger positive values mean lower
%                    water pressure and stronger creep closure.
Pw_over_Pi_zero = Pw_zero ./ Pi;
N_over_Pi_zero = N_zero ./ Pi;

head_float = bed + (rhoi ./ rhow) .* H;
head_zero_minus_float = 0 - head_float;

% Equivalent fN in Ntarget = fN * Pi notation.
fN_equiv = N_over_Pi_zero;

%% Print summary
fprintf('\nSelection mode: %s\n', selection_mode);
fprintf('Shelf nodes: %d\n', sum(is_shelf));
fprintf('Step-5 one-ring grounded nodes: %d\n', sum(step5_grounded_ring));
fprintf('Actual grounded spchead=0 nodes in model: %d\n', sum(actual_grounded_spc0));
fprintf('Selected target nodes: %d\n', sum(target));
fprintf('Overlap target with actual grounded spchead=0: %d\n', sum(target & actual_grounded_spc0));

print_stats('bed [m]', bed(target));
print_stats('thickness H [m]', H(target));
print_stats('head_float for N=0 [m]', head_float(target));
print_stats('head=0 minus head_float [m]', head_zero_minus_float(target));
print_stats('Pw(head=0)/Pi [-]', Pw_over_Pi_zero(target));
print_stats('N(head=0)/Pi [-] = equivalent fN', N_over_Pi_zero(target));

fprintf('\nInterpretation guide:\n');
fprintf('  Pw/Pi = 1, N/Pi = 0: flotation/overburden water pressure.\n');
fprintf('  Pw/Pi < 1, N/Pi > 0: below flotation; larger N/Pi means stronger drainage/closure.\n');
fprintf('  Pw/Pi > 1, N/Pi < 0: above overburden; SHAKTI caps water pressure internally.\n');

%% Save table
if save_csv
    T = table;
    target_ids = find(target);
    T.node_id = target_ids(:);
    T.x = reshape(md.mesh.x(target), [], 1);
    T.y = reshape(md.mesh.y(target), [], 1);
    T.bed_m = bed(target);
    T.thickness_m = H(target);
    T.head_float_m = head_float(target);
    T.head_zero_minus_float_m = head_zero_minus_float(target);
    T.Pw_over_Pi_head0 = Pw_over_Pi_zero(target);
    T.N_over_Pi_head0 = N_over_Pi_zero(target);
    T.equivalent_fN = fN_equiv(target);
    T.actual_spchead_m = spc(target);
    T.is_actual_grounded_spchead0 = actual_grounded_spc0(target);

    csv_file = fullfile(outdir, sprintf('zerohead_flotation_%s.csv', lower(selection_mode)));
    writetable(T, csv_file);
    fprintf('\nSaved CSV: %s\n', csv_file);
end

%% Save quick diagnostic figures
if save_png
    save_node_scatter(md, target, Pw_over_Pi_zero, ...
        'Pw(head=0) / Pi on selected grounded outlet nodes', ...
        fullfile(outdir, sprintf('Pw_over_Pi_%s.png', lower(selection_mode))));

    save_node_scatter(md, target, N_over_Pi_zero, ...
        'N(head=0) / Pi on selected grounded outlet nodes', ...
        fullfile(outdir, sprintf('N_over_Pi_%s.png', lower(selection_mode))));
end

fprintf('\nDone. Output folder: %s\n', outdir);

%% Local functions
function expanded = expand_nodes_by_elements_local(md, seed_nodes, nrings)
    expanded = logical(seed_nodes(:));
    elems = md.mesh.elements;
    for k = 1:nrings
        touched_elements = any(expanded(elems), 2);
        expanded(unique(elems(touched_elements, :))) = true;
    end
end

function print_stats(label, values)
    values = values(:);
    values = values(isfinite(values));
    if isempty(values)
        fprintf('%s: no finite values\n', label);
        return;
    end
    fprintf(['%s: n=%d, min=%.4g, p05=%.4g, median=%.4g, ' ...
        'mean=%.4g, p95=%.4g, max=%.4g\n'], ...
        label, numel(values), min(values), percentile_local(values, 5), ...
        median(values), mean(values), percentile_local(values, 95), max(values));
end

function save_node_scatter(md, target, values, plot_title, output_file)
    fig = figure('Visible', 'off', 'Color', 'w');
    scatter(md.mesh.x(target), md.mesh.y(target), 12, values(target), 'filled');
    axis equal tight;
    grid on;
    colorbar;
    title(plot_title, 'Interpreter', 'none');
    xlabel('x [m]');
    ylabel('y [m]');
    saveas(fig, output_file);
    close(fig);
    fprintf('Saved PNG: %s\n', output_file);
end

function tf = has_member(obj, name)
    if isstruct(obj)
        tf = isfield(obj, name);
    else
        tf = isprop(obj, name);
    end
end

function q = percentile_local(values, pct)
    values = sort(values(:));
    n = numel(values);
    if n == 1
        q = values;
        return;
    end
    pos = 1 + (pct / 100) * (n - 1);
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = values(lo);
    else
        q = values(lo) + (pos - lo) * (values(hi) - values(lo));
    end
end
