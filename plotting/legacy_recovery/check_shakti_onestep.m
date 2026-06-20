%CHECK_SHAKTI_ONESTEP Diagnose legacy Recovery SHAKTI transient results.
%
% Run from the project root:
%   run('plotting/legacy_recovery/check_shakti_onestep.m')
%
% This is kept for the older Recovery workflow, not the current new-mesh run.
clear; clc; close all;

%% Paths and thresholds
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

model_file   = fullfile(paths.models, 'Recovery_SHAKTI_5Steps_300s.mat');
initial_file = fullfile(paths.models, 'Recovery_Hydrology.mat');
outdir       = fullfile(paths.figures, 'Recovery_SHAKTI_5Step_300s_diagnostics');

color_limit_mode = 'minmax';

severe_N_threshold = -1e8;     % Pa
severe_ratio_threshold = 1.5;  % head above overburden ratio
assert(exist(model_file, 'file') == 2, ...
    'Cannot find %s. Run runme.m first to generate the SHAKTI result model.', model_file);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

report_file = fullfile(outdir, 'summary.txt');
if exist(report_file, 'file')
    delete(report_file);
end
diary(report_file);

fprintf('===== SHAKTI transient diagnostics =====\n');
fprintf('Model file: %s\n', model_file);
fprintf('Report dir: %s\n\n', outdir);

%% 璇诲彇妯″瀷鍜岀粨鏋?Sload = load(model_file, 'md');
md = Sload.md;
assert(isfield(md.results, 'TransientSolution') && ~isempty(md.results.TransientSolution), ...
    'No TransientSolution found in %s.', model_file);

% 榛樿璇婃柇鏈€鍚庝竴涓凡淇濆瓨鏃堕棿姝ャ€傝嫢姹傝В鎻愬墠宕╂簝锛岃繖閫氬父灏辨槸宕╂簝鍓?% 鏈€鏈変环鍊肩殑鐘舵€併€?solution_index = numel(md.results.TransientSolution);
S = md.results.TransientSolution(solution_index);
fprintf('Available transient solutions: %d\n', solution_index);
fprintf('Diagnosing solution index: %d\n', solution_index);
fprintf('Available times [yr]:');
for it = 1:solution_index
    fprintf(' %.10g', md.results.TransientSolution(it).time);
end
fprintf('\n');
fprintf('Result time: %.10g yr = %.4f hr\n', S.time, S.time * 8760);
fprintf('Configured dt: %.10g yr = %.4f hr\n', md.timestepping.time_step, md.timestepping.time_step * 8760);
fprintf('Configured final_time: %.10g yr = %.4f hr\n\n', md.timestepping.final_time, md.timestepping.final_time * 8760);
time_label = sprintf('绗?%d 涓凡淇濆瓨鏃堕棿姝ワ紝t = %.4f h (%.10g yr)', ...
    solution_index, S.time * 8760, S.time);
fprintf('Figures diagnose: %s\n\n', time_label);

md0 = [];
if exist(initial_file, 'file') == 2
    S0load = load(initial_file, 'md');
    md0 = S0load.md;
    fprintf('Initial hydrology file loaded: %s\n\n', initial_file);
else
    fprintf('Initial hydrology file not found: %s\n\n', initial_file);
end

fprintf('Mesh: %d vertices, %d elements\n', md.mesh.numberofvertices, md.mesh.numberofelements);
fprintf('spchead nodes: %d\n', sum(~isnan(md.hydrology.spchead)));
fprintf('hydrology.relaxation: %.4g\n', md.hydrology.relaxation);
fprintf('hydrology.storage: min %.4e, max %.4e\n', ...
    min(md.hydrology.storage(:)), max(md.hydrology.storage(:)));
fprintf('gap limits: [%.4e, %.4e] m\n\n', ...
    md.hydrology.gap_height_min, md.hydrology.gap_height_max);

%% 鍙栧嚭涓昏鐗╃悊鍦?rho_i = md.materials.rho_ice;
rho_w = md.materials.rho_freshwater;
g     = md.constants.g;

head = get_result_field(S, 'HydrologyHead');
gap  = get_result_field(S, 'HydrologyGapHeight');
N    = get_result_field(S, 'EffectivePressure');
basal_flux = get_result_field(S, 'HydrologyBasalFlux');
channelization = get_result_field(S, 'DegreeOfChannelization');

% 濡傛灉缁撴灉閲屾病鏈?EffectivePressure锛屽垯鎸夊畾涔変粠姘村ご鍙嶇畻銆?if isempty(N) && ~isempty(head)
    N = rho_i * g * md.geometry.thickness - rho_w * g * (head - md.geometry.base);
end

flotation_head = [];
head_ratio = [];
head_minus_surface = [];
head_minus_flotation = [];

if ~isempty(head)
    flotation_head = md.geometry.base + (rho_i / rho_w) * md.geometry.thickness;
    head_ratio = (head - md.geometry.base) ./ ((rho_i / rho_w) * md.geometry.thickness);
    head_minus_surface = head - md.geometry.surface;
    head_minus_flotation = head - flotation_head;
end

%% 鍦虹粺璁?fprintf('--- Primary SHAKTI result fields ---\n');
summarize_field(md, 'HydrologyHead [m]', head);
summarize_field(md, 'HydrologyGapHeight [m]', gap);
summarize_field(md, 'EffectivePressure [Pa]', N);
summarize_field(md, 'HydrologyBasalFlux', basal_flux);
summarize_field(md, 'DegreeOfChannelization', channelization);

fprintf('\n--- Geometry fields ---\n');
summarize_field(md, 'Thickness H [m]', md.geometry.thickness);
summarize_field(md, 'Bed/base elevation [m]', md.geometry.base);
summarize_field(md, 'Surface elevation [m]', md.geometry.surface);

fprintf('\n--- Derived fields ---\n');
if ~isempty(head)
    summarize_field(md, 'HeadOverburdenRatio [-]', head_ratio);
    summarize_field(md, 'HeadMinusBase [m]', head - md.geometry.base);
    summarize_field(md, 'HeadMinusSurface [m]', head_minus_surface);
    summarize_field(md, 'HeadMinusFloatationHead [m]', head_minus_flotation);
end
if ~isempty(N)
    fprintf('EffectivePressure < 0 nodes: %d\n', sum(N(:) < 0));
    fprintf('EffectivePressure < -1 MPa nodes: %d\n', sum(N(:) < -1e6));
    fprintf('EffectivePressure < -10 MPa nodes: %d\n', sum(N(:) < -1e7));
end
if ~isempty(gap)
    fprintf('Gap at max-limit elements: %d\n', sum(gap(:) >= 0.999 * md.hydrology.gap_height_max));
    fprintf('Gap at min-limit elements: %d\n', sum(gap(:) <= 1.001 * md.hydrology.gap_height_min));
end

%% 寮傚父鐐瑰垎绫?bad_head = false(md.mesh.numberofvertices, 1);
bad_N = false(md.mesh.numberofvertices, 1);
bad_both = false(md.mesh.numberofvertices, 1);
bad_union = false(md.mesh.numberofvertices, 1);

if ~isempty(N) && ~isempty(head_ratio)
    bad_head = head_ratio(:) > severe_ratio_threshold;
    bad_N = N(:) < severe_N_threshold;
    bad_both = bad_head & bad_N;
    bad_union = bad_head | bad_N;

    fprintf('\n--- Severe anomaly criteria ---\n');
    fprintf('Severe N threshold: %.3e Pa\n', severe_N_threshold);
    fprintf('Severe head ratio threshold: %.3f\n', severe_ratio_threshold);
    fprintf('Bad head nodes: %d\n', sum(bad_head));
    fprintf('Bad effective-pressure nodes: %d\n', sum(bad_N));
    fprintf('Bad both nodes: %d\n', sum(bad_both));
    fprintf('Bad union nodes: %d\n', sum(bad_union));

    print_bad_group_summary(md, 'Bad head union group', bad_union);
    print_bad_group_summary(md, 'Bad head-only group', bad_head & ~bad_N);
    print_bad_group_summary(md, 'Bad N-only group', bad_N & ~bad_head);
    print_bad_group_summary(md, 'Bad both group', bad_both);
end

%% 涓庡垵濮嬫按鏂囧満姣旇緝
if ~isempty(md0)
    fprintf('\n--- Changes from initial hydrology model Recovery_Hydrology.mat ---\n');
    if ~isempty(head) && numel(md0.hydrology.head) == numel(head)
        summarize_field(md, 'DeltaHead step1-initial [m]', head - md0.hydrology.head);
    end
    if ~isempty(gap) && numel(md0.hydrology.gap_height) == numel(gap)
        summarize_field(md, 'DeltaGap step1-initial [m]', gap - md0.hydrology.gap_height);
    end
end

%% 鍏抽敭寮傚父浣嶇疆
fprintf('\n--- Top anomaly locations ---\n');
print_top_locations(md, 'largest HydrologyHead', head, 'max', 10);
print_top_locations(md, 'largest HeadOverburdenRatio', head_ratio, 'max', 10);
print_top_locations(md, 'largest HeadMinusSurface', head_minus_surface, 'max', 10);
print_top_locations(md, 'smallest EffectivePressure', N, 'min', 10);
print_top_locations(md, 'largest HydrologyGapHeight', gap, 'max', 10);
if ~isempty(md0) && ~isempty(gap) && numel(md0.hydrology.gap_height) == numel(gap)
    print_top_locations(md, 'largest DeltaGap', gap - md0.hydrology.gap_height, 'max', 10);
end

%% 杈撳嚭寮傚父鐐硅〃
if ~isempty(N) && ~isempty(head_ratio)
    fprintf('\n--- Severe anomaly node tables ---\n');
    write_bad_node_table(md, bad_head, head, N, head_ratio, flotation_head, ...
        fullfile(outdir, 'bad_head_nodes.csv'), 'Bad head nodes');
    write_bad_node_table(md, bad_N, head, N, head_ratio, flotation_head, ...
        fullfile(outdir, 'bad_effective_pressure_nodes.csv'), 'Bad effective-pressure nodes');
    write_bad_node_table(md, bad_both, head, N, head_ratio, flotation_head, ...
        fullfile(outdir, 'bad_both_nodes.csv'), 'Bad both nodes');
    write_bad_node_table(md, bad_union, head, N, head_ratio, flotation_head, ...
        fullfile(outdir, 'bad_union_nodes.csv'), 'Bad union nodes');

    % 鍏煎鏃ф枃浠跺悕锛屼粛鐒惰緭鍑轰竴浠?union 琛ㄣ€?    write_bad_node_table(md, bad_union, head, N, head_ratio, flotation_head, ...
        fullfile(outdir, 'severe_head_pressure_anomalies.csv'), 'Severe union nodes');
end

%% 鐢诲浘
fprintf('\n--- Saving figures ---\n');
make_figure_primary(md, S, head, gap, N, outdir, color_limit_mode, time_label);
make_figure_boundaries(md, head, N, head_ratio, bad_union, outdir, color_limit_mode, time_label);
if ~isempty(md0)
    make_figure_changes(md, md0, head, gap, outdir, color_limit_mode, time_label);
end
if ~isempty(N) && ~isempty(head_ratio)
    make_figure_severe_anomalies(md, head, N, head_ratio, bad_head, bad_N, bad_both, bad_union, outdir, color_limit_mode, time_label);
    make_figure_bad_overlay_physics(md, head, gap, N, head_ratio, bad_union, outdir, color_limit_mode, time_label);
    make_figure_bad_driver_fields(md, S, md0, head, gap, N, head_ratio, flotation_head, bad_union, outdir, color_limit_mode, time_label);
    make_figure_hydraulic_gradient_fields(md, head, flotation_head, bad_union, outdir, color_limit_mode, time_label);
end
fprintf('Figures saved in %s\n', outdir);

diary off;
fprintf('\nDiagnostics complete. Text report: %s\n', report_file);

%% 鏈湴鍑芥暟
function data = get_result_field(S, name)
    if isstruct(S) && isfield(S, name)
        data = S.(name);
    elseif isobject(S) && isprop(S, name)
        data = S.(name);
    else
        data = [];
    end
    if ~isempty(data)
        data = double(data(:));
    end
end

function summarize_field(md, name, data)
    if isempty(data)
        fprintf('%-34s : not available\n', name);
        return;
    end

    data = double(data(:));
    finite = isfinite(data);
    valid = data(finite);
    if isempty(valid)
        fprintf('%-34s : all non-finite (n=%d)\n', name, numel(data));
        return;
    end

    if numel(data) == md.mesh.numberofvertices
        loc = 'node';
    elseif numel(data) == md.mesh.numberofelements
        loc = 'elem';
    else
        loc = 'unknown';
    end

    fprintf('%-34s : %-7s n=%6d finite=%6d min=% .6e p01=% .6e med=% .6e p99=% .6e max=% .6e mean=% .6e\n', ...
        name, loc, numel(data), sum(finite), min(valid), prctile(valid,1), median(valid), ...
        prctile(valid,99), max(valid), mean(valid));
end

function print_bad_group_summary(md, label, mask)
    ids = find(mask(:));
    fprintf('\n%s: %d nodes\n', label, numel(ids));
    if isempty(ids)
        return;
    end

    fprintf('  on boundary: %d\n', sum(md.mesh.vertexonboundary(ids)));
    fprintf('  with spchead: %d\n', sum(~isnan(md.hydrology.spchead(ids))));
    fprintf('  H <= 25.01 m: %d\n', sum(md.geometry.thickness(ids) <= 25.01));
    fprintf('  H <= 50 m: %d\n', sum(md.geometry.thickness(ids) <= 50));
    fprintf('  H min/median/max: %.2f %.2f %.2f m\n', ...
        min(md.geometry.thickness(ids)), median(md.geometry.thickness(ids)), max(md.geometry.thickness(ids)));
    fprintf('  base min/median/max: %.2f %.2f %.2f m\n', ...
        min(md.geometry.base(ids)), median(md.geometry.base(ids)), max(md.geometry.base(ids)));
    fprintf('  surface min/median/max: %.2f %.2f %.2f m\n', ...
        min(md.geometry.surface(ids)), median(md.geometry.surface(ids)), max(md.geometry.surface(ids)));
end

function print_top_locations(md, label, data, mode, nshow)
    if isempty(data)
        fprintf('%s: not available\n', label);
        return;
    end
    data = double(data(:));
    finite = isfinite(data);
    idx = find(finite);
    if isempty(idx)
        fprintf('%s: no finite values\n', label);
        return;
    end
    vals = data(idx);
    if strcmp(mode, 'min')
        [~, order] = sort(vals, 'ascend');
    else
        [~, order] = sort(vals, 'descend');
    end
    idx = idx(order(1:min(nshow, numel(order))));

    fprintf('%s:\n', label);
    for k = 1:numel(idx)
        id = idx(k);
        if numel(data) == md.mesh.numberofvertices
            fprintf(['  node %6d  value=% .6e  x=% .1f  y=% .1f  H=%.2f  ' ...
                     'base=%.2f  surface=%.2f  boundary=%d  spc=%d\n'], ...
                id, data(id), md.mesh.x(id), md.mesh.y(id), ...
                md.geometry.thickness(id), md.geometry.base(id), md.geometry.surface(id), ...
                md.mesh.vertexonboundary(id), ~isnan(md.hydrology.spchead(id)));
        elseif numel(data) == md.mesh.numberofelements
            verts = md.mesh.elements(id,:);
            fprintf('  elem %6d  value=% .6e  x=% .1f  y=% .1f  meanH=%.2f  meanBase=%.2f\n', ...
                id, data(id), mean(md.mesh.x(verts)), mean(md.mesh.y(verts)), ...
                mean(md.geometry.thickness(verts)), mean(md.geometry.base(verts)));
        else
            fprintf('  index %6d  value=% .6e\n', id, data(id));
        end
    end
end

function write_bad_node_table(md, mask, head, N, head_ratio, flotation_head, csv_file, label)
    ids = find(mask(:));
    fprintf('%s: %d nodes -> %s\n', label, numel(ids), csv_file);
    if isempty(ids)
        z = zeros(0,1);
        T = table(z, z, z, z, z, z, z, z, z, z, z, z, z, logical(z), logical(z), ...
            'VariableNames', {'node','x','y','H','base','surface','head','N','head_ratio', ...
            'flotation_head','head_minus_base','head_minus_surface','head_minus_flotation', ...
            'boundary','has_spchead'});
        writetable(T, csv_file);
        return;
    end

    T = table( ...
        ids(:), ...
        md.mesh.x(ids), ...
        md.mesh.y(ids), ...
        md.geometry.thickness(ids), ...
        md.geometry.base(ids), ...
        md.geometry.surface(ids), ...
        head(ids), ...
        N(ids), ...
        head_ratio(ids), ...
        flotation_head(ids), ...
        head(ids) - md.geometry.base(ids), ...
        head(ids) - md.geometry.surface(ids), ...
        head(ids) - flotation_head(ids), ...
        md.mesh.vertexonboundary(ids), ...
        ~isnan(md.hydrology.spchead(ids)), ...
        'VariableNames', {'node','x','y','H','base','surface','head','N','head_ratio', ...
        'flotation_head','head_minus_base','head_minus_surface','head_minus_flotation', ...
        'boundary','has_spchead'} ...
    );

    disp(T);
    writetable(T, csv_file);
end

function make_figure_primary(md, S, head, gap, N, outdir, color_limit_mode, time_label)
    fig = figure('Name', 'SHAKTI 01 鎬昏锛氫富瑕佹按鏂囧満', 'Color', 'w', 'Position', [50 50 1400 950]);
    tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; patchplot(md, head, '姘村ご HydrologyHead [m]', color_limit_mode);
    nexttile; patchplot(md, gap, '绌鸿厰楂樺害 Gap height [m]', color_limit_mode);
    nexttile; patchplot(md, N/1e6, '鏈夋晥鍘嬪姏 N [MPa]', color_limit_mode);
    nexttile; patchplot(md, get_result_field(S, 'HydrologyBasalFlux'), '鍩哄簳姘撮€氶噺 Basal flux', color_limit_mode);
    nexttile; patchplot(md, get_result_field(S, 'DegreeOfChannelization'), '閫氶亾鍖栨寚鏍?Degree of channelization', color_limit_mode);
    nexttile; patchplot(md, double(~isnan(md.hydrology.spchead)), '鍥哄畾姘村ご鑺傜偣 mask', color_limit_mode);

    add_figure_title(fig, ['01 鎬昏锛氫富瑕佹按鏂囧満 - ' time_label]);
    save_figure_pair(fig, outdir, '01_overview_primary_fields.png');
end

function make_figure_boundaries(md, head, N, head_ratio, bad_union, outdir, color_limit_mode, time_label)
    rho_i = md.materials.rho_ice;
    rho_w = md.materials.rho_freshwater;
    flotation_head = md.geometry.base + (rho_i / rho_w) * md.geometry.thickness;
    spc = md.hydrology.spchead;

    fig = figure('Name', 'SHAKTI 02 鍘嬪姏涓庤竟鐣岃瘖鏂?, 'Color', 'w', 'Position', [80 80 1400 900]);
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; patchplot(md, head_ratio, '姘村ご / 鍐拌鍘嬪姏姘村ご', color_limit_mode);
    nexttile; patchplot(md, head - flotation_head, '姘村ご - 鍐拌鍘嬪姏姘村ご [m]', color_limit_mode);
    nexttile; patchplot(md, double(bad_union), '涓ラ噸寮傚父骞堕泦 mask', color_limit_mode);
    nexttile; patchplot(md, spc, '鍥哄畾姘村ご spchead [m]', color_limit_mode);

    add_figure_title(fig, ['02 鍘嬪姏涓庤竟鐣岃瘖鏂?- ' time_label]);
    save_figure_pair(fig, outdir, '02_pressure_boundary_diagnostics.png');
end

function make_figure_severe_anomalies(md, head, N, head_ratio, bad_head, bad_N, bad_both, bad_union, outdir, color_limit_mode, time_label)
    head_ids = find(bad_head);
    N_ids = find(bad_N);
    both_ids = find(bad_both);
    union_ids = find(bad_union);

    fig = figure('Name', 'SHAKTI 04 寮傚父鐐瑰垎绫?, 'Color', 'w', 'Position', [90 90 1500 950]);
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    patchplot(md, head_ratio, '鑳屾櫙锛氭按澶?/ 鍐拌鍘嬪姏姘村ご', color_limit_mode);
    hold on;
    plot(md.mesh.x(head_ids), md.mesh.y(head_ids), 'ro', 'MarkerSize', 4, 'LineWidth', 1);
    title(sprintf('涓ラ噸姘村ご寮傚父鐐?(n=%d)', numel(head_ids)), 'Interpreter', 'none');

    nexttile;
    patchplot(md, N / 1e6, '鑳屾櫙锛氭湁鏁堝帇鍔?[MPa]', color_limit_mode);
    hold on;
    plot(md.mesh.x(N_ids), md.mesh.y(N_ids), 'ko', 'MarkerSize', 5, 'LineWidth', 1.2);
    plot(md.mesh.x(both_ids), md.mesh.y(both_ids), 'ro', 'MarkerSize', 4, 'LineWidth', 1);
    title(sprintf('涓ラ噸璐熸湁鏁堝帇鍔涚偣 榛戝湀(n=%d)锛岄噸鍚堢偣绾㈠湀(n=%d)', numel(N_ids), numel(both_ids)), 'Interpreter', 'none');

    nexttile;
    patchplot(md, double(bad_union), '涓ラ噸寮傚父骞堕泦 mask', color_limit_mode);

    nexttile;
    patchplot(md, md.geometry.surface, '鍐伴潰楂樼▼ surface [m]', color_limit_mode);
    hold on;
    plot(md.mesh.x(union_ids), md.mesh.y(union_ids), 'ro', 'MarkerSize', 4, 'LineWidth', 1);
    title(sprintf('涓ラ噸寮傚父鐐瑰彔鍔犲湪鍐伴潰楂樼▼涓?(n=%d)', numel(union_ids)), 'Interpreter', 'none');

    add_figure_title(fig, ['04 寮傚父鐐瑰垎绫伙細姘村ご寮傚父涓庢瀬绔礋鏈夋晥鍘嬪姏 - ' time_label]);
    save_figure_pair(fig, outdir, '04_anomaly_points.png');
end

function make_figure_bad_overlay_physics(md, head, gap, N, head_ratio, bad_union, outdir, color_limit_mode, time_label)
    % 鎶婁弗閲嶅紓甯哥偣鍙犲姞鍒板熀纭€鐗╃悊鍦轰笂锛岀敤浜庡垽鏂繖浜涚偣涓轰粈涔堝紓甯搞€?    bad_ids = find(bad_union);

    spc_mask = double(isfinite(md.hydrology.spchead(:)));
    head_minus_surface = head - md.geometry.surface(:);

    elem_bad = false(md.mesh.numberofelements, 1);
    if ~isempty(bad_ids)
        elem_bad = any(ismember(md.mesh.elements, bad_ids), 2);
    end
    elem_bad_ids = find(elem_bad);
    elem_x = mean(md.mesh.x(md.mesh.elements(elem_bad_ids,:)), 2);
    elem_y = mean(md.mesh.y(md.mesh.elements(elem_bad_ids,:)), 2);

    fig = figure('Name', 'SHAKTI 05 寮傚父鐐圭墿鐞嗚儗鏅?, ...
        'Color', 'w', 'Position', [40 40 1600 1050]);
    tiledlayout(3,3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    patchplot_with_bad(md, md.geometry.thickness(:), bad_ids, [], [], ...
        '鍐板帤 H [m] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, md.geometry.base(:), bad_ids, [], [], ...
        '鍩哄博楂樼▼ base [m] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, md.geometry.surface(:), bad_ids, [], [], ...
        '鍐伴潰楂樼▼ surface [m] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, head, bad_ids, [], [], ...
        '姘村ご HydrologyHead [m] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, N / 1e6, bad_ids, [], [], ...
        '鏈夋晥鍘嬪姏 N [MPa] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, head_ratio, bad_ids, [], [], ...
        '姘村ご/鍐拌鍘嬪姏姘村ご 姣斿€?+ 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, head_minus_surface, bad_ids, [], [], ...
        '姘村ご - 鍐伴潰 [m] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, gap, [], elem_x, elem_y, ...
        '绌鸿厰楂樺害 gap [m] + 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    nexttile;
    patchplot_with_bad(md, spc_mask, bad_ids, [], [], ...
        '鍥哄畾姘村ご鑺傜偣 mask + 寮傚父鐐?, color_limit_mode);

    add_figure_title(fig, ['05 寮傚父鐐瑰彔鍔犲湪鍩虹鐗╃悊鍦轰笂 - ' time_label]);
    save_figure_pair(fig, outdir, '05_anomaly_physical_context.png');
end

function make_figure_bad_driver_fields(md, S, md0, head, gap, N, head_ratio, flotation_head, bad_union, outdir, color_limit_mode, time_label)
    % 杩欎簺鍥句笓闂ㄧ敤鏉ヨВ閲婂紓甯哥偣涓轰粈涔堝嚭鐜帮細鍒濆鏉′欢銆佺鎺掓按杈圭晫璺濈銆侀€氶噺鍜?gap 鍙樺寲銆?    bad_ids = find(bad_union);
    spc_mask = isfinite(md.hydrology.spchead(:));
    spc_value = md.hydrology.spchead(:);
    spc_distance_km = distance_to_mask_nodes(md, spc_mask) / 1000;
    head_minus_flotation = head - flotation_head;

    initial_head = [];
    initial_head_ratio = [];
    delta_head = [];
    delta_gap = [];
    if ~isempty(md0)
        try
            initial_head_candidate = md0.hydrology.head(:);
        catch
            initial_head_candidate = [];
        end
        if numel(initial_head_candidate) == numel(head)
            initial_head = initial_head_candidate;
            initial_head_ratio = (initial_head - md.geometry.base(:)) ./ ...
                ((md.materials.rho_ice / md.materials.rho_freshwater) * md.geometry.thickness(:));
            delta_head = head - initial_head;
        end
        try
            initial_gap_candidate = md0.hydrology.gap_height(:);
        catch
            initial_gap_candidate = [];
        end
        if numel(initial_gap_candidate) == numel(gap)
            delta_gap = gap - initial_gap_candidate;
        end
    end

    basal_flux = get_result_field(S, 'HydrologyBasalFlux');
    channelization = get_result_field(S, 'DegreeOfChannelization');

    elem_bad = false(md.mesh.numberofelements, 1);
    if ~isempty(bad_ids)
        elem_bad = any(ismember(md.mesh.elements, bad_ids), 2);
    end
    [elem_x, elem_y] = element_centers(md, elem_bad);

    fig = figure('Name', 'SHAKTI 06 寮傚父鐐归┍鍔ㄥ洜绱?, ...
        'Color', 'w', 'Position', [60 60 1700 1050]);
    tiledlayout(3,3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    patchplot_with_bad(md, initial_head, bad_ids, [], [], ...
        '鍒濆姘村ご initial head [m] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, delta_head, bad_ids, [], [], ...
        '姘村ご澧為噺 螖head [m] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, initial_head_ratio, bad_ids, [], [], ...
        '鍒濆姘村ご/鍐拌鍘嬪姏姘村ご 姣斿€?+ 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, spc_distance_km, bad_ids, [], [], ...
        '鍒版渶杩?spchead 璺濈 [km] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, spc_value, bad_ids, [], [], ...
        'spchead 鏁板€?[m] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, head_minus_flotation, bad_ids, [], [], ...
        '姘村ご - 鍐拌鍘嬪姏姘村ご [m] + 寮傚父鐐?, color_limit_mode);

    nexttile;
    patchplot_with_bad(md, basal_flux, [], elem_x, elem_y, ...
        '鍩哄簳姘撮€氶噺 Basal flux + 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    nexttile;
    patchplot_with_bad(md, delta_gap, [], elem_x, elem_y, ...
        '绌鸿厰楂樺害澧為噺 螖gap [m] + 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    nexttile;
    patchplot_with_bad(md, channelization, [], elem_x, elem_y, ...
        '閫氶亾鍖栨寚鏍?+ 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    add_figure_title(fig, ['06 寮傚父鐐归┍鍔ㄥ洜绱狅細鍒濆鍦恒€佽竟鐣岃窛绂汇€侀€氶噺鍜?gap 鍙樺寲 - ' time_label]);
    save_figure_pair(fig, outdir, '06_anomaly_driver_fields.png');
end

function make_figure_hydraulic_gradient_fields(md, head, flotation_head, bad_union, outdir, color_limit_mode, time_label)
    % 姘村ご鏂圭▼鍙楁按鍔涘潯搴︽帶鍒讹紱杩欓噷鐢绘搴﹀ぇ灏忥紝鐪嬪紓甯哥偣鏄惁钀藉湪寮哄潯搴?灏侀棴鍔垮満闄勮繎銆?    bad_ids = find(bad_union);
    elem_bad = false(md.mesh.numberofelements, 1);
    if ~isempty(bad_ids)
        elem_bad = any(ismember(md.mesh.elements, bad_ids), 2);
    end
    [elem_x, elem_y] = element_centers(md, elem_bad);

    head_excess = head - flotation_head;
    grad_head = element_gradient_magnitude(md, head);
    grad_overburden = element_gradient_magnitude(md, flotation_head);
    grad_excess = element_gradient_magnitude(md, head_excess);
    grad_surface = element_gradient_magnitude(md, md.geometry.surface(:));
    grad_base = element_gradient_magnitude(md, md.geometry.base(:));
    grad_thickness = element_gradient_magnitude(md, md.geometry.thickness(:));

    fig = figure('Name', 'SHAKTI 07 姊害涓庣綉鏍煎昂搴︽尟鑽?, ...
        'Color', 'w', 'Position', [80 80 1600 950]);
    tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    patchplot_with_bad(md, grad_head, [], elem_x, elem_y, ...
        '姘村ご姊害 |鈭噃ead| [m/m] + 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    nexttile;
    patchplot_with_bad(md, grad_overburden, [], elem_x, elem_y, ...
        '鍐拌鍘嬪姏姘村ご姊害 |鈭噃over| [m/m] + 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    nexttile;
    patchplot_with_bad(md, grad_excess, [], elem_x, elem_y, ...
        '瓒呭帇姘村ご姊害 |鈭?head-hover)| [m/m] + 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    nexttile;
    patchplot_with_bad(md, grad_surface, [], elem_x, elem_y, ...
        '鍐伴潰鍧″害 |鈭噑urface| [m/m] + 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    nexttile;
    patchplot_with_bad(md, grad_base, [], elem_x, elem_y, ...
        '鍩哄博鍧″害 |鈭嘼ase| [m/m] + 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    nexttile;
    patchplot_with_bad(md, grad_thickness, [], elem_x, elem_y, ...
        '鍐板帤姊害 |鈭嘓| [m/m] + 寮傚父鐩稿叧鍗曞厓', color_limit_mode);

    add_figure_title(fig, ['07 姊害璇婃柇锛氭鏌ユ按澶撮碁鐗囩姸鎸崱 - ' time_label]);
    save_figure_pair(fig, outdir, '07_hydraulic_gradient_fields.png');
end

function make_figure_changes(md, md0, head, gap, outdir, color_limit_mode, time_label)
    fig = figure('Name', 'SHAKTI 03 鐩稿鍒濆鍦哄彉鍖?, 'Color', 'w', 'Position', [110 110 1200 850]);
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    if ~isempty(head) && numel(md0.hydrology.head) == numel(head)
        nexttile; patchplot(md, head - md0.hydrology.head, '姘村ご鍙樺寲 螖head [m]', color_limit_mode);
    else
        nexttile; axis off; title('姘村ご鍙樺寲涓嶅彲鐢?);
    end
    if ~isempty(gap) && numel(md0.hydrology.gap_height) == numel(gap)
        nexttile; patchplot(md, gap - md0.hydrology.gap_height, '绌鸿厰楂樺害鍙樺寲 螖gap [m]', color_limit_mode);
    else
        nexttile; axis off; title('绌鸿厰楂樺害鍙樺寲涓嶅彲鐢?);
    end
    nexttile; patchplot(md, md.geometry.thickness, '鍐板帤 H [m]', color_limit_mode);
    nexttile; patchplot(md, md.geometry.base, '鍩哄博楂樼▼ base [m]', color_limit_mode);

    add_figure_title(fig, ['03 鐩稿 Step 4 鍒濆姘存枃鍦虹殑鍙樺寲 - ' time_label]);
    save_figure_pair(fig, outdir, '03_changes_from_initial.png');
end

function save_figure_pair(fig, outdir, png_name)
    % 鍚屾椂淇濆瓨 PNG 鍜?MATLAB FIG銆侾NG 鐢ㄤ簬蹇€熼瑙堬紝FIG 鐢ㄤ簬鏀惧ぇ鍜屼氦浜掓煡鐪嬨€?    png_file = fullfile(outdir, png_name);
    [~, base_name, ~] = fileparts(png_name);
    fig_file = fullfile(outdir, [base_name '.fig']);
    saveas(fig, png_file);
    savefig(fig, fig_file);
end

function add_figure_title(fig, ttl)
    % 缁欐暣寮犲浘鍔犱腑鏂囨€绘爣棰橈紝骞朵繚鐣欏瓙鍥炬爣棰樸€?    figure(fig);
    try
        sgtitle(ttl, 'Interpreter', 'none');
    catch
        annotation(fig, 'textbox', [0.02 0.965 0.96 0.03], ...
            'String', ttl, 'Interpreter', 'none', ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none');
    end
end

function [elem_x, elem_y] = element_centers(md, elem_mask)
    ids = find(elem_mask(:));
    if isempty(ids)
        elem_x = [];
        elem_y = [];
        return;
    end
    elem_x = mean(md.mesh.x(md.mesh.elements(ids,:)), 2);
    elem_y = mean(md.mesh.y(md.mesh.elements(ids,:)), 2);
end

function grad_mag = element_gradient_magnitude(md, node_data)
    % 瀵逛笁瑙掑舰鍗曞厓璁＄畻鑺傜偣鍦虹殑涓€闃舵搴﹀ぇ灏忥紝鍗曚綅鍙栧喅浜庤緭鍏ュ満銆?    if isempty(node_data) || numel(node_data) ~= md.mesh.numberofvertices
        grad_mag = [];
        return;
    end

    f = double(node_data(:));
    elems = md.mesh.elements;
    x = md.mesh.x(:);
    y = md.mesh.y(:);

    n1 = elems(:,1);
    n2 = elems(:,2);
    n3 = elems(:,3);

    x1 = x(n1); x2 = x(n2); x3 = x(n3);
    y1 = y(n1); y2 = y(n2); y3 = y(n3);
    f1 = f(n1); f2 = f(n2); f3 = f(n3);

    den = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);
    dfdx = ((f2 - f1) .* (y3 - y1) - (f3 - f1) .* (y2 - y1)) ./ den;
    dfdy = ((x2 - x1) .* (f3 - f1) - (x3 - x1) .* (f2 - f1)) ./ den;

    grad_mag = sqrt(dfdx.^2 + dfdy.^2);
    grad_mag(~isfinite(grad_mag)) = NaN;
end

function dmin = distance_to_mask_nodes(md, mask)
    % 閫愬潡璁＄畻姣忎釜鑺傜偣鍒版渶杩戝浐瀹氭按澶磋妭鐐圭殑璺濈锛岄伩鍏嶄竴娆℃€х敓鎴愯繃澶х殑璺濈鐭╅樀銆?    nv = md.mesh.numberofvertices;
    ids = find(mask(:));
    dmin = NaN(nv, 1);
    if isempty(ids)
        return;
    end

    x = md.mesh.x(:);
    y = md.mesh.y(:);
    xb = x(ids)';
    yb = y(ids)';
    chunk = 2000;
    for i1 = 1:chunk:nv
        i2 = min(nv, i1 + chunk - 1);
        dx = bsxfun(@minus, x(i1:i2), xb);
        dy = bsxfun(@minus, y(i1:i2), yb);
        dmin(i1:i2) = sqrt(min(dx.^2 + dy.^2, [], 2));
    end
end

function patchplot(md, data, ttl, color_limit_mode)
    if isempty(data)
        axis off;
        title([ttl ' unavailable'], 'Interpreter', 'none');
        return;
    end

    d = double(data(:));
    nv = md.mesh.numberofvertices;
    ne = md.mesh.numberofelements;

    if numel(d) == nv
        face_color = 'interp';
        cdata = d;
    elseif numel(d) == ne
        face_color = 'flat';
        cdata = d;
    else
        axis off;
        title([ttl ' wrong length'], 'Interpreter', 'none');
        return;
    end

    patch('Faces', md.mesh.elements, ...
          'Vertices', [md.mesh.x, md.mesh.y], ...
          'FaceVertexCData', cdata, ...
          'FaceColor', face_color, ...
          'EdgeColor', 'none');
    axis equal tight; box on;
    valid = d(isfinite(d));
    if ~isempty(valid)
        ttl = sprintf('%s\nmin %.3g, max %.3g', ttl, min(valid), max(valid));
    end
    title(ttl, 'Interpreter', 'none');
    xlabel('X (m)');
    ylabel('Y (m)');
    colorbar;

    if ~isempty(valid) && numel(unique(valid)) > 1
        lo = min(valid);
        hi = max(valid);
        if isfinite(lo) && isfinite(hi) && hi > lo
            clim([lo hi]);
        end
    end
end

function patchplot_with_bad(md, data, bad_node_ids, elem_x, elem_y, ttl, color_limit_mode)
    patchplot(md, data, ttl, color_limit_mode);
    hold on;

    if ~isempty(bad_node_ids)
        plot(md.mesh.x(bad_node_ids), md.mesh.y(bad_node_ids), ...
            'ro', 'MarkerSize', 4, 'LineWidth', 1);
    end

    if ~isempty(elem_x)
        plot(elem_x, elem_y, 'ko', 'MarkerSize', 3, 'LineWidth', 0.8);
    end
end
