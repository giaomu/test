script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(project_root);
paths = shaktiais_paths();

S = load(fullfile(paths.models, 'Recovery_Hydrology.mat'), 'md');
md = S.md;
clear S;

%% ========== 基本信息（纯文本，原样保留） ==========
md   % 整体结构
disp(md.miscellaneous.name);
size(md.mesh.x), md.mesh.numberofvertices, md.mesh.numberofelements

areas = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);
edge_lengths = sqrt(areas*2);
fprintf('边长: min=%.1f, median=%.1f, max=%.1f m\n', ...
    min(edge_lengths), median(edge_lengths), max(edge_lengths));

fprintf('bed:       min=%.1f, max=%.1f, mean=%.1f\n', min(md.geometry.bed), max(md.geometry.bed), mean(md.geometry.bed));
fprintf('surface:   min=%.1f, max=%.1f\n', min(md.geometry.surface), max(md.geometry.surface));
fprintf('thickness: min=%.1f, max=%.1f\n', min(md.geometry.thickness), max(md.geometry.thickness));
fprintf('base:      min=%.1f, max=%.1f\n', min(md.geometry.base), max(md.geometry.base));
resid = md.geometry.thickness - (md.geometry.surface - md.geometry.base);
fprintf('thickness - (surface - base) 最大偏�? %e\n', max(abs(resid)));

n_ice    = sum(md.mask.ice_levelset < 0);
n_no_ice = sum(md.mask.ice_levelset > 0);
fprintf('有冰节点: %d, 无冰节点: %d\n', n_ice, n_no_ice);

fprintf('rheology_n: min=%.1f, max=%.1f\n', min(md.materials.rheology_n), max(md.materials.rheology_n));
fprintf('rheology_B: min=%.3e, max=%.3e\n', min(md.materials.rheology_B), max(md.materials.rheology_B));
fprintf('SSA 开�? %d\n', all(md.flowequation.element_equation==2));

fprintf('friction.coefficient: min=%.2f, max=%.2f, std=%.2f\n', ...
    min(md.friction.coefficient), max(md.friction.coefficient), std(md.friction.coefficient));
assert(min(md.friction.coefficient) >= 0);
assert(max(md.friction.coefficient) <= 4000);
fprintf('iscontrol (应为 0): %d\n', md.inversion.iscontrol);

vel_mod = sqrt(md.initialization.vx.^2 + md.initialization.vy.^2);
vel_obs = md.inversion.vel_obs;
pos = find(md.mask.ice_levelset<0 & vel_obs>10);
rel_err = abs(vel_mod(pos)-vel_obs(pos))./vel_obs(pos);
fprintf('速度相对误差 中位�?%.1f%%, 90%%分位=%.1f%%\n', ...
    median(rel_err)*100, prctile(rel_err,90)*100);


%% ========== Figure 1: 几何与网�?(2x2) ==========
figure('Name','几何与网�?,'Color','w','Position',[50 50 1200 900]);

subplot(2,2,1);
patch('Faces', md.mesh.elements, 'Vertices', [md.mesh.x, md.mesh.y], ...
      'FaceColor', 'none', 'EdgeColor', [0.3 0.3 0.3]);
axis equal tight; box on;
title(sprintf('Mesh (%d nodes, %d elements)', ...
      md.mesh.numberofvertices, md.mesh.numberofelements));
xlabel('X (m)'); ylabel('Y (m)');

subplot(2,2,2);
patchplot(md, md.geometry.bed, 'Bed [m]', 'cmap', 'turbo');

subplot(2,2,3);
patchplot(md, md.geometry.surface, 'Surface [m]', 'cmap', 'turbo');

subplot(2,2,4);
patchplot(md, md.geometry.thickness, 'Thickness [m]', 'cmap', 'turbo');


%% ========== Figure 2: Mask & 流动方程 (1x2) ==========
figure('Name','Mask','Color','w','Position',[50 50 1200 500]);

subplot(1,2,1);
patchplot(md, md.mask.ice_levelset, 'Ice Levelset (<0=ice)', 'cmap', 'turbo');

subplot(1,2,2);
if isprop(md.mask, 'ocean_levelset')
    patchplot(md, md.mask.ocean_levelset, 'Ocean Levelset (<0=ocean)', 'cmap', 'turbo');
end


%% ========== Figure 3: 反演结果 (2x2) ==========
figure('Name','反演结果','Color','w','Position',[50 50 1300 1000]);

subplot(2,2,1);
patchplot(md, md.friction.coefficient, 'Inverted Friction', ...
    'clim', [0 700], 'cmap', 'turbo');

subplot(2,2,2);
patchplot(md, vel_mod, 'Modeled Vel [m/yr]', ...
    'log', true, 'cmap', 'parula', 'cbar_label', 'log_{10} m/yr');

subplot(2,2,3);
patchplot(md, vel_obs, 'Observed Vel [m/yr]', ...
    'log', true, 'cmap', 'parula', 'cbar_label', 'log_{10} m/yr');

subplot(2,2,4);
patchplot(md, log10(abs(vel_mod-vel_obs)+1), 'log_{10}|Residual|', 'cmap', 'hot');


%% ========== 水文字段检查（文本�?==========
fprintf('\nHydrology class: %s\n', class(md.hydrology));
fprintf('head:         min=%.1f, max=%.1f, �?NaN? %d\n', ...
    min(md.hydrology.head), max(md.hydrology.head), any(isnan(md.hydrology.head)));
fprintf('gap_height:   min=%.4f, max=%.4f\n', ...
    min(md.hydrology.gap_height), max(md.hydrology.gap_height));
fprintf('bump_spacing: min=%.2f\n', min(md.hydrology.bump_spacing));
fprintf('bump_height:  min=%.2f, max=%.2f\n', ...
    min(md.hydrology.bump_height), max(md.hydrology.bump_height));
fprintf('reynolds:     min=%.0f, max=%.0f\n', ...
    min(md.hydrology.reynolds), max(md.hydrology.reynolds));
fprintf('englacial_input 总量: %.2e\n', sum(md.hydrology.englacial_input));
fprintf('moulin_input    总量: %.2e\n', sum(md.hydrology.moulin_input));
fprintf('neumannflux     总量: %.2e\n', sum(md.hydrology.neumannflux));
n_outlets = sum(~isnan(md.hydrology.spchead));
fprintf('Dirichlet 出口: %d\n', n_outlets);

% 有效压力
rho_i = md.materials.rho_ice;
rho_w = md.materials.rho_water;
g = md.constants.g;
Neff = rho_i*g*md.geometry.thickness - rho_w*g*(md.hydrology.head - md.geometry.base);
fprintf('有效压力 N: min=%.2e, max=%.2e Pa, N<0 节点=%d\n', ...
    min(Neff), max(Neff), sum(Neff<0));


%% ========== Figure 4: SHAKTI 初始�?(2x2) ==========
figure('Name','SHAKTI 初始�?,'Color','w','Position',[50 50 1300 1000]);

subplot(2,2,1);
patchplot(md, md.hydrology.head, 'Initial Head [m]', 'cmap', 'turbo');

subplot(2,2,2);
patchplot(md, md.hydrology.gap_height, 'Initial Gap Height [m]', 'cmap', 'turbo');

subplot(2,2,3);
outlet_mask = double(~isnan(md.hydrology.spchead));
patchplot(md, outlet_mask, sprintf('Outlets (%d nodes)', n_outlets), ...
    'clim', [0 1], 'cmap', 'gray');

subplot(2,2,4);
patchplot(md, Neff/1e6, 'Initial N [MPa]', 'clim', [-1 10], 'cmap', 'turbo');


%% ========== 耦合开�?& 时间步（文本�?==========
fprintf('\ntransient.ishydrology    : %d\n', md.transient.ishydrology);
fprintf('transient.isstressbalance: %d\n', md.transient.isstressbalance);
fprintf('transient.ismasstransport: %d\n', md.transient.ismasstransport);
fprintf('transient.isthermal      : %d\n', md.transient.isthermal);
fprintf('friction.coupling        : %d\n', md.friction.coupling);
fprintf('inversion.iscontrol      : %d\n', md.inversion.iscontrol);
fprintf('time_step      : %.4f yr = %.2f hr\n', md.timestepping.time_step, md.timestepping.time_step*8760);
fprintf('final_time     : %.4f yr = %.1f days\n', md.timestepping.final_time, md.timestepping.final_time*365);
fprintf('output_frequency: %d\n', md.settings.output_frequency);

% relaxation / storage
fprintf('\nrelaxation: %.3f\n', md.hydrology.relaxation);
fprintf('storage (e_v): %.2e\n', md.hydrology.storage);
if md.hydrology.storage == 0
    fprintf('⚠️ storage=0，PDE 退化为椭圆型\n');
end

% 基底融化水源
fprintf('\ngeothermal flux:  min=%.3f, max=%.3f, mean=%.3f W/m²\n', ...
    min(md.basalforcings.geothermalflux), max(md.basalforcings.geothermalflux), mean(md.basalforcings.geothermalflux));
fprintf('grounded melt:    min=%.2e, max=%.2e (m/yr)\n', ...
    min(md.basalforcings.groundedice_melting_rate), max(md.basalforcings.groundedice_melting_rate));


function patchplot(md, data, ttl, varargin)
% patchplot - �?patch �?ISSM 模型场（顶点场或单元场都行）
% 用法: patchplot(md, data, 'Title')
%       patchplot(md, data, 'Title', 'clim', [0 100])
%       patchplot(md, data, 'Title', 'log', true)
%       patchplot(md, data, 'Title', 'cmap', 'parula')

p = inputParser;
addParameter(p, 'clim', []);
addParameter(p, 'log', false);
addParameter(p, 'cmap', 'parula');
addParameter(p, 'cbar_label', '');
parse(p, varargin{:});

% 数据预处�?
d = double(data(:));
if p.Results.log
    d(d<=0) = NaN;
    d = log10(d);
end

% 判断顶点�?or 单元�?
nv = md.mesh.numberofvertices;
ne = md.mesh.numberofelements;
if numel(d) == nv
    face_color = 'interp';  % 顶点场：线性插�?
    cdata = d;
elseif numel(d) == ne
    face_color = 'flat';    % 单元场：每个三角形一个颜�?
    cdata = d;
else
    error('data 长度既不等于顶点数也不等于单元数');
end

% �?patch
patch('Faces', md.mesh.elements, ...
      'Vertices', [md.mesh.x, md.mesh.y], ...
      'FaceVertexCData', cdata, ...
      'FaceColor', face_color, ...
      'EdgeColor', 'none');

axis equal tight; box on;
title(ttl, 'FontSize', 11, 'FontWeight', 'bold');
xlabel('X (m)'); ylabel('Y (m)');

% colorbar + 范围
cb = colorbar;
if ~isempty(p.Results.cbar_label)
    cb.Label.String = p.Results.cbar_label;
end
if ~isempty(p.Results.clim)
    clim(p.Results.clim);
else
    % 自动去除异常�?
    valid = d(~isnan(d) & isfinite(d));
    if ~isempty(valid)
        lo = prctile(valid, 1);
        hi = prctile(valid, 99);
        if lo < hi, clim([lo hi]); end
    end
end
colormap(gca, p.Results.cmap);
end
