function paths = shaktiais_paths()
%SHAKTIAIS_PATHS Centralized project paths for the SHAKTI AIS workflow.

project_root = fileparts(mfilename('fullpath'));

paths.root = project_root;
paths.issm_root = 'D:\Auniversity\ISSM\ISSMmodel\My_ISSM\ISSM-Windows-MATLAB';
paths.issm_data = fullfile(paths.issm_root, 'examples', 'Data');

paths.inputs = fullfile(project_root, 'inputs');
paths.exp = resolve_folder(fullfile(paths.inputs, 'exp'), project_root, '*.exp');
paths.par = resolve_folder(fullfile(paths.inputs, 'par'), project_root, '*.par');
paths.data = resolve_folder(fullfile(project_root, 'Data'), paths.issm_data, '*.nc');
paths.velocity_file = resolve_file(fullfile(paths.data, 'antarctica_ice_velocity_450m_v2.nc'), ...
    fullfile(paths.issm_data, 'antarctica_ice_velocity_450m_v2.nc'));
paths.newmesh_velocity_file = ...
    'E:\Earthdata\download\NSIDC-0761_1-20260507_071611\antarctica_ice_velocity_2014-2017_450m_v01.1.nc';
paths.bedmachine_file = resolve_file(fullfile(paths.data, 'NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc'), ...
    fullfile(paths.issm_data, 'NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc'));

paths.outputs = fullfile(project_root, 'outputs');
paths.figures = fullfile(paths.outputs, 'figures');
paths.lcurve_figures = fullfile(paths.figures, 'Lcurve');
paths.models = resolve_folder(fullfile(paths.outputs, 'models'), fullfile(project_root, 'Models'), '*.mat');
paths.models_newmesh = fullfile(paths.outputs, 'models_newmesh');
paths.figures_newmesh = fullfile(paths.outputs, 'figures_newmesh');

paths.archive = fullfile(project_root, 'archive');
end

function folder = resolve_folder(clean_folder, legacy_folder, file_pattern)
if folder_has_files(clean_folder, file_pattern) || ~isfolder(legacy_folder)
    folder = clean_folder;
else
    folder = legacy_folder;
end
end

function tf = folder_has_files(folder, file_pattern)
if ~isfolder(folder)
    tf = false;
    return;
end
matches = dir(fullfile(folder, file_pattern));
tf = ~isempty(matches);
end

function file = resolve_file(primary_file, fallback_file)
if isfile(primary_file) || ~isfile(fallback_file)
    file = primary_file;
else
    file = fallback_file;
end
end
