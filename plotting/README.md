# Plotting Scripts

This folder is organized by workflow. Run scripts from the project root with
`run(...)`; each script resolves the project root from its own location.

## New Mesh

Current `runme_newmesh.m` workflow diagnostics:

```matlab
run('plotting/newmesh/check_newmesh_parameterize.m')
run('plotting/newmesh/check_newmesh_inversion.m')
run('plotting/newmesh/compare_geothermal_newmesh.m')
run('plotting/newmesh/check_newmesh_shakti_timeseries.m')
```

Contents:

- `check_newmesh_parameterize.m`: checks Step 2 parameterization fields.
- `check_newmesh_inversion.m`: checks Step 3 inversion fields and misfit.
- `compare_geothermal_newmesh.m`: compares old-model geothermal flux with the new-mesh interpolated field.
- `check_newmesh_shakti_timeseries.m`: checks every saved SHAKTI transient time step.

## Legacy Recovery

Older Recovery mesh and earlier SHAKTI diagnostics:

```matlab
run('plotting/legacy_recovery/check_mat.m')
run('plotting/legacy_recovery/check_shakti_onestep.m')
run('plotting/legacy_recovery/check_bad_nodes_mesh_size.m')
run('plotting/legacy_recovery/check_flake_mesh_regression.m')
run('plotting/legacy_recovery/plot_parameterize_fields.m')
run('plotting/legacy_recovery/plot_model_general.m')
run('plotting/legacy_recovery/plot_recovery_hydrology_bed_lakes.m')
```

## Inversion L-Curve

L-curve and inversion-comparison plotting:

```matlab
run('plotting/inversion_lcurve/plot_lcurve_no_corner.m')
run('plotting/inversion_lcurve/plot_inversion.m')
run('plotting/inversion_lcurve/fig44_inversion_spatial_compare.m')
```

## Docs

- `docs/check_figure_reading_guide.md`: figure-reading notes for older diagnostic outputs.

Project paths are centralized in `../shaktiais_paths.m`. Model outputs are
read from `../outputs/models` and `../outputs/models_newmesh`; generated figures
are saved under `../outputs/figures` and `../outputs/figures_newmesh`.
