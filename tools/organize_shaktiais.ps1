$ErrorActionPreference = 'Stop'

$expectedLeaf = 'SHAKTI_AIS'
$root = (Get-Location).Path
if ((Split-Path -Leaf $root) -ne $expectedLeaf) {
    throw "Run this script from the SHAKTI_AIS project root. Current path: $root"
}

function Move-ChildItemsMerge {
    param(
        [Parameter(Mandatory = $true)]
        [string] $SourceDir,

        [Parameter(Mandatory = $true)]
        [string] $DestinationDir
    )

    if (-not (Test-Path -LiteralPath $SourceDir)) {
        return
    }
    if (-not (Test-Path -LiteralPath $DestinationDir)) {
        New-Item -ItemType Directory -Path $DestinationDir | Out-Null
    }

    Get-ChildItem -LiteralPath $SourceDir -Force | ForEach-Object {
        $target = Join-Path $DestinationDir $_.Name
        if ($_.PSIsContainer -and (Test-Path -LiteralPath $target)) {
            Move-ChildItemsMerge -SourceDir $_.FullName -DestinationDir $target
            if ((Get-ChildItem -LiteralPath $_.FullName -Force | Measure-Object).Count -eq 0) {
                Remove-Item -LiteralPath $_.FullName
            }
        }
        else {
            Move-Item -LiteralPath $_.FullName -Destination $DestinationDir -Force
        }
    }
}

$dirs = @(
    'inputs',
    'inputs\exp',
    'inputs\par',
    'outputs',
    'outputs\models',
    'outputs\figures',
    'archive'
)
foreach ($dir in $dirs) {
    if (-not (Test-Path -LiteralPath $dir)) {
        New-Item -ItemType Directory -Path $dir | Out-Null
    }
}

$expFiles = @(
    'outlet.exp',
    'Recovery.exp',
    'Recovery_all_Noice.exp',
    'Recovery_outlets.exp'
)
foreach ($file in $expFiles) {
    if (Test-Path -LiteralPath $file) {
        Move-Item -LiteralPath $file -Destination 'inputs\exp' -Force
    }
}

$parFiles = @(
    'Recovery.par',
    'Recovery_reset_bed.par'
)
foreach ($file in $parFiles) {
    if (Test-Path -LiteralPath $file) {
        Move-Item -LiteralPath $file -Destination 'inputs\par' -Force
    }
}

if (Test-Path -LiteralPath 'Models') {
    Move-ChildItemsMerge -SourceDir 'Models' -DestinationDir 'outputs\models'
    if ((Get-ChildItem -LiteralPath 'Models' -Force | Measure-Object).Count -eq 0) {
        Remove-Item -LiteralPath 'Models'
    }
}

if (Test-Path -LiteralPath 'Lcurve') {
    if (-not (Test-Path -LiteralPath 'outputs\figures\Lcurve')) {
        New-Item -ItemType Directory -Path 'outputs\figures\Lcurve' | Out-Null
    }
    Move-ChildItemsMerge -SourceDir 'Lcurve' -DestinationDir 'outputs\figures\Lcurve'
    if ((Get-ChildItem -LiteralPath 'Lcurve' -Force | Measure-Object).Count -eq 0) {
        Remove-Item -LiteralPath 'Lcurve'
    }
}

$figureDirs = @(
    'Fig44_inversion_spatial_compare',
    'Fig44_inversion_spatial_compare_fixed',
    'L_curve_models_20000_10_1p000e-05_inversion',
    'L_curve_models_40000_10_1p500e-05_inversion',
    'L_curve_models_50000_10_8p500e-05_inversion',
    'L_curve_models_5000_10_6p444e-06_inversion',
    'L_curve_models_500_10_1p800e-06_inversion',
    'L_curve_models_70000_10_5p500e-05_inversion',
    'Recovery_Inversion_generalplot',
    'Recovery_Parameterize_physics',
    'Recovery_SHAKTI_5Step_300s_diagnostics',
    'Recovery_SHAKTI_5Step_300s_flake_mesh_regression'
)
foreach ($dir in $figureDirs) {
    if (Test-Path -LiteralPath $dir) {
        $dest = Join-Path 'outputs\figures' $dir
        if (-not (Test-Path -LiteralPath $dest)) {
            New-Item -ItemType Directory -Path $dest | Out-Null
        }
        Move-ChildItemsMerge -SourceDir $dir -DestinationDir $dest
        if ((Get-ChildItem -LiteralPath $dir -Force | Measure-Object).Count -eq 0) {
            Remove-Item -LiteralPath $dir
        }
    }
}

$archiveDirs = @(
    'noice_dont_invert',
    'old_geometry',
    '.matlab_prefs_check'
)
foreach ($dir in $archiveDirs) {
    if (Test-Path -LiteralPath $dir) {
        $dest = Join-Path 'archive' $dir
        if (-not (Test-Path -LiteralPath $dest)) {
            New-Item -ItemType Directory -Path $dest | Out-Null
        }
        Move-ChildItemsMerge -SourceDir $dir -DestinationDir $dest
        if ((Get-ChildItem -LiteralPath $dir -Force | Measure-Object).Count -eq 0) {
            Remove-Item -LiteralPath $dir
        }
    }
}

$modelArchiveRoot = 'archive\model_archives'
if (-not (Test-Path -LiteralPath $modelArchiveRoot)) {
    New-Item -ItemType Directory -Path $modelArchiveRoot | Out-Null
}
$modelArchiveDirs = @(
    'noice_dont_invert',
    'old_geometry',
    'old_models'
)
foreach ($dir in $modelArchiveDirs) {
    $source = Join-Path 'outputs\models' $dir
    if (Test-Path -LiteralPath $source) {
        $dest = Join-Path $modelArchiveRoot $dir
        Move-ChildItemsMerge -SourceDir $source -DestinationDir $dest
        if ((Get-ChildItem -LiteralPath $source -Force | Measure-Object).Count -eq 0) {
            Remove-Item -LiteralPath $source
        }
    }
}

if (Test-Path -LiteralPath 'plotting\Recovery_Hydrology_bed_lakes.fig') {
    $hydroFigDir = 'outputs\figures\Recovery_Hydrology_bed_lakes'
    if (-not (Test-Path -LiteralPath $hydroFigDir)) {
        New-Item -ItemType Directory -Path $hydroFigDir | Out-Null
    }
    Move-Item -LiteralPath 'plotting\Recovery_Hydrology_bed_lakes.fig' -Destination $hydroFigDir -Force
}

Write-Host 'SHAKTI_AIS folder organization complete.'
