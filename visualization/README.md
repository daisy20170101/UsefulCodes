# Fault Visualization Scripts

ParaView Python (pvpython) scripts for generating snapshots of fault slip rate data.

## Overview

This directory contains scripts to automatically generate visualizations of fault outputs, specifically:
- **SRs**: Slip Rate Strike component (m/s)
- **SRd**: Slip Rate Dip component (m/s)

The scripts use ParaView's Python interface (pvpython) and can load data directly from **XDMF files** or from ParaView state files.

## Quick Start

```bash
# Generate snapshots from XDMF file (single frame)
pvpython generate_fault_snapshots.py -i fault_output.xdmf

# Generate snapshots for ALL timesteps ⭐ NEW!
pvpython generate_fault_snapshots.py -i fault_output.xdmf --all-timesteps

# Or use the wrapper script
./run_snapshots.sh -i fault_output.xdmf
```

## Prerequisites

- **ParaView** (with pvpython)
- **Input data**: XDMF file with SRs and SRd fields, or
- **State file** (optional): `rakeNew.pvsm` (place in this folder or use original path)

### State File Setup

Copy the ParaView state file to this folder:
```bash
cp /Users/DuoL/Documents/NSHM/Central/Paraviews/rakeNew.pvsm ./visualization/
```

The state file is **optional** when using XDMF input - it's only used for camera/view settings.

### Installing ParaView

Download from: https://www.paraview.org/download/

After installation, `pvpython` should be available in your PATH.

## Scripts

### 1. `generate_fault_snapshots.py`

Generate snapshots of SRs and SRd from XDMF files or ParaView state files.

**Basic Usage:**
```bash
# Load from XDMF file (recommended)
pvpython generate_fault_snapshots.py -i fault_output.xdmf

# Process ALL timesteps in XDMF file ⭐ NEW!
pvpython generate_fault_snapshots.py -i fault_output.xdmf --all-timesteps

# Load from XDMF with custom state file for camera/view settings
pvpython generate_fault_snapshots.py -i fault_output.xdmf -s rakeNew.pvsm

# Process all timesteps with custom state file
pvpython generate_fault_snapshots.py -i fault_output.xdmf -s rakeNew.pvsm --all-timesteps

# Legacy: Load from state file only
pvpython generate_fault_snapshots.py -s /path/to/state.pvsm

# Custom output directory
pvpython generate_fault_snapshots.py -i data.xdmf -o ./my_snapshots

# High resolution output (4K) for all timesteps
pvpython generate_fault_snapshots.py -i data.xdmf -r 3840 2160 --all-timesteps

# Generate only SRs snapshots for all timesteps
pvpython generate_fault_snapshots.py -i data.xdmf --field SRs --all-timesteps

# Generate only SRd snapshot (single frame)
pvpython generate_fault_snapshots.py -i data.xdmf --field SRd
```

**Options:**
- `-i, --input`: Path to XDMF data file (**recommended**)
- `-s, --state-file`: Path to ParaView state file (optional, for camera/view settings)
- `-o, --output-dir`: Output directory for snapshots (default: `./snapshots`)
- `-r, --resolution WIDTH HEIGHT`: Image resolution (default: 1920 1080)
- `-f, --field {SRs,SRd,all}`: Field to visualize (default: all)
- `--output-name`: Output filename prefix (default: fault)
- `--all-timesteps`: **Process all timesteps in data file** ⭐ NEW!

**Output (single frame):**
- `fault_SRs.png`: Snapshot of strike-slip rate component
- `fault_SRd.png`: Snapshot of dip-slip rate component

**Output (with --all-timesteps):**
- `fault_t00000_SRs.png`, `fault_t00000_SRd.png`: Timestep 0
- `fault_t00001_SRs.png`, `fault_t00001_SRd.png`: Timestep 1
- ... (one pair per timestep)

### 2. `batch_fault_snapshots.py`

Batch process multiple data files or timesteps to generate snapshots.

**Usage:**

**Process single XDMF file (all timesteps):**
```bash
# Process all timesteps in a single XDMF file
pvpython batch_fault_snapshots.py -i fault_output.xdmf -o ./snapshots
```

**Process data directory:**
```bash
# Process all XDMF files in directory
pvpython batch_fault_snapshots.py -d /path/to/data -o ./snapshots

# Process specific file pattern
pvpython batch_fault_snapshots.py -d /path/to/data -p "fault_*.vtk"

# Process VTU files
pvpython batch_fault_snapshots.py -d /path/to/data -p "*.vtu"
```

**Process timesteps from state file:**
```bash
# Process all timesteps in state file
pvpython batch_fault_snapshots.py --state-mode -s /path/to/state.pvsm

# With custom resolution
pvpython batch_fault_snapshots.py --state-mode -r 2560 1440
```

**Options:**
- `-i, --input`: Path to single XDMF file (process all timesteps)
- `-d, --data-dir`: Directory containing data files
- `-p, --pattern`: File pattern to match (default: `*.xdmf`)
- `-o, --output-dir`: Output directory (default: `./snapshots_batch`)
- `-r, --resolution WIDTH HEIGHT`: Image resolution (default: 1920 1080)
- `-s, --state-file`: ParaView state file (optional, for camera/view settings)
- `--state-mode`: Process timesteps from state file instead of data files

**Supported file formats:**
- XDMF (`.xdmf`)
- VTK (`.vtk`)
- VTU (`.vtu`)
- PVD (`.pvd`)

**Output:**
For each input file or timestep:
- `<filename>_SRs.png`: Strike-slip rate snapshot
- `<filename>_SRd.png`: Dip-slip rate snapshot

## Examples

### Example 1: Single Snapshot from XDMF File

Generate SRs and SRd snapshots from an XDMF file (recommended workflow):

```bash
cd visualization
pvpython generate_fault_snapshots.py -i /path/to/fault_output.xdmf
```

Output:
```
./snapshots/fault_SRs.png
./snapshots/fault_SRd.png
```

### Example 2: Single Snapshot with State File

Generate snapshots using both XDMF data and state file for camera settings:

```bash
pvpython generate_fault_snapshots.py -i fault_output.xdmf -s rakeNew.pvsm
```

### Example 3: High-Resolution Snapshot

Generate 4K resolution snapshots from XDMF:

```bash
pvpython generate_fault_snapshots.py \
    -i fault_output.xdmf \
    -o ./snapshots_4k \
    -r 3840 2160
```

### Example 4: Process All Timesteps ⭐ NEW!

Generate snapshots for **all timesteps** in an XDMF file:

```bash
# Process all timesteps (recommended method)
pvpython generate_fault_snapshots.py \
    -i fault_output.xdmf \
    --all-timesteps \
    -o ./snapshots_all

# Or with custom resolution
pvpython generate_fault_snapshots.py \
    -i fault_output.xdmf \
    --all-timesteps \
    -r 2560 1440 \
    -o ./snapshots_all_hd
```

This generates numbered snapshots for each timestep:
```
./snapshots_all/fault_t00000_SRs.png
./snapshots_all/fault_t00000_SRd.png
./snapshots_all/fault_t00001_SRs.png
./snapshots_all/fault_t00001_SRd.png
./snapshots_all/fault_t00002_SRs.png
./snapshots_all/fault_t00002_SRd.png
...
```

**Note:** This is the **easiest way** to process all timesteps with `generate_fault_snapshots.py`!

### Example 5: Batch Process Single XDMF File (Alternative Method)

Process all timesteps in a single XDMF file:

```bash
pvpython batch_fault_snapshots.py \
    -i /path/to/fault_output.xdmf \
    -o ./snapshots_batch
```

### Example 5: Batch Process Data Directory

Process all XDMF files in a directory:

```bash
pvpython batch_fault_snapshots.py \
    -d /path/to/simulation/output \
    -p "fault_*.xdmf" \
    -o ./snapshots_batch \
    -r 1920 1080
```

This will generate:
```
./snapshots_batch/fault_001_SRs.png
./snapshots_batch/fault_001_SRd.png
./snapshots_batch/fault_002_SRs.png
./snapshots_batch/fault_002_SRd.png
...
```

### Example 6: Using the Shell Wrapper

The `run_snapshots.sh` script provides an easy-to-use interface:

```bash
# Generate snapshot from XDMF
./run_snapshots.sh -i fault_output.xdmf

# Batch process XDMF with high resolution
./run_snapshots.sh -i fault_output.xdmf -r 3840x2160

# Process directory
./run_snapshots.sh -m batch -d /path/to/data
```

## Workflow Integration

### Creating a State File

1. Open ParaView
2. Load your fault data
3. Set up visualization (color maps, camera angle, etc.)
4. Save state: **File → Save State** → `rakeNew.pvsm`

### Customizing Visualizations

To modify the visualization settings:

1. Edit the state file in ParaView
2. Adjust color maps, camera angles, filters, etc.
3. Save the updated state file
4. Run the scripts with the updated state file

### Automation

For automated processing in a pipeline:

```bash
#!/bin/bash
# Process multiple simulations

for sim_dir in /path/to/simulations/*/; do
    sim_name=$(basename "$sim_dir")
    echo "Processing $sim_name..."

    pvpython batch_fault_snapshots.py \
        -d "$sim_dir" \
        -o "./snapshots_$sim_name" \
        -r 1920 1080
done
```

## Troubleshooting

### pvpython not found

Make sure ParaView is installed and pvpython is in your PATH:

```bash
# Find pvpython location
find /Applications -name pvpython 2>/dev/null  # macOS
which pvpython  # Linux

# Add to PATH (example for macOS)
export PATH="/Applications/ParaView.app/Contents/bin:$PATH"
```

### State file not found

Update the default state file path in the script or specify with `-s`:

```bash
pvpython generate_fault_snapshots.py -s /your/path/to/rakeNew.pvsm
```

### Missing fields (SRs or SRd)

Ensure your data files contain the `SRs` and `SRd` fields. Check available fields:

```python
# In ParaView Python shell
source = GetActiveSource()
print(source.PointData[...])
```

### Memory issues with large datasets

For large datasets, reduce resolution or process files individually:

```bash
# Lower resolution
pvpython batch_fault_snapshots.py -d /path/to/data -r 1280 720

# Process one file at a time
for file in /path/to/data/*.xdmf; do
    pvpython generate_fault_snapshots.py --data-file "$file"
done
```

## Output Format

All snapshots are saved as PNG images with:
- Transparent background option (disabled by default)
- Color bar showing field values
- Configurable resolution
- Field name in color bar title

## Notes

- The scripts preserve the camera angle and view settings from the state file
- Color maps are automatically rescaled to data range for each field
- Both scripts support multiple output resolutions
- Batch processing reports progress and errors for each file

## References

- ParaView Documentation: https://www.paraview.org/documentation/
- ParaView Python Scripting Guide: https://www.paraview.org/paraview-guide/
- pvpython: https://docs.paraview.org/en/latest/ReferenceManual/pythonScripting.html

## Author

Created for NSHM Central fault visualization workflow.

Last updated: 2025-11-05
