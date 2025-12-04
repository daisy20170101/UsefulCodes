# Segmented Fault Meshing Workflow

This workflow separates fault geometry from a `.geo` file into individual segments, meshes each segment independently, and then merges them into a final 3D mesh.

## Overview

The workflow consists of three Python scripts:

1. **`mesh_fault_segments.py`** - Extracts and meshes individual fault segments from `.geo` file
2. **`mesh_ando_from_segments.py`** - Merges pre-meshed segments into 3D mesh with box domain
3. **`run_segmented_meshing.py`** - Convenience script that runs both steps

## Why Use This Approach?

**Advantages:**
- **Independent control**: Each fault segment can be meshed with different parameters
- **Debugging**: Easier to identify problematic segments
- **Reusability**: Individual segment meshes can be reused in different configurations
- **Parallel processing**: Segments can be meshed in parallel (future enhancement)
- **Inspection**: Each segment can be visualized independently before merging

## Quick Start

### Option 1: Run Complete Workflow (Recommended)

```bash
cd /path/to/your/geo/file

python /path/to/MeshPy/run_segmented_meshing.py kaikoura_v11.3hpk2.geo \
    --segment-dir fault_segments \
    --lc-fault 2000 \
    --lc-coarse 10000 \
    --output ando_mesh_merged.msh
```

### Option 2: Run Steps Separately

**Step 1: Extract and mesh individual segments**
```python
from mesh_fault_segments import mesh_individual_fault_segments

segment_files, segment_info = mesh_individual_fault_segments(
    geo_file='kaikoura_v11.3hpk2.geo',
    output_dir='fault_segments',
    lc=2000,  # Mesh size on fault (meters)
    remove_criteria={
        'z_near_zero_tol': 1000.0,  # Remove surfaces with |z| < 1000 m
        'z_below': -14000.0          # Remove surfaces with z < -14000 m
    }
)
```

**Step 2: Merge segments into 3D mesh**
```python
from mesh_ando_from_segments import generate_mesh_from_fault_segments

final_mesh = generate_mesh_from_fault_segments(
    segment_files=segment_files,
    lc_coarse=10000,           # Coarse mesh size (meters)
    box_length=100e3,          # Box padding X (meters)
    box_width=100e3,           # Box padding Y (meters)
    box_depth=100e3,           # Box depth (meters)
    output_mesh='merged.msh',
    use_distance_field=False
)
```

## Output Structure

### Directory Structure
```
your_working_directory/
├── kaikoura_v11.3hpk2.geo          # Input geometry
├── fault_segments/                  # Individual segment meshes
│   ├── fault_segment_6.msh
│   ├── fault_segment_12.msh
│   ├── fault_segment_18.msh
│   └── ... (22 segments total)
└── ando_mesh_merged.msh            # Final 3D mesh
```

### Physical Groups

The final mesh contains these physical groups:

| Group | Name      | Description                    | Elements |
|-------|-----------|--------------------------------|----------|
| 1     | volume    | 3D tetrahedral volume mesh    | ~34,000  |
| 101   | freesurf  | Free surface (top boundary)    | ~1,300   |
| 103   | fault     | All fault segments (embedded)  | ~6,200   |
| 105   | absorb    | Absorbing boundaries (5 sides) | ~3,600   |

## Parameters

### Fault Segment Meshing Parameters

```python
mesh_individual_fault_segments(
    geo_file='file.geo',           # Input .geo file (coordinates in km)
    output_dir='fault_segments',   # Output directory for segment meshes
    lc=2000,                       # Mesh size on fault surfaces (meters)
    remove_criteria={              # Criteria for filtering surfaces
        'z_near_zero_tol': 1000.0, # Remove surfaces with |z| < tolerance
        'z_below': -14000.0        # Remove surfaces with z < value
    }
)
```

### 3D Mesh Generation Parameters

```python
generate_mesh_from_fault_segments(
    segment_files=files,           # List of .msh files for segments
    lc_coarse=10000,              # Coarse mesh size away from faults (m)
    box_length=100e3,             # Extra padding in X direction (m)
    box_width=100e3,              # Extra padding in Y direction (m)
    box_depth=100e3,              # Depth below fault bottom (m)
    output_mesh='output.msh',     # Output mesh filename
    use_distance_field=False,     # Use distance-based refinement
    distance_threshold=5000.0     # Distance for coarsening (m)
)
```

## Command-Line Usage

```bash
python run_segmented_meshing.py <geo_file> [OPTIONS]

Required:
  geo_file              Path to .geo file

Optional:
  --segment-dir DIR     Directory for segment meshes (default: fault_segments)
  --lc-fault SIZE       Fault mesh size in meters (default: 2000)
  --lc-coarse SIZE      Coarse mesh size in meters (default: 10000)
  --box-length SIZE     Box padding X in meters (default: 100000)
  --box-width SIZE      Box padding Y in meters (default: 100000)
  --box-depth SIZE      Box depth in meters (default: 100000)
  --output FILE         Output mesh file (default: ando_fault_mesh_merged.msh)
```

### Example
```bash
python run_segmented_meshing.py kaikoura_v11.3hpk2.geo \
    --segment-dir my_segments \
    --lc-fault 1500 \
    --lc-coarse 8000 \
    --box-depth 80000 \
    --output my_mesh.msh
```

## Workflow Details

### Step 1: Segment Extraction and Meshing

For each surface in the `.geo` file:

1. **Load geometry** - Read .geo file and scale from km to m
2. **Filter surfaces** - Remove surfaces based on z-coordinate criteria
3. **Isolate segment** - Remove all other surfaces, keep only current segment
4. **Set mesh size** - Apply `lc` to all points on the segment
5. **Generate 2D mesh** - Create triangular surface mesh
6. **Save** - Write segment mesh to `fault_segments/fault_segment_<tag>.msh`

Output: Individual `.msh` files for each fault segment

### Step 2: Merging into 3D Mesh

1. **Load segments** - Merge all individual segment `.msh` files
2. **Calculate bounds** - Find combined bounding box of all segments
3. **Create box** - Generate box domain with padding around faults
4. **Embed surfaces** - Use `gmsh.model.mesh.embed()` to embed fault surfaces in volume
5. **Generate 3D mesh** - Create tetrahedral mesh using Frontal/Delaunay/HXT algorithms
6. **Assign groups** - Set physical groups (101, 103, 105, 1)
7. **Write output** - Save final mesh in MSH 2.2 format

## Comparison with Original Approach

| Feature | Original (`mesh_ando.py`) | Segmented Approach |
|---------|---------------------------|-------------------|
| Geometry loading | Load all at once | Load once, process individually |
| Mesh control | Global per surface | Per-segment control possible |
| Debugging | All-or-nothing | Individual segment inspection |
| Reusability | Full re-run needed | Reuse segment meshes |
| Complexity | Single script | Two-step workflow |
| Performance | Faster for simple cases | Better for complex geometries |

## Example Output

```
================================================================================
 SEGMENTED FAULT MESHING WORKFLOW
================================================================================
Input .geo file: kaikoura_v11.3hpk2.geo
Segment directory: fault_segments
Fault mesh size: 2000.0 m
Coarse mesh size: 10000.0 m
Output mesh: ando_mesh_from_segments.msh

================================================================================
 STEP 1: EXTRACTING AND MESHING INDIVIDUAL FAULT SEGMENTS
================================================================================

Processing segment 1/22: Surface 6
  Mesh generated: 331 nodes, 662 elements
  Saved mesh to: fault_segments/fault_segment_6.msh

...

✓ Step 1 complete: 22 fault segments meshed

================================================================================
 STEP 2: MERGING SEGMENTS INTO 3D MESH
================================================================================

Loading fault segment meshes...
Total fault surfaces loaded: 22

Creating box volume...
Box volume created: 1

Embedding 22 fault surface(s) in box volume...
Successfully embedded fault surfaces

Generating 3D mesh...
  ✓ Successfully generated 3D mesh using Frontal

Mesh saved to: ando_mesh_from_segments.msh

================================================================================
 WORKFLOW COMPLETE!
================================================================================
✓ Meshed 22 fault segments
✓ Generated 3D mesh: ando_mesh_from_segments.msh
```

## Verification

After meshing, verify physical groups:

```python
import meshio

mesh = meshio.read('ando_mesh_merged.msh')

print("Physical groups:")
for name, cells in mesh.cells_dict.items():
    print(f"  {name}: {len(cells)} cells")

for tag, name in mesh.field_data.items():
    print(f"  Group {name[0]}: {tag}")
```

Or check element counts:
```bash
python -c "
from mesh_ando_from_segments import generate_mesh_from_fault_segments
import glob
segments = sorted(glob.glob('fault_segments/*.msh'))
print(f'Found {len(segments)} segments')
"
```

## Troubleshooting

### Issue: Physical group already exists
**Solution**: The script now automatically removes old physical groups from `.geo` files

### Issue: No fault surfaces loaded
**Solution**: Check removal criteria - some surfaces may be filtered out

### Issue: Mesh generation fails
**Solution**: Try different meshing algorithms or adjust mesh sizes

### Issue: Missing segment files
**Solution**: Check that Step 1 completed successfully, verify `fault_segments/` directory exists

## Advanced Usage

### Custom mesh size per segment

```python
segment_files = []
for surf_tag, lc_custom in [(6, 1500), (12, 2000), (18, 2500)]:
    # Mesh each segment with custom lc
    files, info = mesh_individual_fault_segments(
        geo_file='file.geo',
        output_dir=f'segment_{surf_tag}',
        lc=lc_custom
    )
    segment_files.extend(files)

# Then merge all segments
final_mesh = generate_mesh_from_fault_segments(segment_files, ...)
```

### Selective segment inclusion

```python
# Mesh all segments
all_files, all_info = mesh_individual_fault_segments(...)

# Only include segments with certain tags
selected_files = [f for f in all_files if any(str(tag) in f for tag in [6, 12, 18])]

# Merge only selected segments
final_mesh = generate_mesh_from_fault_segments(selected_files, ...)
```

## Requirements

- Python 3.7+
- gmsh Python API
- numpy

Install dependencies:
```bash
pip install gmsh numpy
```

## Notes

- Input `.geo` file coordinates must be in **kilometers** (will be converted to meters)
- Output mesh is in **MSH 2.2 format** (SeisSol compatible)
- Embedded surfaces retain their original surface tags
- Lower-dimensional elements (edges, points) are automatically removed from output
