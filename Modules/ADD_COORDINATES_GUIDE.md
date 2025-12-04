# Adding Coordinates from XDMF to Calculate Seismic Moment

If your fault CSV files don't have x, y, z coordinates, you need to extract them from the XDMF file first.

## Problem

Your stress CSV files (e.g., `stress_jp3z_final.csv`) have:
- ✓ ASl (slip)
- ✓ Area
- ✓ fault-tag
- ✗ **Missing**: x, y, z coordinates

But `calculate_seismic_moment_by_fault()` needs coordinates to find closest mu values from `mat3d_fault.csv`.

## Solution

Use `add_coordinates_from_xdmf()` to extract element centroids from the XDMF file.

## Complete Workflow

### In your Jupyter notebook ([cov-analys.ipynb](file:///Users/DuoL/Documents/NSHM/PyScripts/cov-analys.ipynb)):

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from StatsModel.cov_corr_funcs import (add_coordinates_from_xdmf,
                                       calculate_seismic_moment_by_fault,
                                       plot_seismic_moment_summary)

# ============================================================================
# STEP 1: Load fault data
# ============================================================================
modelname = 'jp3z'
ssTable1 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/slab2/stress_{modelname}_final.csv')
ssTable0 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/slab2/stress_{modelname}_t0.csv')

print(f"Loaded {len(ssTable1)} fault elements")

# Calculate derived quantities
ssTable1['sdrop'] = np.sqrt(ssTable0['Td0']**2 + ssTable0['Ts0']**2) - \
                    np.sqrt(ssTable1['Td0']**2 + ssTable1['Ts0']**2)
ssTable1['T0'] = np.sqrt(ssTable0['Td0']**2 + ssTable0['Ts0']**2)
ssTable1['Pf'] = np.abs(ssTable0['Pn0'])
ssTable1['r'] = ssTable1['T0'] / ssTable1['Pf']

# ============================================================================
# STEP 2: Add x, y, z coordinates from XDMF
# ============================================================================
xdmf_file = f'/Users/DuoL/Documents/NSHM/Central/Joint3/data-{modelname}/{modelname}-fault.xdmf'

ssTable1 = add_coordinates_from_xdmf(
    fault_data=ssTable1,
    xdmf_file=xdmf_file,
    method='centroid'  # Use element centroids
)

# Verify coordinates were added
print("\nCoordinates added:")
print(ssTable1[['x', 'y', 'z', 'ASl', 'Area', 'fault-tag']].head())

# ============================================================================
# STEP 3: Filter data (optional)
# ============================================================================
ssTable2 = ssTable1[ssTable1['ASl'].between(0.5, 100)]
print(f"\nFiltered: {len(ssTable1)} → {len(ssTable2)} rows (0.5 < ASl < 100 m)")

# ============================================================================
# STEP 4: Calculate seismic moment by fault segment
# ============================================================================
material_file = '/Users/DuoL/Documents/NSHM/Central/cov/mat3d_fault.csv'

df_moment = calculate_seismic_moment_by_fault(
    fault_data=ssTable2,
    material_file=material_file
)

# ============================================================================
# STEP 5: Display and save results
# ============================================================================
print("\n" + "="*70)
print("SEISMIC MOMENT RESULTS")
print("="*70)
print(df_moment)

# Summary statistics
total_M0 = df_moment['total_M0'].sum()
total_Mw = (2/3) * np.log10(total_M0) - 6.07

print(f"\nTotal seismic moment: {total_M0:.3e} N⋅m")
print(f"Equivalent Mw: {total_Mw:.2f}")

# Save results
df_moment.to_csv(f'moment_by_fault_{modelname}.csv', index=False)
print(f"\n✓ Saved to: moment_by_fault_{modelname}.csv")

# ============================================================================
# STEP 6: Visualize
# ============================================================================
fig, axes = plot_seismic_moment_summary(
    df_moment,
    save_path=f'moment_summary_{modelname}.png'
)
plt.show()
```

## What Each Step Does

### Step 1: Load Data
- Reads your stress CSV files (final and t0)
- Calculates derived quantities (sdrop, T0, Pf, r)

### Step 2: Add Coordinates ⭐ **NEW**
```python
ssTable1 = add_coordinates_from_xdmf(
    fault_data=ssTable1,
    xdmf_file=xdmf_file,
    method='centroid'
)
```

**What it does:**
1. Opens the XDMF file using `seissolxdmf`
2. Reads vertices (node coordinates): `(n_vertices, 3)`
3. Reads connectivity (which vertices form each element): `(n_elements, 3)`
4. For each element:
   - Gets coordinates of its 3 vertices
   - Computes centroid: `(v0 + v1 + v2) / 3`
5. Adds `x`, `y`, `z` columns to your DataFrame

**Output:**
```
======================================================================
ADDING COORDINATES FROM XDMF
======================================================================

Loading XDMF file: /Users/DuoL/Documents/NSHM/Central/Joint3/data-jp3z/jp3z-fault.xdmf
  Reading geometry (vertices)...
  Reading connectivity (elements)...

XDMF data:
  Vertices: 50847
  Elements: 101492
  Vertices per element: 3

Computing element positions using method='centroid'...

Coordinate ranges:
  x: [1234567.00, 1345678.00]
  y: [5678901.00, 5789012.00]
  z: [-50000.00, -5000.00]

✓ Added coordinates to DataFrame
  New columns: ['x', 'y', 'z']
======================================================================
```

### Step 3: Filter (Optional)
- Keep only elements with significant slip
- Common filter: `0.5 < ASl < 100` meters

### Step 4: Calculate Seismic Moment
```python
df_moment = calculate_seismic_moment_by_fault(
    fault_data=ssTable2,
    material_file=material_file
)
```

**What it does:**
1. Loads material properties (x, y, z, mu) from `mat3d_fault.csv`
2. For each fault element:
   - Finds closest material point using KDTree
   - Gets its mu (shear modulus) value
3. Calculates M₀ = mu × ASl × Area for each element
4. Groups by fault-tag and sums
5. Calculates Mw = (2/3) × log₁₀(M₀) - 6.07

### Step 5: Results
**Output DataFrame:**
```
   fault-tag  count      total_M0    Mw  mean_slip  total_area
0         68  44737  2.456e+21  8.45      15.23  8.45e+09
1          3  26861  1.234e+21  8.12      10.45  6.23e+09
2         65  23639  8.765e+20  7.82       6.78  5.12e+09
```

### Step 6: Visualization
Creates 4-panel plot:
1. Total M₀ by fault segment (bar/scatter)
2. Moment magnitude (Mw) by segment
3. M₀ vs slip (bubble size = area, color = mu)
4. Summary statistics table

## File Requirements

### Required Files:

1. **Fault CSV** (from SeisSol output):
   - Location: `/Users/DuoL/Documents/NSHM/Central/cov/slab2/`
   - Files: `stress_{model}_final.csv`, `stress_{model}_t0.csv`
   - Has: ASl, Area, fault-tag, stress components
   - **Missing**: x, y, z (need to add!)

2. **XDMF File** (SeisSol mesh with geometry):
   - Location: `/Users/DuoL/Documents/NSHM/Central/Joint3/data-{model}/`
   - File: `{model}-fault.xdmf`
   - Contains: Vertices and element connectivity
   - **Purpose**: Extract x, y, z coordinates

3. **Material CSV** (3D material properties):
   - Location: `/Users/DuoL/Documents/NSHM/Central/cov/`
   - File: `mat3d_fault.csv`
   - Has: x, y, z, mu, lambda, rho
   - **Purpose**: Provide mu values for moment calculation

## Common Issues

### Issue 1: "Row count mismatch"
```
⚠ WARNING: Row count mismatch!
  fault_data rows: 101492
  XDMF elements: 101421
```

**Cause**: CSV and XDMF from different simulations or files out of sync

**Fix**: Ensure CSV and XDMF are from the same model/simulation

### Issue 2: "seissolxdmf module not found"
```
ImportError: seissolxdmf module is required
```

**Fix**: Install seissolxdmf from SeisSol repository

### Issue 3: "Missing columns in fault_data"
```
ValueError: Missing columns in fault_data: ['x', 'y', 'z']
```

**Cause**: Forgot to run `add_coordinates_from_xdmf()` first

**Fix**: Add Step 2 before Step 4

## Performance

**Typical performance:**
- **Add coordinates**: ~2-5 seconds for 100K elements
- **Calculate moment**: ~3-10 seconds for 50K elements
- **Total workflow**: ~10-20 seconds per model

**Memory usage:**
- ~200-500 MB depending on model size

## Alternative: If Coordinates Already Exist

If your CSV already has x, y, z columns, skip Step 2:

```python
# Load data (already has x, y, z)
ssTable1 = pd.read_csv('stress_model_final.csv')

# Go directly to filtering and moment calculation
ssTable2 = ssTable1[ssTable1['ASl'].between(0.5, 100)]
df_moment = calculate_seismic_moment_by_fault(ssTable2, 'mat3d_fault.csv')
```

## References

- seissolxdmf: https://github.com/SeisSol/SeisSol/tree/master/postprocessing/science/seissolxdmf
- XDMF format: http://www.xdmf.org/
- SeisSol documentation: https://seissol.readthedocs.io/
