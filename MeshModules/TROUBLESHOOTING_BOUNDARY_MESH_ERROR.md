# Troubleshooting: "Could not recover boundary mesh: error 2"

## Problem

Error message:
```
Error   : Could not recover boundary mesh: error 2
```

Additionally: **No Volume physical group** appears in the generated mesh file.

---

## Root Causes

This error occurs when Gmsh cannot properly construct the 3D mesh boundary. Common causes:

### 1. **Missing synchronize() after geometry operations** 🔴
- Geometry changes not propagated to mesh
- Entities have invalid/stale references

### 2. **Physical groups added at wrong stage** 🔴
- Volume physical group must be added **BEFORE** mesh generation
- If added after failed mesh, it won't appear in output

### 3. **Boundary surfaces are not properly closed** 🔴
- Volume boundary has gaps or overlaps
- Embedded surfaces interfere with boundary

### 4. **Embedding issues with geo kernel** 🟡
- `mesh.embed()` not properly synchronized
- Volume entity changed/invalidated after embedding

### 5. **2D mesh failure** 🟡
- 3D mesh depends on successful 2D mesh
- If 2D meshing fails, 3D will also fail

---

## Solutions

### Solution 1: Use mesh_ando_robust.py (Recommended)

The robust version fixes all known issues:

```bash
python mesh_ando_robust.py kaikoura_v11.3hpk2.geo --gui
```

**Key fixes in robust version**:
1. ✅ Synchronize after every geometry operation
2. ✅ Physical groups added BEFORE meshing
3. ✅ Explicit coherence checks
4. ✅ Validate volume exists after embedding
5. ✅ Better error messages

### Solution 2: Run Diagnostics

Check what's wrong with your geometry:

```bash
python diagnose_mesh.py kaikoura_v11.3hpk2.geo
```

This will show:
- Which surfaces are being kept/removed
- Geometric bounds and orientations
- Any validation errors

### Solution 3: Manual Fixes to Original Script

If you must use the original script, apply these fixes:

#### Fix 1: Add synchronize() calls

```python
# After loading geo file
gmsh.merge(geo_file)
gmsh.model.geo.synchronize()  # ← ADD THIS

# After scaling
gmsh.model.geo.dilate(all_entities, 0, 0, 0, 1000, 1000, 1000)
gmsh.model.geo.synchronize()  # ← ADD THIS

# After removing surfaces
gmsh.model.geo.remove(surfaces_to_remove)
gmsh.model.geo.synchronize()  # ← ADD THIS

# After creating volume
box_volume = gmsh.model.geo.addVolume([box_surface_loop])
gmsh.model.geo.synchronize()  # ← ADD THIS

# After embedding
gmsh.model.mesh.embed(2, fault_surfaces, 3, box_volume)
# Do NOT synchronize after embedding!
```

#### Fix 2: Add physical group BEFORE meshing

```python
# WRONG - Physical group added too late
gmsh.model.mesh.generate(3)
gmsh.model.addPhysicalGroup(3, [box_volume], 1)  # ← Too late!

# CORRECT - Add before meshing
gmsh.model.addPhysicalGroup(3, [box_volume], 1)  # ← Before meshing
gmsh.model.setPhysicalName(3, 1, "volume")
gmsh.model.mesh.generate(2)
gmsh.model.mesh.generate(3)
```

#### Fix 3: Verify volume after embedding

```python
# After embedding, verify volume still exists
gmsh.model.mesh.embed(2, fault_surfaces, 3, box_volume)

# Check volumes
volumes = gmsh.model.getEntities(3)
if not volumes:
    raise RuntimeError("Volume disappeared after embedding!")

# Use actual volume tag
actual_volume_tag = volumes[0][1]

# Add physical group with verified tag
gmsh.model.addPhysicalGroup(3, [actual_volume_tag], 1)
```

### Solution 4: Increase Mesh Sizes

If geometry is valid but meshing still fails:

```bash
# Try coarser mesh first
python mesh_ando_robust.py fault.geo --lc 5000 --lc-coarse 20000
```

Large mesh sizes → easier meshing → can identify if problem is geometric or mesh resolution.

### Solution 5: Check Input Geometry

Verify your .geo file is valid:

```bash
# Open in Gmsh GUI
gmsh kaikoura_v11.3hpk2.geo

# Check for:
# - Self-intersecting surfaces
# - Degenerate elements
# - Inconsistent orientations
```

---

## Step-by-Step Debugging

### Step 1: Run Diagnostics

```bash
python diagnose_mesh.py your_fault.geo
```

Check output for:
- ❌ "Error - " messages
- ❌ Surfaces with invalid bounds
- ✅ Number of fault surfaces > 0

### Step 2: Try Robust Version

```bash
python mesh_ando_robust.py your_fault.geo --lc 5000 --gui
```

### Step 3: Check Generated Mesh

```bash
# View in Gmsh
gmsh ando_fault_mesh.msh

# In GUI:
# Tools → Statistics
# Check "Physical entities" section
# Should see: "Volume 1: volume"
```

### Step 4: Verify Physical Groups

```python
import gmsh

gmsh.initialize()
gmsh.open("ando_fault_mesh.msh")

# Get physical groups
groups = gmsh.model.getPhysicalGroups()
print(f"Physical groups: {len(groups)}")

for dim, tag in groups:
    name = gmsh.model.getPhysicalName(dim, tag)
    print(f"  Dim {dim}, Tag {tag}: {name}")

gmsh.finalize()
```

Expected output:
```
Physical groups: 4
  Dim 3, Tag 1: volume          ← Must be present!
  Dim 2, Tag 101: box_top
  Dim 2, Tag 103: fault
  Dim 2, Tag 105: box_sides
```

---

## Common Error Patterns

### Pattern 1: Volume Physical Group Missing

**Symptom**: Mesh file has no volume entities, only surfaces

**Cause**: Physical group added after failed mesh generation

**Fix**: Add physical group before `mesh.generate(3)`

### Pattern 2: "boundary mesh error 2" Only

**Symptom**: Error occurs, but previous versions worked

**Cause**: Geometry changed, surfaces now overlap/intersect

**Fix**: Run diagnostics to identify problematic surfaces

### Pattern 3: Works in GUI, Fails in Script

**Symptom**: Gmsh GUI can mesh it, Python script cannot

**Cause**: Script missing synchronize() calls

**Fix**: Add synchronize() after each geometry operation

### Pattern 4: 2D Mesh OK, 3D Fails

**Symptom**: 2D mesh generates, 3D mesh fails at boundary recovery

**Cause**: Embedded surfaces interfere with volume boundary

**Fix**: Use robust embedding sequence (see mesh_ando_robust.py)

---

## Quick Reference: Synchronize Locations

Always synchronize after:

```python
✅ gmsh.merge(geo_file)
✅ gmsh.model.geo.dilate(...)
✅ gmsh.model.geo.remove(...)
✅ gmsh.model.geo.addVolume(...)
❌ gmsh.model.mesh.embed(...)  # Do NOT synchronize after embed!
```

---

## Files

- **mesh_ando.py** - Original (has issues)
- **mesh_ando_fixed.py** - Fixed embedding, but may still have boundary errors
- **mesh_ando_robust.py** - **Recommended** - Fixes all boundary mesh issues
- **diagnose_mesh.py** - Diagnostic tool

---

## Still Not Working?

If none of the above solutions work:

1. **Share your .geo file** for specific help

2. **Try OCC kernel** instead of geo:
   ```python
   # Convert to OCC (more robust)
   # See mesh_ando_occ.py (TODO: create this)
   ```

3. **Simplify geometry**: Test with just one fault surface first

4. **Check Gmsh version**:
   ```bash
   gmsh --version  # Should be >= 4.8
   ```

5. **Enable debug output**:
   ```python
   gmsh.option.setNumber("General.Verbosity", 99)
   ```

---

## Prevention

To avoid this error in future:

1. ✅ Always synchronize after geometry operations
2. ✅ Add physical groups BEFORE meshing
3. ✅ Verify entities exist after embedding
4. ✅ Use robust version for production
5. ✅ Test with coarse mesh first

---

## Related Issues

- "PLC error" → Similar cause, different error message
- "HXT 3D mesh failed" → May be caused by boundary mesh error
- "Empty volume" → Physical group added too late
