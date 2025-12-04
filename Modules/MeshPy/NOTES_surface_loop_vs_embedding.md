# Surface Loop vs. Embedding for Internal Fault Surfaces

## The Question

Can we replace `gmsh.model.mesh.embed()` with a surface loop that includes both box surfaces and fault surfaces?

```python
# Attempted approach:
all_surfaces = box_surfaces + fault_surfaces
surface_loop = gmsh.model.geo.addSurfaceLoop(all_surfaces)
volume = gmsh.model.geo.addVolume([surface_loop])
```

## The Answer: No

**Internal fault surfaces that don't connect to box boundaries CANNOT be included in a surface loop.**

## Why Not?

### Surface Loop Requirements

A surface loop must form a **closed boundary** for a volume. This means:
- All surfaces must connect to form a watertight shell
- No gaps or "floating" surfaces
- Surfaces define the exterior of the volume

### Problem with Internal Faults

Internal fault surfaces:
- Are embedded **inside** the volume
- Do **NOT** connect to the box boundaries
- Do **NOT** form part of the volume's exterior boundary
- Create "open elements" when included in a surface loop

### What Happens When You Try

When you include internal fault surfaces in a surface loop:

```
Info: 2386 open elements
Error: No elements in volume 1
```

The 3D mesher fails because:
1. The surface loop is not closed (fault surfaces float inside)
2. The volume boundary is ambiguous
3. Cannot determine interior vs exterior

## The Correct Approach: Embedding

### What is Embedding?

```python
gmsh.model.mesh.embed(2, fault_surface_tags, 3, box_volume)
```

This is **NOT** a geometry operation. It's a **mesh constraint** that tells the mesher:
- "When you generate the 3D tetrahedral mesh for this volume..."
- "...make sure the mesh conforms to these internal surfaces"
- "...without splitting the volume geometrically"

### Why Embedding Works

1. **Geometry level**: Volume has a simple closed boundary (6 box surfaces)
2. **Mesh level**: Tetrahedral elements respect internal fault surfaces
3. **Result**: Fault surfaces are represented in the mesh but don't affect volume topology

## Comparison

### Surface Loop (for external boundaries)

```python
# Works for: Volumes bounded by connected surfaces
surfaces = [bottom, top, front, back, left, right]
surface_loop = gmsh.model.geo.addSurfaceLoop(surfaces)
volume = gmsh.model.geo.addVolume([surface_loop])
```

**Use when:** All surfaces form a closed, connected boundary

### Embedding (for internal features)

```python
# Works for: Internal surfaces that don't touch boundaries
box_loop = gmsh.model.geo.addSurfaceLoop(box_surfaces)
volume = gmsh.model.geo.addVolume([box_loop])
gmsh.model.geo.synchronize()
gmsh.model.mesh.embed(2, fault_surfaces, 3, volume)
```

**Use when:** Surfaces are internal features within a volume

## When CAN You Use Surface Loop with Faults?

Surface loop with faults works ONLY if:
1. Fault surfaces connect to the box boundaries
2. Fault surfaces divide the volume into separate regions
3. You want to create multiple volumes (not one volume with internal features)

Example: A fault that cuts completely through the box

```python
# Fault cuts box into two volumes
left_surfaces = [bottom, left, front, back, fault]
right_surfaces = [bottom, right, front, back, fault]

left_loop = gmsh.model.geo.addSurfaceLoop(left_surfaces)
right_loop = gmsh.model.geo.addSurfaceLoop(right_surfaces)

left_volume = gmsh.model.geo.addVolume([left_loop])
right_volume = gmsh.model.geo.addVolume([right_loop])
```

## Conclusion

For the Kaikoura fault segments:
- Fault surfaces are **internal** (don't connect to boundaries)
- Must use **embedding** approach
- Surface loop approach **will fail** with "open elements" error

The embedding approach is not just acceptable—it's the **correct and only** way to handle internal fault surfaces in gmsh.
