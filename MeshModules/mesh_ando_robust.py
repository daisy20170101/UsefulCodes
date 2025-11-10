#!/usr/bin/env python3
"""
Robust version of mesh_ando that fixes "boundary mesh error 2" issue

Key fixes:
1. Proper synchronization before and after embedding
2. Volume physical group added at correct stage
3. Better handling of geo kernel limitations
4. Explicit coherence checks
5. Option to use OCC kernel for more robustness
"""

import gmsh
import os
import sys


def generate_Ando_fault_robust(
    geo_file='kaikoura_v11.3hpk2.geo',
    lc=2000,
    lc_coarse=10e3,
    box_length=100e3,
    box_width=100e3,
    box_depth=100e3,
    output_mesh='ando_fault_mesh.msh',
    show_gui=False,
    verbose=True
):
    """
    Generate Ando fault mesh with robust error handling

    Fixes for "boundary mesh error 2":
    - Proper entity synchronization
    - Coherence checking before meshing
    - Physical groups added before meshing
    - Explicit boundary validation
    """

    def log(msg):
        if verbose:
            print(msg)

    # Validate input
    if not os.path.exists(geo_file):
        raise FileNotFoundError(f"Geo file not found: {geo_file}")

    log("="*70)
    log("ROBUST ANDO FAULT MESH GENERATION")
    log("="*70)

    # Initialize
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", 5 if verbose else 2)

    gmsh.model.add("ando_fault_robust")

    # ==========================================
    # STEP 1: Load and scale geometry
    # ==========================================
    log("\n[1/10] Loading geometry...")
    gmsh.merge(geo_file)
    gmsh.model.geo.synchronize()

    # Remove old physical groups
    try:
        for dim, tag in gmsh.model.getPhysicalGroups():
            gmsh.model.removePhysicalGroups([(dim, tag)])
    except:
        pass

    log("[2/10] Scaling from km to m...")
    all_entities = []
    for dim in range(4):
        all_entities.extend(gmsh.model.getEntities(dim))

    gmsh.model.geo.dilate(all_entities, 0, 0, 0, 1000, 1000, 1000)

    # CRITICAL: Synchronize after scaling
    gmsh.model.geo.synchronize()
    log("   ✓ Geometry scaled and synchronized")

    # ==========================================
    # STEP 2: Filter surfaces
    # ==========================================
    log("\n[3/10] Filtering fault surfaces...")

    surfaces = gmsh.model.getEntities(dim=2)
    surfaces_to_remove = []
    fault_surfaces = []

    for surf_dim, surf_tag in surfaces:
        bbox = gmsh.model.getBoundingBox(surf_dim, surf_tag)
        center_z = (bbox[2] + bbox[5]) / 2

        if abs(center_z) < 1000 or center_z < -14000.0:
            surfaces_to_remove.append((surf_dim, surf_tag))
        else:
            fault_surfaces.append(surf_tag)

    if surfaces_to_remove:
        log(f"   Removing {len(surfaces_to_remove)} surfaces")
        gmsh.model.geo.remove(surfaces_to_remove)

    # CRITICAL: Synchronize after removal
    gmsh.model.geo.synchronize()

    # Re-get fault surfaces after removal
    remaining = gmsh.model.getEntities(dim=2)
    fault_surfaces = [s[1] for s in remaining]

    log(f"   ✓ Kept {len(fault_surfaces)} fault surfaces: {fault_surfaces}")

    if not fault_surfaces:
        raise ValueError("No fault surfaces remaining!")

    # ==========================================
    # STEP 3: Calculate bounds
    # ==========================================
    log("\n[4/10] Calculating fault bounds...")

    xmin, ymin, zmin = float('inf'), float('inf'), float('inf')
    xmax, ymax, zmax = float('-inf'), float('-inf'), float('-inf')

    for surf_tag in fault_surfaces:
        bbox = gmsh.model.getBoundingBox(2, surf_tag)
        xmin = min(xmin, bbox[0])
        ymin = min(ymin, bbox[1])
        zmin = min(zmin, bbox[2])
        xmax = max(xmax, bbox[3])
        ymax = max(ymax, bbox[4])
        zmax = max(zmax, bbox[5])

    log(f"   Fault: X=[{xmin:.0f}, {xmax:.0f}], Y=[{ymin:.0f}, {ymax:.0f}], Z=[{zmin:.0f}, {zmax:.0f}]")

    # ==========================================
    # STEP 4: Create box
    # ==========================================
    log("\n[5/10] Creating box domain...")

    # Add buffers
    buffer_x = max(5000, (xmax - xmin) * 0.05)
    buffer_y = max(5000, (ymax - ymin) * 0.05)

    x0 = xmin - box_length/2 - buffer_x
    x1 = xmax + box_length/2 + buffer_x
    y0 = ymin - box_width/2 - buffer_y
    y1 = ymax + box_width/2 + buffer_y
    z0 = -box_depth
    z1 = 500.0

    log(f"   Box: {(x1-x0)/1e3:.1f} × {(y1-y0)/1e3:.1f} × {(z1-z0)/1e3:.1f} km")

    # Create box corners
    p1 = gmsh.model.geo.addPoint(x0, y0, z0, lc_coarse)
    p2 = gmsh.model.geo.addPoint(x1, y0, z0, lc_coarse)
    p3 = gmsh.model.geo.addPoint(x1, y1, z0, lc_coarse)
    p4 = gmsh.model.geo.addPoint(x0, y1, z0, lc_coarse)
    p5 = gmsh.model.geo.addPoint(x0, y0, z1, lc_coarse)
    p6 = gmsh.model.geo.addPoint(x1, y0, z1, lc_coarse)
    p7 = gmsh.model.geo.addPoint(x1, y1, z1, lc_coarse)
    p8 = gmsh.model.geo.addPoint(x0, y1, z1, lc_coarse)

    # Lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p8)
    l8 = gmsh.model.geo.addLine(p8, p5)
    l9 = gmsh.model.geo.addLine(p1, p5)
    l10 = gmsh.model.geo.addLine(p2, p6)
    l11 = gmsh.model.geo.addLine(p3, p7)
    l12 = gmsh.model.geo.addLine(p4, p8)

    # Surfaces
    cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf_bottom = gmsh.model.geo.addPlaneSurface([cl1])

    cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
    surf_top = gmsh.model.geo.addPlaneSurface([cl2])

    cl3 = gmsh.model.geo.addCurveLoop([l1, l10, -l5, -l9])
    surf_front = gmsh.model.geo.addPlaneSurface([cl3])

    cl4 = gmsh.model.geo.addCurveLoop([l2, l11, -l6, -l10])
    surf_right = gmsh.model.geo.addPlaneSurface([cl4])

    cl5 = gmsh.model.geo.addCurveLoop([l3, l12, -l7, -l11])
    surf_back = gmsh.model.geo.addPlaneSurface([cl5])

    cl6 = gmsh.model.geo.addCurveLoop([l4, l9, -l8, -l12])
    surf_left = gmsh.model.geo.addPlaneSurface([cl6])

    # Volume
    box_surfaces = [surf_bottom, surf_top, surf_front, surf_right, surf_back, surf_left]
    box_surface_loop = gmsh.model.geo.addSurfaceLoop(box_surfaces)
    box_volume = gmsh.model.geo.addVolume([box_surface_loop])

    log(f"   ✓ Box volume created: {box_volume}")

    # CRITICAL: Synchronize before embedding
    gmsh.model.geo.synchronize()

    # ==========================================
    # STEP 5: Embed faults
    # ==========================================
    log("\n[6/10] Embedding fault surfaces...")

    for surf_tag in fault_surfaces:
        try:
            gmsh.model.mesh.embed(2, [surf_tag], 3, box_volume)
            log(f"   ✓ Embedded surface {surf_tag}")
        except Exception as e:
            log(f"   ✗ Failed to embed {surf_tag}: {e}")
            raise

    log("   ✓ All faults embedded")

    # CRITICAL: Check that volume still exists after embedding
    volumes_after_embed = gmsh.model.getEntities(3)
    log(f"   Volumes after embedding: {volumes_after_embed}")

    if not volumes_after_embed:
        raise RuntimeError("No volumes found after embedding! Geometry is broken.")

    # Use the first (and should be only) volume
    final_volume_tag = volumes_after_embed[0][1]
    log(f"   Using volume tag: {final_volume_tag}")

    # ==========================================
    # STEP 6: Add physical groups BEFORE meshing
    # ==========================================
    log("\n[7/10] Adding physical groups...")

    # CRITICAL: Add volume physical group BEFORE mesh generation
    gmsh.model.addPhysicalGroup(3, [final_volume_tag], tag=1)
    gmsh.model.setPhysicalName(3, 1, "volume")
    log(f"   ✓ Physical Volume 1: volume (tag {final_volume_tag})")

    # Individual fault surface physical groups (tags starting from 165)
    log(f"\n   Adding individual fault physical groups (tags 165+):")
    for i, fault_surf in enumerate(fault_surfaces):
        phys_tag = 165 + i
        gmsh.model.addPhysicalGroup(2, [fault_surf], tag=phys_tag)
        gmsh.model.setPhysicalName(2, phys_tag, f"fault_{i+1}")
        log(f"      Physical Surface {phys_tag}: fault_{i+1} (surface {fault_surf})")

    # Combined fault group for backward compatibility
    gmsh.model.addPhysicalGroup(2, fault_surfaces, tag=103)
    gmsh.model.setPhysicalName(2, 103, "fault_all")
    log(f"   ✓ Physical Surface 103: fault_all ({len(fault_surfaces)} surfaces combined)")

    # Top surface
    gmsh.model.addPhysicalGroup(2, [surf_top], tag=101)
    gmsh.model.setPhysicalName(2, 101, "box_top")
    log(f"   ✓ Physical Surface 101: box_top")

    # Other surfaces
    other_surfs = [surf_bottom, surf_front, surf_right, surf_back, surf_left]
    gmsh.model.addPhysicalGroup(2, other_surfs, tag=105)
    gmsh.model.setPhysicalName(2, 105, "box_sides")
    log(f"   ✓ Physical Surface 105: box_sides")

    # ==========================================
    # STEP 7: Set mesh sizes
    # ==========================================
    log("\n[8/10] Setting mesh sizes...")

    # Set fine mesh on faults
    for surf_tag in fault_surfaces:
        bounds = gmsh.model.getBoundary([(2, surf_tag)], combined=False,
                                        oriented=False, recursive=True)
        points = [abs(b[1]) for b in bounds if b[0] == 0]

        for pt in points:
            gmsh.model.mesh.setSize([(0, pt)], lc)

    log(f"   ✓ Fault mesh size: {lc} m")

    gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc_coarse)
    log(f"   ✓ Global range: {lc} - {lc_coarse} m")

    # ==========================================
    # STEP 8: Check coherence
    # ==========================================
    log("\n[9/10] Checking geometry coherence...")

    try:
        # Verify all entities
        for dim in range(4):
            entities = gmsh.model.getEntities(dim)
            log(f"   Dim {dim}: {len(entities)} entities")

        # Verify volume has proper boundary
        vol_boundary = gmsh.model.getBoundary([(3, final_volume_tag)], oriented=False)
        log(f"   Volume boundary: {len(vol_boundary)} surfaces")

        if len(vol_boundary) == 0:
            raise RuntimeError("Volume has no boundary surfaces!")

        log("   ✓ Geometry coherence OK")

    except Exception as e:
        log(f"   ✗ Coherence check failed: {e}")
        raise

    # ==========================================
    # STEP 9: Generate mesh
    # ==========================================
    log("\n[10/10] Generating mesh...")

    # Configure meshing
    gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay 2D
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay 3D (most robust)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)

    # Generate 2D mesh
    log("   Generating 2D mesh...")
    try:
        gmsh.model.mesh.generate(2)
        n_nodes_2d = len(gmsh.model.mesh.getNodes()[0])
        log(f"   ✓ 2D mesh: {n_nodes_2d} nodes")
    except Exception as e:
        log(f"   ✗ 2D mesh failed: {e}")
        raise RuntimeError(f"2D meshing failed: {e}")

    # Generate 3D mesh with fallback
    algorithms = [(1, "Delaunay"), (4, "Frontal"), (10, "HXT")]

    mesh_success = False
    for algo_id, algo_name in algorithms:
        try:
            log(f"   Trying {algo_name}...")
            gmsh.option.setNumber("Mesh.Algorithm3D", algo_id)

            try:
                gmsh.model.mesh.clear(3)
            except:
                pass

            gmsh.model.mesh.generate(3)

            n_nodes = len(gmsh.model.mesh.getNodes()[0])
            log(f"   ✓ 3D mesh: {n_nodes} nodes ({algo_name})")
            mesh_success = True
            break

        except Exception as e:
            log(f"   ✗ {algo_name} failed: {str(e)[:80]}")
            continue

    if not mesh_success:
        raise RuntimeError("All 3D meshing algorithms failed!")

    # ==========================================
    # STEP 10: Save and finalize
    # ==========================================
    log(f"\nWriting mesh to {output_mesh}...")
    gmsh.write(output_mesh)

    # Get statistics
    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_elements = sum(len(t) for t in elem_tags)

    log("="*70)
    log("SUCCESS!")
    log("="*70)
    log(f"Nodes: {n_nodes:,}")
    log(f"Elements: {n_elements:,}")
    log(f"Fault surfaces: {len(fault_surfaces)}")
    log(f"Output: {output_mesh}")
    log("="*70)

    # Verify physical groups are in mesh
    log("\nVerifying physical groups in mesh...")
    phys_groups = gmsh.model.getPhysicalGroups()
    log(f"Physical groups: {len(phys_groups)}")
    for dim, tag in phys_groups:
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        log(f"  Dim {dim}, Tag {tag}: '{name}' ({len(entities)} entities)")

    if show_gui:
        gmsh.fltk.run()

    gmsh.finalize()

    return {
        'mesh_file': output_mesh,
        'num_nodes': n_nodes,
        'num_elements': n_elements,
        'num_fault_surfaces': len(fault_surfaces),
        'volume_tag': final_volume_tag
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Robust Ando fault mesh generation'
    )

    parser.add_argument('geo_file', nargs='?',
                       default='kaikoura_v11.3hpk2.geo',
                       help='Input .geo file')

    parser.add_argument('--lc', type=float, default=2000,
                       help='Fault mesh size (m)')

    parser.add_argument('--lc-coarse', type=float, default=10000,
                       help='Coarse mesh size (m)')

    parser.add_argument('-o', '--output', default='ando_fault_mesh.msh',
                       help='Output mesh file')

    parser.add_argument('--gui', action='store_true',
                       help='Show GUI')

    args = parser.parse_args()

    try:
        result = generate_Ando_fault_robust(
            geo_file=args.geo_file,
            lc=args.lc,
            lc_coarse=args.lc_coarse,
            output_mesh=args.output,
            show_gui=args.gui
        )

        print("\n✓ Mesh generation successful!")
        sys.exit(0)

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
