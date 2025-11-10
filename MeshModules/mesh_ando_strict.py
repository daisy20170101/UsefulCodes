#!/usr/bin/env python3
"""
Strict version with explicit volume validation at every step

This version aggressively checks that the volume exists and is valid
throughout the entire process, and ONLY reports success if 3D mesh
actually contains volume elements.
"""

import gmsh
import os
import sys


def generate_Ando_fault_strict(
    geo_file='kaikoura_v11.3hpk2.geo',
    lc=2000,
    lc_coarse=10e3,
    box_length=100e3,
    box_width=100e3,
    box_depth=100e3,
    output_mesh='ando_fault_mesh.msh',
    show_gui=False
):
    """
    Generate mesh with strict volume validation

    Will fail loudly if volume is not properly created or meshed.
    """

    print("="*70)
    print("STRICT ANDO MESH GENERATION")
    print("="*70)

    if not os.path.exists(geo_file):
        raise FileNotFoundError(f"Geo file not found: {geo_file}")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", 5)

    gmsh.model.add("ando_fault_strict")

    try:
        # ==========================================
        # Load and scale geometry
        # ==========================================
        print("\n[1/12] Loading geometry...")
        gmsh.merge(geo_file)
        gmsh.model.geo.synchronize()

        # Remove old physical groups
        for dim, tag in gmsh.model.getPhysicalGroups():
            gmsh.model.removePhysicalGroups([(dim, tag)])

        print("[2/12] Scaling from km to m...")
        all_entities = []
        for dim in range(4):
            all_entities.extend(gmsh.model.getEntities(dim))

        gmsh.model.geo.dilate(all_entities, 0, 0, 0, 1000, 1000, 1000)
        gmsh.model.geo.synchronize()
        print("   ✓ Scaled")

        # ==========================================
        # Filter surfaces
        # ==========================================
        print("\n[3/12] Filtering surfaces...")
        surfaces = gmsh.model.getEntities(dim=2)
        surfaces_to_remove = []

        for surf_dim, surf_tag in surfaces:
            bbox = gmsh.model.getBoundingBox(surf_dim, surf_tag)
            center_z = (bbox[2] + bbox[5]) / 2

            if abs(center_z) < 1000 or center_z < -14000.0:
                surfaces_to_remove.append((surf_dim, surf_tag))

        if surfaces_to_remove:
            gmsh.model.geo.remove(surfaces_to_remove)

        gmsh.model.geo.synchronize()

        remaining = gmsh.model.getEntities(dim=2)
        fault_surfaces = [s[1] for s in remaining]
        print(f"   ✓ Kept {len(fault_surfaces)} fault surfaces")

        if not fault_surfaces:
            raise ValueError("No fault surfaces remaining!")

        # ==========================================
        # Calculate bounds
        # ==========================================
        print("\n[4/12] Calculating bounds...")
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

        print(f"   Fault extent: {(xmax-xmin)/1e3:.1f} × {(ymax-ymin)/1e3:.1f} × {(zmax-zmin)/1e3:.1f} km")

        # ==========================================
        # Create box - NO EMBEDDING YET
        # ==========================================
        print("\n[5/12] Creating box WITHOUT embedded faults first...")

        buffer_x = max(10000, (xmax - xmin) * 0.1)  # Larger buffer
        buffer_y = max(10000, (ymax - ymin) * 0.1)

        x0 = xmin - box_length/2 - buffer_x
        x1 = xmax + box_length/2 + buffer_x
        y0 = ymin - box_width/2 - buffer_y
        y1 = ymax + box_width/2 + buffer_y
        z0 = -box_depth
        z1 = 1000.0  # Higher top

        print(f"   Box: {(x1-x0)/1e3:.1f} × {(y1-y0)/1e3:.1f} × {(z1-z0)/1e3:.1f} km")

        # Box geometry
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

        print(f"   ✓ Created volume: {box_volume}")

        gmsh.model.geo.synchronize()

        # ==========================================
        # CRITICAL: Verify volume exists BEFORE embedding
        # ==========================================
        print("\n[6/12] Verifying volume exists...")
        volumes_before = gmsh.model.getEntities(3)
        print(f"   Volumes before embedding: {volumes_before}")

        if not volumes_before:
            raise RuntimeError("Volume does not exist before embedding!")

        if len(volumes_before) != 1:
            raise RuntimeError(f"Expected 1 volume, found {len(volumes_before)}")

        vol_tag_before = volumes_before[0][1]
        print(f"   ✓ Volume {vol_tag_before} exists")

        # ==========================================
        # Embed faults
        # ==========================================
        print("\n[7/12] Embedding fault surfaces...")

        for i, surf_tag in enumerate(fault_surfaces):
            print(f"   Embedding surface {surf_tag} ({i+1}/{len(fault_surfaces)})...")
            try:
                gmsh.model.mesh.embed(2, [surf_tag], 3, vol_tag_before)
                print(f"      ✓ Embedded")
            except Exception as e:
                print(f"      ✗ Failed: {e}")
                raise

        # ==========================================
        # CRITICAL: Verify volume AFTER embedding
        # ==========================================
        print("\n[8/12] Verifying volume after embedding...")
        volumes_after = gmsh.model.getEntities(3)
        print(f"   Volumes after embedding: {volumes_after}")

        if not volumes_after:
            raise RuntimeError("CRITICAL: Volume disappeared after embedding!")

        vol_tag_after = volumes_after[0][1]
        print(f"   ✓ Volume still exists: {vol_tag_after}")

        # Check volume boundary
        try:
            boundary = gmsh.model.getBoundary([(3, vol_tag_after)], oriented=False)
            print(f"   Volume boundary: {len(boundary)} surfaces")
            if len(boundary) == 0:
                raise RuntimeError("Volume has no boundary!")
        except Exception as e:
            print(f"   ✗ Cannot get volume boundary: {e}")
            raise

        # ==========================================
        # Add physical groups BEFORE meshing
        # ==========================================
        print("\n[9/12] Adding physical groups...")

        gmsh.model.addPhysicalGroup(3, [vol_tag_after], tag=1)
        gmsh.model.setPhysicalName(3, 1, "volume")
        print(f"   ✓ Physical Volume 1: volume (tag {vol_tag_after})")

        gmsh.model.addPhysicalGroup(2, fault_surfaces, tag=103)
        gmsh.model.setPhysicalName(2, 103, "fault")
        print(f"   ✓ Physical Surface 103: fault ({len(fault_surfaces)} surfaces)")

        gmsh.model.addPhysicalGroup(2, [surf_top], tag=101)
        gmsh.model.setPhysicalName(2, 101, "box_top")

        other_surfs = [surf_bottom, surf_front, surf_right, surf_back, surf_left]
        gmsh.model.addPhysicalGroup(2, other_surfs, tag=105)
        gmsh.model.setPhysicalName(2, 105, "box_sides")

        # ==========================================
        # Set mesh sizes
        # ==========================================
        print("\n[10/12] Setting mesh sizes...")

        for surf_tag in fault_surfaces:
            bounds = gmsh.model.getBoundary([(2, surf_tag)], combined=False,
                                            oriented=False, recursive=True)
            points = [abs(b[1]) for b in bounds if b[0] == 0]
            for pt in points:
                gmsh.model.mesh.setSize([(0, pt)], lc)

        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc_coarse)
        print(f"   ✓ Mesh size: {lc} - {lc_coarse} m")

        # ==========================================
        # Generate 2D mesh
        # ==========================================
        print("\n[11/12] Generating 2D mesh...")

        gmsh.option.setNumber("Mesh.Algorithm", 5)

        try:
            gmsh.model.mesh.generate(2)
            n_nodes_2d = len(gmsh.model.mesh.getNodes()[0])
            print(f"   ✓ 2D mesh: {n_nodes_2d} nodes")
        except Exception as e:
            print(f"   ✗ 2D mesh FAILED: {e}")
            raise RuntimeError(f"2D meshing failed: {e}")

        # ==========================================
        # Generate 3D mesh - TRY ALL ALGORITHMS
        # ==========================================
        print("\n[12/12] Generating 3D mesh...")

        algorithms = [
            (4, "Frontal"),
            (1, "Delaunay"),
            (10, "HXT")
        ]

        mesh_3d_success = False

        for algo_id, algo_name in algorithms:
            print(f"\n   Trying {algo_name}...")

            try:
                gmsh.option.setNumber("Mesh.Algorithm3D", algo_id)

                # Clear previous 3D mesh
                try:
                    gmsh.model.mesh.clear(3)
                except:
                    pass

                # Generate 3D mesh
                gmsh.model.mesh.generate(3)

                # CRITICAL: Check if 3D elements were actually created
                elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)

                if elem_tags and len(elem_tags) > 0:
                    num_3d_elements = sum(len(t) for t in elem_tags)

                    if num_3d_elements > 0:
                        n_nodes = len(gmsh.model.mesh.getNodes()[0])
                        print(f"   ✓ SUCCESS: {n_nodes} nodes, {num_3d_elements} 3D elements")
                        mesh_3d_success = True
                        break
                    else:
                        print(f"   ✗ No 3D elements created (only 2D mesh)")
                        continue
                else:
                    print(f"   ✗ No 3D elements created")
                    continue

            except Exception as e:
                print(f"   ✗ {algo_name} failed: {str(e)[:100]}")
                continue

        if not mesh_3d_success:
            raise RuntimeError(
                "ALL 3D MESHING ALGORITHMS FAILED!\n"
                "Possible causes:\n"
                "1. Fault surfaces interfere with volume boundary\n"
                "2. Geometric intersections/overlaps\n"
                "3. Mesh size too small/large\n"
                "4. Volume is invalid or has holes\n\n"
                "Try:\n"
                "- Increase --lc and --lc-coarse\n"
                "- Check geometry with diagnose_mesh.py\n"
                "- Reduce number of fault surfaces"
            )

        # ==========================================
        # Final validation before saving
        # ==========================================
        print("\n" + "="*70)
        print("FINAL VALIDATION")
        print("="*70)

        # Check volumes
        final_volumes = gmsh.model.getEntities(3)
        print(f"Volumes in model: {len(final_volumes)}")
        if not final_volumes:
            raise RuntimeError("No volumes in final model!")

        # Check physical groups
        phys_groups = gmsh.model.getPhysicalGroups()
        has_volume_phys = any(dim == 3 for dim, tag in phys_groups)
        print(f"Physical groups: {len(phys_groups)}")
        if not has_volume_phys:
            raise RuntimeError("Volume physical group missing!")

        # Check 3D elements
        elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
        num_3d_elements = sum(len(t) for t in elem_tags) if elem_tags else 0
        print(f"3D elements: {num_3d_elements:,}")
        if num_3d_elements == 0:
            raise RuntimeError("No 3D elements in mesh!")

        print("="*70)
        print("✓ ALL CHECKS PASSED")
        print("="*70)

        # ==========================================
        # Save
        # ==========================================
        print(f"\nWriting mesh to {output_mesh}...")
        gmsh.write(output_mesh)

        n_nodes = len(gmsh.model.mesh.getNodes()[0])

        print("\n" + "="*70)
        print("SUCCESS!")
        print("="*70)
        print(f"Output: {output_mesh}")
        print(f"Nodes: {n_nodes:,}")
        print(f"3D Elements: {num_3d_elements:,}")
        print(f"Fault surfaces: {len(fault_surfaces)}")
        print("="*70)

        if show_gui:
            gmsh.fltk.run()

        gmsh.finalize()

        return {
            'mesh_file': output_mesh,
            'num_nodes': n_nodes,
            'num_elements': num_3d_elements,
            'num_fault_surfaces': len(fault_surfaces),
            'volume_tag': vol_tag_after
        }

    except Exception as e:
        print("\n" + "="*70)
        print("FAILED!")
        print("="*70)
        print(f"Error: {e}")
        print("="*70)
        gmsh.finalize()
        raise


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Strict Ando mesh generation')
    parser.add_argument('geo_file', nargs='?', default='kaikoura_v11.3hpk2.geo')
    parser.add_argument('--lc', type=float, default=3000, help='Fault mesh size (m)')
    parser.add_argument('--lc-coarse', type=float, default=15000, help='Coarse mesh size (m)')
    parser.add_argument('--box-length', type=float, default=150000, help='Box length (m)')
    parser.add_argument('--box-width', type=float, default=150000, help='Box width (m)')
    parser.add_argument('--box-depth', type=float, default=120000, help='Box depth (m)')
    parser.add_argument('-o', '--output', default='ando_fault_mesh.msh')
    parser.add_argument('--gui', action='store_true')

    args = parser.parse_args()

    try:
        result = generate_Ando_fault_strict(
            geo_file=args.geo_file,
            lc=args.lc,
            lc_coarse=args.lc_coarse,
            box_length=args.box_length,
            box_width=args.box_width,
            box_depth=args.box_depth,
            output_mesh=args.output,
            show_gui=args.gui
        )

        print("\n✓ Mesh generation complete!")
        sys.exit(0)

    except Exception as e:
        print(f"\n✗ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
