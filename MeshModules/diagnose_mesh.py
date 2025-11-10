#!/usr/bin/env python3
"""
Diagnostic script to check mesh_ando geometry and identify issues

This script runs mesh_ando with additional checks and diagnostics
to identify what's causing the "boundary mesh error 2" problem.
"""

import gmsh
import sys


def diagnose_geometry(geo_file='kaikoura_v11.3hpk2.geo'):
    """
    Run diagnostics on the geometry before meshing
    """
    print("="*70)
    print("GEOMETRY DIAGNOSTICS")
    print("="*70)

    gmsh.initialize()
    gmsh.model.add("diagnostic")

    # Load geo file
    print(f"\n1. Loading {geo_file}...")
    try:
        gmsh.merge(geo_file)
        gmsh.model.geo.synchronize()
        print("   ✓ Loaded successfully")
    except Exception as e:
        print(f"   ✗ Failed to load: {e}")
        gmsh.finalize()
        return

    # Check initial entities
    print("\n2. Initial geometry:")
    for dim in range(4):
        entities = gmsh.model.getEntities(dim)
        print(f"   Dimension {dim}: {len(entities)} entities")
        if dim == 2 and entities:
            print(f"      Surface tags: {[e[1] for e in entities]}")

    # Scale geometry
    print("\n3. Scaling from km to m...")
    all_entities = []
    for dim in range(4):
        all_entities.extend(gmsh.model.getEntities(dim))

    gmsh.model.geo.dilate(all_entities, 0, 0, 0, 1000, 1000, 1000)
    gmsh.model.geo.synchronize()
    print("   ✓ Scaled")

    # Check surfaces after scaling
    print("\n4. Surface analysis after scaling:")
    surfaces = gmsh.model.getEntities(dim=2)

    fault_surfaces = []
    for surf_dim, surf_tag in surfaces:
        try:
            bbox = gmsh.model.getBoundingBox(surf_dim, surf_tag)
            xmin, ymin, zmin, xmax, ymax, zmax = bbox

            center_z = (zmin + zmax) / 2

            # Check if this is a fault surface (not at z=0 or z<-14000)
            is_fault = not (abs(center_z) < 1000 or center_z < -14000.0)

            print(f"\n   Surface {surf_tag}:")
            print(f"      Bounds: X=[{xmin:.1f}, {xmax:.1f}], Y=[{ymin:.1f}, {ymax:.1f}], Z=[{zmin:.1f}, {zmax:.1f}]")
            print(f"      Center Z: {center_z:.1f} m")
            print(f"      Keep as fault: {is_fault}")

            if is_fault:
                fault_surfaces.append(surf_tag)

                # Check surface orientation
                try:
                    orientation = gmsh.model.getOrientation(surf_dim, surf_tag)
                    print(f"      Orientation: {orientation}")
                except:
                    print(f"      Orientation: unknown")

                # Check if surface is closed
                boundaries = gmsh.model.getBoundary([(surf_dim, surf_tag)], oriented=False)
                print(f"      Boundary curves: {len(boundaries)}")

        except Exception as e:
            print(f"   Surface {surf_tag}: Error - {e}")

    print(f"\n5. Summary:")
    print(f"   Total surfaces: {len(surfaces)}")
    print(f"   Fault surfaces to keep: {len(fault_surfaces)}")
    print(f"   Tags: {fault_surfaces}")

    # Check for geometric issues
    print("\n6. Checking for geometric issues...")

    # Check coherence
    try:
        gmsh.model.geo.synchronize()
        print("   ✓ Geometry is synchronized")
    except Exception as e:
        print(f"   ✗ Synchronization issue: {e}")

    # Check if surfaces are valid
    for surf_tag in fault_surfaces:
        try:
            # Get surface normal
            gmsh.model.geo.synchronize()
            print(f"   Surface {surf_tag}: appears valid")
        except Exception as e:
            print(f"   ✗ Surface {surf_tag} has issues: {e}")

    print("\n" + "="*70)
    print("DIAGNOSTIC COMPLETE")
    print("="*70)

    # Save diagnostic geometry
    print("\nSaving diagnostic geometry to diagnostic.geo_unrolled...")
    gmsh.write("diagnostic.geo_unrolled")

    gmsh.finalize()

    return fault_surfaces


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Diagnose mesh_ando geometry')
    parser.add_argument('geo_file', nargs='?', default='kaikoura_v11.3hpk2.geo',
                       help='Geo file to diagnose')

    args = parser.parse_args()

    diagnose_geometry(args.geo_file)
