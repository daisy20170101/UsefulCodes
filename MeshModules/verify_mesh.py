#!/usr/bin/env python3
"""
Verify that a mesh file is valid and contains required physical groups

Usage: python verify_mesh.py ando_fault_mesh.msh
"""

import sys
import os


def verify_mesh(mesh_file):
    """Verify mesh file is valid"""

    print("="*70)
    print("MESH VERIFICATION")
    print("="*70)
    print(f"File: {mesh_file}\n")

    # Check file exists
    if not os.path.exists(mesh_file):
        print(f"✗ File not found: {mesh_file}")
        return False

    file_size = os.path.getsize(mesh_file) / (1024 * 1024)  # MB
    print(f"✓ File exists: {file_size:.2f} MB")

    # Check with Gmsh
    try:
        import gmsh
    except ImportError:
        print("✗ Gmsh Python module not available")
        print("  Install: pip install gmsh")
        return False

    gmsh.initialize()

    try:
        gmsh.open(mesh_file)
        print("✓ File opens successfully")

        # Get entities
        print("\nEntities in mesh:")
        for dim in range(4):
            entities = gmsh.model.getEntities(dim)
            dim_name = ["Points", "Lines", "Surfaces", "Volumes"][dim]
            print(f"  {dim_name}: {len(entities)}")

        # Get physical groups
        print("\nPhysical groups:")
        phys_groups = gmsh.model.getPhysicalGroups()

        if not phys_groups:
            print("  ✗ No physical groups found!")
            gmsh.finalize()
            return False

        has_volume = False
        for dim, tag in phys_groups:
            try:
                name = gmsh.model.getPhysicalName(dim, tag)
            except:
                name = f"unnamed_{tag}"

            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            dim_name = ["Point", "Line", "Surface", "Volume"][dim]

            print(f"  {dim_name} {tag}: '{name}' ({len(entities)} entities)")

            if dim == 3:
                has_volume = True

        # Check mesh statistics
        print("\nMesh statistics:")

        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        num_nodes = len(node_tags)
        print(f"  Nodes: {num_nodes:,}")

        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(3)
        if elem_tags and len(elem_tags) > 0:
            num_elements_3d = sum(len(tags) for tags in elem_tags)
            print(f"  3D elements: {num_elements_3d:,}")
        else:
            print(f"  3D elements: 0 (WARNING)")

        # Critical checks
        print("\nCritical checks:")

        if has_volume:
            print("  ✓ Volume physical group present")
        else:
            print("  ✗ Volume physical group MISSING")

        if num_nodes > 0:
            print(f"  ✓ Mesh contains nodes ({num_nodes:,})")
        else:
            print("  ✗ Mesh has no nodes")

        if elem_tags and len(elem_tags) > 0 and sum(len(tags) for tags in elem_tags) > 0:
            print(f"  ✓ Mesh contains 3D elements")
        else:
            print("  ✗ Mesh has no 3D elements")

        # Overall result
        print("\n" + "="*70)
        if has_volume and num_nodes > 0:
            print("VERIFICATION: PASSED ✓")
            print("="*70)
            gmsh.finalize()
            return True
        else:
            print("VERIFICATION: FAILED ✗")
            print("="*70)
            gmsh.finalize()
            return False

    except Exception as e:
        print(f"✗ Error reading mesh: {e}")
        gmsh.finalize()
        return False


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python verify_mesh.py <mesh_file.msh>")
        sys.exit(1)

    mesh_file = sys.argv[1]
    success = verify_mesh(mesh_file)

    sys.exit(0 if success else 1)
