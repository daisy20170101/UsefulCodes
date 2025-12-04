"""
Generate 3D mesh with fault surfaces by loading pre-meshed fault segments
This approach loads individual .msh files for each fault segment and merges them
"""

import gmsh
import numpy as np
import os


def generate_mesh_from_fault_segments(segment_files,
                                      lc_coarse=10e3,
                                      box_length=100e3,
                                      box_width=100e3,
                                      box_depth=100e3,
                                      output_mesh='ando_fault_mesh_merged.msh',
                                      use_distance_field=True,
                                      distance_threshold=5000.0):
    """
    Generate 3D mesh by loading and merging pre-meshed fault segments

    Parameters:
    -----------
    segment_files : list of str
        List of .msh files for individual fault segments
    lc_coarse : float
        Coarse mesh size away from faults (meters)
    box_length : float
        Extra padding in X direction (meters)
    box_width : float
        Extra padding in Y direction (meters)
    box_depth : float
        Depth below fault (meters)
    output_mesh : str
        Output mesh filename
    use_distance_field : bool
        Use distance-based refinement
    distance_threshold : float
        Distance for mesh coarsening (meters)

    Returns:
    --------
    output_mesh : str
        Path to generated mesh file
    """

    print("=" * 70)
    print("GENERATING 3D MESH FROM FAULT SEGMENTS")
    print("=" * 70)
    print(f"Number of fault segment files: {len(segment_files)}")
    print(f"Coarse mesh size: {lc_coarse} m")
    print(f"Output mesh: {output_mesh}")
    print()

    # Initialize gmsh
    gmsh.initialize()
    gmsh.model.add("ando_fault_merged")

    # Set performance options
    gmsh.option.setNumber("General.NumThreads", 0)
    gmsh.option.setNumber("Mesh.MaxNumThreads3D", 0)

    # Track all fault surface tags
    fault_surface_tags = []
    all_fault_nodes = []

    # Load each fault segment mesh and extract geometry
    print("Loading fault segment meshes...")
    for idx, segment_file in enumerate(segment_files):
        print(f"\n  [{idx + 1}/{len(segment_files)}] Loading: {os.path.basename(segment_file)}")

        if not os.path.exists(segment_file):
            print(f"    ERROR: File not found!")
            continue

        try:
            # Merge the mesh file into current model
            gmsh.merge(segment_file)
            gmsh.model.geo.synchronize()

            # Get surfaces that were just added
            surfaces = gmsh.model.getEntities(dim=2)
            print(f"    Loaded {len(surfaces)} surface(s): {[s[1] for s in surfaces]}")

            # Track the new surface tags
            for surf_dim, surf_tag in surfaces:
                if surf_tag not in fault_surface_tags:
                    fault_surface_tags.append(surf_tag)

                    # Get bounding box
                    bbox = gmsh.model.getBoundingBox(surf_dim, surf_tag)
                    xmin, ymin, zmin, xmax, ymax, zmax = bbox
                    print(f"      Surface {surf_tag}: bbox=[({xmin:.0f}, {ymin:.0f}, {zmin:.0f}) to ({xmax:.0f}, {ymax:.0f}, {zmax:.0f})]")

        except Exception as e:
            print(f"    ERROR loading segment: {e}")
            continue

    if not fault_surface_tags:
        raise ValueError("No fault surfaces loaded! Check segment files.")

    print(f"\nTotal fault surfaces loaded: {len(fault_surface_tags)}")
    print(f"Fault surface tags: {fault_surface_tags}")

    # Calculate combined bounding box for all fault surfaces
    print("\nCalculating combined bounding box...")
    fault_xmin, fault_ymin, fault_zmin = float('inf'), float('inf'), float('inf')
    fault_xmax, fault_ymax, fault_zmax = float('-inf'), float('-inf'), float('-inf')

    for surf_tag in fault_surface_tags:
        bbox = gmsh.model.getBoundingBox(2, surf_tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox
        fault_xmin = min(fault_xmin, xmin)
        fault_ymin = min(fault_ymin, ymin)
        fault_zmin = min(fault_zmin, zmin)
        fault_xmax = max(fault_xmax, xmax)
        fault_ymax = max(fault_ymax, ymax)
        fault_zmax = max(fault_zmax, zmax)

    print(f"Combined fault bounds:")
    print(f"  X: [{fault_xmin:.1f}, {fault_xmax:.1f}] m (extent: {(fault_xmax-fault_xmin)/1e3:.1f} km)")
    print(f"  Y: [{fault_ymin:.1f}, {fault_ymax:.1f}] m (extent: {(fault_ymax-fault_ymin)/1e3:.1f} km)")
    print(f"  Z: [{fault_zmin:.1f}, {fault_zmax:.1f}] m (extent: {(fault_zmax-fault_zmin)/1e3:.1f} km)")

    # Create box domain
    print("\n" + "=" * 70)
    print("CREATING BOX DOMAIN")
    print("=" * 70)

    # Add buffer to prevent geometric intersections
    buffer_x = max(5000.0, (fault_xmax - fault_xmin) * 0.05)
    buffer_y = max(5000.0, (fault_ymax - fault_ymin) * 0.05)
    buffer_z = max(2500.0, (fault_zmax - fault_zmin) * 0.05)

    print(f"Buffer: X={buffer_x/1e3:.1f} km, Y={buffer_y/1e3:.1f} km, Z={buffer_z/1e3:.1f} km")

    # Box bounds
    x0 = fault_xmin - box_length / 2 - buffer_x
    x1 = fault_xmax + box_length / 2 + buffer_x
    y0 = fault_ymin - box_width / 2 - buffer_y
    y1 = fault_ymax + box_width / 2 + buffer_y
    z0 = -box_depth
    z1 = 2000.0

    print(f"Box size: {(x1-x0)/1e3:.1f} × {(y1-y0)/1e3:.1f} × {(z1-z0)/1e3:.1f} km")

    # Create box points
    p1 = gmsh.model.geo.addPoint(x0, y0, z0, lc_coarse)
    p2 = gmsh.model.geo.addPoint(x1, y0, z0, lc_coarse)
    p3 = gmsh.model.geo.addPoint(x1, y1, z0, lc_coarse)
    p4 = gmsh.model.geo.addPoint(x0, y1, z0, lc_coarse)
    p5 = gmsh.model.geo.addPoint(x0, y0, z1, lc_coarse)
    p6 = gmsh.model.geo.addPoint(x1, y0, z1, lc_coarse)
    p7 = gmsh.model.geo.addPoint(x1, y1, z1, lc_coarse)
    p8 = gmsh.model.geo.addPoint(x0, y1, z1, lc_coarse)

    # Create box edges
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

    # Create box surfaces
    cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
    cl3 = gmsh.model.geo.addCurveLoop([l1, l10, -l5, -l9])
    cl4 = gmsh.model.geo.addCurveLoop([l2, l11, -l6, -l10])
    cl5 = gmsh.model.geo.addCurveLoop([l3, l12, -l7, -l11])
    cl6 = gmsh.model.geo.addCurveLoop([l4, l9, -l8, -l12])

    surf_bottom = gmsh.model.geo.addPlaneSurface([cl1])
    surf_top = gmsh.model.geo.addPlaneSurface([cl2])
    surf_front = gmsh.model.geo.addPlaneSurface([cl3])
    surf_right = gmsh.model.geo.addPlaneSurface([cl4])
    surf_back = gmsh.model.geo.addPlaneSurface([cl5])
    surf_left = gmsh.model.geo.addPlaneSurface([cl6])

    print(f"Box surfaces: bottom={surf_bottom}, top={surf_top}, front={surf_front}, right={surf_right}, back={surf_back}, left={surf_left}")

    gmsh.model.geo.synchronize()

    # Create box volume with fault surfaces
    print("\nCreating box volume with fault surfaces...")
    print(f"  Box surfaces: {[surf_bottom, surf_top, surf_front, surf_right, surf_back, surf_left]}")
    print(f"  Fault surfaces: {fault_surface_tags}")

    # IMPORTANT: Internal fault surfaces that don't connect to box boundaries
    # CANNOT be included in a surface loop (they create "open elements").
    # Therefore, we must use the embedding approach.

    # Create box volume from the 6 box surfaces only
    box_surface_loop = gmsh.model.geo.addSurfaceLoop([surf_bottom, surf_top, surf_front, surf_right, surf_back, surf_left])
    box_volume = gmsh.model.geo.addVolume([box_surface_loop])
    print(f"\nBox volume created: {box_volume}")

    gmsh.model.geo.synchronize()

    # Embed fault surfaces in the volume
    # This tells the mesher to conform the 3D mesh to these internal surfaces
    print(f"\nEmbedding {len(fault_surface_tags)} fault surface(s) in box volume...")
    print(f"  Fault surfaces: {fault_surface_tags}")

    gmsh.model.mesh.embed(2, fault_surface_tags, 3, box_volume)
    print(f"✓ Successfully embedded fault surfaces in volume {box_volume}")

    # Verify volume exists
    volumes_after_embed = gmsh.model.getEntities(dim=3)
    print(f"Volumes after embedding: {volumes_after_embed}")

    # Mesh generation
    print("\n" + "=" * 70)
    print("GENERATING MESH")
    print("=" * 70)

    # Set mesh size control
    if use_distance_field:
        print("Using distance-based mesh refinement...")
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", fault_surface_tags)

        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", 2000)  # Fine near faults
        gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_coarse)
        gmsh.model.mesh.field.setNumber(2, "DistMin", distance_threshold)
        gmsh.model.mesh.field.setNumber(2, "DistMax", distance_threshold * 10)

        gmsh.model.mesh.field.setAsBackgroundMesh(2)
    else:
        print("Using simple mesh size control...")
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 2000)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_coarse)

    # Generate 2D mesh
    print("\nGenerating 2D mesh...")
    try:
        gmsh.model.mesh.generate(2)
        num_2d_nodes = len(gmsh.model.mesh.getNodes()[0])
        print(f"  2D mesh generated: {num_2d_nodes} nodes")
    except Exception as e:
        print(f"  Warning: 2D mesh generation issues: {e}")

    # Generate 3D mesh
    print("\nGenerating 3D mesh...")
    algorithms = [
        (4, "Frontal"),
        (1, "Delaunay"),
        (10, "HXT"),
    ]

    mesh_success = False
    for algo_id, algo_name in algorithms:
        try:
            print(f"  Trying {algo_name} algorithm (id={algo_id})...")
            gmsh.option.setNumber("Mesh.Algorithm3D", algo_id)
            try:
                gmsh.model.mesh.clear(3)
            except:
                pass
            gmsh.model.mesh.generate(3)
            print(f"  ✓ Successfully generated 3D mesh using {algo_name}")
            mesh_success = True
            break
        except Exception as e:
            print(f"  ✗ {algo_name} failed: {e}")
            continue

    if not mesh_success:
        raise Exception("All 3D meshing algorithms failed!")

    # Remove duplicate nodes
    print("\nRemoving duplicate nodes...")
    try:
        node_tags_before, _, _ = gmsh.model.mesh.getNodes()
        gmsh.model.mesh.removeDuplicateNodes()
        node_tags_after, _, _ = gmsh.model.mesh.getNodes()
        num_removed = len(node_tags_before) - len(node_tags_after)
        if num_removed > 0:
            print(f"  Removed {num_removed} duplicate nodes")
        else:
            print(f"  No duplicate nodes found")
    except Exception as e:
        print(f"  Warning: {e}")

    # Assign physical groups
    print("\n" + "=" * 70)
    print("ASSIGNING PHYSICAL GROUPS")
    print("=" * 70)

    gmsh.option.setNumber("Mesh.SaveAll", 1)

    # Identify box surfaces vs fault surfaces
    all_surfaces = gmsh.model.getEntities(dim=2)
    box_surfaces = [surf_bottom, surf_top, surf_front, surf_right, surf_back, surf_left]

    freesurf = [surf_top]
    absorb = [s for s in box_surfaces if s != surf_top]

    print(f"\nSurface classification:")
    print(f"  Free surface (top): {freesurf}")
    print(f"  Fault surfaces: {fault_surface_tags}")
    print(f"  Absorbing boundaries: {absorb}")

    # Add physical groups
    gmsh.model.addPhysicalGroup(2, freesurf, 101)
    gmsh.model.setPhysicalName(2, 101, "freesurf")
    print(f"  ✓ Added group 101 (free surface): {len(freesurf)} surfaces")

    gmsh.model.addPhysicalGroup(2, fault_surface_tags, 103)
    gmsh.model.setPhysicalName(2, 103, "fault")
    print(f"  ✓ Added group 103 (fault): {len(fault_surface_tags)} surfaces")

    gmsh.model.addPhysicalGroup(2, absorb, 105)
    gmsh.model.setPhysicalName(2, 105, "absorb")
    print(f"  ✓ Added group 105 (absorbing): {len(absorb)} surfaces")

    all_volumes = gmsh.model.getEntities(dim=3)

    print('all volumes:',all_volumes)

    volume_tag = all_volumes[0][1]
    gmsh.model.addPhysicalGroup(3, [volume_tag], 1)
    gmsh.model.setPhysicalName(3, 1, "volume")
    print(f"  ✓ Added group 1 (volume): {volume_tag}")

    # Mesh statistics
    print("\n" + "=" * 70)
    print("MESH STATISTICS")
    print("=" * 70)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    element_types, element_tags_list, _ = gmsh.model.mesh.getElements()

    print(f"Nodes: {len(node_tags)}")

    total_elements = 0
    for i, elem_type in enumerate(element_types):
        elem_name = gmsh.model.mesh.getElementProperties(elem_type)[0]
        num_elems = len(element_tags_list[i])
        total_elements += num_elems
        print(f"  {elem_name}: {num_elems} elements")

    print(f"Total elements: {total_elements}")

    # Remove lower-dimensional elements
    print("\nRemoving lower-dimensional elements (edges and points)...")
    edges_1d = gmsh.model.getEntities(dim=1)
    if edges_1d:
        for dim, tag in edges_1d:
            try:
                gmsh.model.mesh.clear(dim, tag)
            except:
                pass
        print(f"  Removed 1D edge elements")

    points_0d = gmsh.model.getEntities(dim=0)
    if points_0d:
        for dim, tag in points_0d:
            try:
                gmsh.model.mesh.clear(dim, tag)
            except:
                pass
        print(f"  Removed 0D point elements")

    gmsh.option.setNumber("Mesh.SaveAll", 0)

    # Write mesh
    print("\n" + "=" * 70)
    print("WRITING MESH FILE")
    print("=" * 70)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(output_mesh)
    print(f"Mesh saved to: {output_mesh}")
    print(f"Format: MSH 2.2")

    # Finalize
    gmsh.finalize()

    return output_mesh


if __name__ == "__main__":
    import glob

    # Find all fault segment mesh files
    segment_dir = "fault_segments"
    segment_files = sorted(glob.glob(os.path.join(segment_dir, "fault_segment_*.msh")))

    if not segment_files:
        print(f"ERROR: No fault segment files found in {segment_dir}/")
        print("Run mesh_fault_segments.py first to generate individual segment meshes.")
    else:
        print(f"Found {len(segment_files)} fault segment files")

        output_mesh = generate_mesh_from_fault_segments(
            segment_files=segment_files,
            lc_coarse=10000,
            box_length=100e3,
            box_width=100e3,
            box_depth=100e3,
            output_mesh='ando_fault_mesh_merged.msh',
            use_distance_field=False
        )

        print(f"\n✓ 3D mesh generation complete!")
        print(f"✓ Output: {output_mesh}")
