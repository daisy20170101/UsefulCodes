import gmsh
import os
import sys


def generate_Ando_fault(geo_file='kaikoura_v11.3hpk2.geo', lc=2000, lc_coarse=10e3,
                        distance_threshold=5000.0, box_length=100e3, box_width=100e3,
                        box_depth=100e3, output_mesh='ando_fault_mesh.msh',
                        use_distance_field=False, show_gui=False, verbose=True):
    """
    Generate Ando fault mesh from a .geo file with embedded fault surfaces.

    This function:
    1. Loads fault geometry from .geo file (coordinates in km)
    2. Converts coordinates from km to meters
    3. Removes surfaces based on mass center z-coordinate:
       - Surfaces with mass center at z≈0 (within 1000m tolerance)
       - Surfaces with mass center at z<-14000m
    4. Creates a box domain encompassing the fault with specified padding
    5. Embeds fault surfaces into the volume using gmsh.model.mesh.embed()
    6. Generates 3D tetrahedral mesh with refinement near faults

    Physical Groups:
    - Group 1: Volume (3D tetrahedral mesh)
    - Group 101: Top box surface at z=z_top (includes fault trace)
    - Group 103: Fault surfaces (embedded internal boundaries)
    - Group 105: Other 5 box surfaces (bottom, front, right, back, left)

    Parameters:
    -----------
    geo_file : str
        Path to the .geo file (coordinates should be in km)
    lc : float
        Mesh size on fault surfaces in meters (default: 2000)
    lc_coarse : float
        Coarse mesh size away from fault in meters (default: 10000)
    distance_threshold : float
        Distance from fault where coarsening begins in meters (default: 5000)
    box_length : float
        Extra padding in X direction beyond fault extent in meters (default: 100000)
    box_width : float
        Extra padding in Y direction beyond fault extent in meters (default: 100000)
    box_depth : float
        Depth below fault bottom in meters (default: 100000)
    output_mesh : str
        Output mesh filename (default: 'ando_fault_mesh.msh')
    use_distance_field : bool
        Use distance-based mesh refinement (slower but smoother). Default: False
    show_gui : bool
        Launch Gmsh GUI after meshing (default: False)
    verbose : bool
        Print detailed progress messages (default: True)

    Returns:
    --------
    dict
        Dictionary containing:
        - 'mesh_file': Path to generated mesh file
        - 'num_nodes': Number of mesh nodes
        - 'num_elements': Number of 3D elements
        - 'num_fault_surfaces': Number of embedded fault surfaces
        - 'box_bounds': Bounding box coordinates (xmin, xmax, ymin, ymax, zmin, zmax)

    Raises:
    -------
    FileNotFoundError
        If geo_file doesn't exist
    ValueError
        If no fault surfaces remain after filtering
    Exception
        If mesh generation fails with all algorithms

    Example:
    --------
    >>> result = generate_Ando_fault(
    ...     geo_file='kaikoura_fault.geo',
    ...     lc=2000,
    ...     lc_coarse=10000,
    ...     box_depth=100000,
    ...     output_mesh='kaikoura_mesh.msh'
    ... )
    >>> print(f"Mesh saved with {result['num_nodes']} nodes")
    """

    def log(message):
        """Print message if verbose mode is enabled"""
        if verbose:
            print(message)

    # ========================================
    # 1. Input Validation
    # ========================================
    if not os.path.exists(geo_file):
        raise FileNotFoundError(f"Geo file not found: {geo_file}")

    if lc <= 0 or lc_coarse <= 0:
        raise ValueError("Mesh sizes (lc, lc_coarse) must be positive")

    if lc >= lc_coarse:
        log("Warning: lc >= lc_coarse. Setting lc_coarse = lc * 5")
        lc_coarse = lc * 5

    log("="*70)
    log("ANDO FAULT MESH GENERATION")
    log("="*70)
    log(f"Input file: {geo_file}")
    log(f"Mesh sizes: {lc}m (fault) to {lc_coarse}m (coarse)")
    log(f"Box padding: {box_length/1e3:.1f} km × {box_width/1e3:.1f} km, depth: {box_depth/1e3:.1f} km")
    log(f"Output: {output_mesh}")
    log("="*70)

    # ========================================
    # 2. Initialize Gmsh
    # ========================================
    gmsh.initialize()
    gmsh.model.add("ando_fault")

    # Set performance options
    log("\nConfiguring Gmsh options...")
    gmsh.option.setNumber("General.NumThreads", 0)  # Use all cores
    gmsh.option.setNumber("General.Verbosity", 5 if verbose else 2)
    gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay for 2D
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Start with Delaunay (will try others if needed)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MaxNumThreads3D", 0)

    # ========================================
    # 3. Load and Process Geometry
    # ========================================
    log(f"\nLoading geo file: {geo_file}")
    gmsh.merge(geo_file)
    gmsh.model.geo.synchronize()

    # Remove old physical groups from loaded file
    log("Removing old physical groups from geo file...")
    try:
        all_physical_groups = gmsh.model.getPhysicalGroups()
        if all_physical_groups:
            for dim, tag in all_physical_groups:
                try:
                    gmsh.model.removePhysicalGroups([(dim, tag)])
                except:
                    pass
            log(f"  Removed {len(all_physical_groups)} old physical groups")
    except Exception as e:
        log(f"  Note: {e}")

    # Convert coordinates from km to m
    log("\nConverting coordinates from km to m (scaling by 1000)...")
    all_entities = []
    for dim in range(4):
        entities = gmsh.model.getEntities(dim)
        all_entities.extend(entities)

    gmsh.model.geo.dilate(all_entities, 0, 0, 0, 1000, 1000, 1000)
    gmsh.model.geo.synchronize()
    log("  ✓ Coordinates scaled to meters")

    # ========================================
    # 4. Filter Surfaces
    # ========================================
    log("\n" + "="*70)
    log("FILTERING FAULT SURFACES")
    log("="*70)

    surfaces = gmsh.model.getEntities(dim=2)
    log(f"Total surfaces found: {len(surfaces)}")
    log(f"Surface tags: {[s[1] for s in surfaces]}")

    surfaces_to_remove = []
    surfaces_to_keep = []

    for surf_dim, surf_tag in surfaces:
        try:
            bbox = gmsh.model.getBoundingBox(surf_dim, surf_tag)
            xmin, ymin, zmin, xmax, ymax, zmax = bbox

            center_x = (xmin + xmax) / 2
            center_y = (ymin + ymax) / 2
            center_z = (zmin + zmax) / 2

            log(f"\nSurface {surf_tag}:")
            log(f"  Mass center: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f}) m")

            # Removal criteria
            if abs(center_z) < 1000:  # z ≈ 0 (tolerance 1000m)
                log(f"  → Removing (z ≈ 0)")
                surfaces_to_remove.append((surf_dim, surf_tag))
            elif center_z < -14000.0:
                log(f"  → Removing (z < -14000m)")
                surfaces_to_remove.append((surf_dim, surf_tag))
            else:
                log(f"  → Keeping")
                surfaces_to_keep.append((surf_dim, surf_tag))

        except Exception as e:
            log(f"  Warning: Could not process surface {surf_tag}: {e}")
            surfaces_to_keep.append((surf_dim, surf_tag))

    # Remove filtered surfaces
    if surfaces_to_remove:
        log(f"\nRemoving {len(surfaces_to_remove)} surfaces: {[s[1] for s in surfaces_to_remove]}")
        try:
            gmsh.model.geo.remove(surfaces_to_remove)
            log("  ✓ Surfaces removed")
        except Exception as e:
            log(f"  Warning: {e}")
    else:
        log("\nNo surfaces matched removal criteria")

    gmsh.model.geo.synchronize()

    # Get remaining fault surfaces
    remaining_surfaces = gmsh.model.getEntities(dim=2)
    fault_surface_tags = [s[1] for s in remaining_surfaces]
    log(f"\n✓ After filtering: {len(fault_surface_tags)} fault surface(s) remain")
    log(f"  Tags: {fault_surface_tags}")

    if not fault_surface_tags:
        raise ValueError("No fault surfaces remain after filtering!")

    # ========================================
    # 5. Calculate Fault Bounding Box
    # ========================================
    log("\n" + "="*70)
    log("CALCULATING FAULT BOUNDS")
    log("="*70)

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

    log(f"Fault extent:")
    log(f"  X: [{fault_xmin:.2f}, {fault_xmax:.2f}] m ({(fault_xmax-fault_xmin)/1e3:.1f} km)")
    log(f"  Y: [{fault_ymin:.2f}, {fault_ymax:.2f}] m ({(fault_ymax-fault_ymin)/1e3:.1f} km)")
    log(f"  Z: [{fault_zmin:.2f}, {fault_zmax:.2f}] m ({(fault_zmax-fault_zmin)/1e3:.1f} km)")

    # ========================================
    # 6. Create Box Domain
    # ========================================
    log("\n" + "="*70)
    log("CREATING BOX DOMAIN")
    log("="*70)

    # Add buffer to prevent geometric intersections
    buffer_x = max(5000.0, (fault_xmax - fault_xmin) * 0.05)
    buffer_y = max(5000.0, (fault_ymax - fault_ymin) * 0.05)
    buffer_z = max(2500.0, (fault_zmax - fault_zmin) * 0.05)

    log(f"Buffer: X={buffer_x/1e3:.1f} km, Y={buffer_y/1e3:.1f} km, Z={buffer_z/1e3:.1f} km")

    # Box bounds
    x0 = fault_xmin - box_length / 2 - buffer_x
    x1 = fault_xmax + box_length / 2 + buffer_x
    y0 = fault_ymin - box_width / 2 - buffer_y
    y1 = fault_ymax + box_width / 2 + buffer_y
    z0 = -box_depth
    z1 = 500.0  # Small padding above fault

    actual_box_length = x1 - x0
    actual_box_width = y1 - y0
    actual_box_depth = z1 - z0

    log(f"Box dimensions:")
    log(f"  {actual_box_length/1e3:.1f} km × {actual_box_width/1e3:.1f} km × {actual_box_depth/1e3:.1f} km")
    log(f"  X: [{x0:.2f}, {x1:.2f}]")
    log(f"  Y: [{y0:.2f}, {y1:.2f}]")
    log(f"  Z: [{z0:.2f}, {z1:.2f}]")

    # Create box geometry
    log("\nCreating box geometry...")

    # Box corners
    p1 = gmsh.model.geo.addPoint(x0, y0, z0, lc_coarse)
    p2 = gmsh.model.geo.addPoint(x1, y0, z0, lc_coarse)
    p3 = gmsh.model.geo.addPoint(x1, y1, z0, lc_coarse)
    p4 = gmsh.model.geo.addPoint(x0, y1, z0, lc_coarse)
    p5 = gmsh.model.geo.addPoint(x0, y0, z1, lc_coarse)
    p6 = gmsh.model.geo.addPoint(x1, y0, z1, lc_coarse)
    p7 = gmsh.model.geo.addPoint(x1, y1, z1, lc_coarse)
    p8 = gmsh.model.geo.addPoint(x0, y1, z1, lc_coarse)

    # Box lines
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

    # Box faces
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

    log(f"  ✓ Box surfaces: bottom={surf_bottom}, top={surf_top}, front={surf_front}")
    log(f"                  right={surf_right}, back={surf_back}, left={surf_left}")

    # Create volume
    box_surfaces = [surf_bottom, surf_top, surf_front, surf_right, surf_back, surf_left]
    box_surface_loop = gmsh.model.geo.addSurfaceLoop(box_surfaces)
    box_volume = gmsh.model.geo.addVolume([box_surface_loop])

    log(f"  ✓ Box volume created: {box_volume}")

    gmsh.model.geo.synchronize()

    # ========================================
    # 7. Embed Fault Surfaces
    # ========================================
    log("\n" + "="*70)
    log("EMBEDDING FAULT SURFACES")
    log("="*70)

    log(f"Embedding {len(fault_surface_tags)} fault surface(s) into volume {box_volume}...")

    for i, fault_surf in enumerate(fault_surface_tags):
        try:
            gmsh.model.mesh.embed(2, [fault_surf], 3, box_volume)
            log(f"  ✓ Embedded fault surface {fault_surf}")
        except Exception as e:
            log(f"  ✗ Failed to embed surface {fault_surf}: {e}")
            raise

    log("✓ All fault surfaces embedded successfully")

    # ========================================
    # 8. Set Mesh Sizes
    # ========================================
    log("\n" + "="*70)
    log("CONFIGURING MESH SIZES")
    log("="*70)

    # Set fine mesh on fault surfaces
    log(f"Setting mesh size on fault surfaces: lc = {lc} m")
    total_fault_points = 0

    for surf_tag in fault_surface_tags:
        surf_bounds = gmsh.model.getBoundary([(2, surf_tag)], combined=False,
                                             oriented=False, recursive=True)
        surf_points = [abs(b[1]) for b in surf_bounds if b[0] == 0]
        total_fault_points += len(surf_points)

        for pt in surf_points:
            gmsh.model.mesh.setSize([(0, pt)], lc)

    log(f"  ✓ Set mesh size on {total_fault_points} fault points")

    # Configure mesh refinement strategy
    if use_distance_field:
        log(f"\nUsing distance-based mesh refinement:")
        log(f"  Near fault: {lc} m")
        log(f"  Transition: {distance_threshold} m to {distance_threshold*10} m")
        log(f"  Far field: {lc_coarse} m")

        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", fault_surface_tags)

        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", lc)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_coarse)
        gmsh.model.mesh.field.setNumber(2, "DistMin", distance_threshold)
        gmsh.model.mesh.field.setNumber(2, "DistMax", distance_threshold * 10)

        gmsh.model.mesh.field.setAsBackgroundMesh(2)
    else:
        log(f"\nUsing simple mesh size control:")
        log(f"  Min: {lc} m")
        log(f"  Max: {lc_coarse} m")

        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc_coarse)

    # ========================================
    # 9. Add Physical Groups
    # ========================================
    log("\n" + "="*70)
    log("ADDING PHYSICAL GROUPS")
    log("="*70)

    # Volume
    gmsh.model.addPhysicalGroup(3, [box_volume], 1)
    gmsh.model.setPhysicalName(3, 1, "volume")
    log(f"  ✓ Physical group 1: Volume (tag {box_volume})")

    # Individual fault surface physical groups (tags starting from 165)
    log(f"\n  Adding individual fault physical groups (tags 165+):")
    for i, fault_surf in enumerate(fault_surface_tags):
        phys_tag = 165 + i
        gmsh.model.addPhysicalGroup(2, [fault_surf], tag=phys_tag)
        gmsh.model.setPhysicalName(2, phys_tag, f"fault_{i+1}")
        log(f"    Physical Surface {phys_tag}: fault_{i+1} (surface {fault_surf})")

    # Combined fault group for backward compatibility
    gmsh.model.addPhysicalGroup(2, fault_surface_tags, 103)
    gmsh.model.setPhysicalName(2, 103, "fault_all")
    log(f"  ✓ Physical group 103: fault_all ({len(fault_surface_tags)} surfaces combined)")

    # Top surface
    gmsh.model.addPhysicalGroup(2, [surf_top], 101)
    gmsh.model.setPhysicalName(2, 101, "box_top")
    log(f"  ✓ Physical group 101: Top surface")

    # Other box surfaces
    other_surfaces = [surf_bottom, surf_front, surf_right, surf_back, surf_left]
    gmsh.model.addPhysicalGroup(2, other_surfaces, 105)
    gmsh.model.setPhysicalName(2, 105, "box_sides")
    log(f"  ✓ Physical group 105: Box sides ({len(other_surfaces)} surfaces)")

    # ========================================
    # 10. Generate Mesh
    # ========================================
    log("\n" + "="*70)
    log("GENERATING MESH")
    log("="*70)

    # Generate 2D mesh
    log("Generating 2D mesh...")
    try:
        gmsh.model.mesh.generate(2)
        num_2d_nodes = len(gmsh.model.mesh.getNodes()[0])
        log(f"  ✓ 2D mesh generated ({num_2d_nodes} nodes)")
    except Exception as e:
        log(f"  Warning: 2D mesh generation issue: {e}")

    # Generate 3D mesh with fallback algorithms
    log("\nGenerating 3D mesh (this may take several minutes)...")

    algorithms = [
        (1, "Delaunay"),
        (4, "Frontal"),
        (10, "HXT")
    ]

    mesh_success = False
    last_error = None
    successful_algo = None

    for algo_id, algo_name in algorithms:
        try:
            log(f"\n  Trying {algo_name} algorithm (id={algo_id})...")
            gmsh.option.setNumber("Mesh.Algorithm3D", algo_id)

            # Clear previous 3D mesh attempts
            try:
                gmsh.model.mesh.clear(3)
            except:
                pass

            # Generate 3D mesh
            gmsh.model.mesh.generate(3)

            log(f"  ✓ Successfully generated 3D mesh using {algo_name}")
            mesh_success = True
            successful_algo = algo_name
            break

        except Exception as e:
            error_msg = str(e)
            log(f"  ✗ {algo_name} failed: {error_msg[:100]}")
            last_error = error_msg
            continue

    if not mesh_success:
        error_msg = f"\nAll 3D meshing algorithms failed!\n"
        error_msg += f"Last error: {last_error}\n\n"
        error_msg += "Suggestions:\n"
        error_msg += f"  1. Increase mesh sizes: lc={lc*2}, lc_coarse={lc_coarse*2}\n"
        error_msg += f"  2. Reduce box size or check fault geometry\n"
        error_msg += f"  3. Check that fault surfaces are valid (not self-intersecting)\n"
        raise Exception(error_msg)

    # ========================================
    # 11. Get Statistics and Save
    # ========================================
    log("\n" + "="*70)
    log("MESH STATISTICS")
    log("="*70)

    num_nodes = len(gmsh.model.mesh.getNodes()[0])
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(3)
    num_elements = sum(len(tags) for tags in elem_tags) if elem_tags else 0

    log(f"Algorithm used: {successful_algo}")
    log(f"Nodes: {num_nodes:,}")
    log(f"3D elements: {num_elements:,}")
    log(f"Embedded fault surfaces: {len(fault_surface_tags)}")

    # Write mesh
    log(f"\nWriting mesh to: {output_mesh}")
    gmsh.write(output_mesh)
    log("  ✓ Mesh file saved")

    # Prepare result dictionary
    result = {
        'mesh_file': output_mesh,
        'num_nodes': num_nodes,
        'num_elements': num_elements,
        'num_fault_surfaces': len(fault_surface_tags),
        'box_bounds': (x0, x1, y0, y1, z0, z1),
        'algorithm': successful_algo
    }

    # Show GUI if requested
    if show_gui:
        log("\nLaunching Gmsh GUI...")
        gmsh.fltk.run()

    # Finalize
    gmsh.finalize()

    log("\n" + "="*70)
    log("SUCCESS!")
    log("="*70)

    return result


if __name__ == "__main__":
    """
    Example usage and command-line interface
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate 3D mesh with embedded fault surfaces from .geo file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python mesh_ando_fixed.py kaikoura_fault.geo

  # Custom mesh sizes
  python mesh_ando_fixed.py fault.geo --lc 1000 --lc-coarse 5000

  # Custom box dimensions
  python mesh_ando_fixed.py fault.geo --box-length 150000 --box-depth 80000

  # With distance field and GUI
  python mesh_ando_fixed.py fault.geo --distance-field --gui
        """
    )

    parser.add_argument('geo_file', type=str,
                       help='Input .geo file (coordinates in km)')

    parser.add_argument('--lc', type=float, default=2000,
                       help='Mesh size on fault surface in meters (default: 2000)')

    parser.add_argument('--lc-coarse', type=float, default=10000,
                       help='Coarse mesh size in meters (default: 10000)')

    parser.add_argument('--distance-threshold', type=float, default=5000,
                       help='Distance from fault for coarsening in meters (default: 5000)')

    parser.add_argument('--box-length', type=float, default=100000,
                       help='Box padding in X direction in meters (default: 100000)')

    parser.add_argument('--box-width', type=float, default=100000,
                       help='Box padding in Y direction in meters (default: 100000)')

    parser.add_argument('--box-depth', type=float, default=100000,
                       help='Box depth in meters (default: 100000)')

    parser.add_argument('-o', '--output', type=str, default='ando_fault_mesh.msh',
                       help='Output mesh file (default: ando_fault_mesh.msh)')

    parser.add_argument('--distance-field', action='store_true',
                       help='Use distance-based mesh refinement (slower)')

    parser.add_argument('--gui', action='store_true',
                       help='Launch Gmsh GUI after meshing')

    parser.add_argument('--quiet', action='store_true',
                       help='Suppress detailed output')

    args = parser.parse_args()

    # Run mesh generation
    try:
        result = generate_Ando_fault(
            geo_file=args.geo_file,
            lc=args.lc,
            lc_coarse=args.lc_coarse,
            distance_threshold=args.distance_threshold,
            box_length=args.box_length,
            box_width=args.box_width,
            box_depth=args.box_depth,
            output_mesh=args.output,
            use_distance_field=args.distance_field,
            show_gui=args.gui,
            verbose=not args.quiet
        )

        print("\n" + "="*70)
        print("MESH GENERATION COMPLETE")
        print("="*70)
        print(f"Output file: {result['mesh_file']}")
        print(f"Nodes: {result['num_nodes']:,}")
        print(f"Elements: {result['num_elements']:,}")
        print(f"Fault surfaces: {result['num_fault_surfaces']}")
        print(f"Algorithm: {result['algorithm']}")
        print("="*70)

        sys.exit(0)

    except Exception as e:
        print(f"\n{'='*70}")
        print("ERROR")
        print("="*70)
        print(f"{e}")
        print("="*70)
        sys.exit(1)
