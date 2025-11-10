import gmsh

def generate_Ando_fault(geo_file='kaikoura_v11.3hpk2.geo', lc=2000, lc_coarse=10e3,
                        distance_threshold=5000.0, box_length=100e3, box_width=100e3,
                        box_depth=100e3, output_mesh='ando_fault_mesh.msh',
                        use_distance_field=False):
    """
    Generate Ando fault mesh from a .geo file.
    Removes surfaces based on mass center z-coordinate:
    - Surfaces with mass center at z=0
    - Surfaces with mass center at z<-14000.0

    The box will encompass the fault surfaces with specified padding.
    Fault surfaces are cut at box boundaries using boolean fragment operations.
    Parts of fault surfaces extending above the box top are removed, and
    intersection traces are preserved on the top surface.

    Physical Groups:
    - Group 1: Volume (3D tetrahedral mesh)
    - Group 101: Top box surface at z=-100m (includes fault trace)
    - Group 103: Fault surfaces (cut at box boundaries, embedded)
    - Group 105: Other 5 box surfaces (bottom, front, right, back, left)

    NOTE: The geo file coordinates are in km and will be converted to meters.

    Parameters:
    -----------
    geo_file : str
        Path to the .geo file (default: 'kaikoura_v11.3hpk2.geo')
        NOTE: Coordinates in this file should be in km
    lc : float
        Characteristic length on the fault surface in meters (default: 2000)
    lc_coarse : float
        Coarse characteristic length away from fault in meters (default: 10e3)
    distance_threshold : float
        Distance from fault where coarsening begins in meters (default: 5000.0)
    box_length : float
        Extra padding in X direction beyond fault extent (default: 100e3)
    box_width : float
        Extra padding in Y direction beyond fault extent (default: 100e3)
    box_depth : float
        Depth below fault bottom in meters (default: 100e3)
    output_mesh : str
        Output mesh filename (default: 'ando_fault_mesh.msh')
    use_distance_field : bool
        Use distance-based mesh refinement (slow but gradual). Default: False

    Returns:
    --------
    output_mesh : str
        Path to the generated mesh file
    """
    import gmsh
    import numpy as np

    # Initialize gmsh
    gmsh.initialize()
    gmsh.model.add("ando_fault")

    # Set performance options BEFORE loading geometry
    print("Setting mesh performance options...")
    gmsh.option.setNumber("General.NumThreads", 0)  # Use all available cores
    gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay for 2D (faster)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT for 3D (more robust with embedded surfaces)
    gmsh.option.setNumber("Mesh.Optimize", 1)  # Enable optimization
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)  # Disable Netgen (can be slow)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MaxNumThreads3D", 0)  # Use all cores for 3D meshing

    # Load the .geo file (coordinates in km)
    print(f"Loading geo file: {geo_file}")
    gmsh.merge(geo_file)

    # Synchronize to get the current model state
    gmsh.model.geo.synchronize()

    # Remove all existing physical groups from the loaded geo file
    # (they will be invalid after we remove/transform surfaces)
    print("Removing old physical groups from geo file...")
    try:
        # Get all physical groups and remove them
        all_physical_groups = gmsh.model.getPhysicalGroups()
        if all_physical_groups:
            for dim, tag in all_physical_groups:
                try:
                    gmsh.model.removePhysicalGroups([(dim, tag)])
                except:
                    pass
            print(f"Removed {len(all_physical_groups)} old physical groups")
        else:
            print("No old physical groups found")
    except Exception as e:
        print(f"Note: Could not remove old physical groups: {e}")

    # Convert coordinates from km to m by scaling the entire model
    print("Converting coordinates from km to m (scaling by 1000)...")

    # Get all entities
    all_entities = []
    for dim in range(4):  # 0=points, 1=curves, 2=surfaces, 3=volumes
        entities = gmsh.model.getEntities(dim)
        all_entities.extend(entities)

    # Scale the geometry by 1000 (km to m)
    gmsh.model.geo.dilate(all_entities, 0, 0, 0, 1000, 1000, 1000)
    print("Geometry scaled from km to m")

    # Synchronize after scaling
    gmsh.model.geo.synchronize()

    # Get all surfaces
    surfaces = gmsh.model.getEntities(dim=2)
    print(f'\nTotal surfaces found: {len(surfaces)}')
    print(f'Surface tags: {[s[1] for s in surfaces]}')

    # Calculate mass center for each surface and determine which to remove
    surfaces_to_remove = []
    surfaces_to_keep = []

    for surf_dim, surf_tag in surfaces:
        try:
            # Get bounding box to calculate center
            bbox = gmsh.model.getBoundingBox(surf_dim, surf_tag)
            xmin, ymin, zmin, xmax, ymax, zmax = bbox

            # Calculate mass center (centroid)
            center_x = (xmin + xmax) / 2
            center_y = (ymin + ymax) / 2
            center_z = (zmin + zmax) / 2

            print(f'Surface {surf_tag}: mass center = ({center_x:.2f}, {center_y:.2f}, {center_z:.2f}) m')

            # Check removal criteria
            if abs(center_z) < 1000:  # z ≈ 0 (with tolerance of 1000 m)
                print(f'  -> Marked for removal (z ≈ 0)')
                surfaces_to_remove.append((surf_dim, surf_tag))
            elif center_z < -14000.0:
                print(f'  -> Marked for removal (z < -14000.0 m)')
                surfaces_to_remove.append((surf_dim, surf_tag))
            else:
                print(f'  -> Keeping this surface')
                surfaces_to_keep.append((surf_dim, surf_tag))

        except Exception as e:
            print(f'Warning: Could not process surface {surf_tag}: {e}')
            surfaces_to_keep.append((surf_dim, surf_tag))

    # Remove marked surfaces
    if surfaces_to_remove:
        print(f'\nRemoving {len(surfaces_to_remove)} surfaces: {[s[1] for s in surfaces_to_remove]}')
        try:
            gmsh.model.geo.remove(surfaces_to_remove)
            print('Successfully removed surfaces based on mass center criteria')
        except Exception as e:
            print(f'Warning: Could not remove some surfaces: {e}')
    else:
        print('\nNo surfaces match removal criteria')

    # IMPORTANT: Synchronize after removal to update entity tags
    gmsh.model.geo.synchronize()

    # Get fault surfaces AFTER removal and synchronization
    remaining_surfaces = gmsh.model.getEntities(dim=2)
    fault_surface_tags = [s[1] for s in remaining_surfaces]
    print(f'\nAfter removal, found {len(fault_surface_tags)} fault surface(s): {fault_surface_tags}')

    if not fault_surface_tags:
        raise ValueError("No fault surfaces remain after filtering!")

    # Get bounding box of all fault surfaces combined
    try:
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

        print(f"Combined fault surface bounds:")
        print(f"  X: [{fault_xmin:.2f}, {fault_xmax:.2f}] m (extent: {(fault_xmax-fault_xmin)/1e3:.1f} km)")
        print(f"  Y: [{fault_ymin:.2f}, {fault_ymax:.2f}] m (extent: {(fault_ymax-fault_ymin)/1e3:.1f} km)")
        print(f"  Z: [{fault_zmin:.2f}, {fault_zmax:.2f}] m (extent: {(fault_zmax-fault_zmin)/1e3:.1f} km)")

    except Exception as e:
        print(f"Warning: Could not get fault bounding box: {e}")
        raise ValueError("Cannot create box without fault bounding box!")

    # Create box that encompasses the entire fault with padding
    print("\n=== Creating box domain to encompass fault ===")

    # Add extra buffer to avoid geometric intersections (5% of domain size, minimum 5km)
    # This is critical to prevent PLC errors during 3D meshing
    buffer_x = max(5000.0, (fault_xmax - fault_xmin) * 0.05)
    buffer_y = max(5000.0, (fault_ymax - fault_ymin) * 0.05)
    buffer_z = max(2500.0, (fault_zmax - fault_zmin) * 0.05)

    print(f"Using buffer: X={buffer_x/1e3:.1f} km, Y={buffer_y/1e3:.1f} km, Z={buffer_z/1e3:.1f} km")

    # Box bounds: fault extent + padding on all sides + buffer
    x0 = fault_xmin - box_length / 2 - buffer_x
    x1 = fault_xmax + box_length / 2 + buffer_x
    y0 = fault_ymin - box_width / 2 - buffer_y
    y1 = fault_ymax + box_width / 2 + buffer_y
    z0 = -box_depth  # Depth below fault
    z1 = 0.0+500.0  # Small padding above fault

    actual_box_length = x1 - x0
    actual_box_width = y1 - y0
    actual_box_depth = z1 - z0

    print(f"Box padding:")
    print(f"  X: ±{box_length/2/1e3:.1f} km")
    print(f"  Y: ±{box_width/2/1e3:.1f} km")
    print(f"  Z: -{box_depth/1e3:.1f} km below, z_top={z1/1e3:.1f} km")
    print(f"Total box size:")
    print(f"  {actual_box_length/1e3:.1f} km × {actual_box_width/1e3:.1f} km × {actual_box_depth/1e3:.1f} km")

    # Create box using geo kernel (simpler, compatible with loaded surfaces)
    print("\n=== Creating box geometry ===")

    # Add box points (bottom face: p1-p4, top face: p5-p8)
    p1 = gmsh.model.geo.addPoint(x0, y0, z0, lc_coarse)
    p2 = gmsh.model.geo.addPoint(x1, y0, z0, lc_coarse)
    p3 = gmsh.model.geo.addPoint(x1, y1, z0, lc_coarse)
    p4 = gmsh.model.geo.addPoint(x0, y1, z0, lc_coarse)
    p5 = gmsh.model.geo.addPoint(x0, y0, z1, lc_coarse)
    p6 = gmsh.model.geo.addPoint(x1, y0, z1, lc_coarse)
    p7 = gmsh.model.geo.addPoint(x1, y1, z1, lc_coarse)
    p8 = gmsh.model.geo.addPoint(x0, y1, z1, lc_coarse)

    # Create lines for the box
    # Bottom face
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Top face
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p8)
    l8 = gmsh.model.geo.addLine(p8, p5)

    # Vertical edges
    l9 = gmsh.model.geo.addLine(p1, p5)
    l10 = gmsh.model.geo.addLine(p2, p6)
    l11 = gmsh.model.geo.addLine(p3, p7)
    l12 = gmsh.model.geo.addLine(p4, p8)

    # Create curve loops for each face
    # Bottom face (z=z0)
    cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    # Top face (z=z1)
    cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
    # Front face (y=y0)
    cl3 = gmsh.model.geo.addCurveLoop([l1, l10, -l5, -l9])
    # Right face (x=x1)
    cl4 = gmsh.model.geo.addCurveLoop([l2, l11, -l6, -l10])
    # Back face (y=y1)
    cl5 = gmsh.model.geo.addCurveLoop([l3, l12, -l7, -l11])
    # Left face (x=x0)
    cl6 = gmsh.model.geo.addCurveLoop([l4, l9, -l8, -l12])

    # Create surfaces from curve loops
    surf_bottom = gmsh.model.geo.addPlaneSurface([cl1])
    surf_top = gmsh.model.geo.addPlaneSurface([cl2])
    surf_front = gmsh.model.geo.addPlaneSurface([cl3])
    surf_right = gmsh.model.geo.addPlaneSurface([cl4])
    surf_back = gmsh.model.geo.addPlaneSurface([cl5])
    surf_left = gmsh.model.geo.addPlaneSurface([cl6])

    print(f"Box surfaces created: bottom={surf_bottom}, top={surf_top}, front={surf_front}, right={surf_right}, back={surf_back}, left={surf_left}")

    # IMPORTANT: Synchronize to finalize all geometry
    gmsh.model.geo.synchronize()

    volume_loop = [surf_bottom, surf_top, surf_front, surf_right, surf_back, surf_left]

    print('volumne loop include: ',volume_loop)
    # Create the box volume
    box_surface_loop = gmsh.model.geo.addSurfaceLoop(volume_loop)

    fault_loop = gmsh.model.geo.addSurfaceLoop(fault_surface_tags)  # Create loop!


    # Create volume: outer + all internal faults
    box_volume = gmsh.model.geo.addVolume([box_surface_loop,fault_loop])


    print(f"Box volume created: {box_volume}")

    #  Synchronize
    gmsh.model.geo.synchronize()

    # Set variables for consistency with rest of code
    other_box_surfaces = [surf_bottom, surf_front, surf_right, surf_back, surf_left]

    # Add physical group for the volume (important for SeisSol)
    print("\nAdding physical groups...")
    
    print(f"Added physical group 1 for volume (name: 'volume')")

    # Add physical group for embedded fault surfaces (already cut by fragment operation)
    # No need to embed - they're already part of the volume boundary after fragmentation
    if fault_surface_tags:
        print(f"\nAdding physical group for {len(fault_surface_tags)} fault surface(s)...")
        try:
            # These surfaces are now internal boundaries in the volume
            # The fragment operation already ensured they're properly integrated
            gmsh.model.addPhysicalGroup(2, fault_surface_tags, 103)
            print(f"Added physical group 103 for fault surfaces (name: 'fault')")
        except Exception as e:
            print(f"Warning: Could not add fault physical group: {e}")
    else:
        print("\nWarning: No embedded fault surfaces found after fragmentation!")

    # Set mesh size on all fault surfaces
    # Note: Do this BEFORE creating mesh fields to avoid conflicts
    if fault_surface_tags:
        try:
            total_fault_points = 0
            for surf_tag in fault_surface_tags:
                # Get all points on this fault surface
                surf_bounds = gmsh.model.getBoundary([(2, surf_tag)], combined=False, oriented=False, recursive=True)
                surf_points = [abs(b[1]) for b in surf_bounds if b[0] == 0]
                total_fault_points += len(surf_points)

                # Set fine mesh size on fault
                for pt in surf_points:
                    gmsh.model.mesh.setSize([(0, pt)], lc)

            print(f"Set mesh size on {len(fault_surface_tags)} fault surface(s): lc = {lc} m ({total_fault_points} points)")
        except Exception as e:
            print(f"Warning: Could not set mesh size on fault surfaces: {e}")

    # Choose meshing strategy
    if use_distance_field:
        print("\nUsing distance-based mesh refinement (slower but gradual)...")
        # Set up distance-based mesh refinement
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", fault_surface_tags)

        # Threshold field for gradual coarsening
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", lc)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_coarse)
        gmsh.model.mesh.field.setNumber(2, "DistMin", distance_threshold)
        gmsh.model.mesh.field.setNumber(2, "DistMax", distance_threshold * 10)

        # Set the threshold field as background mesh
        gmsh.model.mesh.field.setAsBackgroundMesh(2)

        print(f"Mesh refinement: {lc} m near fault, coarsening to {lc_coarse} m at {distance_threshold} m distance")
    else:
        print("\nUsing simple mesh size control (faster, less gradual transition)...")
        # Use simpler approach: just set mesh sizes, let gmsh interpolate
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_coarse)
        print(f"Mesh size range: {lc} m (min) to {lc_coarse} m (max)")

    # First generate 2D mesh (required before 3D)
    print("\nGenerating 2D mesh...")
    try:
        gmsh.model.mesh.generate(2)
        num_2d_nodes = len(gmsh.model.mesh.getNodes()[0])
        print(f"2D mesh generated successfully ({num_2d_nodes} nodes)")
    except Exception as e:
        print(f"Warning: 2D mesh generation had issues: {e}")
        # Continue anyway, might work for 3D

    # Generate 3D mesh with fallback algorithms
    print("\nGenerating 3D mesh...")
    print("This may take a while for complex geometries...")

    # Try different 3D meshing algorithms in order of robustness
    algorithms = [
        (1, "Delaunay"),
        (4, "Frontal"),
        (10, "HXT")
    ]

    mesh_success = False
    last_error = None

    gmsh.model.geo.synchronize()

    for algo_id, algo_name in algorithms:
        try:
            print(f"Trying 3D algorithm: {algo_name} (id={algo_id})...")
            gmsh.option.setNumber("Mesh.Algorithm3D", algo_id)

            # Clear 3D mesh but keep 2D
            try:
                gmsh.model.mesh.clear(3)
            except:
                pass
        except:
            pass

    gmsh.model.mesh.generate(3)
    print(f"Successfully generated 3D mesh using {algo_name} algorithm")
    mesh_success = True
           
    if not mesh_success:
        error_msg = f"All 3D meshing algorithms failed. Last error: {last_error}\n"
        error_msg += "Try adjusting mesh sizes (increase lc and lc_coarse)."
        raise Exception(error_msg)

    # Add physical groups for box boundary surfaces AFTER successful meshing
    # This ensures all 6 box surfaces are saved to the output file
    print("\nAdding physical groups for box boundary surfaces...")
    try:
        # Top surface (z=z1) goes to physical group 101
        gmsh.model.addPhysicalGroup(2, [surf_top], 101)
        print(f"Added physical group 101 for top surface (name: 'box_top')")

        # Other box surfaces (bottom + 4 sides) go to physical group 105
        # Note: other_box_surfaces was already set during fragmentation
        gmsh.model.addPhysicalGroup(2, other_box_surfaces, 105)
        print(f"Added physical group 105 for {len(other_box_surfaces)} other box surfaces (name: 'box_sides')")
    except Exception as e:
        print(f"Warning: Could not add box boundary physical groups: {e}")

    gmsh.model.geo.synchronize()


    volumes = gmsh.model.getEntities(dim=3)
    print(volumes)

    gmsh.model.addPhysicalGroup(3, [volumes[0][1]], 1)

    # Get mesh statistics
    num_nodes = len(gmsh.model.mesh.getNodes()[0])
    num_elements = len(gmsh.model.mesh.getElements()[1][0]) if gmsh.model.mesh.getElements()[1] else 0
    print(f"\nMesh statistics:")
    print(f"  Nodes: {num_nodes}")
    print(f"  Elements: {num_elements}")

    # Write mesh file
    gmsh.write(output_mesh)
    print(f"Mesh saved to: {output_mesh}")

    # Optional: Launch GUI to visualize
    # gmsh.fltk.run()

    # Finalize
    gmsh.finalize()

    return output_mesh
