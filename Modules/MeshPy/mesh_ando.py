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

    The box will encompass the entire fault surface with specified padding.

    Physical Groups:
    - Group 1: Volume (3D tetrahedral mesh)
    - Group 101: Top box surface (free surface)
    - Group 103: Fault surfaces (embedded, from geo file)
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
    # gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay for 2D (faster)
    # gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT for 3D (more robust with embedded surfaces)
    # gmsh.option.setNumber("Mesh.Optimize", 1)  # Enable optimization
    # gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)  # Disable Netgen (can be slow)
    # gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
    # gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    # gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
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
            print(f"Found {len(all_physical_groups)} old physical groups to remove")
            for dim, tag in all_physical_groups:
                try:
                    # Correct API: gmsh.model.removePhysicalGroups takes a list of tuples
                    gmsh.model.removePhysicalGroups([(dim, tag)])
                    print(f"  Removed physical group: dim={dim}, tag={tag}")
                except Exception as e:
                    print(f"  Warning: Could not remove group (dim={dim}, tag={tag}): {e}")

            # Verify removal
            remaining_groups = gmsh.model.getPhysicalGroups()
            if remaining_groups:
                print(f"  WARNING: {len(remaining_groups)} physical groups still remain!")
            else:
                print(f"  Successfully removed all old physical groups")
        else:
            print("  No old physical groups found")
    except Exception as e:
        print(f"  Note: Could not remove old physical groups: {e}")

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

    # Remove specific fault surfaces (PK surfaces)
    # surfaces_to_remove_pk = [(2, 54)]
    # try:
    #     gmsh.model.geo.remove(surfaces_to_remove_pk)
    #     print(f'Removed PK surfaces: {[s[1] for s in surfaces_to_remove_pk]}')
    # except Exception as e:
    #     print(f'Warning: Could not remove PK surfaces: {e}')

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
    z0 =  - box_depth   # Depth below fault
    z1 = 0.0  # Small padding above fault

    actual_box_length = x1 - x0
    actual_box_width = y1 - y0
    actual_box_depth = z1 - z0

    print(f"Box padding:")
    print(f"  X: ±{box_length/2/1e3:.1f} km")
    print(f"  Y: ±{box_width/2/1e3:.1f} km")
    print(f"  Z: -{box_depth/1e3:.1f} km below, +{box_depth/4/1e3:.1f} km above")
    print(f"Total box size:")
    print(f"  {actual_box_length/1e3:.1f} km × {actual_box_width/1e3:.1f} km × {actual_box_depth/1e3:.1f} km")

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

    # Create box volume first (without fault surfaces)
    print(f"\nCreating box volume...")
    box_surface_loop = gmsh.model.geo.addSurfaceLoop([surf_bottom, surf_top, surf_front, surf_right, surf_back, surf_left])
    box_volume = gmsh.model.geo.addVolume([box_surface_loop])
    print(f"Box volume created: {box_volume}")

    # Synchronize before embedding
    gmsh.model.geo.synchronize()

    # EMBED the fault surfaces in the volume
    # This is the correct approach for internal surfaces that don't connect to boundaries
    print(f"\nEmbedding {len(fault_surface_tags)} fault surface(s) in box volume...")
    print(f"  Fault surfaces: {fault_surface_tags}")

    gmsh.model.mesh.embed(2, fault_surface_tags, 3, box_volume)
    # print(f"Successfully embedded fault surfaces in volume {box_volume}")

    # Verify volume exists
    volumes_after_embed = gmsh.model.getEntities(dim=3)
    print(f"Volumes after embedding: {volumes_after_embed}")


    # # Add physical group for all fault 
    # gmsh.model.addPhysicalGroup(2, fault_surface_tags, 103)
    # # gmsh.model.setPhysicalName(2, 103, "fault")
    # print(f"Added physical group 103 for {len(fault_surface_tags)} fault surface(s) (name: 'fault')")


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

                # Set fine mesh size on fault - use the regular mesh API after sync
                for pt in surf_points:
                    try:
                        gmsh.model.mesh.setSize([(0, pt)], lc)
                    except:
                        # If mesh API doesn't work, try geo API
                        gmsh.model.geo.mesh.setSize([(0, pt)], lc)

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
        (4, "Frontal"),
        (1, "Delaunay"),
        (10, 'HXT'),
    ]

    mesh_success = False
    last_error = None
    for algo_id, algo_name in algorithms:
        try:
            print(f"Trying 3D algorithm: {algo_name} (id={algo_id})...")
            gmsh.option.setNumber("Mesh.Algorithm3D", algo_id)

            # Clear 3D mesh but keep 2D
            try:
                gmsh.model.mesh.clear(3)
            except:
                pass

            gmsh.model.mesh.generate(3)
            print(f"Successfully generated 3D mesh using {algo_name} algorithm")
            mesh_success = True
            break
        except Exception as e:
            print(f"Algorithm {algo_name} failed: {e}")
            last_error = e
            continue

    if not mesh_success:
        error_msg = f"All 3D meshing algorithms failed. Last error: {last_error}\n"
        error_msg += "Try adjusting mesh sizes (increase lc and lc_coarse)."
        raise Exception(error_msg)

    # Remove duplicate nodes after meshing
    print("\nRemoving duplicate nodes...")
    try:
        # Get current node count before removal
        node_tags_before, _, _ = gmsh.model.mesh.getNodes()
        num_before = len(node_tags_before)

        # Remove duplicate nodes (tolerance in model units - meters)
        gmsh.model.mesh.removeDuplicateNodes()

        # Get node count after removal
        node_tags_after, _, _ = gmsh.model.mesh.getNodes()
        num_after = len(node_tags_after)

        num_removed = num_before - num_after
        if num_removed > 0:
            print(f"  Removed {num_removed} duplicate nodes")
            print(f"  Nodes before: {num_before}")
            print(f"  Nodes after: {num_after}")
        else:
            print(f"  No duplicate nodes found")
    except Exception as e:
        print(f"  Warning: Could not remove duplicate nodes: {e}")

    print("\n✓ 3D mesh generation complete")

    # Add physical groups for box boundary surfaces AFTER successful meshing
    # This ensures all 6 box surfaces are saved to the output file
    # print("\nAdding physical groups for box boundary surfaces...")
    # try:
    #     # Top surface (z=z1) goes to physical group 101
    #     gmsh.model.addPhysicalGroup(2, [surf_top], 101)
    #     # gmsh.model.setPhysicalName(2, 101, "box_top")
    #     print(f"Added physical group 101 for top surface (name: 'box_top')")

    #     # Other box surfaces (bottom + 4 sides) go to physical group 105
    #     other_box_surfaces = [surf_bottom, surf_front, surf_right, surf_back, surf_left]
    #     gmsh.model.addPhysicalGroup(2, other_box_surfaces, 105)
    #     # gmsh.model.setPhysicalName(2, 105, "box_sides")
    #     print(f"Added physical group 105 for {len(other_box_surfaces)} other box surfaces (name: 'box_sides')")
    # except Exception as e:
    #     print(f"Warning: Could not add box boundary physical groups: {e}")

    gmsh.model.geo.synchronize()

    # IMPORTANT: Set SaveAll=1 temporarily to ensure physical groups are recognized
    # We'll set it back to 0 before writing to exclude lower-dimensional elements
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    print("Set Mesh.SaveAll=1 (temporary - for physical group setup)")

    all_surfaces = gmsh.model.getEntities(dim=2)

   

    faults_end = []
    freesurf = []
    absorb = []

    xc = 0.5*(fault_xmax+fault_xmin)
    yc = 0.5*(fault_ymax+fault_ymin)
    zc = 0.5*(fault_zmax+fault_zmin)

    # Known box surface tags
    box_surfaces = [surf_bottom, surf_top, surf_front, surf_right, surf_back, surf_left]

    for surface in all_surfaces:
        try:
            surf_tag = surface[1]
            bbox = gmsh.model.getBoundingBox(surface[0], surf_tag)
            xmin, ymin, zmin, xmax, ymax, zmax = bbox

            x0 = 0.5*(xmin + xmax)
            y0 = 0.5*(ymin + ymax)
            z0 = 0.5*(zmin + zmax)

            # Check if this is a box surface
            if surf_tag in box_surfaces:
                # Top surface (free surface)
                if surf_tag == surf_top:
                    freesurf.append(surf_tag)
                    print(f'  Surface {surf_tag}: Box top -> Group 101 (free surface)')
                # Other box surfaces (absorbing boundaries)
                else:
                    absorb.append(surf_tag)
                    print(f'  Surface {surf_tag}: Box side/bottom -> Group 105 (absorbing)')
            else:
                # This is a fault surface (embedded)
                faults_end.append(surf_tag)
                print(f'  Surface {surf_tag}: Fault surface -> Group 103')

        except Exception as e:
            print(f'Error processing surface {surface[1]}: {e}')

    # DO NOT synchronize here - it can invalidate the surface lists!
    # gmsh.model.geo.synchronize()

    print('\n=== Physical Group Assignment ===')
    print(f'Group 101 (Free Surface): {freesurf}')
    print(f'Group 103 (Fault Surfaces): {faults_end}')
    print(f'Group 105 (Absorbing Boundaries): {absorb}')


    # Get all volumes after meshing
    print('\n=== Checking for Volumes ===')
    all_volumes = gmsh.model.getEntities(dim=3)
    print(f'Found {len(all_volumes)} volume(s): {all_volumes}')

    if not all_volumes:
        print("ERROR: No volumes found after meshing!")
        print("This indicates the volume was not properly created or was lost.")
        print("Checking mesh elements...")

        # Check what was actually meshed
        element_types, _, _ = gmsh.model.mesh.getElements()
        print(f"Element types in mesh: {element_types}")
        if 4 not in element_types:
            print("  No tetrahedral elements (type 4) found!")
            print("  The 3D mesh was not generated.")

        raise ValueError("No volumes found in the mesh! Check volume creation and meshing.")

    print('\n=== Summary ===')
    print(f'Volumes: {all_volumes}')
    print(f'Free surface (Group 101): {freesurf}')
    print(f'Fault surfaces (Group 103): {faults_end}')
    print(f'Absorbing boundaries (Group 105): {absorb}')

    # Add physical groups - DO THIS BEFORE ANY SYNCHRONIZE
    # IMPORTANT: For embedded surfaces to appear in the output, we must set SaveAll=1
    # This is already set earlier at line 446

    # First, add all surface groups (including embedded fault surfaces)
    print('\n=== Adding Surface Physical Groups ===')

    # Free surface (group 101)
    if freesurf:
        gmsh.model.addPhysicalGroup(2, freesurf, 101)
        gmsh.model.setPhysicalName(2, 101, "freesurf")
        print(f'  ✓ Added group 101 (free surface): {freesurf}')
    else:
        print('  WARNING: No free surface found!')

    # Fault surfaces (group 103) - use original fault_surface_tags
    if fault_surface_tags:
        print(f'  Adding fault surfaces to group 103...')
        print(f'    Fault surface tags: {fault_surface_tags}')
        gmsh.model.addPhysicalGroup(2, fault_surface_tags, 103)
        gmsh.model.setPhysicalName(2, 103, "fault")
        print(f'  ✓ Added group 103 (fault): {fault_surface_tags}')
    else:
        print('  WARNING: No fault surfaces found!')

    # Absorbing boundaries (group 105)
    if absorb:
        gmsh.model.addPhysicalGroup(2, absorb, 105)
        gmsh.model.setPhysicalName(2, 105, "absorb")
        print(f'  ✓ Added group 105 (absorbing): {absorb}')
    else:
        print('  WARNING: No absorbing boundaries found!')

    # Finally, add volume group 1
    print('\n=== Adding Volume Physical Group ===')

    # Use the volume tag that was created (box_volume should match all_volumes[0][1])
    volume_tag = all_volumes[0][1]
    print(f'Volume tag from meshing: {volume_tag}')
    print(f'Original box_volume tag: {box_volume}')

    gmsh.model.addPhysicalGroup(3, [volume_tag], 1)
    gmsh.model.setPhysicalName(3, 1, "volume")
    print(f'  ✓ Added group 1 (volume): {volume_tag}')

    # Get mesh statistics and check for issues
    print(f"\n=== Mesh Statistics and Diagnostics ===")

    # Get all nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    num_nodes = len(node_tags)

    # Reshape coordinates to (N, 3) for analysis
    coords = node_coords.reshape(-1, 3)

    print(f"  Nodes: {num_nodes}")

    # Check for duplicate nodes by comparing coordinates
    unique_coords = np.unique(coords, axis=0)
    num_unique = len(unique_coords)
    num_duplicates = num_nodes - num_unique

    if num_duplicates > 0:
        print(f"  WARNING: Found {num_duplicates} duplicate node coordinates!")
        print(f"  Unique positions: {num_unique}")
    else:
        print(f"  ✓ No duplicate nodes detected")

    # Get element statistics
    element_types, element_tags_list, _ = gmsh.model.mesh.getElements()

    total_elements = 0
    for i, elem_type in enumerate(element_types):
        elem_name = gmsh.model.mesh.getElementProperties(elem_type)[0]
        num_elems = len(element_tags_list[i])
        total_elements += num_elems
        print(f"  {elem_name}: {num_elems} elements")

    print(f"  Total elements: {total_elements}")

    # Check mesh quality for 3D elements (tetrahedra = type 4)
    if 4 in element_types:
        idx = element_types.tolist().index(4)
        tet_tags = element_tags_list[idx]

        print(f"\n  Tetrahedral mesh quality:")
        print(f"    Number of tetrahedra: {len(tet_tags)}")

        # Get qualities for all tetrahedra
        try:
            qualities = gmsh.model.mesh.getElementQualities(tet_tags, "minSJ")
            if len(qualities) > 0:
                min_qual = np.min(qualities)
                max_qual = np.max(qualities)
                avg_qual = np.mean(qualities)
                print(f"    Quality (min scaled Jacobian):")
                print(f"      Min: {min_qual:.4f}")
                print(f"      Max: {max_qual:.4f}")
                print(f"      Avg: {avg_qual:.4f}")

                # Count poor quality elements
                poor_quality = np.sum(qualities < 0.1)
                if poor_quality > 0:
                    print(f"    WARNING: {poor_quality} elements with quality < 0.1")
        except:
            print(f"    (Quality calculation not available)")

    # Check mesh size distribution
    print(f"\n  Coordinate ranges:")
    print(f"    X: [{coords[:, 0].min()/1e3:.2f}, {coords[:, 0].max()/1e3:.2f}] km")
    print(f"    Y: [{coords[:, 1].min()/1e3:.2f}, {coords[:, 1].max()/1e3:.2f}] km")
    print(f"    Z: [{coords[:, 2].min()/1e3:.2f}, {coords[:, 2].max()/1e3:.2f}] km")

    # Estimate average element size for 3D mesh
    domain_volume = (coords[:, 0].max() - coords[:, 0].min()) * \
                    (coords[:, 1].max() - coords[:, 1].min()) * \
                    (coords[:, 2].max() - coords[:, 2].min())

    if total_elements > 0:
        avg_elem_volume = domain_volume / total_elements
        avg_elem_size = avg_elem_volume ** (1/3)
        print(f"\n  Domain volume: {domain_volume/1e9:.2f} km³")
        print(f"  Average element size: {avg_elem_size:.1f} m")
        print(f"  Expected element size (lc): {lc:.1f} m")

        if avg_elem_size < lc * 0.5:
            print(f"  WARNING: Average element size much smaller than lc!")
            print(f"           Mesh may be over-refined.")

    # Check for extremely small elements that could indicate issues
    try:
        # Sample some element volumes for tetrahedra
        if 4 in element_types:
            idx = element_types.tolist().index(4)
            sample_tags = element_tags_list[idx][:min(1000, len(element_tags_list[idx]))]

            sample_volumes = []
            for tag in sample_tags:
                vol = gmsh.model.mesh.getElementQualities([tag], "volume")
                if len(vol) > 0:
                    sample_volumes.append(vol[0])

            if sample_volumes:
                min_vol = np.min(sample_volumes)
                max_vol = np.max(sample_volumes)
                avg_vol = np.mean(sample_volumes)

                min_size = min_vol ** (1/3)
                max_size = max_vol ** (1/3)

                print(f"\n  Element size range (sampled):")
                print(f"    Min: {min_size:.1f} m")
                print(f"    Max: {max_size:.1f} m")
                print(f"    Avg: {avg_vol**(1/3):.1f} m")

                if max_size / min_size > 100:
                    print(f"  WARNING: Large size variation (ratio: {max_size/min_size:.1f})")
                    print(f"           This may indicate meshing issues.")
    except Exception as e:
        print(f"\n  (Element size analysis not available: {e})")

    # Verify physical groups before saving
    physical_groups = gmsh.model.getPhysicalGroups()
    print(f"\n=== Physical Groups in Model ===")
    for dim, tag in physical_groups:
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        group_type = {0: 'Point', 1: 'Curve', 2: 'Surface', 3: 'Volume'}[dim]
        print(f"  {group_type} Group {tag} ({name}): {len(entities)} entities -> {entities}")

    # Set mesh format options
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # Use 2.2 for SeisSol compatibility

    # Remove lower-dimensional mesh elements (0D points, 1D edges)
    # Keep only 2D surface elements and 3D volume elements
    print("\nRemoving lower-dimensional elements (points and edges)...")
    try:
        # Remove 1D elements (edges)
        edges_1d = gmsh.model.getEntities(dim=1)
        if edges_1d:
            print(f"  Found {len(edges_1d)} edge entities")
            for dim, tag in edges_1d:
                try:
                    gmsh.model.mesh.clear(dim, tag)
                except:
                    pass
            print(f"  Removed 1D edge elements")

        # Remove 0D elements (points)
        points_0d = gmsh.model.getEntities(dim=0)
        if points_0d:
            print(f"  Found {len(points_0d)} point entities")
            for dim, tag in points_0d:
                try:
                    gmsh.model.mesh.clear(dim, tag)
                except:
                    pass
            print(f"  Removed 0D point elements")

        # Alternative: Use SaveAll=0 and only save elements in physical groups
        # This automatically excludes 0D and 1D elements not in physical groups
        gmsh.option.setNumber("Mesh.SaveAll", 0)
        print("  Set Mesh.SaveAll=0 (only save elements in physical groups)")

    except Exception as e:
        print(f"  Warning: Could not remove lower-dimensional elements: {e}")

    # Write mesh file
    print(f"\n=== Writing Mesh File ===")
    gmsh.write(output_mesh)
    print(f"Mesh saved to: {output_mesh}")
    print(f"Format: MSH 2.2 - 3D volume mesh with surfaces only (no edges/points)")

    # print(fault_surface_tags)

    # Optional: Launch GUI to visualize
    # gmsh.fltk.run()

    # Finalize
    gmsh.finalize()

    return output_mesh
