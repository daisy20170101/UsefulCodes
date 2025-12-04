from MeshFunc.locate_topo_data import locate_modify_topo
import gmsh
import numpy as np


def generate_mesh_fault(MeshFile, stl_fault, stl_box, xc, yc, zmax, zmin, L, B, ncfile):
    """
    Generate a 3D mesh from a watertight box STL and a fault STL.
    xc, yc, zmin: center and the bottom of the domain
    L, B: length and width of the domain
    """

    # Check if gmsh is already initialized
    try:
        gmsh.initialize()
    except:
        # If already initialized, just continue
        pass

    gmsh.clear()
    gmsh.model.add("KaikouraFault")
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay

    # --- Step 1: merge the fault STL (surface only)
    gmsh.merge(stl_fault)

    gmsh.model.geo.synchronize()

    # Collect all current surfaces (faults)
    fault_surfaces = gmsh.model.getEntities(dim=2)

    fault_faces = [tag for (dim, tag) in fault_surfaces]

    # --- Step 2: merge the box STL (watertight volume)
    gmsh.merge(stl_box)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.geo.synchronize()

    # Convert imported STL surfaces into a closed volume

    box_surfaces = gmsh.model.getEntities(dim=2)

    gmsh.model.geo.synchronize()

    box_surface_tags = [s[1] for s in box_surfaces]

    loop = gmsh.model.geo.addSurfaceLoop(box_surface_tags)
    v1 = gmsh.model.geo.addVolume([loop])
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(3, [v1], 1)

    # --- Step 3: distance-based refinement near faults
    print("Stage 2: Applying distance-based refinement...")
    resolution = 500.0  # base element size
    r = 1000.0  # distance scaling

    # distance = gmsh.model.mesh.field.add("Distance")
    # gmsh.model.mesh.field.setNumbers(distance, "FacesList", fault_faces)
    # gmsh.model.mesh.field.setNumber(distance, "Sampling", 100)

    # threshold1 = gmsh.model.mesh.field.add("Threshold")
    # gmsh.model.mesh.field.setNumber(threshold1, "IField", distance)
    # gmsh.model.mesh.field.setNumber(threshold1, "LcMin", 1.0 * resolution)
    # gmsh.model.mesh.field.setNumber(threshold1, "LcMax", 20.0 * resolution)
    # gmsh.model.mesh.field.setNumber(threshold1, "DistMin", 1.0 * r)
    # gmsh.model.mesh.field.setNumber(threshold1, "DistMax", 3 * r)

    # minimum = gmsh.model.mesh.field.add("Min")
    # gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [threshold1])
    # gmsh.model.mesh.field.setAsBackgroundMesh(minimum)

    # --- Step 4: generate mesh

    volumes = gmsh.model.getEntities(dim=3)

    print("volumes:", volumes)

    gmsh.model.mesh.generate(3)
    gmsh.model.geo.synchronize()

    # --- Step 5: classify surfaces for BCs
    faces = gmsh.model.getEntities(dim=2)

    print("All faces:", faces)

    absorb_end = []
    fault_group = []

    for surface in faces:

        nodeTags, nodeCoord, _ = gmsh.model.mesh.getNodes(2, surface[1])
        vxyz = nodeCoord.reshape((-1, 3))
        com = vxyz.mean(axis=0)

        # print(com)

        if np.allclose(com, [xc, yc, zmax], atol=1e3):

            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], 101)

            print("Free surface:", surface[1])

        elif np.sqrt((com[0] - xc) ** 2 + (com[1] - yc) ** 2) > 0.75 * (L + B) / 4.0:
            absorb_end.append(surface[1])

        elif np.allclose(com, [xc, yc, zmin], atol=1e3):
            absorb_end.append(surface[1])

        else:
            # print("fault:", surface[1])
            fault_group.append(surface[1])

    print("Absorbing BC:", absorb_end)
    print("fault: ", fault_group)

    # if absorb_end:
    gmsh.model.addPhysicalGroup(2, absorb_end, 105)

    # if fault_group:
    gmsh.model.addPhysicalGroup(2, fault_group, 103)

    gmsh.model.geo.synchronize()

    # gmsh.option.setNumber("Mesh.MshFileVersion", 4)
    # gmsh.option.setNumber("Mesh.SaveAll", 1)

    # --- Step 6: save mesh
    gmsh.write(MeshFile + ".msh2")
    gmsh.finalize()


def split_stl_surfaces(input_stl_file, output_dir="./surfaces/"):
    """
    Load an STL file and save each surface (group of connected triangles) as an independent STL file.
    This function preserves the original surface grouping from the STL file.

    Parameters:
    -----------
    input_stl_file : str
        Path to the input STL file
    output_dir : str, optional
        Directory to save individual surface STL files (default: "./surfaces/")

    Returns:
    --------
    list
        List of output file paths for each surface
    """
    import os

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Check if gmsh is already initialized
    try:
        gmsh.initialize()
    except:
        # If already initialized, just continue
        pass

    try:
        # Clear any existing models
        gmsh.clear()

        # Create a new model
        gmsh.model.add("stl_surface_splitter")

        # Load the STL file
        print(f"Loading STL file: {input_stl_file}")
        gmsh.merge(input_stl_file)
        gmsh.model.geo.synchronize()

        # Get all surfaces from the imported STL
        surfaces = gmsh.model.getEntities(dim=2)
        print(f"Found {len(surfaces)} surfaces in STL file")

        output_files = []

        # Process each surface
        for i, (dim, surface_tag) in enumerate(surfaces):
            surface_id = i + 1

            # Create a new model for this surface
            surface_model = f"surface_{surface_id}"
            gmsh.model.add(surface_model)
            gmsh.model.setCurrent(surface_model)

            try:
                # Switch back to original model to get surface data
                gmsh.model.setCurrent("stl_surface_splitter")

                # Get all triangular elements belonging to this surface
                element_types, element_tags, element_node_tags = (
                    gmsh.model.mesh.getElements(dim=2, tag=surface_tag)
                )

                if len(element_types) == 0:
                    print(f"Warning: Surface {surface_id} has no elements, skipping...")
                    gmsh.model.setCurrent(surface_model)
                    gmsh.model.remove()
                    continue

                # Get all nodes for this surface
                node_tags, node_coords, _ = gmsh.model.mesh.getNodes(
                    dim=2, tag=surface_tag
                )

                if len(node_tags) == 0:
                    print(f"Warning: Surface {surface_id} has no nodes, skipping...")
                    gmsh.model.setCurrent(surface_model)
                    gmsh.model.remove()
                    continue

                # Switch back to surface model
                gmsh.model.setCurrent(surface_model)

                # Create mapping from node tags to coordinates
                node_coords_reshaped = node_coords.reshape(-1, 3)
                node_map = {
                    int(tag): coord
                    for tag, coord in zip(node_tags, node_coords_reshaped)
                }

                # Add nodes to the new model
                node_id_map = {}  # Map from old node tag to new node tag
                for old_tag, coord in zip(node_tags, node_coords_reshaped):
                    new_tag = gmsh.model.geo.addPoint(coord[0], coord[1], coord[2])
                    node_id_map[int(old_tag)] = new_tag

                # Add triangular elements
                triangle_count = 0
                for elem_type, elem_tags, elem_nodes in zip(
                    element_types, element_tags, element_node_tags
                ):
                    if elem_type == 2:  # Triangle elements
                        triangles = elem_nodes.reshape(-1, 3)

                        for triangle in triangles:
                            triangle_count += 1
                            try:
                                # Map old node IDs to new node IDs
                                new_nodes = [
                                    node_id_map[int(node_id)] for node_id in triangle
                                ]

                                # Create lines
                                l1 = gmsh.model.geo.addLine(new_nodes[0], new_nodes[1])
                                l2 = gmsh.model.geo.addLine(new_nodes[1], new_nodes[2])
                                l3 = gmsh.model.geo.addLine(new_nodes[2], new_nodes[0])

                                # Create surface
                                loop = gmsh.model.geo.addCurveLoop([l1, l2, l3])
                                surf = gmsh.model.geo.addPlaneSurface([loop])

                            except Exception as e:
                                print(
                                    f"Warning: Could not create triangle {triangle_count} in surface {surface_id}: {e}"
                                )
                                continue

                print(f"Surface {surface_id}: Created {triangle_count} triangles")

                # Synchronize and generate mesh
                gmsh.model.geo.synchronize()
                gmsh.model.mesh.generate(2)

                # Save surface
                output_file = os.path.join(output_dir, f"surface_{surface_id:03d}.stl")
                gmsh.write(output_file)
                output_files.append(output_file)
                print(f"Saved surface {surface_id} to: {output_file}")

            except Exception as e:
                print(f"Error processing surface {surface_id}: {e}")

            # Remove surface model to free memory
            try:
                gmsh.model.setCurrent(surface_model)
                gmsh.model.remove()
            except:
                pass

        print(f"Successfully split {len(output_files)} surfaces from {input_stl_file}")
        return output_files

    except Exception as e:
        print(f"Error processing STL file: {e}")
        return []

    finally:
        # Don't finalize here since gmsh might be used elsewhere
        pass


def diagnose_stl_surfaces(input_stl_file):
    """
    Diagnose potential issues with STL surfaces that might cause missing elements.

    Parameters:
    -----------
    input_stl_file : str
        Path to the input STL file

    Returns:
    --------
    dict
        Dictionary containing diagnostic information
    """
    import numpy as np

    # Check if gmsh is already initialized
    try:
        gmsh.initialize()
    except:
        pass

    try:
        gmsh.clear()
        gmsh.model.add("stl_diagnostics")

        print(f"Diagnosing STL file: {input_stl_file}")
        gmsh.merge(input_stl_file)
        gmsh.model.geo.synchronize()

        # Get surfaces
        surfaces = gmsh.model.getEntities(dim=2)
        print(f"Total surfaces found: {len(surfaces)}")

        diagnostics = {
            "total_surfaces": len(surfaces),
            "surface_info": [],
            "potential_issues": [],
        }

        all_nodes = set()
        all_triangles = []

        for i, (dim, surface_tag) in enumerate(surfaces):
            surface_id = i + 1

            try:
                # Get elements and nodes for this surface
                element_types, element_tags, element_node_tags = (
                    gmsh.model.mesh.getElements(dim=2, tag=surface_tag)
                )
                node_tags, node_coords, _ = gmsh.model.mesh.getNodes(
                    dim=2, tag=surface_tag
                )

                if len(element_types) == 0:
                    diagnostics["potential_issues"].append(
                        f"Surface {surface_id}: No elements found"
                    )
                    continue

                triangle_count = 0
                for elem_type, elem_tags, elem_nodes in zip(
                    element_types, element_tags, element_node_tags
                ):
                    if elem_type == 2:  # Triangle elements
                        triangles = elem_nodes.reshape(-1, 3)
                        triangle_count += len(triangles)
                        all_triangles.extend(triangles.tolist())

                # Collect all node IDs for this surface
                surface_nodes = set(int(tag) for tag in node_tags)
                all_nodes.update(surface_nodes)

                # Calculate surface bounds
                if len(node_coords) > 0:
                    coords = node_coords.reshape(-1, 3)
                    bounds = {
                        "x_min": float(np.min(coords[:, 0])),
                        "x_max": float(np.max(coords[:, 0])),
                        "y_min": float(np.min(coords[:, 1])),
                        "y_max": float(np.max(coords[:, 1])),
                        "z_min": float(np.min(coords[:, 2])),
                        "z_max": float(np.max(coords[:, 2])),
                    }
                else:
                    bounds = None

                surface_info = {
                    "surface_id": surface_id,
                    "surface_tag": surface_tag,
                    "triangle_count": triangle_count,
                    "node_count": len(node_tags),
                    "bounds": bounds,
                }

                diagnostics["surface_info"].append(surface_info)
                print(
                    f"Surface {surface_id}: {triangle_count} triangles, {len(node_tags)} nodes"
                )

            except Exception as e:
                diagnostics["potential_issues"].append(
                    f"Surface {surface_id}: Error reading data - {e}"
                )

        # Check for shared nodes between surfaces (potential intersections)
        shared_nodes = {}
        for surface_info in diagnostics["surface_info"]:
            try:
                node_tags, _, _ = gmsh.model.mesh.getNodes(
                    dim=2, tag=surface_info["surface_tag"]
                )
                for node_tag in node_tags:
                    node_id = int(node_tag)
                    if node_id not in shared_nodes:
                        shared_nodes[node_id] = []
                    shared_nodes[node_id].append(surface_info["surface_id"])
            except:
                continue

        # Find nodes shared by multiple surfaces
        intersection_nodes = {
            node_id: surfaces
            for node_id, surfaces in shared_nodes.items()
            if len(surfaces) > 1
        }

        if intersection_nodes:
            diagnostics["potential_issues"].append(
                f"Found {len(intersection_nodes)} nodes shared between surfaces (potential intersections)"
            )

            # Group by surface pairs
            surface_pairs = {}
            for node_id, surface_list in intersection_nodes.items():
                for i in range(len(surface_list)):
                    for j in range(i + 1, len(surface_list)):
                        pair = tuple(sorted([surface_list[i], surface_list[j]]))
                        if pair not in surface_pairs:
                            surface_pairs[pair] = 0
                        surface_pairs[pair] += 1

            print("\nSurface intersection analysis:")
            for (surf1, surf2), count in sorted(surface_pairs.items()):
                print(f"  Surfaces {surf1}-{surf2}: {count} shared nodes")
                if count > 10:  # Significant intersection
                    diagnostics["potential_issues"].append(
                        f"Significant intersection between surfaces {surf1} and {surf2}: {count} shared nodes"
                    )

        # Check for potential overlapping triangles
        triangle_set = set()
        duplicate_triangles = 0
        for triangle in all_triangles:
            triangle_sorted = tuple(sorted(triangle))
            if triangle_sorted in triangle_set:
                duplicate_triangles += 1
            triangle_set.add(triangle_sorted)

        if duplicate_triangles > 0:
            diagnostics["potential_issues"].append(
                f"Found {duplicate_triangles} potentially duplicate triangles"
            )

        diagnostics["total_triangles"] = len(all_triangles)
        diagnostics["total_nodes"] = len(all_nodes)

        print(f"\nSummary:")
        print(f"  Total triangles: {diagnostics['total_triangles']}")
        print(f"  Total nodes: {diagnostics['total_nodes']}")
        print(f"  Potential issues found: {len(diagnostics['potential_issues'])}")

        if diagnostics["potential_issues"]:
            print("\nPotential issues:")
            for issue in diagnostics["potential_issues"]:
                print(f"  - {issue}")

        return diagnostics

    except Exception as e:
        print(f"Error during diagnosis: {e}")
        return {"error": str(e)}

    finally:
        pass


def split_stl_surfaces_robust(
    input_stl_file, output_dir="./surfaces/", fix_intersections=True
):
    """
    Robust version of split_stl_surfaces with intersection handling.

    Parameters:
    -----------
    input_stl_file : str
        Path to the input STL file
    output_dir : str, optional
        Directory to save individual surface STL files
    fix_intersections : bool, optional
        Whether to attempt fixing intersection issues

    Returns:
    --------
    list
        List of output file paths for each surface
    """
    import os

    # First, diagnose the file
    print("Running diagnostics...")
    diagnostics = diagnose_stl_surfaces(input_stl_file)

    if "error" in diagnostics:
        print(f"Diagnostic failed: {diagnostics['error']}")
        return []

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    try:
        gmsh.clear()
        gmsh.model.add("stl_robust_splitter")

        # Load with healing options if requested
        if fix_intersections:
            gmsh.option.setNumber("Mesh.ToleranceInitialDelaunay", 1e-12)
            gmsh.option.setNumber("Mesh.ToleranceEdgeLength", 1e-12)

        gmsh.merge(input_stl_file)
        gmsh.model.geo.synchronize()

        # Get surfaces
        surfaces = gmsh.model.getEntities(dim=2)
        output_files = []

        for i, (dim, surface_tag) in enumerate(surfaces):
            surface_id = i + 1

            try:
                # Create new model for this surface
                surface_model = f"surface_{surface_id}"
                gmsh.model.add(surface_model)
                gmsh.model.setCurrent(surface_model)

                # Switch back to get data
                gmsh.model.setCurrent("stl_robust_splitter")

                # Get surface data with error handling
                try:
                    element_types, element_tags, element_node_tags = (
                        gmsh.model.mesh.getElements(dim=2, tag=surface_tag)
                    )
                    node_tags, node_coords, _ = gmsh.model.mesh.getNodes(
                        dim=2, tag=surface_tag
                    )
                except Exception as e:
                    print(f"Warning: Could not read surface {surface_id} data: {e}")
                    gmsh.model.setCurrent(surface_model)
                    gmsh.model.remove()
                    continue

                if len(element_types) == 0 or len(node_tags) == 0:
                    print(f"Warning: Surface {surface_id} has no data, skipping...")
                    gmsh.model.setCurrent(surface_model)
                    gmsh.model.remove()
                    continue

                # Switch to surface model
                gmsh.model.setCurrent(surface_model)

                # Add nodes with duplicate checking
                node_coords_reshaped = node_coords.reshape(-1, 3)
                node_id_map = {}

                for old_tag, coord in zip(node_tags, node_coords_reshaped):
                    try:
                        new_tag = gmsh.model.geo.addPoint(coord[0], coord[1], coord[2])
                        node_id_map[int(old_tag)] = new_tag
                    except Exception as e:
                        print(
                            f"Warning: Could not add point for surface {surface_id}: {e}"
                        )
                        continue

                # Add triangular elements with validation
                triangle_count = 0
                successful_triangles = 0

                for elem_type, elem_tags, elem_nodes in zip(
                    element_types, element_tags, element_node_tags
                ):
                    if elem_type == 2:  # Triangle elements
                        triangles = elem_nodes.reshape(-1, 3)

                        for triangle in triangles:
                            triangle_count += 1
                            try:
                                # Check if all nodes exist
                                new_nodes = []
                                for node_id in triangle:
                                    if int(node_id) in node_id_map:
                                        new_nodes.append(node_id_map[int(node_id)])
                                    else:
                                        raise KeyError(f"Node {node_id} not found")

                                if len(new_nodes) == 3:
                                    # Create lines with duplicate checking
                                    l1 = gmsh.model.geo.addLine(
                                        new_nodes[0], new_nodes[1]
                                    )
                                    l2 = gmsh.model.geo.addLine(
                                        new_nodes[1], new_nodes[2]
                                    )
                                    l3 = gmsh.model.geo.addLine(
                                        new_nodes[2], new_nodes[0]
                                    )

                                    # Create surface
                                    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3])
                                    surf = gmsh.model.geo.addPlaneSurface([loop])
                                    successful_triangles += 1

                            except Exception as e:
                                # Don't print every failed triangle to avoid spam
                                if triangle_count <= 10 or triangle_count % 100 == 0:
                                    print(
                                        f"Warning: Could not create triangle {triangle_count} in surface {surface_id}: {e}"
                                    )
                                continue

                print(
                    f"Surface {surface_id}: {successful_triangles}/{triangle_count} triangles successfully created"
                )

                if successful_triangles == 0:
                    print(
                        f"Warning: No triangles created for surface {surface_id}, skipping save"
                    )
                    gmsh.model.remove()
                    continue

                # Synchronize and generate mesh
                gmsh.model.geo.synchronize()
                gmsh.model.mesh.generate(2)

                # Save surface
                output_file = os.path.join(output_dir, f"surface_{surface_id:03d}.stl")
                gmsh.write(output_file)
                output_files.append(output_file)
                print(f"Saved surface {surface_id} to: {output_file}")

            except Exception as e:
                print(f"Error processing surface {surface_id}: {e}")

            # Cleanup
            try:
                gmsh.model.setCurrent(surface_model)
                gmsh.model.remove()
            except:
                pass

        print(
            f"Successfully created {len(output_files)} surface files from {diagnostics['total_surfaces']} surfaces"
        )
        return output_files

    except Exception as e:
        print(f"Error in robust surface splitting: {e}")
        return []

    finally:
        pass
