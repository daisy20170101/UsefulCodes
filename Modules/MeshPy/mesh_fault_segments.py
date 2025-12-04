"""
Extract individual fault segments from .geo file and mesh them separately
Each segment is saved as a separate .msh file for later merging
"""

import gmsh
import numpy as np
import os


def mesh_individual_fault_segments(geo_file='kaikoura_v11.3hpk2.geo',
                                   output_dir='fault_segments',
                                   lc=2000,
                                   remove_criteria=None):
    """
    Extract and mesh each fault surface from a .geo file separately

    Parameters:
    -----------
    geo_file : str
        Path to the .geo file
    output_dir : str
        Directory to save individual fault segment meshes
    lc : float
        Characteristic length for fault surface mesh (in meters)
    remove_criteria : dict, optional
        Criteria for removing surfaces:
        - 'z_near_zero_tol': Remove surfaces with |z| < tolerance (default: 1000 m)
        - 'z_below': Remove surfaces with z < value (default: -14000 m)

    Returns:
    --------
    segment_files : list
        List of paths to generated .msh files for each fault segment
    segment_info : list
        List of dictionaries containing info about each segment:
        - 'file': path to .msh file
        - 'tag': original surface tag
        - 'bbox': bounding box (xmin, ymin, zmin, xmax, ymax, zmax)
        - 'center': mass center (x, y, z)
    """

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Set default removal criteria
    if remove_criteria is None:
        remove_criteria = {
            'z_near_zero_tol': 1000.0,  # Remove surfaces with |z| < 1000 m
            'z_below': -14000.0          # Remove surfaces with z < -14000 m
        }

    print("=" * 70)
    print("MESHING INDIVIDUAL FAULT SEGMENTS")
    print("=" * 70)
    print(f"Input .geo file: {geo_file}")
    print(f"Output directory: {output_dir}")
    print(f"Mesh size (lc): {lc} m")
    print(f"Removal criteria: {remove_criteria}")
    print()

    # Initialize gmsh
    gmsh.initialize()
    gmsh.model.add("fault_segments")

    # Load the .geo file
    print(f"Loading geo file: {geo_file}")
    gmsh.merge(geo_file)
    gmsh.model.geo.synchronize()

    # Remove old physical groups
    try:
        all_physical_groups = gmsh.model.getPhysicalGroups()
        if all_physical_groups:
            for dim, tag in all_physical_groups:
                gmsh.model.removePhysicalGroups([(dim, tag)])
            print(f"Removed {len(all_physical_groups)} old physical groups")
    except Exception as e:
        print(f"Note: Could not remove old physical groups: {e}")

    # Convert coordinates from km to m (scaling by 1000)
    print("\nConverting coordinates from km to m (scaling by 1000)...")
    all_entities = []
    for dim in range(4):
        entities = gmsh.model.getEntities(dim)
        all_entities.extend(entities)

    gmsh.model.geo.dilate(all_entities, 0, 0, 0, 1000, 1000, 1000)
    gmsh.model.geo.synchronize()
    print("Geometry scaled from km to m")

    # Get all surfaces
    surfaces = gmsh.model.getEntities(dim=2)
    print(f"\nTotal surfaces found: {len(surfaces)}")

    # Filter surfaces based on removal criteria
    surfaces_to_keep = []
    surfaces_to_remove = []

    for surf_dim, surf_tag in surfaces:
        try:
            bbox = gmsh.model.getBoundingBox(surf_dim, surf_tag)
            xmin, ymin, zmin, xmax, ymax, zmax = bbox

            # Calculate mass center
            center_x = (xmin + xmax) / 2
            center_y = (ymin + ymax) / 2
            center_z = (zmin + zmax) / 2

            # Check removal criteria
            remove = False
            if abs(center_z) < remove_criteria['z_near_zero_tol']:
                remove = True
                reason = f"z ≈ 0 (|z| < {remove_criteria['z_near_zero_tol']} m)"
            elif center_z < remove_criteria['z_below']:
                remove = True
                reason = f"z < {remove_criteria['z_below']} m"

            if remove:
                surfaces_to_remove.append((surf_dim, surf_tag))
                print(f"  Surface {surf_tag}: REMOVE - {reason}")
            else:
                surfaces_to_keep.append({
                    'tag': surf_tag,
                    'bbox': bbox,
                    'center': (center_x, center_y, center_z)
                })
                print(f"  Surface {surf_tag}: KEEP - center=({center_x:.1f}, {center_y:.1f}, {center_z:.1f}) m")

        except Exception as e:
            print(f"  Warning: Could not process surface {surf_tag}: {e}")

    print(f"\nSurfaces to keep: {len(surfaces_to_keep)}")
    print(f"Surfaces to remove: {len(surfaces_to_remove)}")

    # Finalize this gmsh session
    gmsh.finalize()

    # Now mesh each surface individually
    segment_info = []
    segment_files = []

    for idx, surf_info in enumerate(surfaces_to_keep):
        surf_tag = surf_info['tag']

        print(f"\n{'=' * 70}")
        print(f"Processing segment {idx + 1}/{len(surfaces_to_keep)}: Surface {surf_tag}")
        print(f"{'=' * 70}")

        # Start a new gmsh session for this surface
        gmsh.initialize()
        gmsh.model.add(f"segment_{surf_tag}")

        # Load the original geo file
        gmsh.merge(geo_file)
        gmsh.model.geo.synchronize()

        # Remove all physical groups from this instance
        try:
            all_physical_groups_temp = gmsh.model.getPhysicalGroups()
            if all_physical_groups_temp:
                for dim_pg, tag_pg in all_physical_groups_temp:
                    try:
                        gmsh.model.removePhysicalGroups([(dim_pg, tag_pg)])
                    except:
                        pass
        except:
            pass

        # Scale coordinates
        all_entities_temp = []
        for dim in range(4):
            entities_temp = gmsh.model.getEntities(dim)
            all_entities_temp.extend(entities_temp)
        gmsh.model.geo.dilate(all_entities_temp, 0, 0, 0, 1000, 1000, 1000)
        gmsh.model.geo.synchronize()

        # Remove all surfaces except the current one
        all_surfaces_temp = gmsh.model.getEntities(dim=2)
        surfaces_to_remove_temp = [(2, s[1]) for s in all_surfaces_temp if s[1] != surf_tag]

        if surfaces_to_remove_temp:
            try:
                gmsh.model.geo.remove(surfaces_to_remove_temp)
                gmsh.model.geo.synchronize()
                print(f"  Removed {len(surfaces_to_remove_temp)} other surfaces")
            except Exception as e:
                print(f"  Warning: Could not remove other surfaces: {e}")

        # Verify only our surface remains
        remaining = gmsh.model.getEntities(dim=2)
        print(f"  Remaining surfaces: {[s[1] for s in remaining]}")

        if len(remaining) != 1 or remaining[0][1] != surf_tag:
            print(f"  ERROR: Surface {surf_tag} not found or multiple surfaces remain!")
            gmsh.finalize()
            continue

        # Set mesh size on this surface
        try:
            surf_bounds = gmsh.model.getBoundary([(2, surf_tag)], combined=False, oriented=False, recursive=True)
            surf_points = [abs(b[1]) for b in surf_bounds if b[0] == 0]

            for pt in surf_points:
                try:
                    gmsh.model.mesh.setSize([(0, pt)], lc)
                except:
                    pass

            print(f"  Set mesh size lc={lc} m on {len(surf_points)} points")
        except Exception as e:
            print(f"  Warning: Could not set mesh size: {e}")

        # Add physical group for this surface
        gmsh.model.addPhysicalGroup(2, [surf_tag], 1)
        gmsh.model.setPhysicalName(2, 1, f"fault_segment_{surf_tag}")
        print(f"  Added physical group 1 for surface {surf_tag}")

        # Generate 2D mesh
        print(f"  Generating 2D mesh...")
        try:
            gmsh.model.mesh.generate(2)

            # Get mesh statistics
            node_tags, _, _ = gmsh.model.mesh.getNodes()
            element_types, element_tags_list, _ = gmsh.model.mesh.getElements(dim=2)

            num_nodes = len(node_tags)
            num_elements = sum(len(tags) for tags in element_tags_list)

            print(f"  Mesh generated: {num_nodes} nodes, {num_elements} elements")
        except Exception as e:
            print(f"  ERROR: Mesh generation failed: {e}")
            gmsh.finalize()
            continue

        # Save mesh file
        output_file = os.path.join(output_dir, f"fault_segment_{surf_tag}.msh")
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.SaveAll", 0)  # Only save elements in physical groups

        try:
            gmsh.write(output_file)
            print(f"  Saved mesh to: {output_file}")

            segment_files.append(output_file)
            segment_info.append({
                'file': output_file,
                'tag': surf_tag,
                'bbox': surf_info['bbox'],
                'center': surf_info['center'],
                'num_nodes': num_nodes,
                'num_elements': num_elements
            })
        except Exception as e:
            print(f"  ERROR: Could not save mesh: {e}")

        # Finalize this gmsh session
        gmsh.finalize()

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print(f"Total fault segments meshed: {len(segment_files)}")
    print(f"Output directory: {output_dir}")
    print()

    for info in segment_info:
        print(f"  {os.path.basename(info['file'])}")
        print(f"    Tag: {info['tag']}")
        print(f"    Nodes: {info['num_nodes']}, Elements: {info['num_elements']}")
        print(f"    Center: ({info['center'][0]:.1f}, {info['center'][1]:.1f}, {info['center'][2]:.1f}) m")

    return segment_files, segment_info


if __name__ == "__main__":
    # Example usage
    segment_files, segment_info = mesh_individual_fault_segments(
        geo_file='kaikoura_v11.3hpk2.geo',
        output_dir='fault_segments',
        lc=2000,
        remove_criteria={
            'z_near_zero_tol': 1000.0,
            'z_below': -14000.0
        }
    )

    print(f"\n✓ Successfully meshed {len(segment_files)} fault segments")
