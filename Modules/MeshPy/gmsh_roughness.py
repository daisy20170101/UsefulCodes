"""
Functions to add surface roughness to fault interfaces in gmsh.

This module provides functions to add realistic roughness to fault surfaces
by shifting points in the normal direction based on a random patch distribution.
"""

import numpy as np
import gmsh
from typing import List, Tuple, Optional


def generate_random_patches(
    bounds: Tuple[float, float, float, float, float, float],
    n_patches: int,
    radius_range: Tuple[float, float],
    random_seed: Optional[int] = None
) -> List[Tuple[float, float, float, float]]:
    """
    Generate random circular patches with varying radii in 3D space.

    Parameters:
    -----------
    bounds : tuple
        (x_min, x_max, y_min, y_max, z_min, z_max) bounds for patch centers
    n_patches : int
        Number of patches to generate
    radius_range : tuple
        (min_radius, max_radius) in same units as bounds
    random_seed : int, optional
        Random seed for reproducibility

    Returns:
    --------
    list of tuples
        List of (center_x, center_y, center_z, radius) for each patch
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    min_radius, max_radius = radius_range

    patches = []
    for _ in range(n_patches):
        center_x = np.random.uniform(x_min, x_max)
        center_y = np.random.uniform(y_min, y_max)
        center_z = np.random.uniform(z_min, z_max)
        radius = np.random.uniform(min_radius, max_radius)
        patches.append((center_x, center_y, center_z, radius))

    return patches


def calculate_roughness_at_point(
    point: Tuple[float, float, float],
    patches: List[Tuple[float, float, float, float]],
    amplitude_factor: float = 0.1,
    blend_exponent: float = 2.0
) -> float:
    """
    Calculate roughness height at a point based on nearby patches.

    The roughness amplitude at a point is proportional to the radius of nearby
    patches, with Gaussian-like decay with distance from patch centers.

    Parameters:
    -----------
    point : tuple
        (x, y, z) coordinates of the point
    patches : list of tuples
        List of (center_x, center_y, center_z, radius) for each patch
    amplitude_factor : float, optional
        Factor controlling roughness amplitude (amplitude = factor * radius)
        Default: 0.1 (roughness is 10% of patch radius)
    blend_exponent : float, optional
        Exponent for distance-based blending (higher = sharper transitions)
        Default: 2.0

    Returns:
    --------
    float
        Roughness height (positive = outward, negative = inward)
    """
    x, y, z = point
    total_roughness = 0.0

    for center_x, center_y, center_z, radius in patches:
        # Calculate 3D Euclidean distance from point to patch center
        distance = np.sqrt((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)

        # Gaussian-like influence function
        # Influence decays with distance, limited to ~2*radius
        influence = np.exp(-0.5 * (distance / radius)**blend_exponent)

        # Amplitude proportional to patch radius
        amplitude = amplitude_factor * radius

        # Random phase for this patch (alternating positive/negative)
        # Use hash of patch center for deterministic "random" phase
        phase = 1.0 if (hash((center_x, center_y, center_z)) % 2) == 0 else -1.0

        # Add contribution from this patch
        total_roughness += phase * amplitude * influence

    return total_roughness


def compute_mesh_normals(
    surface_tag: int,
    gmsh_model=None
) -> dict:
    """
    Compute normal vectors at all nodes on a meshed surface.

    Uses element geometry to compute normals at each node by averaging
    normals of all elements sharing that node.

    Parameters:
    -----------
    surface_tag : int
        Tag of the surface
    gmsh_model : gmsh model, optional
        gmsh model to use (uses default if None)

    Returns:
    --------
    dict
        Dictionary mapping node_tag -> normal vector [nx, ny, nz]
    """
    if gmsh_model is None:
        gmsh_model = gmsh.model

    # Get all elements on the surface
    element_types, element_tags, node_tags_per_elem = gmsh_model.mesh.getElements(dim=2, tag=surface_tag)

    # Get all nodes
    all_node_tags, node_coords, _ = gmsh_model.mesh.getNodes(dim=2, tag=surface_tag)
    node_coords = node_coords.reshape(-1, 3)

    # Create dictionary for node coordinates
    node_coord_dict = {int(tag): coords for tag, coords in zip(all_node_tags, node_coords)}

    # Dictionary to accumulate normals for each node
    node_normals = {int(tag): np.zeros(3) for tag in all_node_tags}

    # Process each element type
    for elem_type, elem_tags, node_tags in zip(element_types, element_tags, node_tags_per_elem):
        # Assuming triangular elements (type 2)
        if elem_type == 2:  # 3-node triangle
            nodes_per_elem = 3
            node_tags = node_tags.reshape(-1, nodes_per_elem)

            # Compute normal for each triangle
            for elem_nodes in node_tags:
                # Get coordinates of the three vertices
                p0 = node_coord_dict[int(elem_nodes[0])]
                p1 = node_coord_dict[int(elem_nodes[1])]
                p2 = node_coord_dict[int(elem_nodes[2])]

                # Compute edge vectors
                v1 = p1 - p0
                v2 = p2 - p0

                # Compute normal using cross product
                normal = np.cross(v1, v2)
                norm_length = np.linalg.norm(normal)

                if norm_length > 1e-10:
                    normal = normal / norm_length

                    # Add this normal to all three nodes
                    for node in elem_nodes:
                        node_normals[int(node)] += normal

    # Normalize accumulated normals
    for node_tag in node_normals:
        normal = node_normals[node_tag]
        norm_length = np.linalg.norm(normal)
        if norm_length > 1e-10:
            node_normals[node_tag] = normal / norm_length
        else:
            # Fallback to z-direction if no valid normal
            node_normals[node_tag] = np.array([0.0, 0.0, 1.0])

    return node_normals


def add_roughness_to_surface(
    surface_tag: int,
    n_patches: int = 50,
    radius_range: Tuple[float, float] = (1000.0, 5000.0),
    amplitude_factor: float = 0.1,
    blend_exponent: float = 2.0,
    normal_direction: Optional[np.ndarray] = None,
    random_seed: Optional[int] = None,
    gmsh_model=None
) -> None:
    """
    Add roughness to a gmsh surface by shifting points along normals.

    This function:
    1. Gets all points on the surface
    2. Generates random patches
    3. Calculates roughness amplitude at each point
    4. Shifts points along surface normal

    Parameters:
    -----------
    surface_tag : int
        gmsh tag of the surface to roughen
    n_patches : int, optional
        Number of roughness patches. Default: 50
    radius_range : tuple, optional
        (min_radius, max_radius) for patches in model units.
        Default: (1000, 5000) meters
    amplitude_factor : float, optional
        Roughness amplitude as fraction of patch radius.
        Default: 0.1 (10% of radius)
    blend_exponent : float, optional
        Blending exponent for smooth transitions. Default: 2.0
    normal_direction : np.ndarray, optional
        Global normal direction [nx, ny, nz]. If None, computed per-point.
        Default: None
    random_seed : int, optional
        Random seed for reproducibility
    gmsh_model : gmsh model, optional
        gmsh model to use (uses default if None)

    Returns:
    --------
    None
        Modifies the gmsh model in place

    Example:
    --------
    >>> import gmsh
    >>> gmsh.initialize()
    >>> gmsh.model.add("fault")
    >>> # ... create surface with tag 1 ...
    >>> add_roughness_to_surface(
    ...     surface_tag=1,
    ...     n_patches=30,
    ...     radius_range=(2000.0, 8000.0),
    ...     amplitude_factor=0.15,
    ...     random_seed=42
    ... )
    >>> gmsh.write("rough_fault.msh")
    """
    if gmsh_model is None:
        gmsh_model = gmsh.model

    # Try to get points from meshed surface first
    is_meshed = False
    try:
        # For meshed surfaces, get all nodes on the surface
        node_tags, node_coords, _ = gmsh_model.mesh.getNodes(dim=2, tag=surface_tag)

        if len(node_tags) > 0:
            # Reshape coordinates from flat array to (n_nodes, 3)
            node_coords = node_coords.reshape(-1, 3)
            points_coords = [(tag, coords) for tag, coords in zip(node_tags, node_coords)]
            is_meshed = True
            print(f"Processing {len(points_coords)} mesh nodes on surface {surface_tag}")
        else:
            raise ValueError("No mesh nodes found")
    except:
        # Fallback: try to get CAD points from boundary
        dim_tags = gmsh_model.getBoundary([(2, surface_tag)], combined=False, oriented=False, recursive=True)
        point_tags = [tag for dim, tag in dim_tags if dim == 0]

        if len(point_tags) == 0:
            raise ValueError(f"Surface {surface_tag} has no points. Please mesh the surface first or check surface tag.")

        points_coords = []
        for point_tag in point_tags:
            coords = gmsh_model.getValue(0, point_tag, [])
            points_coords.append((point_tag, coords))

        print(f"Processing {len(points_coords)} CAD points on surface {surface_tag}")

    if len(points_coords) == 0:
        raise ValueError(f"Could not extract any points from surface {surface_tag}")

    # Calculate bounds for patch generation
    all_coords = np.array([coords for _, coords in points_coords])
    x_min, y_min, z_min = all_coords.min(axis=0)
    x_max, y_max, z_max = all_coords.max(axis=0)

    bounds = (x_min, x_max, y_min, y_max, z_min, z_max)

    # Generate random patches
    patches = generate_random_patches(bounds, n_patches, radius_range, random_seed)
    print(f"Generated {len(patches)} roughness patches in 3D space")

    # Compute normals
    if normal_direction is None:
        if is_meshed:
            # Compute per-node normals from mesh geometry
            print("Computing surface normals from mesh geometry...")
            node_normals = compute_mesh_normals(surface_tag, gmsh_model)
            use_per_node_normals = True
        else:
            # For CAD geometry, use global normal
            normal_direction = np.array([0.0, 0.0, 1.0])
            print("Using default normal direction: [0, 0, 1]")
            use_per_node_normals = False
    else:
        # Normalize provided direction
        normal_direction = np.array(normal_direction)
        normal_direction = normal_direction / np.linalg.norm(normal_direction)
        print(f"Using provided normal direction: {normal_direction}")
        use_per_node_normals = False

    # Apply roughness to each point
    max_roughness = 0.0
    for point_tag, coords in points_coords:
        # Calculate roughness height at this point
        roughness = calculate_roughness_at_point(
            coords, patches, amplitude_factor, blend_exponent
        )
        max_roughness = max(max_roughness, abs(roughness))

        # Get normal for this point
        if use_per_node_normals:
            point_normal = node_normals[int(point_tag)]
        else:
            point_normal = normal_direction

        # Shift point along normal
        new_coords = coords + roughness * point_normal

        # Update point position in gmsh (different method for meshed vs CAD)
        if is_meshed:
            # For mesh nodes, use setNode with empty parametric coords
            gmsh_model.mesh.setNode(int(point_tag), [new_coords[0], new_coords[1], new_coords[2]], [])
        else:
            # For CAD points, use setCoordinates
            gmsh_model.setCoordinates(int(point_tag), new_coords[0], new_coords[1], new_coords[2])

    print(f"Applied roughness: max amplitude = {max_roughness:.2f} model units")
    print("Surface roughness applied successfully!")


def add_anisotropic_roughness(
    surface_tag: int,
    n_patches: int = 50,
    radius_range_along: Tuple[float, float] = (5000.0, 15000.0),
    radius_range_across: Tuple[float, float] = (1000.0, 3000.0),
    strike_angle: float = 0.0,
    amplitude_factor: float = 0.1,
    normal_direction: Optional[np.ndarray] = None,
    random_seed: Optional[int] = None,
    gmsh_model=None
) -> None:
    """
    Add anisotropic roughness with elliptical patches aligned with fault strike.

    Useful for realistic fault roughness where along-strike correlation length
    is typically larger than down-dip correlation length.

    Parameters:
    -----------
    surface_tag : int
        gmsh tag of the surface to roughen
    n_patches : int, optional
        Number of roughness patches. Default: 50
    radius_range_along : tuple, optional
        (min, max) radius along strike direction (meters). Default: (5000, 15000)
    radius_range_across : tuple, optional
        (min, max) radius across strike direction (meters). Default: (1000, 3000)
    strike_angle : float, optional
        Strike angle in degrees (0 = North). Default: 0.0
    amplitude_factor : float, optional
        Roughness amplitude factor. Default: 0.1
    normal_direction : np.ndarray, optional
        Surface normal direction. Default: None (auto-detect)
    random_seed : int, optional
        Random seed for reproducibility
    gmsh_model : gmsh model, optional
        gmsh model to use

    Example:
    --------
    >>> # Add roughness with 10 km along-strike, 2 km down-dip correlation
    >>> add_anisotropic_roughness(
    ...     surface_tag=1,
    ...     radius_range_along=(8000, 12000),
    ...     radius_range_across=(1500, 2500),
    ...     strike_angle=45.0,  # NE-SW strike
    ...     amplitude_factor=0.2
    ... )
    """
    if gmsh_model is None:
        gmsh_model = gmsh.model

    # Try to get points from meshed surface first
    is_meshed = False
    try:
        # For meshed surfaces, get all nodes on the surface
        node_tags, node_coords, _ = gmsh_model.mesh.getNodes(dim=2, tag=surface_tag)

        if len(node_tags) > 0:
            # Reshape coordinates from flat array to (n_nodes, 3)
            node_coords = node_coords.reshape(-1, 3)
            points_coords = [(tag, coords) for tag, coords in zip(node_tags, node_coords)]
            is_meshed = True
            print(f"Found {len(points_coords)} mesh nodes on surface {surface_tag}")
        else:
            raise ValueError("No mesh nodes found")
    except:
        # Fallback: try to get CAD points from boundary
        dim_tags = gmsh_model.getBoundary([(2, surface_tag)], combined=False, oriented=False, recursive=True)
        point_tags = [tag for dim, tag in dim_tags if dim == 0]

        if len(point_tags) == 0:
            raise ValueError(f"Surface {surface_tag} has no points. Please mesh the surface first or check surface tag.")

        points_coords = []
        for point_tag in point_tags:
            coords = gmsh_model.getValue(0, point_tag, [])
            points_coords.append((point_tag, coords))

        print(f"Found {len(points_coords)} CAD points on surface {surface_tag}")

    if len(points_coords) == 0:
        raise ValueError(f"Could not extract any points from surface {surface_tag}")

    all_coords = np.array([coords for _, coords in points_coords])
    x_min, y_min, z_min = all_coords.min(axis=0)
    x_max, y_max, z_max = all_coords.max(axis=0)

    # Generate patches with anisotropic radii in 3D
    if random_seed is not None:
        np.random.seed(random_seed)

    strike_rad = np.deg2rad(strike_angle)
    cos_strike = np.cos(strike_rad)
    sin_strike = np.sin(strike_rad)

    # Rotation matrix for strike alignment in 2D (horizontal plane)
    rotation = np.array([
        [cos_strike, -sin_strike],
        [sin_strike, cos_strike]
    ])

    patches = []
    for _ in range(n_patches):
        center_x = np.random.uniform(x_min, x_max)
        center_y = np.random.uniform(y_min, y_max)
        center_z = np.random.uniform(z_min, z_max)
        radius_along = np.random.uniform(*radius_range_along)
        radius_across = np.random.uniform(*radius_range_across)

        # Store 3D center, both radii and rotation info
        patches.append((center_x, center_y, center_z, radius_along, radius_across, rotation))

    print(f"Generated {len(patches)} anisotropic patches in 3D space")

    # Compute normals
    if normal_direction is None:
        if is_meshed:
            # Compute per-node normals from mesh geometry
            print("Computing surface normals from mesh geometry...")
            node_normals = compute_mesh_normals(surface_tag, gmsh_model)
            use_per_node_normals = True
        else:
            # For CAD geometry, use global normal
            normal_direction = np.array([0.0, 0.0, 1.0])
            print("Using default normal direction: [0, 0, 1]")
            use_per_node_normals = False
    else:
        # Normalize provided direction
        normal_direction = np.array(normal_direction)
        normal_direction = normal_direction / np.linalg.norm(normal_direction)
        print(f"Using provided normal direction: {normal_direction}")
        use_per_node_normals = False

    # Apply anisotropic roughness
    max_roughness = 0.0
    for point_tag, coords in points_coords:
        x, y, z = coords
        total_roughness = 0.0

        for center_x, center_y, center_z, r_along, r_across, rot in patches:
            # Calculate 3D distance components
            dx = x - center_x
            dy = y - center_y
            dz = z - center_z

            # Transform horizontal components to patch-local coordinates (strike-aligned)
            local_xy = rot.T @ np.array([dx, dy])

            # Elliptical distance in 3D: anisotropic in horizontal plane, isotropic in vertical
            # Use geometric mean for vertical direction
            r_vertical = np.sqrt(r_along * r_across)
            distance = np.sqrt((local_xy[0] / r_along)**2 + (local_xy[1] / r_across)**2 + (dz / r_vertical)**2)

            # Influence function
            influence = np.exp(-0.5 * distance**2)

            # Use geometric mean of radii for amplitude
            amplitude = amplitude_factor * np.sqrt(r_along * r_across)
            phase = 1.0 if (hash((center_x, center_y, center_z)) % 2) == 0 else -1.0

            total_roughness += phase * amplitude * influence

        max_roughness = max(max_roughness, abs(total_roughness))

        # Get normal for this point
        if use_per_node_normals:
            point_normal = node_normals[int(point_tag)]
        else:
            point_normal = normal_direction

        # Apply shift
        new_coords = coords + total_roughness * point_normal

        # Update point position in gmsh (different method for meshed vs CAD)
        if is_meshed:
            # For mesh nodes, use setNode with empty parametric coords
            gmsh_model.mesh.setNode(int(point_tag), [new_coords[0], new_coords[1], new_coords[2]], [])
        else:
            # For CAD points, use setCoordinates
            gmsh_model.setCoordinates(int(point_tag), new_coords[0], new_coords[1], new_coords[2])

    print(f"Applied anisotropic roughness: max amplitude = {max_roughness:.2f}")
    print("Anisotropic roughness applied successfully!")


if __name__ == "__main__":
    # Example usage
    print("Example: Adding roughness to a gmsh surface")
    print("=" * 60)
    print("""
    import gmsh
    from DyFltFunc.gmsh_roughness import add_roughness_to_surface

    # Initialize gmsh and create/load model
    gmsh.initialize()
    gmsh.model.add("fault_with_roughness")

    # ... create or load your surface geometry ...
    # Let's say your surface has tag = 1

    # Add isotropic roughness
    add_roughness_to_surface(
        surface_tag=1,
        n_patches=50,
        radius_range=(2000.0, 8000.0),  # 2-8 km patches
        amplitude_factor=0.15,           # 15% of patch radius
        random_seed=42
    )

    # Or add anisotropic roughness (realistic for faults)
    from DyFltFunc.gmsh_roughness import add_anisotropic_roughness

    add_anisotropic_roughness(
        surface_tag=1,
        n_patches=30,
        radius_range_along=(8000, 15000),   # 8-15 km along strike
        radius_range_across=(1500, 3000),   # 1.5-3 km down-dip
        strike_angle=45.0,                  # NE-SW strike
        amplitude_factor=0.2
    )

    # Generate mesh and save
    gmsh.model.mesh.generate(2)
    gmsh.write("rough_fault.msh")
    gmsh.finalize()
    """)
