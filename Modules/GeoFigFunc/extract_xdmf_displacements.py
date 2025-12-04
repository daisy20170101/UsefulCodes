"""
Extract displacement components (u1, u2, u3) from XDMF surface files at selected positions
Uses seissolxdmf to read SeisSol XDMF files and scipy for spatial interpolation
"""

import numpy as np
import seissolxdmf
from scipy.interpolate import griddata, LinearNDInterpolator
from scipy.spatial import cKDTree


def load_xdmf_surface(xdmf_file, timestep=-1):
    """
    Load XDMF surface file using seissolxdmf and extract mesh data

    Parameters:
    -----------
    xdmf_file : str
        Path to XDMF file (e.g., 'surface.xdmf')
    timestep : int
        Which timestep to load (-1 for last, 0 for first)

    Returns:
    --------
    points : ndarray
        Mesh node coordinates (N, 3) [x, y, z]
    point_data : dict
        Dictionary of point data arrays (e.g., {'u1': array, 'u2': array, 'u3': array})
    sx : seissolxdmf object
        SeisSol XDMF reader object for further access
    connect : ndarray
        Cell connectivity array
    """
    print(f"Loading XDMF file: {xdmf_file}")

    # Read the mesh using seissolxdmf
    sx = seissolxdmf.seissolxdmf(xdmf_file)

    # Get geometry (vertex coordinates)
    geometry = sx.ReadGeometry()
    points = np.array(geometry)

    # Get connectivity
    connect = sx.ReadConnect()

    print(f"  Number of vertices: {len(points)}")
    print(f"  Number of cells: {len(connect)}")
    print(f"  Number of timesteps: {sx.ndt}")

    # Get available data fields
    data_names = sx.ReadAvailableDataFields()
    print(f"  Available data fields: {data_names}")

    # Handle timestep index
    if timestep == -1:
        timestep = sx.ndt - 1
    elif timestep >= sx.ndt:
        raise ValueError(f"Timestep {timestep} exceeds available timesteps (0-{sx.ndt-1})")

    print(f"  Loading timestep: {timestep}")

    # Read all available data at the specified timestep
    point_data = {}
    for field_name in data_names:
        data = sx.ReadData(field_name, timestep)
        point_data[field_name] = np.array(data)
        if len(data.shape) == 1:
            print(f"    {field_name}: shape ({len(data)},)")
        else:
            print(f"    {field_name}: shape {data.shape}")

    return points, point_data, sx, connect


def extract_displacements_at_positions(xdmf_file,
                                       query_positions,
                                       displacement_fields=None,
                                       timestep=-1,
                                       method='linear',
                                       search_radius=None):
    """
    Extract displacement components (u1, u2, u3) at specified (x, y) positions

    Parameters:
    -----------
    xdmf_file : str
        Path to XDMF surface file
    query_positions : array-like
        Query positions as (N, 2) array of [x, y] coordinates
        or (N, 3) array of [x, y, z] coordinates
    displacement_fields : list or str, optional
        Field names for displacement components. Can be:
        - List of 3 field names: ['u1', 'u2', 'u3']
        - Single field name containing vector: 'u'
        - None: auto-detect (tries 'u', 'u1', 'u2', 'u3')
    timestep : int
        Timestep to extract (-1 for last)
    method : str
        Interpolation method: 'linear', 'nearest', 'cubic'
    search_radius : float, optional
        Maximum search radius for nearest neighbor (in same units as coordinates)
        If None, uses global interpolation

    Returns:
    --------
    results : dict
        Dictionary containing:
        - 'positions': Query positions (N, 2) or (N, 3)
        - 'u1': u1 displacement component at each position
        - 'u2': u2 displacement component at each position
        - 'u3': u3 displacement component at each position
        - 'displacement_magnitude': sqrt(u1^2 + u2^2 + u3^2)
        - 'mesh_points': Nearest mesh point for each query position
        - 'distances': Distance to nearest mesh point
    """

    # Load XDMF data
    points, point_data, sx, connect = load_xdmf_surface(xdmf_file, timestep)

    n_vertices = len(points)
    n_cells = len(connect)

    # Auto-detect displacement fields if not specified
    if displacement_fields is None:
        # Check for common field names
        available_fields = list(point_data.keys())

        if 'u' in available_fields:
            displacement_fields = 'u'
        elif 'u1' in available_fields and 'u2' in available_fields and 'u3' in available_fields:
            displacement_fields = ['u1', 'u2', 'u3']
        elif 'displacement' in available_fields:
            displacement_fields = 'displacement'
        else:
            raise ValueError(f"Cannot auto-detect displacement fields. Available: {available_fields}")

        print(f"  Auto-detected displacement fields: {displacement_fields}")

    # Extract displacement components
    if isinstance(displacement_fields, str):
        # Single field containing vector
        if displacement_fields not in point_data:
            raise ValueError(f"Displacement field '{displacement_fields}' not found. "
                           f"Available fields: {list(point_data.keys())}")

        displacements = point_data[displacement_fields]

        # Ensure displacements have 3 components
        if len(displacements.shape) == 1:
            # Scalar field - assume it's a single component
            raise ValueError(f"Field '{displacement_fields}' is scalar. Need vector field or separate u1/u2/u3 fields.")
        elif displacements.shape[1] != 3:
            raise ValueError(f"Expected 3 displacement components, got {displacements.shape[1]}")

    elif isinstance(displacement_fields, list):
        # Separate fields for u1, u2, u3
        if len(displacement_fields) != 3:
            raise ValueError(f"Expected 3 field names, got {len(displacement_fields)}")

        for field in displacement_fields:
            if field not in point_data:
                raise ValueError(f"Field '{field}' not found. Available: {list(point_data.keys())}")

        # Stack into (N, 3) array
        displacements = np.column_stack([
            point_data[displacement_fields[0]],
            point_data[displacement_fields[1]],
            point_data[displacement_fields[2]]
        ])
    else:
        raise ValueError(f"displacement_fields must be str or list, got {type(displacement_fields)}")

    print(f"\nDisplacement fields: {displacement_fields}")
    print(f"  Displacement data shape: {displacements.shape}")

    # Check if data is cell-based or vertex-based
    if len(displacements) == n_cells:
        print(f"  Data is cell-based (converting to vertices)")
        # Convert cell data to vertex data by averaging
        vertex_displacements = np.zeros((n_vertices, 3))
        vertex_counts = np.zeros(n_vertices)

        for i, cell in enumerate(connect):
            for vertex_idx in cell:
                vertex_displacements[vertex_idx] += displacements[i]
                vertex_counts[vertex_idx] += 1

        # Average the contributions
        for i in range(n_vertices):
            if vertex_counts[i] > 0:
                vertex_displacements[i] /= vertex_counts[i]

        displacements = vertex_displacements
        print(f"  Converted to vertex data shape: {displacements.shape}")

    elif len(displacements) == n_vertices:
        print(f"  Data is vertex-based")
    else:
        raise ValueError(f"Displacement data size ({len(displacements)}) doesn't match "
                        f"vertices ({n_vertices}) or cells ({n_cells})")

    print(f"  Range u1: [{displacements[:, 0].min():.6e}, {displacements[:, 0].max():.6e}]")
    print(f"  Range u2: [{displacements[:, 1].min():.6e}, {displacements[:, 1].max():.6e}]")
    print(f"  Range u3: [{displacements[:, 2].min():.6e}, {displacements[:, 2].max():.6e}]")

    # Convert query positions to numpy array
    query_positions = np.asarray(query_positions)

    # Handle 2D (x, y) or 3D (x, y, z) query positions
    if query_positions.ndim == 1:
        query_positions = query_positions.reshape(1, -1)

    n_queries = len(query_positions)
    query_dim = query_positions.shape[1]

    print(f"\nQuery positions: {n_queries} points")
    print(f"  Dimension: {query_dim}D")

    # Extract x, y coordinates from mesh points for interpolation
    mesh_xy = points[:, :2]  # Use only x, y for 2D interpolation

    # For 2D query positions, use x, y
    if query_dim == 2:
        query_xy = query_positions
    else:
        query_xy = query_positions[:, :2]

    # Initialize result arrays
    u1_interp = np.zeros(n_queries)
    u2_interp = np.zeros(n_queries)
    u3_interp = np.zeros(n_queries)

    print(f"\nInterpolation method: {method}")

    if search_radius is not None:
        # Use local interpolation with search radius
        print(f"  Search radius: {search_radius}")

        # Build KD-tree for fast nearest neighbor search
        tree = cKDTree(mesh_xy)

        for i, qpos in enumerate(query_xy):
            # Find points within search radius
            indices = tree.query_ball_point(qpos, search_radius)

            if len(indices) == 0:
                print(f"  Warning: No points found within radius for query {i}")
                u1_interp[i] = np.nan
                u2_interp[i] = np.nan
                u3_interp[i] = np.nan
                continue

            # Use local points for interpolation
            local_points = mesh_xy[indices]
            local_u1 = displacements[indices, 0]
            local_u2 = displacements[indices, 1]
            local_u3 = displacements[indices, 2]

            if method == 'nearest' or len(indices) == 1:
                # Find nearest point
                local_tree = cKDTree(local_points)
                dist, idx = local_tree.query(qpos)
                u1_interp[i] = local_u1[idx]
                u2_interp[i] = local_u2[idx]
                u3_interp[i] = local_u3[idx]
            else:
                # Linear interpolation
                try:
                    u1_interp[i] = griddata(local_points, local_u1, qpos, method='linear')
                    u2_interp[i] = griddata(local_points, local_u2, qpos, method='linear')
                    u3_interp[i] = griddata(local_points, local_u3, qpos, method='linear')
                except:
                    # Fallback to nearest if linear fails
                    local_tree = cKDTree(local_points)
                    dist, idx = local_tree.query(qpos)
                    u1_interp[i] = local_u1[idx]
                    u2_interp[i] = local_u2[idx]
                    u3_interp[i] = local_u3[idx]

    else:
        # Global interpolation
        print("  Using global interpolation")

        if method == 'linear':
            # Use LinearNDInterpolator for better performance
            interp_u1 = LinearNDInterpolator(mesh_xy, displacements[:, 0])
            interp_u2 = LinearNDInterpolator(mesh_xy, displacements[:, 1])
            interp_u3 = LinearNDInterpolator(mesh_xy, displacements[:, 2])

            u1_interp = interp_u1(query_xy)
            u2_interp = interp_u2(query_xy)
            u3_interp = interp_u3(query_xy)

        elif method == 'nearest':
            # Use KD-tree for nearest neighbor
            tree = cKDTree(mesh_xy)
            distances, indices = tree.query(query_xy)

            u1_interp = displacements[indices, 0]
            u2_interp = displacements[indices, 1]
            u3_interp = displacements[indices, 2]

        else:
            # Use griddata for cubic or other methods
            u1_interp = griddata(mesh_xy, displacements[:, 0], query_xy, method=method)
            u2_interp = griddata(mesh_xy, displacements[:, 1], query_xy, method=method)
            u3_interp = griddata(mesh_xy, displacements[:, 2], query_xy, method=method)

    # Find nearest mesh points and distances
    tree = cKDTree(mesh_xy)
    distances, nearest_indices = tree.query(query_xy)
    nearest_points = points[nearest_indices]

    # Calculate displacement magnitude
    displacement_magnitude = np.sqrt(u1_interp**2 + u2_interp**2 + u3_interp**2)

    # Compile results
    results = {
        'positions': query_positions,
        'u1': u1_interp,
        'u2': u2_interp,
        'u3': u3_interp,
        'displacement_magnitude': displacement_magnitude,
        'nearest_mesh_points': nearest_points,
        'distances_to_mesh': distances
    }

    # Print summary
    print(f"\nExtraction complete:")
    print(f"  Valid results: {np.sum(~np.isnan(u1_interp))} / {n_queries}")
    print(f"  Average distance to mesh: {np.nanmean(distances):.6e}")

    return results


def extract_displacements_grid(xdmf_file,
                               x_range,
                               y_range,
                               nx=50,
                               ny=50,
                               displacement_fields=None,
                               timestep=-1,
                               method='linear'):
    """
    Extract displacements on a regular grid for visualization

    Parameters:
    -----------
    xdmf_file : str
        Path to XDMF surface file
    x_range : tuple
        (x_min, x_max) for grid
    y_range : tuple
        (y_min, y_max) for grid
    nx, ny : int
        Number of grid points in x and y directions
    displacement_fields : list or str, optional
        Field names for displacement components (auto-detect if None)
    timestep : int
        Timestep to extract (-1 for last)
    method : str
        Interpolation method

    Returns:
    --------
    results : dict
        Dictionary containing:
        - 'x_grid': 2D array of x coordinates
        - 'y_grid': 2D array of y coordinates
        - 'u1_grid': 2D array of u1 displacements
        - 'u2_grid': 2D array of u2 displacements
        - 'u3_grid': 2D array of u3 displacements
        - 'magnitude_grid': 2D array of displacement magnitudes
    """

    # Create regular grid
    x = np.linspace(x_range[0], x_range[1], nx)
    y = np.linspace(y_range[0], y_range[1], ny)
    x_grid, y_grid = np.meshgrid(x, y)

    # Flatten grid for query
    query_positions = np.column_stack([x_grid.ravel(), y_grid.ravel()])

    # Extract displacements
    results = extract_displacements_at_positions(
        xdmf_file,
        query_positions,
        displacement_fields=displacement_fields,
        timestep=timestep,
        method=method
    )

    # Reshape results to grid
    grid_results = {
        'x_grid': x_grid,
        'y_grid': y_grid,
        'u1_grid': results['u1'].reshape(ny, nx),
        'u2_grid': results['u2'].reshape(ny, nx),
        'u3_grid': results['u3'].reshape(ny, nx),
        'magnitude_grid': results['displacement_magnitude'].reshape(ny, nx)
    }

    return grid_results


# ============================================================================
# Example usage
# ============================================================================

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Example 1: Extract displacements at specific positions
    print("=" * 70)
    print("Example 1: Extract displacements at specific positions")
    print("=" * 70)

    xdmf_file = '/path/to/your/surface.xdmf'

    # Define query positions (x, y)
    query_positions = np.array([
        [0.0, 0.0],
        [1000.0, 0.0],
        [0.0, 1000.0],
        [1000.0, 1000.0],
        [500.0, 500.0]
    ])

    # Extract displacements (will auto-detect field names)
    results = extract_displacements_at_positions(
        xdmf_file,
        query_positions,
        timestep=-1,  # Last timestep
        method='linear'
    )

    # Or explicitly specify field names:
    # results = extract_displacements_at_positions(
    #     xdmf_file,
    #     query_positions,
    #     displacement_fields=['u1', 'u2', 'u3'],  # or 'u' for vector field
    #     timestep=-1,
    #     method='linear'
    # )

    # Print results
    print("\nResults:")
    print(f"{'X':>12} {'Y':>12} {'U1':>12} {'U2':>12} {'U3':>12} {'|U|':>12} {'Dist':>12}")
    print("-" * 84)
    for i in range(len(query_positions)):
        print(f"{query_positions[i, 0]:12.2f} {query_positions[i, 1]:12.2f} "
              f"{results['u1'][i]:12.6e} {results['u2'][i]:12.6e} "
              f"{results['u3'][i]:12.6e} {results['displacement_magnitude'][i]:12.6e} "
              f"{results['distances_to_mesh'][i]:12.6e}")


    # Example 2: Extract displacements on a regular grid
    print("\n" + "=" * 70)
    print("Example 2: Extract displacements on regular grid and plot")
    print("=" * 70)

    grid_results = extract_displacements_grid(
        xdmf_file,
        x_range=(0, 10000),
        y_range=(0, 10000),
        nx=100,
        ny=100,
        timestep=-1,
        method='linear'
    )

    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # U1 component
    im1 = axes[0, 0].pcolormesh(grid_results['x_grid'], grid_results['y_grid'],
                                 grid_results['u1_grid'], shading='auto', cmap='RdBu_r')
    axes[0, 0].set_title('U1 Displacement')
    axes[0, 0].set_xlabel('X')
    axes[0, 0].set_ylabel('Y')
    plt.colorbar(im1, ax=axes[0, 0], label='U1 (m)')

    # U2 component
    im2 = axes[0, 1].pcolormesh(grid_results['x_grid'], grid_results['y_grid'],
                                 grid_results['u2_grid'], shading='auto', cmap='RdBu_r')
    axes[0, 1].set_title('U2 Displacement')
    axes[0, 1].set_xlabel('X')
    axes[0, 1].set_ylabel('Y')
    plt.colorbar(im2, ax=axes[0, 1], label='U2 (m)')

    # U3 component
    im3 = axes[1, 0].pcolormesh(grid_results['x_grid'], grid_results['y_grid'],
                                 grid_results['u3_grid'], shading='auto', cmap='RdBu_r')
    axes[1, 0].set_title('U3 Displacement')
    axes[1, 0].set_xlabel('X')
    axes[1, 0].set_ylabel('Y')
    plt.colorbar(im3, ax=axes[1, 0], label='U3 (m)')

    # Magnitude
    im4 = axes[1, 1].pcolormesh(grid_results['x_grid'], grid_results['y_grid'],
                                 grid_results['magnitude_grid'], shading='auto', cmap='viridis')
    axes[1, 1].set_title('Displacement Magnitude')
    axes[1, 1].set_xlabel('X')
    axes[1, 1].set_ylabel('Y')
    plt.colorbar(im4, ax=axes[1, 1], label='|U| (m)')

    plt.tight_layout()
    plt.savefig('displacement_components.png', dpi=300, bbox_inches='tight')
    print("\nSaved plot to: displacement_components.png")

    plt.show()
