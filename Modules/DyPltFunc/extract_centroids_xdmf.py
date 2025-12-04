"""
Function to extract centroid positions of triangular elements from XDMF file.
Filters elements at specific depth (z=-10 km) with lateral spacing (every 10 km).
"""

import numpy as np
import seissolxdmf


def extract_triangle_centroids(
    xdmf_file, target_depth=-10.0, lateral_spacing=30.0, depth_tolerance=0.05
):
    """
    Extract centroid positions of triangular elements from an XDMF file.

    Filters elements at a specific depth with regular lateral spacing.

    Parameters
    ----------
    xdmf_file : str
        Path to the XDMF file
    target_depth : float, optional
        Target depth in km (default: -10.0 km, negative = below surface)
    lateral_spacing : float, optional
        Desired lateral spacing between selected centroids in km (default: 10.0)
    depth_tolerance : float, optional
        Tolerance for depth filtering in km (default: 0.5)

    Returns
    -------
    dict
        Dictionary containing:
        - 'centroids': (N, 3) array of centroid coordinates [x, y, z]
        - 'element_ids': indices of selected elements
        - 'all_centroids': (M, 3) array of all centroids at target depth
        - 'n_total': total number of elements at target depth
        - 'n_selected': number of selected elements after spacing filter

    Examples
    --------
    >>> results = extract_triangle_centroids('fault_mesh.xdmf',
    ...                                       target_depth=-10.0,
    ...                                       lateral_spacing=10.0)
    >>> print(f"Selected {results['n_selected']} elements")
    >>> centroids = results['centroids']
    """

    # Load XDMF file using seissolxdmf module
    print(f"Loading XDMF file: {xdmf_file}")
    sx = seissolxdmf.seissolxdmf(xdmf_file)

    # Read geometry (vertices)
    vertices = sx.ReadGeometry()  # Shape: (n_vertices, 3) - [x, y, z]

    # Read connectivity (cells/elements)
    cells = sx.ReadConnect()  # Shape: (n_cells, 3) - triangle vertex indices

    # fault_tag = sx.ReadData('fault-tag',idt=1)

    print(f"\nMesh statistics:")
    print(f"  Total vertices: {len(vertices)}")
    print(f"  Total triangular elements: {len(cells)}")

    # Calculate centroids for all triangular elements
    # Centroid = average of three vertices
    v0 = vertices[cells[:, 0]]  # First vertex of each triangle
    v1 = vertices[cells[:, 1]]  # Second vertex
    v2 = vertices[cells[:, 2]]  # Third vertex

    all_centroids = (v0 + v1 + v2) / 3.0  # Shape: (n_cells, 3)

    # Convert depth to km if vertices are in meters
    # Assuming z-coordinates are in meters, convert to km
    if np.abs(all_centroids[:, 2]).max() > 100:  # Likely in meters
        all_centroids[:, 2] /= 1000.0
        vertices[:, 2] /= 1000.0
        print(f"  Converted z-coordinates from meters to kilometers")

    print(
        f"  Depth range: {all_centroids[:, 2].min():.2f} to {all_centroids[:, 2].max():.2f} km"
    )

    # Filter elements at target depth
    depth_mask = np.abs(all_centroids[:, 2] - target_depth) < depth_tolerance
    centroids_at_depth = all_centroids[depth_mask]
    element_ids_at_depth = np.where(depth_mask)[0]

    print(f"\nFiltering at depth z = {target_depth} ± {depth_tolerance} km:")
    print(f"  Elements at target depth: {len(centroids_at_depth)}")

    if len(centroids_at_depth) == 0:
        print("  WARNING: No elements found at target depth!")
        return {
            "centroids": np.array([]),
            "element_ids": np.array([]),
            "all_centroids": np.array([]),
            "n_total": 0,
            "n_selected": 0,
        }

    # Apply lateral spacing filter using greedy selection
    selected_centroids, selected_ids = apply_lateral_spacing(
        centroids_at_depth, element_ids_at_depth, lateral_spacing
    )

    print(f"\nApplying lateral spacing of {lateral_spacing} km:")
    print(f"  Selected elements: {len(selected_centroids)}")
    print(f"  Spatial extent:")
    print(
        f"    X: {selected_centroids[:, 0].min():.2f} to {selected_centroids[:, 0].max():.2f} km"
    )
    print(
        f"    Y: {selected_centroids[:, 1].min():.2f} to {selected_centroids[:, 1].max():.2f} km"
    )

    return {
        "centroids": selected_centroids,
        "element_ids": selected_ids,
        "all_centroids": centroids_at_depth,
        "n_total": len(centroids_at_depth),
        "n_selected": len(selected_centroids),
        # "fault-tag": fault_tag[selected_ids],
    }


def apply_lateral_spacing(centroids, element_ids, spacing):
    """
    Apply lateral spacing filter to select one element per spatial bin.

    Creates a 2D grid with bins of size 'spacing' and selects one element
    (closest to bin center) from each bin that contains elements.

    Parameters
    ----------
    centroids : ndarray
        (N, 3) array of centroid coordinates
    element_ids : ndarray
        (N,) array of element indices
    spacing : float
        Lateral spacing (bin size) in km

    Returns
    -------
    selected_centroids : ndarray
        (M, 3) array of selected centroids
    selected_ids : ndarray
        (M,) array of corresponding element IDs
    """

    # Extract only x-y coordinates for lateral binning
    xy_coords = centroids[:, :2].copy()

    # Get spatial extent
    x_min, y_min = xy_coords.min(axis=0)
    x_max, y_max = xy_coords.max(axis=0)

    # Create bins
    x_bins = np.arange(x_min, x_max + spacing, spacing)
    y_bins = np.arange(y_min, y_max + spacing, spacing)

    print(f"  Creating grid: {len(x_bins)-1} x {len(y_bins)-1} bins")

    # Assign each point to a bin
    x_indices = np.digitize(xy_coords[:, 0], x_bins) - 1
    y_indices = np.digitize(xy_coords[:, 1], y_bins) - 1

    # Create dictionary to store points in each bin
    bins_dict = {}
    for i in range(len(centroids)):
        bin_key = (x_indices[i], y_indices[i])
        if bin_key not in bins_dict:
            bins_dict[bin_key] = []
        bins_dict[bin_key].append(i)

    # Select one element per bin (the one closest to bin center)
    selected_indices = []

    for (x_idx, y_idx), indices in bins_dict.items():
        if len(indices) == 0:
            continue

        # Calculate bin center
        bin_center_x = x_bins[x_idx] + spacing / 2
        bin_center_y = y_bins[y_idx] + spacing / 2
        bin_center = np.array([bin_center_x, bin_center_y])

        # Find element closest to bin center
        distances = np.linalg.norm(xy_coords[indices] - bin_center, axis=1)
        closest_idx = indices[np.argmin(distances)]
        selected_indices.append(closest_idx)

    selected_indices = np.array(selected_indices)

    print(f"  Populated bins: {len(bins_dict)}")
    print(f"  Selected elements: {len(selected_indices)}")

    return centroids[selected_indices], element_ids[selected_indices]


def plot_centroids(results, save_path=None):
    """
    Create visualization of selected centroids.

    Parameters
    ----------
    results : dict
        Output from extract_triangle_centroids()
    save_path : str, optional
        Path to save figure
    """
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(14, 5))

    # 2D plan view
    ax1 = fig.add_subplot(121)

    if len(results["all_centroids"]) > 0:
        ax1.scatter(
            results["all_centroids"][:, 0],
            results["all_centroids"][:, 1],
            c="lightgray",
            s=10,
            alpha=0.3,
            label="All elements at depth",
        )

    if len(results["centroids"]) > 0:
        ax1.scatter(
            results["centroids"][:, 0],
            results["centroids"][:, 1],
            c="red",
            s=80,
            marker="o",
            edgecolors="darkred",
            linewidth=1.5,
            alpha=0.8,
            label="Selected centroids",
            zorder=3,
        )

    ax1.set_xlabel("X (km)", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Y (km)", fontsize=12, fontweight="bold")
    ax1.set_title(
        "Plan View - Selected Centroids", fontsize=13, fontweight="bold", pad=15
    )
    ax1.legend(frameon=True, shadow=True, fontsize=10)
    ax1.grid(True, alpha=0.3, linestyle=":")
    ax1.set_aspect("equal", adjustable="box")

    # 3D view
    ax2 = fig.add_subplot(122, projection="3d")

    if len(results["centroids"]) > 0:
        ax2.scatter(
            results["centroids"][:, 0],
            results["centroids"][:, 1],
            results["centroids"][:, 2],
            c="red",
            s=100,
            marker="o",
            edgecolors="darkred",
            linewidth=1.5,
            alpha=0.8,
        )

    ax2.set_xlabel("X (km)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Y (km)", fontsize=11, fontweight="bold")
    ax2.set_zlabel("Z (km)", fontsize=11, fontweight="bold")
    ax2.set_title("3D View", fontsize=13, fontweight="bold", pad=15)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"\nFigure saved to: {save_path}")

    plt.show()


# Example usage
if __name__ == "__main__":
    # Example call
    xdmf_file = "fault_mesh.xdmf"

    results = extract_triangle_centroids(
        xdmf_file,
        target_depth=-10.0,  # 10 km depth
        lateral_spacing=10.0,  # 10 km spacing
        depth_tolerance=0.5,  # ±0.5 km tolerance
    )

    # Plot results
    plot_centroids(results, save_path="centroids_z10km.png")

    # Access results
    print(f"\nCentroid coordinates shape: {results['centroids'].shape}")
    print(f"First 5 centroids:\n{results['centroids'][:5]}")
    print(f"Element IDs: {results['element_ids'][:10]}...")
