"""
Basin receiver placement utilities.

This module provides functions to generate receiver locations in the Wellington Basin,
including receivers at the basin bottom with adaptive spacing.
"""

import numpy as np


def find_basin_bottom_depth(
    grdvel, grdz, basin_threshold=800.0, bedrock_threshold=1340.0
):
    """
    Find the actual basin bottom by detecting the sediment-bedrock interface.

    This function handles geological structures where basin sediments may be
    buried below water layers or other low-velocity zones. It looks for the
    deepest continuous low-velocity zone that represents basin sediments.

    Parameters:
    -----------
    grdvel : ndarray
        3D velocity array (Nz, Ny, Nx)
    grdz : ndarray
        Z coordinate array (depth values)
    basin_threshold : float, default=1000.0
        Velocity threshold below which material is considered basin sediment
    bedrock_threshold : float, default=1500.0
        Velocity threshold above which material is considered bedrock

    Returns:
    --------
    ndarray
        2D array of basin bottom depths (Ny, Nx)
    """
    _, Ny, Nx = grdvel.shape
    basin_bottom = np.zeros((Ny, Nx))

    for i in range(Nx):
        for j in range(Ny):
            velocity_profile = grdvel[:, j, i]

            # Find the actual basin bottom by looking for sediment-bedrock transition
            basin_bottom_depth = find_sediment_bedrock_interface(
                velocity_profile, grdz, basin_threshold, bedrock_threshold
            )

            basin_bottom[j, i] = basin_bottom_depth

    return basin_bottom


def find_sediment_bedrock_interface(
    velocity_profile, grdz, basin_threshold=800.0, bedrock_threshold=1500.0
):
    """
    Find the sediment-bedrock interface in a 1D velocity profile.

    This looks for the deepest point where we transition from basin sediments
    (< basin_threshold) to bedrock (> bedrock_threshold).

    Parameters:
    -----------
    velocity_profile : ndarray
        1D velocity profile with depth
    grdz : ndarray
        Depth coordinates
    basin_threshold : float
        Velocity threshold for basin sediment
    bedrock_threshold : float
        Velocity threshold for bedrock

    Returns:
    --------
    float
        Depth of basin bottom
    """
    # Find all basin sediment points
    basin_mask = velocity_profile < basin_threshold
    bedrock_mask = velocity_profile > bedrock_threshold

    if not np.any(basin_mask):
        # No basin sediment found
        return grdz[0]  # Return surface

    # Look for transitions from basin to bedrock
    basin_bottom_idx = 0

    # Method 1: Find the deepest continuous basin sequence before hitting bedrock
    basin_indices = np.where(basin_mask)[0]

    if len(basin_indices) == 0:
        return grdz[0]

    # Find the deepest basin point
    deepest_basin_idx = basin_indices[-1]

    # Check if there's a clear transition to bedrock below this point
    if deepest_basin_idx < len(velocity_profile) - 1:
        # Look ahead to see if velocities increase to bedrock levels
        remaining_profile = velocity_profile[deepest_basin_idx + 1 :]
        if len(remaining_profile) > 0:
            # Check if we hit bedrock velocities within the next few layers
            next_few_layers = remaining_profile[: min(3, len(remaining_profile))]
            if np.any(next_few_layers > bedrock_threshold):
                # Clear transition to bedrock found
                basin_bottom_idx = deepest_basin_idx
            else:
                # No clear bedrock transition, use deepest basin point
                basin_bottom_idx = deepest_basin_idx
        else:
            basin_bottom_idx = deepest_basin_idx
    else:
        # At the bottom of the model
        basin_bottom_idx = deepest_basin_idx

    # Alternative method: Look for the steepest velocity gradient
    # This can help identify the sediment-bedrock boundary more precisely
    if len(velocity_profile) > 3:
        velocity_gradient = np.diff(velocity_profile)
        # Find the largest positive gradient (biggest velocity increase with depth)
        # within the basin sediment zone

        # Only consider gradients within or just below basin sediment zones
        gradient_search_zone = []
        for idx in basin_indices:
            if idx < len(velocity_gradient):
                gradient_search_zone.append(idx)
            if idx + 1 < len(velocity_gradient):
                gradient_search_zone.append(idx + 1)

        if gradient_search_zone:
            search_zone = list(set(gradient_search_zone))  # Remove duplicates
            search_zone.sort()

            max_gradient_idx = None
            max_gradient = 0

            for idx in search_zone:
                if (
                    velocity_gradient[idx] > max_gradient
                    and velocity_gradient[idx] > 200
                ):  # Significant increase
                    max_gradient = velocity_gradient[idx]
                    max_gradient_idx = idx

            # If we found a significant velocity increase, use that as basin bottom
            if max_gradient_idx is not None:
                # The basin bottom is at the depth before the big velocity jump
                basin_bottom_idx = max_gradient_idx

    return grdz[basin_bottom_idx]


def generate_basin_bottom_receivers(
    grdx, grdy, grdz, grdvel, spacing=100.0, depth_offset=10.0, basin_threshold=800.0
):
    """
    Generate receiver locations at the bottom of the basin with specified spacing.

    Parameters:
    -----------
    grdx, grdy : ndarray
        X and Y coordinate arrays
    grdz : ndarray
        Z coordinate array (depth values)
    grdvel : ndarray
        3D velocity array (Nz, Ny, Nx)
    spacing : float, default=100.0
        Desired spacing between receivers in meters
    depth_offset : float, default=10.0
        Depth below basin bottom to place receivers (meters)
    basin_threshold : float, default=1000.0
        Velocity threshold to identify basin sediment

    Returns:
    --------
    tuple
        (receiver_coords, basin_depths) where:
        - receiver_coords: ndarray of shape (N, 3) with Z, X, Y coordinates
        - basin_depths: ndarray of basin bottom depths at receiver locations
    """
    # Find basin bottom depths
    basin_bottom = find_basin_bottom_depth(grdvel, grdz, basin_threshold)

    # Create grid with desired spacing
    x_min, x_max = grdx.min(), grdx.max()
    y_min, y_max = grdy.min(), grdy.max()

    # Generate receiver grid points
    x_receivers = np.arange(x_min, x_max + spacing, spacing)
    y_receivers = np.arange(y_min, y_max + spacing, spacing)

    receiver_coords = []
    basin_depths = []

    for x_rec in x_receivers:
        for y_rec in y_receivers:
            # Find nearest grid point in velocity model
            i_idx = np.argmin(np.abs(grdx - x_rec))
            j_idx = np.argmin(np.abs(grdy - y_rec))

            # Get basin bottom depth at this location
            bottom_depth = basin_bottom[j_idx, i_idx]

            # Check if there's a basin anywhere in the velocity column (not just at surface)
            has_basin = np.any(grdvel[:, j_idx, i_idx] < basin_threshold)
            if has_basin:
                # Place receiver depth_offset below basin bottom
                receiver_depth = bottom_depth - depth_offset

                receiver_coords.append([receiver_depth, grdx[i_idx], grdy[j_idx]])
                basin_depths.append(bottom_depth)

    return np.array(receiver_coords), np.array(basin_depths)


def generate_basin_profile_receivers(
    grdx,
    grdy,
    grdz,
    grdvel,
    profile_y=None,
    spacing=100.0,
    depth_offset=10.0,
    basin_threshold=800.0,
):
    """
    Generate receivers along a specific Y profile across the basin.

    Parameters:
    -----------
    grdx, grdy : ndarray
        X and Y coordinate arrays
    grdz : ndarray
        Z coordinate array
    grdvel : ndarray
        3D velocity array (Nz, Ny, Nx)
    profile_y : float, optional
        Y coordinate for profile. If None, uses middle of domain
    spacing : float, default=50.0
        Spacing between receivers along profile
    depth_offset : float, default=10.0
        Depth below basin bottom to place receivers
    basin_threshold : float, default=1000.0
        Velocity threshold to identify basin sediment

    Returns:
    --------
    ndarray
        Array of receiver coordinates (N, 3) with Z, X, Y coordinates
    """
    if profile_y is None:
        profile_y = (grdy.min() + grdy.max()) / 2

    # Find nearest Y index
    j_idx = np.argmin(np.abs(grdy - profile_y))
    actual_y = grdy[j_idx]

    # Generate X coordinates with desired spacing
    x_min, x_max = grdx.min(), grdx.max()
    x_profile = np.arange(x_min, x_max + spacing, spacing)

    receiver_coords = []

    for x_rec in x_profile:
        # Find nearest X index
        i_idx = np.argmin(np.abs(grdx - x_rec))

        # Check if this location has basin sediment anywhere in the column
        if np.any(grdvel[:, j_idx, i_idx] < basin_threshold):
            # Find basin bottom depth
            basin_mask = grdvel[:, j_idx, i_idx] < basin_threshold
            if np.any(basin_mask):
                basin_indices = np.where(basin_mask)[0]
                deepest_basin_idx = basin_indices[-1]
                bottom_depth = grdz[deepest_basin_idx]

                # Place receiver below basin bottom
                receiver_depth = bottom_depth - depth_offset
                receiver_coords.append([receiver_depth, grdx[i_idx], actual_y])

    return np.array(receiver_coords)


def save_receivers_xyz(receiver_coords, filename, labels=None):
    """
    Save receiver coordinates to text file in XYZ format.

    Parameters:
    -----------
    receiver_coords : ndarray
        Array of receiver coordinates (N, 3)
    filename : str
        Output filename
    labels : list, optional
        List of receiver labels/names
    """
    if labels is None:
        labels = [f"REC_{i:04d}" for i in range(len(receiver_coords))]

    with open(filename, "w") as f:
        f.write("# X Y Z Label\n")
        for i, (x, y, z) in enumerate(receiver_coords):
            f.write(f"{x:.1f} {y:.1f} {z:.1f} {labels[i]}\n")

    print(f"Saved {len(receiver_coords)} receivers to {filename}")


def visualize_receivers_2d(
    grdx, grdy, grdvel, receiver_coords, depth_slice=0, basin_threshold=800.0
):
    """
    Create a 2D visualization of receiver locations on basin map.

    Parameters:
    -----------
    grdx, grdy : ndarray
        X and Y coordinate arrays
    grdvel : ndarray
        3D velocity array
    receiver_coords : ndarray
        Receiver coordinates (N, 3)
    depth_slice : int, default=0
        Which depth slice to show for basin structure
    basin_threshold : float, default=500.0
        Velocity threshold for basin visualization

    Returns:
    --------
    tuple
        (fig, ax) matplotlib objects
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 8))

    # Show basin structure (surface velocities)
    extent = [grdx.min(), grdx.max(), grdy.min(), grdy.max()]
    im = ax.imshow(
        grdvel[depth_slice, :, :],
        origin="lower",
        extent=extent,
        cmap="viridis",
        vmax=basin_threshold,
    )

    # Plot receivers
    if len(receiver_coords) > 0:
        ax.scatter(
            receiver_coords[:, 0],
            receiver_coords[:, 1],
            c="red",
            s=20,
            marker="^",
            label=f"{len(receiver_coords)} Basin Receivers",
        )

    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_title("Basin Bottom Receivers")
    ax.legend()

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Velocity (m/s)")

    return fig, ax
