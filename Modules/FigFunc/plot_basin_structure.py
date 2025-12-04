"""
Basin structure plotting utilities for geological visualization with receiver locations.

This module provides functions to create geological maps and cross-sections
showing basin structure and receiver placements.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def plot_geological_map_with_receivers(
    grdx,
    grdy,
    grdvel,
    receiver_coords,
    depth_slice=0,
    basin_threshold=800.0,
    figsize=(12, 10),
    save_path=None,
):
    """
    Create a geological map visualization with receiver locations.

    Parameters:
    -----------
    grdx, grdy : ndarray
        X and Y coordinate arrays
    grdvel : ndarray
        3D velocity array
    receiver_coords : ndarray
        Receiver coordinates (N, 3) with Z, X, Y format
    depth_slice : int, default=0
        Which depth slice to show for geological structure
    basin_threshold : float, default=800.0
        Velocity threshold for basin/bedrock distinction
    figsize : tuple, default=(12, 10)
        Figure size
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    tuple
        (fig, ax) matplotlib objects
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Create extent for the map
    extent = [grdx.min(), grdx.max(), grdy.min(), grdy.max()]

    # Get velocity slice for geological visualization
    vel_slice = grdvel[depth_slice, :, :]

    # Create custom colormap for geological interpretation
    # Low velocities (sediments) in warm colors, high velocities (bedrock) in cool colors
    colors = [
        "#8B4513",
        "#DAA520",
        "#228B22",
        "#4682B4",
        "#191970",
    ]  # Brown to dark blue
    n_bins = 100
    cmap = mcolors.LinearSegmentedColormap.from_list("geological", colors, N=n_bins)

    # Plot velocity structure
    im = ax.imshow(
        vel_slice, origin="lower", extent=extent, cmap=cmap, aspect="equal", alpha=0.8
    )

    # Add contour lines for better geological interpretation
    contour_levels = [400, 800, 1200, 1600, 2000, 2500, 3000]
    contours = ax.contour(
        vel_slice,
        levels=contour_levels,
        extent=extent,
        colors="black",
        alpha=0.3,
        linewidths=0.5,
    )
    ax.clabel(contours, inline=True, fontsize=8, fmt="%d m/s")

    # Highlight basin boundary
    basin_contour = ax.contour(
        vel_slice,
        levels=[basin_threshold],
        extent=extent,
        colors="red",
        linewidths=2,
        alpha=0.7,
    )
    ax.clabel(basin_contour, inline=True, fontsize=10, fmt="Basin boundary (%d m/s)")

    # Plot receivers if provided
    if len(receiver_coords) > 0:
        # Extract X, Y coordinates (remember format is Z, X, Y)
        recv_x = receiver_coords[:, 1]  # X coordinates
        recv_y = receiver_coords[:, 2]  # Y coordinates
        recv_z = receiver_coords[:, 0]  # Z coordinates (depths)

        # Plot receivers with color coding by depth
        scatter = ax.scatter(
            recv_x,
            recv_y,
            c=recv_z,
            s=30,
            marker="^",
            cmap="plasma_r",
            edgecolors="white",
            linewidths=0.5,
            alpha=0.9,
            label=f"{len(receiver_coords)} Basin Receivers",
        )

        # Add colorbar for receiver depths
        cbar_recv = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)
        cbar_recv.set_label("Receiver Depth (m)", rotation=270, labelpad=15)

    # Formatting
    ax.set_xlabel("Easting (m)", fontsize=12)
    ax.set_ylabel("Northing (m)", fontsize=12)
    ax.set_title(
        f"Geological Map with Basin Receivers\n(Depth Slice: {depth_slice}, Velocity at surface)",
        fontsize=14,
        fontweight="bold",
    )

    # Add grid
    ax.grid(True, alpha=0.3, linestyle="--")

    # Add legend
    if len(receiver_coords) > 0:
        ax.legend(loc="upper right", framealpha=0.9)

    # Add colorbar for velocity structure
    cbar_vel = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar_vel.set_label("Velocity (m/s)", rotation=270, labelpad=15)

    # Add geological interpretation text
    textstr = f"""Geological Interpretation:
    • Brown/Yellow: Basin sediments (< {basin_threshold} m/s)
    • Green/Blue: Bedrock (> {basin_threshold} m/s)
    • Red line: Basin-bedrock boundary
    • Triangles: Receiver locations"""

    props = dict(boxstyle="round", facecolor="wheat", alpha=0.8)
    ax.text(
        0.02,
        0.98,
        textstr,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=props,
    )

    plt.tight_layout()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Geological map saved to: {save_path}")

    return fig, ax


def plot_basin_cross_section(
    grdx,
    grdy,
    grdz,
    grdvel,
    receiver_coords,
    profile_y=None,
    basin_threshold=800.0,
    figsize=(14, 8),
    save_path=None,
):
    """
    Create a cross-section view of the basin with receiver locations.

    Parameters:
    -----------
    grdx, grdy, grdz : ndarray
        X, Y, Z coordinate arrays
    grdvel : ndarray
        3D velocity array
    receiver_coords : ndarray
        Receiver coordinates (N, 3) with Z, X, Y format
    profile_y : float, optional
        Y coordinate for cross-section. If None, uses middle of domain
    basin_threshold : float, default=800.0
        Velocity threshold for basin visualization
    figsize : tuple, default=(14, 8)
        Figure size
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    tuple
        (fig, ax) matplotlib objects
    """
    if profile_y is None:
        profile_y = (grdy.min() + grdy.max()) / 2

    # Find nearest Y index
    j_idx = np.argmin(np.abs(grdy - profile_y))
    actual_y = grdy[j_idx]

    fig, ax = plt.subplots(figsize=figsize)

    # Extract velocity profile at this Y location
    vel_profile = grdvel[:, j_idx, :]  # Shape: (Nz, Nx)

    # Create meshgrid for plotting
    X, Z = np.meshgrid(grdx, grdz)

    # Create custom colormap
    colors = ["#8B4513", "#DAA520", "#228B22", "#4682B4", "#191970"]
    cmap = mcolors.LinearSegmentedColormap.from_list("geological", colors, N=100)

    # Plot velocity cross-section
    im = ax.pcolormesh(X, Z, vel_profile, cmap=cmap, shading="auto")

    # Add velocity contours
    contour_levels = [400, 800, 1200, 1600, 2000, 2500, 3000]
    contours = ax.contour(
        X,
        Z,
        vel_profile,
        levels=contour_levels,
        colors="black",
        alpha=0.5,
        linewidths=0.8,
    )
    ax.clabel(contours, inline=True, fontsize=8, fmt="%d")

    # Highlight basin boundary
    basin_contour = ax.contour(
        X,
        Z,
        vel_profile,
        levels=[basin_threshold],
        colors="red",
        linewidths=3,
        alpha=0.8,
    )

    # Plot receivers in this cross-section
    if len(receiver_coords) > 0:
        # Filter receivers near this Y profile (within 100m)
        recv_mask = np.abs(receiver_coords[:, 2] - actual_y) < 100
        if np.any(recv_mask):
            recv_near = receiver_coords[recv_mask]
            ax.scatter(
                recv_near[:, 1],
                recv_near[:, 0],
                c="yellow",
                s=50,
                marker="v",
                edgecolors="red",
                linewidths=1,
                label=f"Receivers near Y={actual_y:.0f}m",
            )

    # Formatting
    ax.set_xlabel("Easting (m)", fontsize=12)
    ax.set_ylabel("Depth (m)", fontsize=12)
    ax.set_title(
        f"Basin Cross-Section at Y = {actual_y:.0f} m\nVelocity Structure and Receiver Placement",
        fontsize=14,
        fontweight="bold",
    )

    # Set Y axis to show depth correctly (depth increases downward)
    # ax.set_ylim(grdz.max(), grdz.min())

    # Add grid
    ax.grid(True, alpha=0.3, linestyle="--")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("shear velocity (m/s)", rotation=270, labelpad=15)

    # Add legend if receivers are plotted
    if len(receiver_coords) > 0:
        ax.legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()
    ax.invert_yaxis()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Cross-section saved to: {save_path}")

    return fig, ax


def plot_3d_basin_structure(
    grdx,
    grdy,
    grdz,
    grdvel,
    receiver_coords=None,
    basin_threshold=800.0,
    subsample=2,
    figsize=(12, 9),
    save_path=None,
):
    """
    Create a 3D visualization of basin structure with optional receiver locations.

    Parameters:
    -----------
    grdx, grdy, grdz : ndarray
        X, Y, Z coordinate arrays
    grdvel : ndarray
        3D velocity array
    receiver_coords : ndarray, optional
        Receiver coordinates (N, 3) with Z, X, Y format
    basin_threshold : float, default=800.0
        Velocity threshold for basin visualization
    subsample : int, default=2
        Subsampling factor to reduce data size for 3D plotting
    figsize : tuple, default=(12, 9)
        Figure size
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    tuple
        (fig, ax) matplotlib objects
    """
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    # Subsample data for performance
    x_sub = grdx[::subsample]
    y_sub = grdy[::subsample]
    z_sub = grdz[::subsample]
    vel_sub = grdvel[::subsample, ::subsample, ::subsample]

    # Create 3D meshgrid
    X, Y, Z = np.meshgrid(x_sub, y_sub, z_sub, indexing="ij")

    # Find basin sediment locations
    basin_mask = vel_sub < basin_threshold

    # Plot basin sediments as scatter points
    if np.any(basin_mask):
        basin_x = X[basin_mask]
        basin_y = Y[basin_mask]
        basin_z = Z[basin_mask]
        basin_vel = vel_sub[basin_mask]

        scatter = ax.scatter(
            basin_x,
            basin_y,
            basin_z,
            c=basin_vel,
            cmap="Reds_r",
            s=10,
            alpha=0.6,
            label="Basin Sediments",
        )

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=20)
        cbar.set_label("Velocity (m/s)", rotation=270, labelpad=15)

    # Plot receivers if provided
    if receiver_coords is not None and len(receiver_coords) > 0:
        recv_x = receiver_coords[:, 1]
        recv_y = receiver_coords[:, 2]
        recv_z = receiver_coords[:, 0]

        ax.scatter(
            recv_x,
            recv_y,
            recv_z,
            c="blue",
            s=50,
            marker="^",
            alpha=0.8,
            label="Receivers",
        )

    # Formatting
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_zlabel("Depth (m)")
    ax.set_title(
        f"3D Basin Structure\n(Basin threshold: {basin_threshold} m/s)",
        fontsize=14,
        fontweight="bold",
    )

    # Invert Z axis (depth increases downward)
    ax.invert_zaxis()

    # Add legend
    ax.legend()

    plt.tight_layout()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"3D plot saved to: {save_path}")

    return fig, ax


def plot_velocity_histogram(
    grdvel, basin_threshold=800.0, figsize=(10, 6), save_path=None
):
    """
    Plot histogram of velocity distribution to help understand geological structure.

    Parameters:
    -----------
    grdvel : ndarray
        3D velocity array
    basin_threshold : float, default=800.0
        Velocity threshold for basin/bedrock distinction
    figsize : tuple, default=(10, 6)
        Figure size
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    tuple
        (fig, ax) matplotlib objects
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Flatten velocity array
    velocities = grdvel.flatten()

    # Remove any NaN or infinite values
    velocities = velocities[np.isfinite(velocities)]

    # Create histogram
    n_bins = 50
    counts, bins, patches = ax.hist(
        velocities, bins=n_bins, alpha=0.7, edgecolor="black", linewidth=0.5
    )

    # Color bars based on basin threshold
    for i, (patch, bin_val) in enumerate(zip(patches, bins[:-1])):
        if bin_val < basin_threshold:
            patch.set_facecolor("#8B4513")  # Brown for basin sediments
        else:
            patch.set_facecolor("#4682B4")  # Blue for bedrock

    # Add vertical line at threshold
    ax.axvline(
        basin_threshold,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Basin threshold ({basin_threshold} m/s)",
    )

    # Formatting
    ax.set_xlabel("Velocity (m/s)", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title("Velocity Distribution in Basin Model", fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Add statistics text
    basin_mask = velocities < basin_threshold
    basin_percent = np.sum(basin_mask) / len(velocities) * 100

    stats_text = f"""Statistics:
    Total points: {len(velocities):,}
    Basin sediments: {np.sum(basin_mask):,} ({basin_percent:.1f}%)
    Bedrock: {np.sum(~basin_mask):,} ({100-basin_percent:.1f}%)

    Velocity range: {velocities.min():.0f} - {velocities.max():.0f} m/s
    Mean velocity: {velocities.mean():.0f} m/s
    Median velocity: {np.median(velocities):.0f} m/s"""

    props = dict(boxstyle="round", facecolor="wheat", alpha=0.8)
    ax.text(
        0.98,
        0.98,
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=props,
    )

    plt.tight_layout()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Histogram saved to: {save_path}")

    return fig, ax


def plot_current_cell_cross_section(
    grdx,
    grdy,
    grdz,
    grdvel,
    receiver_coords=None,
    profile_y=None,
    basin_threshold=800.0,
    figsize=(15, 10),
    save_path=None,
):
    """
    Plot basin cross-section using current Jupyter cell data with enhanced visualization.

    Parameters:
    -----------
    grdx, grdy, grdz : ndarray
        X, Y, Z coordinate arrays from current cell
    grdvel : ndarray
        3D velocity array from current cell
    receiver_coords : ndarray, optional
        Receiver coordinates (N, 3) with Z, X, Y format
    profile_y : float, optional
        Y coordinate for cross-section. If None, uses middle of domain
    basin_threshold : float, default=800.0
        Velocity threshold for basin visualization
    figsize : tuple, default=(15, 10)
        Figure size
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    tuple
        (fig, (ax1, ax2)) matplotlib objects - main plot and depth profile
    """
    if profile_y is None:
        profile_y = (grdy.min() + grdy.max()) / 2

    # Find nearest Y index
    j_idx = np.argmin(np.abs(grdy - profile_y))
    actual_y = grdy[j_idx]

    # Create subplot layout
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, height_ratios=[3, 1])

    # Main cross-section plot (top)
    vel_profile = grdvel[:, j_idx, :]  # Shape: (Nz, Nx)
    X, Z = np.meshgrid(grdx, grdz)

    # Custom colormap for geological interpretation
    colors = ["#8B4513", "#DAA520", "#228B22", "#4682B4", "#191970"]
    cmap = mcolors.LinearSegmentedColormap.from_list("geological", colors, N=100)

    # Plot velocity cross-section
    im = ax1.pcolormesh(X, Z, vel_profile, cmap=cmap, shading="auto")

    # Add velocity contours
    contour_levels = [200, 400, 600, 800, 1000, 1200, 1600, 2000, 2500, 3000]
    contours = ax1.contour(
        X,
        Z,
        vel_profile,
        levels=contour_levels,
        colors="white",
        alpha=0.6,
        linewidths=0.8,
    )
    ax1.clabel(contours, inline=True, fontsize=8, fmt="%d", colors="white")

    # Highlight basin boundary with thick red line
    basin_contour = ax1.contour(
        X,
        Z,
        vel_profile,
        levels=[basin_threshold],
        colors="red",
        linewidths=3,
        alpha=0.9,
    )

    # Plot receivers in cross-section
    if receiver_coords is not None and len(receiver_coords) > 0:
        # Filter receivers near this Y profile (within 200m)
        recv_mask = np.abs(receiver_coords[:, 2] - actual_y) < 200
        if np.any(recv_mask):
            recv_near = receiver_coords[recv_mask]
            ax1.scatter(
                recv_near[:, 1],
                recv_near[:, 0],
                c="yellow",
                s=60,
                marker="v",
                edgecolors="red",
                linewidths=1.5,
                alpha=0.9,
                label=f"Receivers (±200m from Y={actual_y:.0f})",
                zorder=10,
            )

    # Format main plot
    ax1.set_ylabel("Depth (m)", fontsize=12)
    ax1.set_title(
        f"Wellington Basin Cross-Section at Y = {actual_y:.0f} m\nVelocity Structure and Receiver Placement",
        fontsize=14,
        fontweight="bold",
    )
    ax1.invert_yaxis()
    ax1.grid(True, alpha=0.3, linestyle="--")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax1, shrink=0.8)
    cbar.set_label("S-wave Velocity (m/s)", rotation=270, labelpad=15)

    # Add legend
    if receiver_coords is not None and len(receiver_coords) > 0:
        ax1.legend(loc="upper right", framealpha=0.9)

    # Bottom plot: Velocity vs depth profile at center of domain
    center_x_idx = len(grdx) // 2
    center_velocity_profile = vel_profile[:, center_x_idx]

    ax2.plot(
        center_velocity_profile,
        grdz,
        "b-",
        linewidth=2,
        label="Center velocity profile",
    )
    ax2.axvline(
        basin_threshold,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Basin threshold ({basin_threshold} m/s)",
    )

    # Highlight basin layers
    basin_mask = center_velocity_profile < basin_threshold
    if np.any(basin_mask):
        basin_depths = grdz[basin_mask]
        basin_vels = center_velocity_profile[basin_mask]
        ax2.fill_betweenx(
            basin_depths,
            0,
            basin_vels,
            alpha=0.3,
            color="brown",
            label="Basin sediments",
        )

    ax2.set_xlabel("Velocity (m/s)", fontsize=12)
    ax2.set_ylabel("Depth (m)", fontsize=12)
    ax2.set_title(f"Velocity Profile at X = {grdx[center_x_idx]:.0f} m", fontsize=12)
    ax2.invert_yaxis()
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="lower right")

    # Add geological interpretation text
    textstr = f"""Basin Analysis Summary:
    • Profile location: Y = {actual_y:.0f} m (Northing)
    • Basin threshold: {basin_threshold} m/s
    • Cross-section shows velocity structure from surface to {grdz.min():.0f} m depth
    • Red line marks sediment-bedrock boundary
    • Yellow triangles show receiver locations within ±200m"""

    props = dict(boxstyle="round", facecolor="lightblue", alpha=0.8)
    fig.text(
        0.02,
        0.02,
        textstr,
        fontsize=9,
        verticalalignment="bottom",
        bbox=props,
        transform=fig.transFigure,
    )

    plt.tight_layout()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Cross-section saved to: {save_path}")

    return fig, (ax1, ax2)


def plot_basin_summary_dashboard(
    grdx,
    grdy,
    grdz,
    grdvel,
    receiver_coords,
    basin_threshold=800.0,
    figsize=(16, 12),
    save_path=None,
):
    """
    Create a comprehensive dashboard showing multiple views of basin structure.

    Parameters:
    -----------
    grdx, grdy, grdz : ndarray
        Coordinate arrays from current cell
    grdvel : ndarray
        3D velocity array from current cell
    receiver_coords : ndarray
        Receiver coordinates (N, 3) with Z, X, Y format
    basin_threshold : float, default=800.0
        Velocity threshold for basin visualization
    figsize : tuple, default=(16, 12)
        Figure size
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    tuple
        (fig, axes) matplotlib objects
    """
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1, 0.8])

    # Top left: Surface geological map
    ax1 = fig.add_subplot(gs[0, 0])
    surface_vel = grdvel[0, :, :]
    extent = [grdx.min(), grdx.max(), grdy.min(), grdy.max()]

    colors = ["#8B4513", "#DAA520", "#228B22", "#4682B4", "#191970"]
    cmap = mcolors.LinearSegmentedColormap.from_list("geological", colors, N=100)

    im1 = ax1.imshow(
        surface_vel, origin="lower", extent=extent, cmap=cmap, aspect="equal"
    )

    # Plot receivers on map
    if len(receiver_coords) > 0:
        ax1.scatter(
            receiver_coords[:, 1],
            receiver_coords[:, 2],
            c=receiver_coords[:, 0],
            s=20,
            marker="^",
            cmap="plasma_r",
            edgecolors="white",
            linewidths=0.3,
        )

    ax1.set_title("Surface Velocity Map", fontweight="bold")
    ax1.set_xlabel("Easting (m)")
    ax1.set_ylabel("Northing (m)")
    plt.colorbar(im1, ax=ax1, shrink=0.8, label="Velocity (m/s)")

    # Top middle: Cross-section
    ax2 = fig.add_subplot(gs[0, 1])
    profile_y = (grdy.min() + grdy.max()) / 2
    j_idx = np.argmin(np.abs(grdy - profile_y))
    vel_profile = grdvel[:, j_idx, :]
    X, Z = np.meshgrid(grdx, grdz)

    im2 = ax2.pcolormesh(X, Z, vel_profile, cmap=cmap, shading="auto")
    ax2.contour(X, Z, vel_profile, levels=[basin_threshold], colors="red", linewidths=2)

    # Plot receivers in cross-section
    if len(receiver_coords) > 0:
        recv_mask = np.abs(receiver_coords[:, 2] - grdy[j_idx]) < 100
        if np.any(recv_mask):
            recv_near = receiver_coords[recv_mask]
            ax2.scatter(recv_near[:, 1], recv_near[:, 0], c="yellow", s=40, marker="v")

    ax2.set_title("Cross-Section View", fontweight="bold")
    ax2.set_xlabel("Easting (m)")
    ax2.set_ylabel("Depth (m)")
    ax2.invert_yaxis()
    plt.colorbar(im2, ax=ax2, shrink=0.8, label="Velocity (m/s)")

    # Top right: Velocity histogram
    ax3 = fig.add_subplot(gs[0, 2])
    velocities = grdvel.flatten()
    velocities = velocities[np.isfinite(velocities)]

    ax3.hist(velocities, bins=30, alpha=0.7, edgecolor="black", color="skyblue")
    ax3.axvline(basin_threshold, color="red", linestyle="--", linewidth=2)
    ax3.set_title("Velocity Distribution", fontweight="bold")
    ax3.set_xlabel("Velocity (m/s)")
    ax3.set_ylabel("Frequency")
    ax3.grid(True, alpha=0.3)

    # Bottom: Receiver statistics and information
    ax4 = fig.add_subplot(gs[1, :])
    ax4.axis("off")

    # Calculate statistics
    if len(receiver_coords) > 0:
        recv_stats = f"""
RECEIVER GENERATION SUMMARY
═══════════════════════════════════════════════════════════════════

Generated Receivers: {len(receiver_coords):,}
Coordinate Format: [Depth, Easting, Northing] (Z, X, Y)

Depth Statistics:
  • Min depth: {receiver_coords[:, 0].min():.1f} m
  • Max depth: {receiver_coords[:, 0].max():.1f} m
  • Mean depth: {receiver_coords[:, 0].mean():.1f} m

Spatial Coverage:
  • X range: {receiver_coords[:, 1].min():.0f} to {receiver_coords[:, 1].max():.0f} m
  • Y range: {receiver_coords[:, 2].min():.0f} to {receiver_coords[:, 2].max():.0f} m

Basin Parameters:
  • Basin threshold: {basin_threshold} m/s
  • Velocity range: {grdvel.min():.0f} - {grdvel.max():.0f} m/s
  • Basin sediment percentage: {(grdvel < basin_threshold).sum() / grdvel.size * 100:.1f}%

Model Domain:
  • Grid size: {grdvel.shape[2]} × {grdvel.shape[1]} × {grdvel.shape[0]} (X × Y × Z)
  • X extent: {grdx.min():.0f} to {grdx.max():.0f} m ({(grdx.max()-grdx.min()):.0f} m wide)
  • Y extent: {grdy.min():.0f} to {grdy.max():.0f} m ({(grdy.max()-grdy.min()):.0f} m wide)
  • Z extent: {grdz.max():.0f} to {grdz.min():.0f} m ({abs(grdz.min()-grdz.max()):.0f} m deep)
        """
    else:
        recv_stats = "No receivers generated - check basin threshold and velocity model"

    ax4.text(
        0.05,
        0.95,
        recv_stats,
        transform=ax4.transAxes,
        fontsize=10,
        verticalalignment="top",
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8),
    )

    plt.tight_layout()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Dashboard saved to: {save_path}")

    return fig, (ax1, ax2, ax3, ax4)
