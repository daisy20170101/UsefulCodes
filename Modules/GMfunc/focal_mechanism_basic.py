#!/usr/bin/env python3
"""
Basic focal mechanism plotting without ObsPy dependencies
Uses matplotlib and basic trigonometry to create beachball plots
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
import matplotlib.patches as patches


def strike_dip_rake_to_planes(strike, dip, rake):
    """
    Convert strike, dip, rake to nodal plane orientations
    This is a simplified calculation for visualization
    """
    # Convert to radians
    strike_rad = np.radians(strike)
    dip_rad = np.radians(dip)
    rake_rad = np.radians(rake)

    # Calculate auxiliary plane (simplified)
    # In reality, this involves more complex spherical trigonometry
    aux_strike = strike + 90
    if aux_strike >= 360:
        aux_strike -= 360

    aux_dip = 90 - dip
    aux_rake = 180 - rake

    return (strike, dip, rake), (aux_strike, aux_dip, aux_rake)


def plot_beachball(
    strike,
    dip,
    rake,
    ax=None,
    center=(0, 0),
    radius=1,
    compressive_color="black",
    extensive_color="white",
):
    """
    Plot a focal mechanism beachball with proper spherical projection
    Following Aki & Richards (1980) convention

    Parameters:
    -----------
    strike : float
        Strike angle in degrees (0-360), measured clockwise from North
    dip : float
        Dip angle in degrees (0-90)
    rake : float
        Rake angle in degrees (-180 to 180)
    ax : matplotlib axis
        Axis to plot on
    center : tuple
        Center coordinates (x, y)
    radius : float
        Radius of the beachball
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    # Convert to radians
    phi = np.radians(strike)  # strike
    delta = np.radians(dip)   # dip
    lam = np.radians(rake)     # rake

    # Calculate fault plane normal vector (pointing up from lower hemisphere)
    # Using Aki & Richards (1980) convention
    # Coordinate system: x=North, y=East, z=Down
    n = np.array([
        -np.sin(delta) * np.sin(phi),   # North component
        np.sin(delta) * np.cos(phi),    # East component
        -np.cos(delta)                   # Down component
    ])

    # Calculate slip vector on fault plane
    # Coordinate system: x=North, y=East, z=Down
    # Slip vector = rake direction along the fault plane
    # Using Aki & Richards formulation
    s = np.array([
        np.cos(lam) * np.cos(phi) + np.sin(lam) * np.cos(delta) * np.sin(phi),   # North
        np.cos(lam) * np.sin(phi) - np.sin(lam) * np.cos(delta) * np.cos(phi),   # East
        -np.sin(lam) * np.sin(delta)                                               # Down
    ])

    # Calculate P and T axes (principal axes)
    # P-axis (compression) and T-axis (tension)
    p_axis = (n + s) / np.linalg.norm(n + s)
    t_axis = (n - s) / np.linalg.norm(n - s)

    # Create a grid for the beachball
    npoints = 300
    x_grid = np.linspace(-1, 1, npoints)
    y_grid = np.linspace(-1, 1, npoints)
    X, Y = np.meshgrid(x_grid, y_grid)

    # Mask points outside the unit circle
    R = np.sqrt(X**2 + Y**2)
    mask = R <= 1.0

    # Convert from 2D projection back to 3D sphere (equal-area Schmidt projection)
    # For lower hemisphere projection
    r = np.sqrt(X**2 + Y**2)
    r_safe = np.where(mask, r, 1.0)

    # Equal-area (Schmidt) projection - most common for focal mechanisms
    # Convert (x,y) to point on lower hemisphere
    theta_point = np.arctan2(Y, X)  # azimuth on projection

    # For equal-area projection: r = sqrt(2) * sin(incidence/2)
    # where incidence is angle from downward vertical
    incidence = 2 * np.arcsin(r_safe / np.sqrt(2))
    incidence = np.where(incidence > np.pi/2, np.pi/2, incidence)

    # Convert to Cartesian coordinates on unit sphere (lower hemisphere)
    # x = North, y = East, z = Down
    # theta_point is azimuth measured from +X axis (right/East in plot) counter-clockwise
    # We need: x=North (up in plot), y=East (right in plot)
    # So: North = sin(incidence) * cos(theta_point), East = sin(incidence) * sin(theta_point)
    x_sphere = np.sin(incidence) * np.cos(theta_point)   # North component (up in plot)
    y_sphere = np.sin(incidence) * np.sin(theta_point)   # East component (right in plot)
    z_sphere = -np.cos(incidence)                         # Down component

    # Stack into array for easier computation
    points = np.stack([x_sphere, y_sphere, z_sphere], axis=-1)

    # Determine compression vs extension by testing which side of fault plane
    # Dot product with normal gives sign
    # For each point, we need to determine if it's in compression or extension quadrant
    # This is done by checking the sign of the dot product with both nodal planes

    dot_n = x_sphere * n[0] + y_sphere * n[1] + z_sphere * n[2]
    dot_s = x_sphere * s[0] + y_sphere * s[1] + z_sphere * s[2]

    # Compression occurs where both dot products have the same sign
    compression = ((dot_n * dot_s) > 0) & mask

    # Fill the background with white (extension)
    background = Circle(center, radius, facecolor=extensive_color,
                       edgecolor='none', zorder=1)
    ax.add_patch(background)

    # Plot compression quadrants in black
    ax.contourf(X*radius + center[0], Y*radius + center[1],
                compression.astype(float),
                levels=[0.5, 1.5], colors=[compressive_color],
                antialiased=True, zorder=2)

    # Draw nodal planes (where dot product is zero)
    ax.contour(X*radius + center[0], Y*radius + center[1],
               dot_n, levels=[0], colors='black', linewidths=0.5, zorder=3)
    ax.contour(X*radius + center[0], Y*radius + center[1],
               dot_s, levels=[0], colors='black', linewidths=0.5, zorder=3)

    # Draw outer circle
    circle = Circle(center, radius, fill=False, edgecolor="black",
                   linewidth=0.5, zorder=4)
    ax.add_patch(circle)

    # Set equal aspect ratio and limits
    ax.set_xlim(center[0] - radius * 1.2, center[0] + radius * 1.2)
    ax.set_ylim(center[1] - radius * 1.2, center[1] + radius * 1.2)
    ax.set_aspect("equal")

    return ax


def plot_multiple_mechanisms():
    """
    Create a figure showing different types of focal mechanisms
    """
    mechanisms = [
        (0, 90, 0, "Pure Strike-slip\n(0°, 90°, 0°)"),
        (0, 90, 180, "Pure Strike-slip\n(0°, 90°, 180°)"),
        (0, 60, -90, "Normal fault\n(0°, 60°, -90°)"),
        (0, 30, 90, "Thrust fault\n(0°, 30°, 90°)"),
        (45, 70, -45, "Oblique normal\n(45°, 70°, -45°)"),
        (120, 45, 135, "Oblique thrust\n(120°, 45°, 135°)"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Focal Mechanisms: Strike, Dip, Rake", fontsize=16, fontweight="bold")

    for i, (strike, dip, rake, title) in enumerate(mechanisms):
        row = i // 3
        col = i % 3
        ax = axes[row, col]

        plot_beachball(strike, dip, rake, ax, center=(0, 0), radius=1)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])

        # Add compass directions
        ax.text(0, 1.3, "N", ha="center", va="center", fontsize=10, fontweight="bold")
        ax.text(1.3, 0, "E", ha="center", va="center", fontsize=10, fontweight="bold")

    plt.tight_layout()
    plt.savefig("focal_mechanisms_basic.png", dpi=300, bbox_inches="tight")
    plt.show()

    return fig


def plot_mechanism_map():
    """
    Plot focal mechanisms on a simple map
    """
    # Sample earthquake data for New Zealand region
    earthquakes = [
        # (lon, lat, strike, dip, rake, magnitude, name)
        (174.0, -41.5, 45, 60, -90, 6.2, "Wellington"),
        (172.5, -43.5, 90, 45, 180, 5.8, "Canterbury"),
        (176.5, -38.5, 120, 30, 45, 5.5, "Bay of Plenty"),
        (171.0, -44.0, 30, 75, -120, 6.0, "West Coast"),
        (178.0, -37.0, 200, 85, 10, 5.2, "East Cape"),
    ]

    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot coastline (simplified New Zealand outline)
    # This is just for illustration - use real coastline data for production
    ax.plot(
        [170, 178, 178, 170, 170],
        [-46, -46, -34, -34, -46],
        "k-",
        linewidth=1,
        alpha=0.3,
    )

    for lon, lat, strike, dip, rake, mag, name in earthquakes:
        # Scale mechanism size by magnitude
        radius = mag * 0.15

        # Plot the focal mechanism
        plot_beachball(strike, dip, rake, ax, center=(lon, lat), radius=radius)

        # Add earthquake info
        ax.text(
            lon + 0.3,
            lat + 0.3,
            f"{name}\nM{mag}",
            fontsize=9,
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="gray"
            ),
        )

    ax.set_xlim(169, 179)
    ax.set_ylim(-46.5, -36)
    ax.set_xlabel("Longitude (°E)", fontsize=12)
    ax.set_ylabel("Latitude (°S)", fontsize=12)
    ax.set_title(
        "Focal Mechanisms - New Zealand Region", fontsize=14, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")

    # Add legend
    legend_elements = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="black",
            markersize=10,
            label="Compressive quadrants",
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="white",
            markeredgecolor="black",
            markersize=10,
            label="Extensive quadrants",
        ),
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    plt.tight_layout()
    plt.savefig("focal_mechanisms_map_basic.png", dpi=300, bbox_inches="tight")
    plt.show()

    return fig


def plot_cross_section(
    earthquakes,
    section_start,
    section_end,
    section_width=50,
    max_depth=200,
    title="Focal Mechanisms - Vertical Cross-Section",
):
    """
    Plot focal mechanisms on a vertical cross-section

    Parameters:
    -----------
    earthquakes : list
        List of tuples (lon, lat, depth, strike, dip, rake, magnitude, name)
    section_start : tuple
        (lon, lat) coordinates of section start point
    section_end : tuple
        (lon, lat) coordinates of section end point
    section_width : float
        Width of the cross-section in km (earthquakes within this distance)
    max_depth : float
        Maximum depth to display in km
    title : str
        Title for the plot
    """
    import math

    fig, (ax_map, ax_section) = plt.subplots(2, 1, figsize=(14, 12))

    # Calculate section parameters
    start_lon, start_lat = section_start
    end_lon, end_lat = section_end

    # Calculate section length and azimuth
    section_length = (
        math.sqrt((end_lon - start_lon) ** 2 + (end_lat - start_lat) ** 2) * 111
    )  # Approximate km
    section_azimuth = math.atan2(end_lon - start_lon, end_lat - start_lat)

    # === TOP PANEL: Map view ===
    ax_map.set_title(f"{title} - Map View", fontsize=14, fontweight="bold")

    # Plot section line
    ax_map.plot(
        [start_lon, end_lon],
        [start_lat, end_lat],
        "r-",
        linewidth=3,
        label="Cross-section line",
    )

    # Add section width bounds (simplified)
    width_offset = section_width / 111  # Convert km to degrees (approximate)
    perp_azimuth = section_azimuth + math.pi / 2

    # Plot section bounds
    bound1_start = (
        start_lon + width_offset * math.sin(perp_azimuth),
        start_lat + width_offset * math.cos(perp_azimuth),
    )
    bound1_end = (
        end_lon + width_offset * math.sin(perp_azimuth),
        end_lat + width_offset * math.cos(perp_azimuth),
    )
    bound2_start = (
        start_lon - width_offset * math.sin(perp_azimuth),
        start_lat - width_offset * math.cos(perp_azimuth),
    )
    bound2_end = (
        end_lon - width_offset * math.sin(perp_azimuth),
        end_lat - width_offset * math.cos(perp_azimuth),
    )

    ax_map.plot(
        [bound1_start[0], bound1_end[0]],
        [bound1_start[1], bound1_end[1]],
        "r--",
        alpha=0.5,
    )
    ax_map.plot(
        [bound2_start[0], bound2_end[0]],
        [bound2_start[1], bound2_end[1]],
        "r--",
        alpha=0.5,
    )

    # Plot earthquakes on map
    selected_earthquakes = []

    for eq in earthquakes:
        lon, lat, depth, strike, dip, rake, mag, name = eq

        # Calculate distance from earthquake to section line
        # Simplified distance calculation
        dist_to_line = abs(
            (end_lon - start_lon) * (start_lat - lat)
            - (start_lon - lon) * (end_lat - start_lat)
        ) / math.sqrt((end_lon - start_lon) ** 2 + (end_lat - start_lat) ** 2)
        dist_to_line *= 111  # Convert to km

        # Color code by depth
        if depth <= 50:
            color = "red"
        elif depth <= 100:
            color = "orange"
        else:
            color = "blue"

        # Plot earthquake location
        ax_map.scatter(lon, lat, c=color, s=mag * 20, alpha=0.7, edgecolors="black")

        # Select earthquakes within section width
        if dist_to_line <= section_width:
            # Calculate distance along section
            along_section = (
                (lon - start_lon) * (end_lon - start_lon)
                + (lat - start_lat) * (end_lat - start_lat)
            ) / ((end_lon - start_lon) ** 2 + (end_lat - start_lat) ** 2)
            distance_along = along_section * section_length

            selected_earthquakes.append(
                (distance_along, depth, strike, dip, rake, mag, name, color)
            )

    # Set map extent
    all_lons = [eq[0] for eq in earthquakes] + [start_lon, end_lon]
    all_lats = [eq[1] for eq in earthquakes] + [start_lat, end_lat]
    margin = 0.5
    ax_map.set_xlim(min(all_lons) - margin, max(all_lons) + margin)
    ax_map.set_ylim(min(all_lats) - margin, max(all_lats) + margin)
    ax_map.set_xlabel("Longitude (°E)")
    ax_map.set_ylabel("Latitude (°S)")
    ax_map.grid(True, alpha=0.3)
    ax_map.legend()

    # Add depth legend for map
    legend_elements = [
        plt.scatter([], [], c="red", s=50, label="0-50 km depth"),
        plt.scatter([], [], c="orange", s=50, label="50-100 km depth"),
        plt.scatter([], [], c="blue", s=50, label=">100 km depth"),
    ]
    ax_map.legend(handles=legend_elements, loc="upper right")

    # === BOTTOM PANEL: Cross-section ===
    ax_section.set_title(
        f"{title} - Vertical Cross-Section", fontsize=14, fontweight="bold"
    )

    # Plot focal mechanisms on cross-section
    for dist_along, depth, strike, dip, rake, mag, name, color in selected_earthquakes:
        # Scale mechanism size by magnitude
        radius = mag * 2  # Adjust scaling as needed

        # Plot focal mechanism at (distance, depth) coordinates
        plot_beachball(
            strike, dip, rake, ax_section, center=(dist_along, depth), radius=radius
        )

        # Add earthquake label
        if mag >= 5.0:  # Only label larger earthquakes to avoid clutter
            ax_section.text(
                dist_along + radius + 2,
                depth - radius,
                f"{name}\nM{mag}",
                fontsize=8,
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    alpha=0.8,
                    edgecolor="gray",
                ),
            )

    # Set cross-section axes
    ax_section.set_xlim(0, section_length)
    ax_section.set_ylim(max_depth, 0)  # Depth increases downward
    ax_section.set_xlabel("Distance along section (km)")
    ax_section.set_ylabel("Depth (km)")
    ax_section.grid(True, alpha=0.3)

    # Add section endpoints labels
    ax_section.text(
        0,
        -10,
        f"A ({start_lat:.1f}°S, {start_lon:.1f}°E)",
        ha="left",
        va="top",
        fontweight="bold",
    )
    ax_section.text(
        section_length,
        -10,
        f"B ({end_lat:.1f}°S, {end_lon:.1f}°E)",
        ha="right",
        va="top",
        fontweight="bold",
    )

    plt.tight_layout()
    plt.savefig("focal_mechanisms_cross_section.png", dpi=300, bbox_inches="tight")
    plt.show()

    return fig


def example_subduction_zone_cross_section():
    """
    Example showing focal mechanisms in a subduction zone cross-section
    """
    # Sample subduction zone earthquake data (Hikurangi subduction zone, NZ)
    # (lon, lat, depth, strike, dip, rake, magnitude, name)
    subduction_earthquakes = [
        # Shallow interface events
        (178.0, -37.5, 15, 200, 15, 90, 6.8, "Interface1"),
        (177.8, -37.8, 18, 205, 18, 85, 6.2, "Interface2"),
        (177.5, -38.2, 22, 210, 20, 95, 5.9, "Interface3"),
        # Intermediate depth events
        (177.3, -38.5, 35, 215, 25, 90, 5.5, "Interface4"),
        (177.0, -38.8, 45, 220, 30, 85, 5.8, "Interface5"),
        # Intraslab events (normal faulting within slab)
        (176.8, -39.1, 80, 30, 60, -90, 6.1, "Intraslab1"),
        (176.5, -39.4, 95, 45, 70, -85, 5.7, "Intraslab2"),
        (176.2, -39.7, 120, 60, 65, -95, 5.4, "Intraslab3"),
        # Deep intraslab events
        (175.8, -40.2, 150, 75, 75, -80, 5.2, "Deep1"),
        (175.5, -40.5, 180, 90, 80, -90, 5.0, "Deep2"),
        # Upper plate events (crustal)
        (178.5, -37.0, 8, 45, 85, 0, 5.3, "Crustal1"),
        (178.3, -37.3, 12, 60, 80, 10, 4.8, "Crustal2"),
    ]

    # Define cross-section from offshore to onshore (perpendicular to trench)
    section_start = (179.0, -36.5)  # Offshore start
    section_end = (175.0, -41.0)  # Onshore end

    print("\n=== Subduction Zone Cross-Section Example ===")

    fig = plot_cross_section(
        subduction_earthquakes,
        section_start,
        section_end,
        section_width=75,  # 75 km width
        max_depth=200,
        title="Hikurangi Subduction Zone - Focal Mechanisms",
    )

    print("✓ Subduction zone cross-section created")
    print("This shows:")
    print("  • Interface thrust events (shallow)")
    print("  • Intraslab normal events (intermediate/deep)")
    print("  • Upper plate crustal events")

    return fig


def main():
    """
    Run focal mechanism plotting examples
    """
    print("Basic Focal Mechanism Plotting")
    print("=" * 50)
    print("Creating examples of focal mechanism plots...")

    # Plot different mechanism types
    plot_multiple_mechanisms()
    print("✓ Multiple mechanism types plotted")

    # Plot mechanisms on map
    plot_mechanism_map()
    print("✓ Map with focal mechanisms created")

    # Plot cross-section example
    example_subduction_zone_cross_section()

    print("\n" + "=" * 50)
    print("Understanding Focal Mechanisms:")
    print("• Black areas: Compressive quadrants")
    print("• White areas: Extensive quadrants")
    print("• Lines: Nodal planes (possible fault orientations)")
    print("\nParameters:")
    print("• Strike: Azimuth of fault line (0-360°)")
    print("• Dip: Inclination from horizontal (0-90°)")
    print("• Rake: Slip direction on fault (-180° to +180°)")
    print("  - Normal faulting: -90° ± 30°")
    print("  - Strike-slip: 0° or ±180° ± 30°")
    print("  - Thrust faulting: +90° ± 30°")

    print("\nNote: These are simplified beachballs.")
    print("For accurate focal mechanisms, use ObsPy when available:")
    print("conda install obspy")


if __name__ == "__main__":
    main()
