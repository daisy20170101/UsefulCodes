"""
Plot 3D seismological inversion setup with fault geometry and seismic stations
Shows fault polygons from VTK files and synthetic station distribution on free surface
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pyvista as pv

def read_vtk_fault(vtk_file):
    """
    Read VTK file and extract fault geometry

    Parameters:
    -----------
    vtk_file : str
        Path to VTK file

    Returns:
    --------
    vertices : ndarray
        Vertex coordinates (N, 3)
    faces : ndarray
        Face connectivity
    """
    mesh = pv.read(vtk_file)

    # Get vertices
    vertices = mesh.points

    # Get faces
    faces = mesh.faces.reshape(-1, 4)[:, 1:]  # Skip the count column

    return vertices, faces


def generate_station_grid(lon_range, lat_range, n_stations=20, z=0):
    """
    Generate synthetic seismic station locations on free surface

    Parameters:
    -----------
    lon_range : tuple
        (min_lon, max_lon)
    lat_range : tuple
        (min_lat, max_lat)
    n_stations : int
        Number of stations (will be squared for grid)
    z : float
        Elevation of free surface (default: 0)

    Returns:
    --------
    stations : ndarray
        Station coordinates (N, 3)
    """
    n_side = int(np.sqrt(n_stations))
    lons = np.linspace(lon_range[0], lon_range[1], n_side)
    lats = np.linspace(lat_range[0], lat_range[1], n_side)

    lon_grid, lat_grid = np.meshgrid(lons, lats)

    stations = np.column_stack([
        lon_grid.ravel(),
        lat_grid.ravel(),
        np.full(n_side**2, z)
    ])

    return stations


def plot_seismo_inversion_3d(vtk_files,
                              output_file='seismo_inversion_3d.png',
                              station_params=None,
                              fault_colors=None,
                              fault_labels=None,
                              fault_alpha=0.7,
                              station_color='red',
                              station_size=50,
                              station_marker='^',
                              view_elev=20,
                              view_azim=45,
                              figsize=(14, 10),
                              add_axes_labels=True,
                              show_grid=True):
    """
    Create 3D visualization of seismological inversion setup

    Parameters:
    -----------
    vtk_files : list of str
        List of VTK file paths for fault geometries
    output_file : str
        Output filename for the figure
    station_params : dict, optional
        Parameters for station generation:
        - 'lon_range': (min, max) longitude
        - 'lat_range': (min, max) latitude
        - 'n_stations': number of stations
        - 'z': elevation of free surface
        If None, auto-generates based on fault extent
    fault_colors : list, optional
        List of colors for each fault (default: color cycle)
    fault_labels : list, optional
        List of labels for each fault
    fault_alpha : float
        Transparency of fault polygons (0-1)
    station_color : str
        Color for station markers
    station_size : int
        Size of station markers
    station_marker : str
        Marker style for stations ('^', 'v', 'o', 's', etc.)
    view_elev : float
        Elevation angle for 3D view (degrees)
    view_azim : float
        Azimuth angle for 3D view (degrees)
    figsize : tuple
        Figure size (width, height)
    add_axes_labels : bool
        Whether to add axis labels
    show_grid : bool
        Whether to show grid

    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """

    print("Creating 3D seismological inversion visualization...")

    # Create figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    # Default colors if not provided
    if fault_colors is None:
        fault_colors = plt.cm.tab10(np.linspace(0, 1, len(vtk_files)))

    # Default labels if not provided
    if fault_labels is None:
        fault_labels = [f"Fault {i+1}" for i in range(len(vtk_files))]

    # Track extent for auto-generating station grid
    all_vertices = []

    # ========================================================================
    # Plot fault geometries
    # ========================================================================
    print(f"  Loading and plotting {len(vtk_files)} fault geometries...")

    for i, vtk_file in enumerate(vtk_files):
        print(f"    Processing: {vtk_file}")

        try:
            vertices, faces = read_vtk_fault(vtk_file)
            all_vertices.append(vertices)

            # Create polygon collection for 3D plotting
            polygons = []
            for face in faces:
                polygon = vertices[face]
                polygons.append(polygon)

            # Create 3D polygon collection
            poly_collection = Poly3DCollection(
                polygons,
                alpha=fault_alpha,
                facecolor=fault_colors[i],
                edgecolor='black',
                linewidth=0.5,
                label=fault_labels[i]
            )

            ax.add_collection3d(poly_collection)

            print(f"      ✓ Loaded {len(vertices)} vertices, {len(faces)} faces")

        except Exception as e:
            print(f"      ✗ Error loading {vtk_file}: {e}")

    # ========================================================================
    # Generate and plot seismic stations
    # ========================================================================
    if all_vertices:
        all_verts = np.vstack(all_vertices)

        # Auto-generate station parameters if not provided
        if station_params is None:
            lon_min, lon_max = all_verts[:, 0].min(), all_verts[:, 0].max()
            lat_min, lat_max = all_verts[:, 1].min(), all_verts[:, 1].max()

            # Expand range slightly for better coverage
            lon_margin = (lon_max - lon_min) * 0.2
            lat_margin = (lat_max - lat_min) * 0.2

            station_params = {
                'lon_range': (lon_min - lon_margin, lon_max + lon_margin),
                'lat_range': (lat_min - lat_margin, lat_max + lat_margin),
                'n_stations': 25,
                'z': 0  # Free surface
            }

        print(f"  Generating synthetic seismic stations...")
        print(f"    Lon range: {station_params['lon_range']}")
        print(f"    Lat range: {station_params['lat_range']}")
        print(f"    Number of stations: {station_params['n_stations']}")

        stations = generate_station_grid(
            station_params['lon_range'],
            station_params['lat_range'],
            station_params['n_stations'],
            station_params.get('z', 0)
        )

        # Plot stations on free surface
        ax.scatter(
            stations[:, 0],
            stations[:, 1],
            stations[:, 2],
            c=station_color,
            marker=station_marker,
            s=station_size,
            edgecolor='black',
            linewidth=1,
            label='Seismic Stations',
            zorder=100
        )

        print(f"    ✓ Plotted {len(stations)} stations")

        # Plot free surface as transparent plane
        lon_surf = np.linspace(station_params['lon_range'][0],
                               station_params['lon_range'][1], 20)
        lat_surf = np.linspace(station_params['lat_range'][0],
                               station_params['lat_range'][1], 20)
        lon_surf_grid, lat_surf_grid = np.meshgrid(lon_surf, lat_surf)
        z_surf = np.zeros_like(lon_surf_grid)

        ax.plot_surface(
            lon_surf_grid,
            lat_surf_grid,
            z_surf,
            alpha=0.1,
            color='gray',
            rstride=2,
            cstride=2,
            linewidth=0.5,
            edgecolor='gray'
        )

        # Set axis limits based on data
        ax.set_xlim(station_params['lon_range'])
        ax.set_ylim(station_params['lat_range'])

        z_min = all_verts[:, 2].min()
        z_max = max(all_verts[:, 2].max(), station_params.get('z', 0))
        z_margin = (z_max - z_min) * 0.1
        ax.set_zlim(z_min - z_margin, z_max + z_margin)

    # ========================================================================
    # Configure plot appearance
    # ========================================================================
    if add_axes_labels:
        ax.set_xlabel('Longitude (°)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Latitude (°)', fontsize=12, fontweight='bold')
        ax.set_zlabel('Depth (km)', fontsize=12, fontweight='bold')

    # Set view angle
    ax.view_init(elev=view_elev, azim=view_azim)

    # Grid
    if show_grid:
        ax.grid(True, alpha=0.3)

    # Legend
    ax.legend(loc='upper left', fontsize=10, framealpha=0.9)

    # Title
    ax.set_title('Seismological Inversion Setup\nFault Geometry and Station Distribution',
                 fontsize=14, fontweight='bold', pad=20)

    # ========================================================================
    # Save figure
    # ========================================================================
    print(f"\nSaving figure to: {output_file}")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print("Done!")

    return fig, ax


def plot_seismo_inversion_multiple_views(vtk_files,
                                          output_file='seismo_inversion_multi.png',
                                          station_params=None,
                                          fault_colors=None,
                                          fault_labels=None):
    """
    Create multiple view angles of the seismological inversion setup

    Parameters:
    -----------
    vtk_files : list of str
        List of VTK file paths
    output_file : str
        Output filename
    station_params : dict, optional
        Station generation parameters
    fault_colors : list, optional
        Colors for each fault
    fault_labels : list, optional
        Labels for each fault

    Returns:
    --------
    fig : matplotlib figure object
    """

    print("Creating multi-view seismological inversion visualization...")

    # Create figure with subplots
    fig = plt.figure(figsize=(18, 12))

    # Define view angles (elev, azim)
    views = [
        (20, 45, 'View 1: Southeast'),
        (20, 135, 'View 2: Southwest'),
        (60, 45, 'View 3: Top-Southeast'),
        (10, 0, 'View 4: South')
    ]

    for idx, (elev, azim, title) in enumerate(views):
        ax = fig.add_subplot(2, 2, idx+1, projection='3d')

        # Load and plot faults
        if fault_colors is None:
            colors = plt.cm.tab10(np.linspace(0, 1, len(vtk_files)))
        else:
            colors = fault_colors

        if fault_labels is None:
            labels = [f"Fault {i+1}" for i in range(len(vtk_files))]
        else:
            labels = fault_labels

        all_vertices = []

        for i, vtk_file in enumerate(vtk_files):
            try:
                vertices, faces = read_vtk_fault(vtk_file)
                all_vertices.append(vertices)

                polygons = []
                for face in faces:
                    polygon = vertices[face]
                    polygons.append(polygon)

                poly_collection = Poly3DCollection(
                    polygons,
                    alpha=0.7,
                    facecolor=colors[i],
                    edgecolor='black',
                    linewidth=0.5,
                    label=labels[i] if idx == 0 else None
                )
                ax.add_collection3d(poly_collection)

            except Exception as e:
                print(f"    Warning: {e}")

        # Generate and plot stations
        if all_vertices:
            all_verts = np.vstack(all_vertices)

            if station_params is None:
                lon_min, lon_max = all_verts[:, 0].min(), all_verts[:, 0].max()
                lat_min, lat_max = all_verts[:, 1].min(), all_verts[:, 1].max()
                lon_margin = (lon_max - lon_min) * 0.2
                lat_margin = (lat_max - lat_min) * 0.2

                params = {
                    'lon_range': (lon_min - lon_margin, lon_max + lon_margin),
                    'lat_range': (lat_min - lat_margin, lat_max + lat_margin),
                    'n_stations': 25,
                    'z': 0
                }
            else:
                params = station_params

            stations = generate_station_grid(
                params['lon_range'],
                params['lat_range'],
                params['n_stations'],
                params.get('z', 0)
            )

            ax.scatter(
                stations[:, 0],
                stations[:, 1],
                stations[:, 2],
                c='red',
                marker='^',
                s=30,
                edgecolor='black',
                linewidth=0.5,
                label='Stations' if idx == 0 else None,
                zorder=100
            )

            # Set limits
            ax.set_xlim(params['lon_range'])
            ax.set_ylim(params['lat_range'])
            z_min = all_verts[:, 2].min()
            z_max = max(all_verts[:, 2].max(), 0)
            ax.set_zlim(z_min * 1.1, z_max * 1.1)

        # Configure appearance
        ax.set_xlabel('Lon (°)', fontsize=10)
        ax.set_ylabel('Lat (°)', fontsize=10)
        ax.set_zlabel('Depth (km)', fontsize=10)
        ax.view_init(elev=elev, azim=azim)
        ax.grid(True, alpha=0.3)
        ax.set_title(title, fontsize=12, fontweight='bold')

        if idx == 0:
            ax.legend(loc='upper left', fontsize=8, framealpha=0.9)

    plt.suptitle('Seismological Inversion Setup - Multiple Views',
                 fontsize=16, fontweight='bold', y=0.98)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved to: {output_file}")

    return fig


# ============================================================================
# Example usage
# ============================================================================

if __name__ == "__main__":

    # Example 1: Single view with specific VTK files
    print("="*70)
    print("Example 1: Single 3D view")
    print("="*70 + "\n")

    vtk_files = [
        '/Users/DuoL/Documents/NSHM/Central/Paraviews/slab_hik_edge.vtk',
        '/Users/DuoL/Documents/NSHM/Central/Paraviews/welfault_edge.vtk',
        '/Users/DuoL/Documents/NSHM/Central/Paraviews/wairarapa_edge.vtk'
    ]

    fault_labels = [
        'Hikurangi Slab',
        'Wellington Fault',
        'Wairarapa Fault'
    ]

    fault_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, orange, green

    fig, ax = plot_seismo_inversion_3d(
        vtk_files=vtk_files,
        output_file='seismo_inversion_3d.png',
        fault_colors=fault_colors,
        fault_labels=fault_labels,
        fault_alpha=0.7,
        station_color='red',
        station_size=60,
        station_marker='^',
        view_elev=25,
        view_azim=50
    )

    plt.show()

    # Example 2: Multiple views
    print("\n" + "="*70)
    print("Example 2: Multiple views")
    print("="*70 + "\n")

    fig_multi = plot_seismo_inversion_multiple_views(
        vtk_files=vtk_files,
        output_file='seismo_inversion_multi.png',
        fault_colors=fault_colors,
        fault_labels=fault_labels
    )

    plt.show()

    print("\n" + "="*70)
    print("Complete!")
    print("="*70)


"""
Additional usage examples:

# Custom station grid
station_params = {
    'lon_range': (174.5, 175.5),
    'lat_range': (-41.5, -40.5),
    'n_stations': 36,
    'z': 0
}

fig, ax = plot_seismo_inversion_3d(
    vtk_files=['fault1.vtk', 'fault2.vtk'],
    station_params=station_params,
    output_file='custom_setup.png'
)

# Single fault with many stations
fig, ax = plot_seismo_inversion_3d(
    vtk_files=['slab_hik_edge.vtk'],
    fault_labels=['Hikurangi Subduction Interface'],
    fault_colors=['steelblue'],
    station_params={'n_stations': 49},  # 7x7 grid
    station_marker='v',  # Inverted triangle
    view_elev=30,
    view_azim=60
)
"""
