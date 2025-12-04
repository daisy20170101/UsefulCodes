"""
Plot New Zealand geological map with plate boundaries using Matplotlib and Cartopy
Supports separate colormaps for land topography and ocean bathymetry
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch
from matplotlib import patheffects

def plot_nz_geological_map(output_file='nz_geological_map.png',
                           show_earthquakes=False,
                           eq_startdate='2014-01-01',
                           eq_minmag=5.0,
                           figsize=(12, 14),
                           topography_file=None,
                           land_cmap='terrain',
                           ocean_cmap='Blues_r',
                           land_vmin=None,
                           land_vmax=None,
                           ocean_vmin=None,
                           ocean_vmax=None,
                           add_hillshade=True,
                           lon_range=None,
                           lat_range=None):
    """
    Create a geological map of New Zealand with plate boundaries
    Supports separate colormaps for land and ocean

    Parameters:
    -----------
    output_file : str
        Output filename (png, pdf, etc.)
    show_earthquakes : bool
        Whether to plot recent earthquakes
    eq_startdate : str
        Start date for earthquake data
    eq_minmag : float
        Minimum magnitude for earthquakes
    figsize : tuple
        Figure size (width, height) in inches
    topography_file : str, optional
        Path to topography data file (.grd, .nc, .tif, etc.)
    land_cmap : str
        Colormap for land topography (elevation > 0)
        Options: 'terrain', 'gist_earth', 'YlOrBr', 'copper', 'Greens', etc.
    ocean_cmap : str
        Colormap for ocean bathymetry (elevation < 0)
        Options: 'Blues_r', 'viridis', 'ocean', 'Blues', etc.
    land_vmin : float, optional
        Minimum elevation for land colormap (default: 0)
    land_vmax : float, optional
        Maximum elevation for land colormap (default: auto)
    ocean_vmin : float, optional
        Minimum depth for ocean colormap (default: auto)
    ocean_vmax : float, optional
        Maximum depth for ocean colormap (default: 0)
    add_hillshade : bool
        Whether to add hillshade effect (default: True)
    lon_range : tuple, optional
        Longitude range (min_lon, max_lon). If None, uses default NZ extent [164, 180]
    lat_range : tuple, optional
        Latitude range (min_lat, max_lat). If None, uses default NZ extent [-48, -34]
    """

    print("Creating New Zealand geological map with Matplotlib/Cartopy...")

    # Define New Zealand region
    if lon_range is None:
        lon_range = [164, 180]
    if lat_range is None:
        lat_range = [-48, -34]

    extent = [lon_range[0], lon_range[1], lat_range[0], lat_range[1]]  # [west, east, south, north]

    print(f"  Map extent: Lon [{extent[0]}, {extent[1]}], Lat [{extent[2]}, {extent[3]}]")

    # Create figure and axis
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    # ========================================================================
    # Add topography with separate land/ocean colormaps
    # ========================================================================
    if topography_file:
        print(f"  Loading topography from: {topography_file}")

        try:
            import xarray as xr
            from matplotlib.colors import LightSource

            # Load topography data
            topo_data = xr.open_dataset(topography_file)

            # Get variable name
            var_names = list(topo_data.data_vars)
            if len(var_names) == 1:
                topo_var = var_names[0]
            else:
                for name in ['z', 'elevation', 'Band1', 'topo', 'altitude']:
                    if name in var_names:
                        topo_var = name
                        break
                else:
                    topo_var = var_names[0]

            print(f"    Using variable: {topo_var}")

            # Get coordinates
            topo = topo_data[topo_var]
            coord_names = list(topo.coords)
            lon_names = [c for c in coord_names if c.lower() in ['lon', 'longitude', 'x']]
            lat_names = [c for c in coord_names if c.lower() in ['lat', 'latitude', 'y']]

            if lon_names and lat_names:
                lon_coord = lon_names[0]
                lat_coord = lat_names[0]

                lons = topo[lon_coord].values
                lats = topo[lat_coord].values
                elevation = topo.values

                # Squeeze if needed
                if elevation.ndim > 2:
                    elevation = np.squeeze(elevation)

                print(f"    Elevation range: {np.nanmin(elevation):.0f} to {np.nanmax(elevation):.0f} m")

                # Create meshgrid if needed
                if lons.ndim == 1 and lats.ndim == 1:
                    lons, lats = np.meshgrid(lons, lats)

                # Separate land and ocean data
                land_mask = elevation >= 0
                ocean_mask = elevation < 0

                # Set default ranges if not provided
                if land_vmin is None:
                    land_vmin = 0
                if land_vmax is None:
                    land_vmax = np.nanpercentile(elevation[land_mask], 98) if np.any(land_mask) else 3000

                if ocean_vmin is None:
                    ocean_vmin = np.nanpercentile(elevation[ocean_mask], 2) if np.any(ocean_mask) else -6000
                if ocean_vmax is None:
                    ocean_vmax = 0

                print(f"    Land colormap range: {land_vmin:.0f} to {land_vmax:.0f} m")
                print(f"    Ocean colormap range: {ocean_vmin:.0f} to {ocean_vmax:.0f} m")

                # Create hillshade if requested
                if add_hillshade:
                    ls = LightSource(azdeg=315, altdeg=45)

                    # Apply hillshade to land
                    if np.any(land_mask):
                        land_data = np.where(land_mask, elevation, np.nan)
                        land_shaded = ls.shade(land_data, cmap=plt.get_cmap(land_cmap),
                                             vmin=land_vmin, vmax=land_vmax,
                                             vert_exag=1, blend_mode='soft')

                        ax.pcolormesh(lons, lats, land_data,
                                     cmap=land_cmap,
                                     vmin=land_vmin, vmax=land_vmax,
                                     transform=ccrs.PlateCarree(),
                                     shading='auto', zorder=1, alpha=0.9)

                    # Apply hillshade to ocean
                    if np.any(ocean_mask):
                        ocean_data = np.where(ocean_mask, elevation, np.nan)
                        ocean_shaded = ls.shade(ocean_data, cmap=plt.get_cmap(ocean_cmap),
                                              vmin=ocean_vmin, vmax=ocean_vmax,
                                              vert_exag=0.5, blend_mode='soft')

                        ax.pcolormesh(lons, lats, ocean_data,
                                     cmap=ocean_cmap,
                                     vmin=ocean_vmin, vmax=ocean_vmax,
                                     transform=ccrs.PlateCarree(),
                                     shading='auto', zorder=1, alpha=0.9)
                else:
                    # Plot without hillshade
                    if np.any(land_mask):
                        land_data = np.where(land_mask, elevation, np.nan)
                        im_land = ax.pcolormesh(lons, lats, land_data,
                                               cmap=land_cmap,
                                               vmin=land_vmin, vmax=land_vmax,
                                               transform=ccrs.PlateCarree(),
                                               shading='auto', zorder=1, alpha=0.8)

                    if np.any(ocean_mask):
                        ocean_data = np.where(ocean_mask, elevation, np.nan)
                        im_ocean = ax.pcolormesh(lons, lats, ocean_data,
                                                cmap=ocean_cmap,
                                                vmin=ocean_vmin, vmax=ocean_vmax,
                                                transform=ccrs.PlateCarree(),
                                                shading='auto', zorder=1, alpha=0.8)

                # Add colorbars for both land and ocean
                if np.any(land_mask):
                    # Create land colorbar
                    sm_land = plt.cm.ScalarMappable(
                        cmap=land_cmap,
                        norm=mcolors.Normalize(vmin=land_vmin, vmax=land_vmax)
                    )
                    sm_land.set_array([])
                    cbar_land = plt.colorbar(sm_land, ax=ax, orientation='vertical',
                                            pad=0.08, shrink=0.4, aspect=15,
                                            anchor=(0, 0.6))
                    cbar_land.set_label('Land Elevation (m)', fontsize=10, fontweight='bold')

                if np.any(ocean_mask):
                    # Create ocean colorbar
                    sm_ocean = plt.cm.ScalarMappable(
                        cmap=ocean_cmap,
                        norm=mcolors.Normalize(vmin=ocean_vmin, vmax=ocean_vmax)
                    )
                    sm_ocean.set_array([])
                    cbar_ocean = plt.colorbar(sm_ocean, ax=ax, orientation='vertical',
                                             pad=0.08, shrink=0.4, aspect=15,
                                             anchor=(0, 0.0))
                    cbar_ocean.set_label('Ocean Depth (m)', fontsize=10, fontweight='bold')

                print("    ✓ Topography plotted successfully with dual colormaps")

            else:
                print("    Warning: Could not identify lon/lat coordinates")
                topography_file = None

        except ImportError:
            print("    Warning: xarray not installed. Install with: pip install xarray")
            topography_file = None
        except Exception as e:
            print(f"    Warning: Could not load topography: {e}")
            topography_file = None

    # ========================================================================
    # Add features (if topography not loaded)
    # ========================================================================
    if not topography_file:
        print("  Adding base map features...")
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.5)
        ax.add_feature(cfeature.LAND, facecolor='wheat', edgecolor='black', linewidth=0.5)

    # Always add coastlines and borders
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8, zorder=4)
    ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5, alpha=0.5, zorder=4)

    # Add gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    # ========================================================================
    # Plot plate boundaries and faults (optional - currently disabled)
    # ========================================================================
    # To enable specific faults, uncomment the relevant sections below

    # # Alpine Fault
    # alpine_fault = np.array([
    #     [166.5, -44.0], [167.0, -43.5], [167.5, -43.0], [168.0, -42.5],
    #     [168.5, -42.0], [169.0, -41.5], [169.5, -41.0], [170.0, -40.5], [170.5, -40.0]
    # ])
    # ax.plot(alpine_fault[:, 0], alpine_fault[:, 1], color='darkred', linewidth=3,
    #         linestyle='-', transform=ccrs.PlateCarree(), label='Alpine Fault', zorder=5)

    # # Hikurangi Subduction Zone
    # hikurangi = np.array([
    #     [178.5, -37.5], [178.0, -38.0], [177.5, -38.5], [177.0, -39.0],
    #     [176.5, -39.5], [176.0, -40.0], [175.5, -40.5], [175.0, -41.0],
    #     [174.5, -41.5], [174.0, -42.0]
    # ])
    # ax.plot(hikurangi[:, 0], hikurangi[:, 1], color='darkorange', linewidth=3,
    #         linestyle='-', transform=ccrs.PlateCarree(), label='Hikurangi Subduction Zone', zorder=5)

    # # Puysegur Trench
    # puysegur = np.array([[166.0, -45.5], [165.5, -46.0], [165.0, -46.5]])
    # ax.plot(puysegur[:, 0], puysegur[:, 1], color='orange', linewidth=3,
    #         linestyle='-', transform=ccrs.PlateCarree(), label='Puysegur Trench', zorder=5)

    # # Hope Fault
    # hope_fault = np.array([[172.0, -42.5], [172.5, -42.4], [173.0, -42.3], [173.5, -42.2]])
    # ax.plot(hope_fault[:, 0], hope_fault[:, 1], color='red', linewidth=2,
    #         linestyle='--', transform=ccrs.PlateCarree(), label='Major Faults', zorder=5)

    # # Wellington Fault
    # wellington_fault = np.array([[174.8, -41.3], [174.85, -41.2], [174.9, -41.1]])
    # ax.plot(wellington_fault[:, 0], wellington_fault[:, 1], color='red', linewidth=2,
    #         linestyle='--', transform=ccrs.PlateCarree(), zorder=5)

    # ========================================================================
    # Add major cities
    # ========================================================================
    print("  Adding major cities...")

    cities = pd.DataFrame({
        'name': ['Auckland', 'Wellington', 'Christchurch', 'Dunedin', 'Hamilton'],
        'lon': [174.76, 174.78, 172.64, 170.50, 175.28],
        'lat': [-36.85, -41.29, -43.53, -45.87, -37.78]
    })

    ax.scatter(cities['lon'], cities['lat'], s=100, c='white', edgecolor='black',
               linewidth=1.5, transform=ccrs.PlateCarree(), zorder=10)

    for idx, city in cities.iterrows():
        txt = ax.text(city['lon'] + 0.3, city['lat'] + 0.2, city['name'],
                      fontsize=10, fontweight='bold', transform=ccrs.PlateCarree(),
                      zorder=11, bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                          alpha=0.7, edgecolor='black'))

    # ========================================================================
    # Optionally plot earthquakes
    # ========================================================================
    if show_earthquakes:
        print(f"  Downloading earthquakes M≥{eq_minmag} since {eq_startdate}...")
        try:
            import requests
            from datetime import datetime

            url = "https://earthquake.usgs.gov/fdsnws/event/1/query"
            params = {
                'format': 'geojson',
                'starttime': eq_startdate,
                'endtime': datetime.now().strftime('%Y-%m-%d'),
                'minmagnitude': eq_minmag,
                'minlatitude': extent[2], 'maxlatitude': extent[3],
                'minlongitude': extent[0], 'maxlongitude': extent[1]
            }

            response = requests.get(url, params=params, timeout=60)
            if response.status_code == 200:
                data = response.json()
                features = data.get('features', [])
                if features:
                    eq_lons = [f['geometry']['coordinates'][0] for f in features]
                    eq_lats = [f['geometry']['coordinates'][1] for f in features]
                    eq_mags = [f['properties']['mag'] for f in features]
                    eq_sizes = [10 * (2 ** m) for m in eq_mags]

                    print(f"    Plotting {len(features)} earthquakes...")
                    ax.scatter(eq_lons, eq_lats, s=eq_sizes, c='yellow',
                               edgecolor='black', linewidth=0.5, alpha=0.6,
                               transform=ccrs.PlateCarree(),
                               label=f'Earthquakes M≥{eq_minmag}', zorder=8)
        except Exception as e:
            print(f"    Warning: {e}")

    # ========================================================================
    # Add plate motion arrows (optional - currently disabled)
    # ========================================================================
    # To enable plate motion vectors, uncomment the section below

    # print("  Adding plate motion vectors...")
    # arrow1 = FancyArrowPatch((168, -42), (169, -41), arrowstyle='->',
    #                          mutation_scale=30, linewidth=3, color='blue',
    #                          transform=ccrs.PlateCarree(), zorder=6)
    # ax.add_patch(arrow1)
    # ax.text(168.5, -41.5, 'Pacific Plate\n~40 mm/yr', fontsize=11,
    #         fontweight='bold', color='blue', transform=ccrs.PlateCarree(),
    #         zorder=11, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # arrow2 = FancyArrowPatch((177, -39), (176, -39), arrowstyle='->',
    #                          mutation_scale=30, linewidth=3, color='orange',
    #                          transform=ccrs.PlateCarree(), zorder=6)
    # ax.add_patch(arrow2)
    # ax.text(177.5, -38.5, 'Subduction\n~45 mm/yr', fontsize=11,
    #         fontweight='bold', color='orange', transform=ccrs.PlateCarree(),
    #         zorder=11, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # ========================================================================
    # Add title and legend
    # ========================================================================
    ax.set_title('New Zealand Tectonic Setting\nPacific-Australian Plate Boundary',
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(loc='upper left', fontsize=10, framealpha=0.9, edgecolor='black')

    # Add scale bar
    scale_lon = 165
    scale_lat = -47
    scale_length = 2
    ax.plot([scale_lon, scale_lon + scale_length], [scale_lat, scale_lat],
            'k-', linewidth=3, transform=ccrs.PlateCarree(), zorder=10)
    ax.text(scale_lon + scale_length/2, scale_lat - 0.3, '~200 km',
            fontsize=10, ha='center', transform=ccrs.PlateCarree(), zorder=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # Save figure
    print(f"\nSaving figure to: {output_file}")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print("Done!")

    return fig, ax


# ============================================================================
# Example usage
# ============================================================================

if __name__ == "__main__":
    import sys

    # Example: NZ map with dual colormaps
    print("="*70)
    print("Creating New Zealand map with dual colormaps")
    print("="*70 + "\n")

    # Example with topography file (update path as needed)
    topography_file = None  # Set to your .nc or .grd file path

    if len(sys.argv) > 1:
        topography_file = sys.argv[1]
        print(f"Using topography file: {topography_file}\n")

    fig, ax = plot_nz_geological_map(
        output_file='nz_dual_cmap.png',
        topography_file=topography_file,
        land_cmap='terrain',      # Land: brown/green terrain
        ocean_cmap='Blues_r',      # Ocean: blue (reversed)
        land_vmin=0,
        land_vmax=3000,
        ocean_vmin=-6000,
        ocean_vmax=0,
        add_hillshade=True,
        show_earthquakes=True,
        eq_startdate='2020-01-01',
        eq_minmag=5.5
    )

    plt.show()

    print("\n" + "="*70)
    print("Complete!")
    print("="*70)

"""
Usage examples:

# Without topography data (basic features):
python plot_nz_geology_matplotlib.py

# With topography data:
python plot_nz_geology_matplotlib.py /path/to/nz_topo.nc

# Programmatic use:
from plot_nz_geology_matplotlib import plot_nz_geological_map

fig, ax = plot_nz_geological_map(
    output_file='my_nz_map.png',
    topography_file='nz_topo.nc',
    land_cmap='YlOrBr',        # Yellow-Orange-Brown for land
    ocean_cmap='ocean',         # Ocean colormap for bathymetry
    add_hillshade=True
)
"""
