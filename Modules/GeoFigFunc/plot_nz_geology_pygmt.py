"""
Plot New Zealand geological map with plate boundaries using PyGMT
Highlights the Pacific-Australian plate boundary
"""

import pygmt
import pandas as pd
import numpy as np

def plot_nz_geological_map(output_file='nz_geological_map.png',
                           show_earthquakes=False,
                           eq_startdate='2014-01-01',
                           eq_minmag=5.0):
    """
    Create a geological map of New Zealand with plate boundaries

    Parameters:
    -----------
    output_file : str
        Output filename (png, pdf, etc.)
    show_earthquakes : bool
        Whether to plot recent earthquakes
    eq_startdate : str
        Start date for earthquake data (if show_earthquakes=True)
    eq_minmag : float
        Minimum magnitude for earthquakes
    """

    print("Creating New Zealand geological map with PyGMT...")

    # Define New Zealand region
    # [west, east, south, north]
    region = [164, 180, -48, -34]

    # Create figure
    fig = pygmt.Figure()

    # ========================================================================
    # Base map with topography/bathymetry
    # ========================================================================
    print("  Plotting topography and bathymetry...")

    fig.basemap(
        region=region,
        projection="M15c",  # Mercator projection, 15 cm width
        frame=["WSne", "xa4f2", "ya4f2"]  # Frame with lat/lon labels
    )

    # Plot Earth relief (topography/bathymetry)
    # Use high resolution for New Zealand
    fig.grdimage(
        "@earth_relief_30s",  # 30 arc-second resolution (~1 km)
        region=region,
        cmap="geo",  # Geological colormap (brown land, blue ocean)
        shading=True,
        transparency=20
    )

    # Add colorbar for elevation
    fig.colorbar(
        frame=["x+lElevation", "y+lm"],
        position="JMR+o1c/0c+w8c/0.5c"
    )

    # ========================================================================
    # Plot coastlines and borders
    # ========================================================================
    print("  Plotting coastlines...")

    fig.coast(
        shorelines="1/0.5p,black",
        borders="1/0.3p,gray50",
        area_thresh=100,  # Don't plot small islands < 100 km²
        resolution="h"  # High resolution
    )

    # ========================================================================
    # Plot plate boundaries
    # ========================================================================
    print("  Plotting plate boundaries...")

    # Download and plot plate boundaries from Nuvel1 model
    # Bird (2003) plate boundary dataset
    fig.plot(
        data="@PB2002_boundaries.txt",
        region=region,
        pen="2p,red",
        label="Plate Boundary"
    )

    # Alternative: Plot specific plate boundaries manually
    # Alpine Fault (simplified)
    alpine_fault = np.array([
        [166.5, -44.0],
        [167.0, -43.5],
        [167.5, -43.0],
        [168.0, -42.5],
        [168.5, -42.0],
        [169.0, -41.5],
        [169.5, -41.0],
        [170.0, -40.5],
        [170.5, -40.0],
    ])

    fig.plot(
        x=alpine_fault[:, 0],
        y=alpine_fault[:, 1],
        pen="3p,darkred,-",
        label="Alpine Fault"
    )

    # Hikurangi Subduction Zone (simplified)
    hikurangi = np.array([
        [178.5, -37.5],
        [178.0, -38.0],
        [177.5, -38.5],
        [177.0, -39.0],
        [176.5, -39.5],
        [176.0, -40.0],
        [175.5, -40.5],
        [175.0, -41.0],
        [174.5, -41.5],
        [174.0, -42.0],
    ])

    fig.plot(
        x=hikurangi[:, 0],
        y=hikurangi[:, 1],
        pen="3p,orange,-",
        label="Hikurangi Subduction"
    )

    # Puysegur Trench (southwestern subduction)
    puysegur = np.array([
        [166.0, -45.5],
        [165.5, -46.0],
        [165.0, -46.5],
    ])

    fig.plot(
        x=puysegur[:, 0],
        y=puysegur[:, 1],
        pen="3p,orange,-",
        label="Puysegur Trench"
    )

    # ========================================================================
    # Plot major faults
    # ========================================================================
    print("  Plotting major faults...")

    # Wellington Fault
    wellington_fault = np.array([
        [174.8, -41.3],
        [174.85, -41.2],
        [174.9, -41.1],
    ])

    fig.plot(
        x=wellington_fault[:, 0],
        y=wellington_fault[:, 1],
        pen="1.5p,red,--",
        label="Major Faults"
    )

    # Hope Fault
    hope_fault = np.array([
        [172.0, -42.5],
        [172.5, -42.4],
        [173.0, -42.3],
        [173.5, -42.2],
    ])

    fig.plot(
        x=hope_fault[:, 0],
        y=hope_fault[:, 1],
        pen="1.5p,red,--"
    )

    # ========================================================================
    # Add major cities
    # ========================================================================
    print("  Adding major cities...")

    cities = pd.DataFrame({
        'name': ['Auckland', 'Wellington', 'Christchurch', 'Dunedin', 'Hamilton'],
        'lon': [174.76, 174.78, 172.64, 170.50, 175.28],
        'lat': [-36.85, -41.29, -43.53, -45.87, -37.78]
    })

    fig.plot(
        x=cities['lon'],
        y=cities['lat'],
        style="c0.3c",
        fill="white",
        pen="1p,black"
    )

    # Add city labels
    for idx, city in cities.iterrows():
        fig.text(
            x=city['lon'],
            y=city['lat'],
            text=city['name'],
            font="10p,Helvetica-Bold,black",
            justify="BL",
            offset="0.2c/0.2c",
            fill="white@30",
            pen="0.5p,black"
        )

    # ========================================================================
    # Optionally plot earthquakes
    # ========================================================================
    if show_earthquakes:
        print(f"  Downloading earthquakes M≥{eq_minmag} since {eq_startdate}...")

        try:
            # Download earthquake data from USGS
            import requests
            from datetime import datetime

            url = "https://earthquake.usgs.gov/fdsnws/event/1/query"
            params = {
                'format': 'geojson',
                'starttime': eq_startdate,
                'endtime': datetime.now().strftime('%Y-%m-%d'),
                'minmagnitude': eq_minmag,
                'minlatitude': region[2],
                'maxlatitude': region[3],
                'minlongitude': region[0],
                'maxlongitude': region[1]
            }

            response = requests.get(url, params=params, timeout=60)

            if response.status_code == 200:
                data = response.json()
                features = data.get('features', [])

                if features:
                    eq_lons = []
                    eq_lats = []
                    eq_mags = []

                    for feat in features:
                        coords = feat['geometry']['coordinates']
                        mag = feat['properties']['mag']
                        eq_lons.append(coords[0])
                        eq_lats.append(coords[1])
                        eq_mags.append(mag)

                    # Convert magnitude to symbol size
                    # Size = 0.05 * 2^mag (exponential scaling)
                    eq_sizes = [0.05 * (2 ** m) for m in eq_mags]

                    print(f"    Plotting {len(features)} earthquakes...")

                    # Plot earthquakes
                    fig.plot(
                        x=eq_lons,
                        y=eq_lats,
                        style="cc",
                        size=eq_sizes,
                        fill="yellow",
                        pen="0.5p,black",
                        transparency=30,
                        label=f"Earthquakes M≥{eq_minmag}"
                    )
                else:
                    print("    No earthquakes found in region")
        except Exception as e:
            print(f"    Warning: Could not download earthquake data: {e}")

    # ========================================================================
    # Add tectonic motion arrows
    # ========================================================================
    print("  Adding plate motion vectors...")

    # Pacific plate motion relative to Australian plate
    # Approximate convergence rates

    # Western side (Alpine Fault): oblique convergence ~35-40 mm/yr
    fig.plot(
        x=[168],
        y=[-42],
        style="v0.5c+e+a45",  # Vector arrow
        direction=[[45], [3]],  # Azimuth and length
        fill="blue",
        pen="1.5p,blue"
    )

    fig.text(
        x=168.5,
        y=-41.5,
        text="Pacific Plate",
        font="12p,Helvetica-Bold,blue",
        justify="BL"
    )

    # Eastern side (Hikurangi): subduction ~40-50 mm/yr
    fig.plot(
        x=[177],
        y=[-39],
        style="v0.5c+e+a270",  # Vector arrow pointing west
        direction=[[270], [3]],
        fill="orange",
        pen="1.5p,orange"
    )

    fig.text(
        x=177.5,
        y=-38.5,
        text="Subduction",
        font="12p,Helvetica-Bold,orange",
        justify="BL"
    )

    # ========================================================================
    # Add legend
    # ========================================================================
    print("  Adding legend...")

    fig.legend(
        position="JTL+jTL+o0.5c",
        box="+gwhite+p1p,black"
    )

    # ========================================================================
    # Add title and annotations
    # ========================================================================
    fig.text(
        x=172,
        y=-33.5,
        text="New Zealand Tectonic Setting",
        font="20p,Helvetica-Bold,black",
        justify="CM"
    )

    fig.text(
        x=172,
        y=-34,
        text="Pacific-Australian Plate Boundary",
        font="14p,Helvetica,gray30",
        justify="CM"
    )

    # Add scale bar
    fig.basemap(
        map_scale="jBL+w200k+o1c/1c+f+l"
    )

    # ========================================================================
    # Save figure
    # ========================================================================
    print(f"\nSaving figure to: {output_file}")
    fig.savefig(output_file, dpi=300)

    print("Done!")

    return fig


def plot_detailed_kaikoura_region(output_file='kaikoura_region_detailed.png'):
    """
    Create a detailed map of the Kaikoura region
    Focused on the 2016 M7.8 Kaikoura earthquake area
    """

    print("Creating detailed Kaikoura region map...")

    # Kaikoura region
    region = [171, 175, -43, -41]

    fig = pygmt.Figure()

    # Base map
    fig.basemap(
        region=region,
        projection="M20c",
        frame=["WSne+tKaikoura Region - 2016 M7.8 Earthquake", "xa1f0.5", "ya1f0.5"]
    )

    # High-resolution topography
    fig.grdimage(
        "@earth_relief_15s",  # 15 arc-second (~500m)
        region=region,
        cmap="geo",
        shading=True
    )

    fig.colorbar(
        frame=["x+lElevation", "y+lm"],
        position="JMR+o1c/0c+w10c/0.5c"
    )

    # Coastline
    fig.coast(
        shorelines="1/0.5p,black",
        resolution="f"  # Full resolution
    )

    # Plot major faults in the region
    hope_fault = np.array([
        [172.0, -42.5],
        [172.5, -42.4],
        [173.0, -42.3],
        [173.5, -42.2],
    ])

    fig.plot(
        x=hope_fault[:, 0],
        y=hope_fault[:, 1],
        pen="3p,red,-",
        label="Hope Fault"
    )

    # Kaikoura epicenter (2016-11-13, M7.8)
    epicenter_lon = 172.74
    epicenter_lat = -42.69

    fig.plot(
        x=[epicenter_lon],
        y=[epicenter_lat],
        style="a0.6c",  # Star
        fill="red",
        pen="1.5p,black",
        label="M7.8 Kaikoura 2016"
    )

    # Add Kaikoura town
    fig.plot(
        x=[173.68],
        y=[-42.40],
        style="c0.3c",
        fill="yellow",
        pen="1p,black"
    )

    fig.text(
        x=173.68,
        y=-42.35,
        text="Kaikoura",
        font="12p,Helvetica-Bold,black",
        justify="BC"
    )

    # Legend
    fig.legend(
        position="JTR+jTR+o0.5c",
        box="+gwhite+p1p,black"
    )

    # Scale bar
    fig.basemap(
        map_scale="jBL+w50k+o1c/1c+f+l"
    )

    # Save
    print(f"Saving to: {output_file}")
    fig.savefig(output_file, dpi=300)

    print("Done!")

    return fig


# ============================================================================
# Example usage
# ============================================================================

if __name__ == "__main__":

    # Plot full New Zealand map
    print("="*70)
    print("Plot 1: Full New Zealand geological map")
    print("="*70 + "\n")

    fig1 = plot_nz_geological_map(
        output_file='nz_geological_map.png',
        show_earthquakes=True,
        eq_startdate='2020-01-01',
        eq_minmag=5.5
    )

    # Plot detailed Kaikoura region
    print("\n" + "="*70)
    print("Plot 2: Detailed Kaikoura region")
    print("="*70 + "\n")

    fig2 = plot_detailed_kaikoura_region(
        output_file='kaikoura_region_detailed.png'
    )

    print("\n" + "="*70)
    print("All plots completed!")
    print("="*70)


# ============================================================================
# Additional examples
# ============================================================================
"""
# Just the map without earthquakes
fig = plot_nz_geological_map(
    output_file='nz_geology_only.png',
    show_earthquakes=False
)

# With more recent earthquakes
fig = plot_nz_geological_map(
    output_file='nz_recent_earthquakes.png',
    show_earthquakes=True,
    eq_startdate='2023-01-01',
    eq_minmag=4.0
)

# Save as PDF instead
fig = plot_nz_geological_map(
    output_file='nz_geological_map.pdf',
    show_earthquakes=True
)
"""
