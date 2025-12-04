"""
Download USGS earthquake catalog data for major events
Includes magnitude, seismic moment, location, depth, and derived parameters
"""

import requests
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import json

def calculate_seismic_moment(magnitude, magnitude_type='mww'):
    """
    Calculate seismic moment from moment magnitude
    M0 = 10^(1.5 * Mw + 9.1) in N⋅m

    Parameters:
    -----------
    magnitude : float
        Moment magnitude
    magnitude_type : str
        Type of magnitude

    Returns:
    --------
    float : Seismic moment in N⋅m
    """
    if magnitude_type.lower() in ['mww', 'mw', 'mwc', 'mwb', 'mwr']:
        # Moment magnitude relationship
        M0 = 10 ** (1.5 * magnitude + 9.1)
        return M0
    else:
        # Approximate conversion for other magnitude types
        M0 = 10 ** (1.5 * magnitude + 9.1)
        return M0


def download_usgs_catalog(starttime=None, endtime=None, minmagnitude=6.0,
                          maxmagnitude=None, output_file=None):
    """
    Download earthquake catalog from USGS using their API

    Parameters:
    -----------
    starttime : str or datetime
        Start time (default: 10 years ago)
        Format: 'YYYY-MM-DD' or datetime object
    endtime : str or datetime
        End time (default: today)
        Format: 'YYYY-MM-DD' or datetime object
    minmagnitude : float
        Minimum magnitude (default: 6.0 for major earthquakes)
    maxmagnitude : float
        Maximum magnitude (default: None)
    output_file : str
        Output CSV filename (default: 'usgs_earthquake_catalog.csv')

    Returns:
    --------
    pandas.DataFrame : Earthquake catalog
    """

    # Set default time range (last 10 years)
    if endtime is None:
        endtime = datetime.now()
    elif isinstance(endtime, str):
        endtime = datetime.strptime(endtime, '%Y-%m-%d')

    if starttime is None:
        starttime = endtime - timedelta(days=365*10)
    elif isinstance(starttime, str):
        starttime = datetime.strptime(starttime, '%Y-%m-%d')

    # Format dates for USGS API
    starttime_str = starttime.strftime('%Y-%m-%d')
    endtime_str = endtime.strftime('%Y-%m-%d')

    # Default output filename
    if output_file is None:
        output_file = f'usgs_earthquake_catalog_{starttime_str}_to_{endtime_str}.csv'

    print(f"{'='*70}")
    print(f"Downloading USGS Earthquake Catalog")
    print(f"{'='*70}")
    print(f"Time range: {starttime_str} to {endtime_str}")
    print(f"Minimum magnitude: {minmagnitude}")
    if maxmagnitude:
        print(f"Maximum magnitude: {maxmagnitude}")
    print(f"Output file: {output_file}\n")

    # Build USGS API query
    base_url = "https://earthquake.usgs.gov/fdsnws/event/1/query"

    params = {
        'format': 'geojson',
        'starttime': starttime_str,
        'endtime': endtime_str,
        'minmagnitude': minmagnitude,
        'orderby': 'time-asc'
    }

    if maxmagnitude:
        params['maxmagnitude'] = maxmagnitude

    print("Querying USGS API...")

    try:
        response = requests.get(base_url, params=params, timeout=120)

        if response.status_code == 200:
            data = response.json()

            # Extract earthquake data
            earthquakes = []
            features = data.get('features', [])

            print(f"Found {len(features)} earthquakes\n")
            print(f"Processing earthquake data...")

            for feature in features:
                props = feature['properties']
                geom = feature['geometry']

                # Get all available magnitude types
                mag = props.get('mag')
                mag_type = props.get('magType', 'unknown')

                # Calculate seismic moment
                seismic_moment = calculate_seismic_moment(mag, mag_type) if mag else None

                # Convert time to datetime
                time_ms = props.get('time')
                time_dt = datetime.fromtimestamp(time_ms / 1000) if time_ms else None

                eq_data = {
                    'event_id': feature.get('id'),
                    'time': time_dt,
                    'latitude': geom['coordinates'][1] if geom else None,
                    'longitude': geom['coordinates'][0] if geom else None,
                    'depth_km': geom['coordinates'][2] if geom else None,
                    'magnitude': mag,
                    'magnitude_type': mag_type,
                    'seismic_moment_Nm': seismic_moment,
                    'place': props.get('place'),
                    'event_type': props.get('type'),
                    'status': props.get('status'),
                    'tsunami': props.get('tsunami', 0),
                    'significance': props.get('sig'),
                    'felt_reports': props.get('felt'),
                    'cdi': props.get('cdi'),  # Community Decimal Intensity
                    'mmi': props.get('mmi'),  # Modified Mercalli Intensity
                    'alert': props.get('alert'),
                    'network': props.get('net'),
                    'code': props.get('code'),
                    'ids': props.get('ids'),
                    'sources': props.get('sources'),
                    'types': props.get('types'),
                    'nst': props.get('nst'),  # Number of stations
                    'dmin': props.get('dmin'),  # Minimum distance to station
                    'gap': props.get('gap'),  # Azimuthal gap
                    'rms': props.get('rms'),  # RMS travel time residual
                    'url': props.get('url'),
                    'detail_url': props.get('detail')
                }

                earthquakes.append(eq_data)

            # Create DataFrame
            df = pd.DataFrame(earthquakes)

            # Add derived columns
            if not df.empty:
                # Calculate log10(M0) for easier comparison
                df['log10_M0'] = df['seismic_moment_Nm'].apply(
                    lambda x: np.log10(x) if x and pd.notna(x) else None
                )

                # Add year and month columns
                df['year'] = df['time'].dt.year
                df['month'] = df['time'].dt.month

                # Sort by time
                df = df.sort_values('time', ascending=False)

            # Save to CSV
            df.to_csv(output_file, index=False)

            print(f"\n{'='*70}")
            print(f"Download Complete!")
            print(f"{'='*70}")
            print(f"Total earthquakes: {len(df)}")
            print(f"Saved to: {output_file}\n")

            # Print summary statistics
            if not df.empty:
                print("Magnitude distribution:")
                print(df['magnitude'].describe())
                print(f"\nMagnitude range: {df['magnitude'].min():.1f} - {df['magnitude'].max():.1f}")
                print(f"Events with tsunami flag: {df['tsunami'].sum()}")

                # Print largest events
                print(f"\n{'='*70}")
                print(f"Top 10 Largest Earthquakes:")
                print(f"{'='*70}")
                top10 = df.nlargest(10, 'magnitude')[['time', 'magnitude', 'place', 'depth_km']]
                for idx, row in top10.iterrows():
                    print(f"M{row['magnitude']:.1f} - {row['time'].strftime('%Y-%m-%d')} - {row['place']}")

            return df

        else:
            print(f"Error: HTTP {response.status_code}")
            print(response.text)
            return None

    except Exception as e:
        print(f"Error downloading data: {e}")
        return None


def download_detailed_event_info(event_ids, output_dir='.'):
    """
    Download detailed information for specific events
    This includes finite-fault products with duration information

    Parameters:
    -----------
    event_ids : list
        List of USGS event IDs
    output_dir : str
        Output directory
    """
    from pathlib import Path

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    detailed_data = []

    print(f"\nDownloading detailed information for {len(event_ids)} events...\n")

    for event_id in event_ids:
        try:
            url = f"https://earthquake.usgs.gov/earthquakes/feed/v1.0/detail/{event_id}.geojson"
            response = requests.get(url, timeout=30)

            if response.status_code == 200:
                data = response.json()
                props = data['properties']

                # Check for finite-fault product
                products = props.get('products', {})
                finite_fault = products.get('finite-fault', [])

                event_data = {
                    'event_id': event_id,
                    'magnitude': props.get('mag'),
                    'has_finite_fault': len(finite_fault) > 0
                }

                # Extract duration from finite-fault product if available
                if finite_fault:
                    ff_props = finite_fault[0].get('properties', {})
                    event_data['rupture_duration_s'] = ff_props.get('rupture-duration')
                    event_data['rupture_velocity_km_s'] = ff_props.get('rupture-velocity')
                    event_data['model_dip'] = ff_props.get('model-dip')
                    event_data['model_length_km'] = ff_props.get('model-length')
                    event_data['model_strike'] = ff_props.get('model-strike')
                    event_data['model_top_km'] = ff_props.get('model-top')
                    event_data['model_width_km'] = ff_props.get('model-width')
                    event_data['maximum_slip_m'] = ff_props.get('maximum-slip')

                detailed_data.append(event_data)
                print(f"✓ {event_id}: M{props.get('mag'):.1f}")

            else:
                print(f"✗ {event_id}: HTTP {response.status_code}")

        except Exception as e:
            print(f"✗ {event_id}: {e}")

    # Create DataFrame
    df_detailed = pd.DataFrame(detailed_data)

    # Save to CSV
    output_file = output_path / 'detailed_event_info.csv'
    df_detailed.to_csv(output_file, index=False)

    print(f"\nSaved detailed data to: {output_file}")
    print(f"Events with finite-fault models: {df_detailed['has_finite_fault'].sum()}/{len(df_detailed)}")

    return df_detailed


# ============================================================================
# Example usage
# ============================================================================

if __name__ == "__main__":

    # Download catalog for last 10 years, magnitude >= 6.5
    print("Downloading earthquake catalog...\n")

    df = download_usgs_catalog(
        starttime='2014-01-01',
        endtime=None,  # Today
        minmagnitude=6.5,
        output_file='usgs_major_earthquakes_2014_2024.csv'
    )

    # Optional: Download detailed information for largest events
    if df is not None and not df.empty:
        print("\n" + "="*70)
        print("Downloading detailed information for M≥7.5 events...")
        print("="*70 + "\n")

        large_events = df[df['magnitude'] >= 7.5]['event_id'].tolist()

        if large_events:
            df_detailed = download_detailed_event_info(
                large_events,
                output_dir='./detailed_events'
            )

            # Merge with main catalog
            df_merged = df.merge(df_detailed, on='event_id', how='left', suffixes=('', '_detailed'))
            df_merged.to_csv('usgs_earthquakes_with_duration.csv', index=False)

            print(f"\nMerged catalog saved to: usgs_earthquakes_with_duration.csv")


# ============================================================================
# Alternative: Download specific time ranges or regions
# ============================================================================
"""
# Last 5 years, M≥7.0
df = download_usgs_catalog(
    starttime='2019-01-01',
    minmagnitude=7.0,
    output_file='usgs_m7_2019_2024.csv'
)

# All M≥8.0 since 2000
df = download_usgs_catalog(
    starttime='2000-01-01',
    minmagnitude=8.0,
    output_file='usgs_m8_2000_2024.csv'
)

# Specific date range
df = download_usgs_catalog(
    starttime='2010-01-01',
    endtime='2020-12-31',
    minmagnitude=6.5,
    output_file='usgs_earthquakes_2010s.csv'
)
"""
