#!/usr/bin/env python3
"""
Download seismic moment release data for major earthquakes (M 7.5-9.0)
Uses USGS Earthquake Catalog API
"""

import requests
import pandas as pd
import numpy as np
from datetime import datetime
import json

def magnitude_to_moment(magnitude):
    """
    Convert moment magnitude to seismic moment (N·m)
    Using the relationship: Mw = (2/3) * log10(M0) - 10.7
    where M0 is in dyne·cm, converted to N·m

    Parameters:
    -----------
    magnitude : float
        Moment magnitude (Mw)

    Returns:
    --------
    float : Seismic moment in N·m
    """
    # M0 in dyne·cm
    M0_dyne_cm = 10 ** (1.5 * magnitude + 16.05)
    # Convert to N·m (1 dyne·cm = 1e-7 N·m)
    M0_Nm = M0_dyne_cm * 1e-7
    return M0_Nm

def download_usgs_earthquake_data(start_date='1900-01-01', end_date=None,
                                  min_mag=7.5, max_mag=9.0):
    """
    Download earthquake data from USGS earthquake catalog

    Parameters:
    -----------
    start_date : str
        Start date in format 'YYYY-MM-DD'
    end_date : str
        End date in format 'YYYY-MM-DD' (default: today)
    min_mag : float
        Minimum magnitude
    max_mag : float
        Maximum magnitude

    Returns:
    --------
    pandas.DataFrame : Earthquake data
    """
    if end_date is None:
        end_date = datetime.now().strftime('%Y-%m-%d')

    print(f"Downloading earthquake data from USGS...")
    print(f"Date range: {start_date} to {end_date}")
    print(f"Magnitude range: {min_mag} to {max_mag}")

    # USGS API endpoint
    base_url = "https://earthquake.usgs.gov/fdsnws/event/1/query"

    params = {
        'format': 'geojson',
        'starttime': start_date,
        'endtime': end_date,
        'minmagnitude': min_mag,
        'maxmagnitude': max_mag,
        'orderby': 'time'
    }

    # Add headers to avoid 403 errors
    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }

    try:
        response = requests.get(base_url, params=params, headers=headers, timeout=120)
        response.raise_for_status()
        data = response.json()

        print(f"Successfully downloaded data for {data['metadata']['count']} earthquakes")

        # Extract relevant information
        earthquakes = []
        for feature in data['features']:
            props = feature['properties']
            coords = feature['geometry']['coordinates']

            eq_data = {
                'time': pd.to_datetime(props['time'], unit='ms'),
                'latitude': coords[1],
                'longitude': coords[0],
                'depth_km': coords[2],
                'magnitude': props['mag'],
                'magnitude_type': props.get('magType', 'unknown'),
                'place': props.get('place', 'unknown'),
                'event_id': feature['id']
            }

            # Calculate seismic moment
            eq_data['seismic_moment_Nm'] = magnitude_to_moment(props['mag'])

            earthquakes.append(eq_data)

        df = pd.DataFrame(earthquakes)
        df = df.sort_values('time').reset_index(drop=True)

        return df

    except requests.exceptions.RequestException as e:
        print(f"Error downloading data: {e}")
        return None

def calculate_moment_release_rate(df, time_window_years=1):
    """
    Calculate cumulative moment release and moment release rate

    Parameters:
    -----------
    df : pandas.DataFrame
        Earthquake dataframe with 'time' and 'seismic_moment_Nm' columns
    time_window_years : float
        Time window for rate calculation (in years)

    Returns:
    --------
    pandas.DataFrame : Dataframe with cumulative and rate information
    """
    df = df.copy()

    # Calculate cumulative moment
    df['cumulative_moment_Nm'] = df['seismic_moment_Nm'].cumsum()

    # Calculate time in decimal years from start
    start_time = df['time'].iloc[0]
    df['years_since_start'] = (df['time'] - start_time).dt.total_seconds() / (365.25 * 24 * 3600)

    # Calculate moment release rate (moving window)
    time_window_days = time_window_years * 365.25
    df['moment_release_rate_Nm_per_year'] = np.nan

    for i in range(len(df)):
        current_time = df['time'].iloc[i]
        window_start = current_time - pd.Timedelta(days=time_window_days)

        # Get events in the time window
        mask = (df['time'] >= window_start) & (df['time'] <= current_time)
        moment_in_window = df.loc[mask, 'seismic_moment_Nm'].sum()

        # Calculate rate (moment per year)
        actual_window_years = (current_time - max(window_start, start_time)).total_seconds() / (365.25 * 24 * 3600)
        if actual_window_years > 0:
            df.loc[df.index[i], 'moment_release_rate_Nm_per_year'] = moment_in_window / actual_window_years

    return df

def save_data(df, filename='major_earthquakes_moment_data.csv'):
    """
    Save earthquake data to CSV file

    Parameters:
    -----------
    df : pandas.DataFrame
        Earthquake dataframe
    filename : str
        Output filename
    """
    df.to_csv(filename, index=False)
    print(f"\nData saved to: {filename}")
    print(f"Total earthquakes: {len(df)}")
    print(f"Total seismic moment released: {df['seismic_moment_Nm'].sum():.3e} N·m")
    print(f"Date range: {df['time'].min()} to {df['time'].max()}")

def print_summary_statistics(df):
    """
    Print summary statistics of the earthquake data
    """
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)

    print(f"\nTotal number of earthquakes: {len(df)}")
    print(f"\nMagnitude statistics:")
    print(f"  Min: {df['magnitude'].min():.2f}")
    print(f"  Max: {df['magnitude'].max():.2f}")
    print(f"  Mean: {df['magnitude'].mean():.2f}")
    print(f"  Median: {df['magnitude'].median():.2f}")

    print(f"\nSeismic moment statistics:")
    print(f"  Total moment released: {df['seismic_moment_Nm'].sum():.3e} N·m")
    print(f"  Mean moment per event: {df['seismic_moment_Nm'].mean():.3e} N·m")
    print(f"  Median moment per event: {df['seismic_moment_Nm'].median():.3e} N·m")

    print(f"\nLargest events:")
    largest = df.nlargest(10, 'magnitude')[['time', 'magnitude', 'place', 'seismic_moment_Nm']]
    for idx, row in largest.iterrows():
        print(f"  M{row['magnitude']:.1f} - {row['time'].strftime('%Y-%m-%d')} - {row['place']}")
        print(f"         Moment: {row['seismic_moment_Nm']:.3e} N·m")

    print("="*70)

def main():
    """
    Main function to download and process earthquake data
    """
    print("="*70)
    print("SEISMIC MOMENT RELEASE DATA DOWNLOADER")
    print("Major Earthquakes (M 7.5 - 9.0)")
    print("="*70)

    # Download data
    df = download_usgs_earthquake_data(
        start_date='1900-01-01',
        min_mag=7.5,
        max_mag=9.0
    )

    if df is None or len(df) == 0:
        print("No data downloaded. Exiting.")
        return

    # Calculate moment release rates
    print("\nCalculating moment release rates...")
    df = calculate_moment_release_rate(df, time_window_years=1)

    # Print summary
    print_summary_statistics(df)

    # Save to file
    save_data(df, 'major_earthquakes_moment_data.csv')

    print("\n" + "="*70)
    print("DONE!")
    print("="*70)

    # Additional output: save a simplified version
    df_simple = df[['time', 'latitude', 'longitude', 'depth_km', 'magnitude',
                    'seismic_moment_Nm', 'cumulative_moment_Nm',
                    'moment_release_rate_Nm_per_year', 'place']]
    df_simple.to_csv('major_earthquakes_moment_data_simple.csv', index=False)
    print(f"Simplified data saved to: major_earthquakes_moment_data_simple.csv")

if __name__ == '__main__':
    main()
