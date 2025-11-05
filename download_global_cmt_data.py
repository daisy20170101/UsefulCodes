#!/usr/bin/env python3
"""
Download seismic moment release data for major earthquakes (M 7.5-9.0)
Uses Global CMT Catalog
"""

import requests
import pandas as pd
import numpy as np
from datetime import datetime
import re
import io

def parse_ndk_line(line):
    """Parse a line from NDK format (used by Global CMT)"""
    # This is a simplified parser for NDK format
    # Full format details: http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/allorder.ndk_explained
    return line

def download_global_cmt_data(start_year=1976, end_year=None, min_mag=7.5, max_mag=9.0):
    """
    Download earthquake data from Global CMT catalog

    Parameters:
    -----------
    start_year : int
        Start year
    end_year : int
        End year (default: current year)
    min_mag : float
        Minimum magnitude
    max_mag : float
        Maximum magnitude

    Returns:
    --------
    pandas.DataFrame : Earthquake data
    """
    if end_year is None:
        end_year = datetime.now().year

    print(f"Downloading earthquake data from Global CMT Catalog...")
    print(f"Year range: {start_year} to {end_year}")
    print(f"Magnitude range: {min_mag} to {max_mag}")
    print(f"Note: Global CMT catalog starts from 1976")

    # Global CMT web catalog URL
    base_url = "https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form"

    # Prepare form data
    form_data = {
        'itype': 'ymd',  # Search by year-month-day
        'yr': str(start_year),
        'mo': '1',
        'day': '1',
        'otype': 'ymd',
        'oyr': str(end_year),
        'omo': '12',
        'oday': '31',
        'jyr': '1976',
        'jday': '1',
        'ojyr': str(end_year),
        'ojday': '365',
        'nday': '1',
        'lmw': str(min_mag),
        'umw': str(max_mag),
        'lms': '0',
        'ums': '10',
        'lmb': '0',
        'umb': '10',
        'llat': '-90',
        'ulat': '90',
        'llon': '-180',
        'ulon': '180',
        'lhd': '0',
        'uhd': '1000',
        'lts': '-9999',
        'uts': '9999',
        'lpe1': '0',
        'upe1': '90',
        'lpe2': '0',
        'upe2': '90',
        'list': 'Comprehensive'
    }

    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36'
    }

    try:
        print("Sending request to Global CMT...")
        response = requests.post(base_url, data=form_data, headers=headers, timeout=120)
        response.raise_for_status()

        # Parse the response
        content = response.text

        # The response should contain NDK format data
        # Parse it line by line
        earthquakes = []
        lines = content.split('\n')

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            # NDK format has 5 lines per event
            if line and not line.startswith('<') and not line.startswith('Error'):
                # Line 1: Hypocenter line
                if len(line) > 60:
                    try:
                        # Parse basic info from first line
                        # Format: catalog, date, time, latitude, longitude, depth, magnitudes, location
                        parts = line.split()

                        # Extract date/time
                        if len(parts) > 5:
                            # Look for year-month-day pattern
                            date_str = parts[0] if '/' in parts[0] else None

                            if i + 4 < len(lines):
                                # Read next 4 lines for complete event info
                                line2 = lines[i + 1] if i + 1 < len(lines) else ""
                                line3 = lines[i + 2] if i + 2 < len(lines) else ""
                                line4 = lines[i + 3] if i + 3 < len(lines) else ""
                                line5 = lines[i + 4] if i + 4 < len(lines) else ""

                                # Parse moment tensor from line 4 and 5
                                # This is a simplified parsing
                                eq_data = parse_ndk_event(line, line2, line3, line4, line5)
                                if eq_data and eq_data['magnitude'] >= min_mag and eq_data['magnitude'] <= max_mag:
                                    earthquakes.append(eq_data)

                                i += 5
                                continue
                    except Exception as e:
                        pass

            i += 1

        if len(earthquakes) > 0:
            print(f"Successfully parsed {len(earthquakes)} earthquakes")
            df = pd.DataFrame(earthquakes)
            df = df.sort_values('time').reset_index(drop=True)
            return df
        else:
            print("No earthquakes found in response. Trying alternative method...")
            return None

    except requests.exceptions.RequestException as e:
        print(f"Error downloading data: {e}")
        return None

def parse_ndk_event(line1, line2, line3, line4, line5):
    """
    Parse a 5-line NDK format event
    Simplified parser for key information
    """
    try:
        # Line 1: Basic hypocenter info
        # Format is complex, but key fields are at fixed positions
        if len(line1) < 80:
            return None

        # Extract date and time (positions 5-26)
        date_time_str = line1[5:26].strip()

        # Extract coordinates and depth
        lat_str = line1[27:33].strip()
        lon_str = line1[34:41].strip()
        depth_str = line1[42:47].strip()

        # Line 2: CMT info including magnitude
        if len(line2) < 60:
            return None

        # Event name contains date info
        event_name = line2[0:16].strip()

        # Extract moment magnitude and other mags from line 2
        # Position 48-51 typically has Mw value

        # Line 4 contains moment tensor exponents and values
        # Position 0-2 has exponent for moment
        exp_str = line4[0:2].strip() if len(line4) > 2 else "0"

        # Moment tensor values are in dyne-cm
        # Parse exponent
        try:
            exponent = int(exp_str)
        except:
            exponent = 0

        # Create earthquake dictionary
        # Extract magnitude from event name (usually in format: MMDDHHMMCmagnitude)
        mag_match = re.search(r'([0-9]\.[0-9])', line1[50:80])
        magnitude = float(mag_match.group(1)) if mag_match else 0.0

        # Parse date from event name (format: YYMMDDHHMMSS)
        if len(event_name) >= 12:
            year = int("20" + event_name[0:2]) if int(event_name[0:2]) < 50 else int("19" + event_name[0:2])
            month = int(event_name[2:4])
            day = int(event_name[4:6])
            hour = int(event_name[6:8])
            minute = int(event_name[8:10])
            eq_time = datetime(year, month, day, hour, minute)
        else:
            eq_time = datetime.now()

        # Calculate scalar seismic moment from exponent
        # M0 = 10^exponent dyne-cm
        M0_dyne_cm = 10 ** (exponent)
        M0_Nm = M0_dyne_cm * 1e-7

        eq_data = {
            'time': eq_time,
            'latitude': float(lat_str) if lat_str else 0.0,
            'longitude': float(lon_str) if lon_str else 0.0,
            'depth_km': float(depth_str) if depth_str else 0.0,
            'magnitude': magnitude,
            'seismic_moment_Nm': M0_Nm,
            'event_name': event_name
        }

        return eq_data

    except Exception as e:
        return None

def download_usgs_recent(min_mag=7.5, max_mag=9.0, years_back=30):
    """
    Download recent earthquake data from USGS (smaller date range)
    """
    end_date = datetime.now()
    start_date = end_date.replace(year=end_date.year - years_back)

    print(f"Downloading recent earthquake data from USGS...")
    print(f"Date range: {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}")
    print(f"Magnitude range: {min_mag} to {max_mag}")

    base_url = "https://earthquake.usgs.gov/fdsnws/event/1/query"

    params = {
        'format': 'csv',
        'starttime': start_date.strftime('%Y-%m-%d'),
        'endtime': end_date.strftime('%Y-%m-%d'),
        'minmagnitude': min_mag,
        'maxmagnitude': max_mag,
        'orderby': 'time'
    }

    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36'
    }

    try:
        response = requests.get(base_url, params=params, headers=headers, timeout=120)
        response.raise_for_status()

        # Parse CSV response
        df = pd.read_csv(io.StringIO(response.text))

        print(f"Successfully downloaded data for {len(df)} earthquakes")

        # Extract relevant columns and rename
        df_clean = pd.DataFrame({
            'time': pd.to_datetime(df['time']),
            'latitude': df['latitude'],
            'longitude': df['longitude'],
            'depth_km': df['depth'],
            'magnitude': df['mag'],
            'magnitude_type': df['magType'],
            'place': df['place'],
            'event_id': df['id']
        })

        # Calculate seismic moment
        df_clean['seismic_moment_Nm'] = magnitude_to_moment(df_clean['magnitude'])

        return df_clean

    except Exception as e:
        print(f"Error downloading data: {e}")
        return None

def magnitude_to_moment(magnitude):
    """Convert magnitude to seismic moment (N·m)"""
    # Handle both scalar and array inputs
    M0_dyne_cm = 10 ** (1.5 * magnitude + 16.05)
    M0_Nm = M0_dyne_cm * 1e-7
    return M0_Nm

def calculate_moment_release_rate(df, time_window_years=1):
    """Calculate cumulative moment release and moment release rate"""
    df = df.copy()

    # Sort by time
    df = df.sort_values('time').reset_index(drop=True)

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

def print_summary_statistics(df):
    """Print summary statistics"""
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

    print(f"\nLargest events:")
    largest = df.nlargest(min(10, len(df)), 'magnitude')[['time', 'magnitude', 'seismic_moment_Nm']]
    for idx, row in largest.iterrows():
        date_str = row['time'].strftime('%Y-%m-%d')
        print(f"  M{row['magnitude']:.1f} - {date_str} - Moment: {row['seismic_moment_Nm']:.3e} N·m")

    print("="*70)

def main():
    """Main function"""
    print("="*70)
    print("SEISMIC MOMENT RELEASE DATA DOWNLOADER")
    print("Major Earthquakes (M 7.5 - 9.0)")
    print("="*70)

    # Try USGS first (recent data, usually more reliable)
    print("\nAttempting to download from USGS (last 30 years)...")
    df = download_usgs_recent(min_mag=7.5, max_mag=9.0, years_back=30)

    if df is None or len(df) == 0:
        print("\nUSGS download failed. Trying Global CMT catalog...")
        df = download_global_cmt_data(start_year=1976, min_mag=7.5, max_mag=9.0)

    if df is None or len(df) == 0:
        print("\nERROR: Could not download data from any source.")
        print("Please check your internet connection and try again.")
        return

    # Calculate moment release rates
    print("\nCalculating moment release rates...")
    df = calculate_moment_release_rate(df, time_window_years=1)

    # Print summary
    print_summary_statistics(df)

    # Save to file
    output_file = 'major_earthquakes_moment_data.csv'
    df.to_csv(output_file, index=False)
    print(f"\nData saved to: {output_file}")
    print(f"Date range: {df['time'].min()} to {df['time'].max()}")

    print("\n" + "="*70)
    print("DONE!")
    print("="*70)

if __name__ == '__main__':
    main()
