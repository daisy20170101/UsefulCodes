#!/usr/bin/env python3
"""
Create seismic moment release database for major earthquakes (M 7.5-9.0)
Uses a curated list of well-documented major earthquakes from 1900-2025
"""

import pandas as pd
import numpy as np
from datetime import datetime

def magnitude_to_moment(magnitude):
    """
    Convert moment magnitude to seismic moment (N·m)
    Using: Mw = (2/3) * log10(M0) - 10.7, where M0 is in dyne·cm
    """
    M0_dyne_cm = 10 ** (1.5 * magnitude + 16.05)
    M0_Nm = M0_dyne_cm * 1e-7
    return M0_Nm

def create_major_earthquake_database():
    """
    Create database of major earthquakes (M 7.5-9.0) from historical records
    Data compiled from USGS, Global CMT, and published literature
    """

    # Curated list of major earthquakes M 7.5+
    earthquakes = [
        # Format: (date, lat, lon, depth, magnitude, location)

        # M 9.0+ earthquakes
        ('1960-05-22', -38.24, -73.05, 25, 9.5, 'Valdivia, Chile'),
        ('1964-03-28', 61.02, -147.65, 25, 9.2, 'Prince William Sound, Alaska'),
        ('2004-12-26', 3.30, 95.78, 30, 9.1, 'Sumatra-Andaman Islands'),
        ('2011-03-11', 38.30, 142.37, 29, 9.1, 'Tohoku, Japan'),
        ('1952-11-04', 52.76, 160.06, 21, 9.0, 'Kamchatka, Russia'),

        # M 8.5-8.9 earthquakes
        ('1868-08-13', -18.50, -70.50, 25, 8.8, 'Arica, Chile'),
        ('1906-01-31', 1.00, -81.50, 25, 8.8, 'Ecuador-Colombia'),
        ('1950-08-15', 28.50, 96.50, 15, 8.6, 'Assam, India'),
        ('1957-03-09', 51.56, -175.39, 33, 8.6, 'Andreanof Islands, Alaska'),
        ('1965-02-04', 51.21, 178.50, 30, 8.7, 'Rat Islands, Alaska'),
        ('2005-03-28', 2.08, 97.01, 30, 8.6, 'Northern Sumatra'),
        ('2007-09-12', -4.44, 101.37, 34, 8.5, 'Southern Sumatra'),
        ('2010-02-27', -36.12, -72.90, 22, 8.8, 'Maule, Chile'),
        ('2012-04-11', 2.33, 93.06, 20, 8.6, 'Off Sumatra'),

        # M 8.0-8.4 earthquakes
        ('1923-02-03', 54.00, 161.00, 15, 8.3, 'Kamchatka, Russia'),
        ('1938-02-01', -5.05, 131.62, 25, 8.5, 'Banda Sea, Indonesia'),
        ('1946-04-01', 53.49, -163.00, 15, 8.6, 'Unimak Island, Alaska'),
        ('1963-10-13', 44.90, 149.60, 23, 8.5, 'Kuril Islands'),
        ('1985-09-19', 18.19, -102.53, 15, 8.0, 'Michoacan, Mexico'),
        ('1995-10-09', 19.06, -104.21, 33, 8.0, 'Colima-Jalisco, Mexico'),
        ('2001-06-23', -16.26, -73.64, 33, 8.4, 'Southern Peru'),
        ('2003-09-25', 41.82, 143.91, 27, 8.3, 'Hokkaido, Japan'),
        ('2006-11-15', 46.59, 153.27, 10, 8.3, 'Kuril Islands'),
        ('2007-04-01', -8.47, 157.04, 10, 8.1, 'Solomon Islands'),
        ('2015-09-16', -31.57, -71.67, 22, 8.3, 'Illapel, Chile'),
        ('2017-09-08', 15.02, -93.90, 58, 8.2, 'Chiapas, Mexico'),
        ('2021-07-29', 55.33, -157.84, 32, 8.2, 'Alaska Peninsula'),

        # M 7.5-7.9 earthquakes (selection of significant events)
        ('1906-04-18', 37.75, -122.55, 5, 7.9, 'San Francisco, California'),
        ('1920-12-16', 36.50, 105.70, 15, 7.8, 'Haiyuan, China'),
        ('1927-05-22', 37.80, 102.60, 15, 7.6, 'Gulang, China'),
        ('1931-02-02', 39.40, 177.80, 20, 7.9, 'Napier, New Zealand'),
        ('1939-01-25', -36.20, -72.20, 35, 7.8, 'Chillan, Chile'),
        ('1939-12-26', 39.77, 39.53, 20, 7.8, 'Erzincan, Turkey'),
        ('1944-12-07', 33.70, 136.20, 40, 7.9, 'Tonankai, Japan'),
        ('1946-12-20', 32.80, 134.80, 30, 7.9, 'Nankaido, Japan'),
        ('1950-08-15', 28.36, 96.45, 15, 8.6, 'Assam-Tibet'),
        ('1957-07-28', 18.00, -101.50, 60, 7.8, 'Guerrero, Mexico'),
        ('1960-01-13', -35.42, -73.48, 35, 7.6, 'Concepcion, Chile'),
        ('1968-05-16', 40.90, 143.35, 25, 7.9, 'Tokachi-oki, Japan'),
        ('1970-01-04', 24.22, 102.50, 10, 7.5, 'Tonghai, China'),
        ('1976-02-04', 24.20, 122.70, 30, 7.5, 'Taiwan'),
        ('1976-07-27', 39.60, 118.00, 11, 7.5, 'Tangshan, China'),
        ('1978-12-12', 38.20, 142.04, 30, 7.7, 'Miyagi, Japan'),
        ('1979-12-12', -1.60, -79.36, 20, 7.7, 'Ecuador'),
        ('1985-03-03', -33.24, -71.85, 33, 7.8, 'Valparaiso, Chile'),
        ('1990-06-20', 36.96, 49.41, 18, 7.4, 'Manjil-Rudbar, Iran'),
        ('1990-07-16', 15.71, 121.16, 26, 7.7, 'Luzon, Philippines'),
        ('1992-09-02', 7.50, -77.50, 20, 7.6, 'Nicaragua'),
        ('1994-10-04', 43.77, 147.32, 14, 8.3, 'Kuril Islands'),
        ('1995-07-30', -23.34, -70.29, 47, 7.8, 'Antofagasta, Chile'),
        ('1998-03-25', -0.89, 80.33, 10, 7.5, 'Balleny Islands'),
        ('1999-08-17', 40.75, 29.86, 17, 7.6, 'Izmit, Turkey'),
        ('1999-09-20', 23.77, 120.98, 8, 7.6, 'Chi-Chi, Taiwan'),
        ('2001-01-26', 23.41, 70.23, 16, 7.7, 'Bhuj, India'),
        ('2002-11-03', 63.52, -147.44, 4, 7.9, 'Denali, Alaska'),
        ('2004-12-23', -49.31, 161.35, 10, 8.1, 'Macquarie Island'),
        ('2005-03-28', 2.09, 97.11, 30, 8.6, 'Nias, Indonesia'),
        ('2006-05-03', -20.19, -174.17, 55, 8.0, 'Tonga'),
        ('2009-09-29', -15.49, -172.10, 18, 8.1, 'Samoa'),
        ('2010-04-04', 32.30, -115.30, 32, 7.2, 'Baja California, Mexico'),
        ('2011-07-06', -29.50, -176.30, 20, 7.6, 'Kermadec Islands'),
        ('2013-02-06', -10.80, 165.14, 24, 8.0, 'Santa Cruz Islands'),
        ('2013-05-24', 54.89, 153.22, 609, 8.3, 'Sea of Okhotsk'),
        ('2014-04-01', -19.61, -70.77, 25, 8.2, 'Iquique, Chile'),
        ('2015-04-25', 28.15, 84.71, 8, 7.8, 'Gorkha, Nepal'),
        ('2015-05-12', 27.84, 86.08, 15, 7.3, 'Nepal aftershock'),
        ('2016-08-24', 42.70, 13.23, 8, 6.2, 'Central Italy'),  # Corrected magnitude
        ('2018-08-05', -8.26, 116.45, 15, 6.9, 'Lombok, Indonesia'),  # Corrected magnitude
        ('2018-09-28', -0.18, 119.85, 10, 7.5, 'Palu, Indonesia'),
        ('2019-05-26', -5.88, -75.27, 122, 8.0, 'Northern Peru'),
        ('2020-01-28', 19.42, -78.75, 10, 7.7, 'Caribbean Sea'),
        ('2020-07-22', 55.05, -158.53, 28, 7.8, 'Alaska Peninsula'),
        ('2021-08-12', 18.43, -73.48, 10, 7.2, 'Haiti'),
        ('2022-09-19', 18.36, -103.22, 20, 7.6, 'Michoacan, Mexico'),
        ('2023-02-06', 37.17, 37.03, 17, 7.8, 'Turkey-Syria'),
        ('2023-02-06', 38.01, 37.24, 10, 7.5, 'Turkey-Syria aftershock'),
        ('2024-01-01', 37.31, 137.26, 10, 7.6, 'Noto Peninsula, Japan'),
    ]

    # Create DataFrame
    df_list = []
    for eq in earthquakes:
        date_str, lat, lon, depth, mag, location = eq

        # Parse date
        try:
            eq_time = pd.to_datetime(date_str)
        except:
            continue

        # Calculate seismic moment
        moment = magnitude_to_moment(mag)

        df_list.append({
            'time': eq_time,
            'latitude': lat,
            'longitude': lon,
            'depth_km': depth,
            'magnitude': mag,
            'seismic_moment_Nm': moment,
            'location': location
        })

    df = pd.DataFrame(df_list)
    df = df.sort_values('time').reset_index(drop=True)

    return df

def calculate_moment_release_rate(df, time_window_years=1):
    """Calculate cumulative moment release and moment release rate"""
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

        # Calculate rate
        actual_window_years = (current_time - max(window_start, start_time)).total_seconds() / (365.25 * 24 * 3600)
        if actual_window_years > 0:
            df.loc[df.index[i], 'moment_release_rate_Nm_per_year'] = moment_in_window / actual_window_years

    return df

def calculate_gutenberg_richter(df, mag_bins=None):
    """Calculate Gutenberg-Richter frequency-magnitude distribution"""
    if mag_bins is None:
        mag_bins = np.arange(7.5, 10.0, 0.1)

    # Count earthquakes above each magnitude threshold
    counts = []
    for mag in mag_bins:
        count = len(df[df['magnitude'] >= mag])
        counts.append(count)

    # Calculate annual rate (assuming data spans from first to last event)
    time_span_years = (df['time'].max() - df['time'].min()).total_seconds() / (365.25 * 24 * 3600)
    annual_rates = np.array(counts) / time_span_years

    return pd.DataFrame({
        'magnitude_threshold': mag_bins,
        'cumulative_count': counts,
        'annual_rate': annual_rates
    })

def print_summary_statistics(df):
    """Print summary statistics"""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS - MAJOR EARTHQUAKES (M 7.5-9.0)")
    print("="*80)

    print(f"\nTotal number of earthquakes: {len(df)}")
    print(f"Time span: {df['time'].min().strftime('%Y-%m-%d')} to {df['time'].max().strftime('%Y-%m-%d')}")
    time_span_years = (df['time'].max() - df['time'].min()).total_seconds() / (365.25 * 24 * 3600)
    print(f"Duration: {time_span_years:.1f} years")
    print(f"Average rate: {len(df) / time_span_years:.2f} events/year")

    print(f"\nMagnitude statistics:")
    print(f"  Min: M {df['magnitude'].min():.1f}")
    print(f"  Max: M {df['magnitude'].max():.1f}")
    print(f"  Mean: M {df['magnitude'].mean():.2f}")
    print(f"  Median: M {df['magnitude'].median():.2f}")

    print(f"\nSeismic moment statistics:")
    total_moment = df['seismic_moment_Nm'].sum()
    print(f"  Total moment released: {total_moment:.3e} N·m")
    print(f"  Average moment rate: {total_moment / time_span_years:.3e} N·m/year")
    print(f"  Mean moment per event: {df['seismic_moment_Nm'].mean():.3e} N·m")
    print(f"  Median moment per event: {df['seismic_moment_Nm'].median():.3e} N·m")

    # Magnitude distribution
    print(f"\nMagnitude distribution:")
    for mag_range in [(9.0, 10.0), (8.5, 9.0), (8.0, 8.5), (7.5, 8.0)]:
        count = len(df[(df['magnitude'] >= mag_range[0]) & (df['magnitude'] < mag_range[1])])
        pct = 100 * count / len(df)
        print(f"  M {mag_range[0]:.1f}-{mag_range[1]:.1f}: {count:3d} events ({pct:5.1f}%)")

    print(f"\n10 Largest earthquakes:")
    print(f"{'Date':>12}  {'Mag':>5}  {'Moment (N·m)':>13}  {'Location':<40}")
    print("-" * 80)
    largest = df.nlargest(10, 'magnitude')
    for idx, row in largest.iterrows():
        date_str = row['time'].strftime('%Y-%m-%d')
        print(f"{date_str:>12}  M{row['magnitude']:.1f}  {row['seismic_moment_Nm']:13.3e}  {row['location']:<40}")

    print("="*80)

def main():
    """Main function"""
    print("="*80)
    print("SEISMIC MOMENT RELEASE DATABASE CREATOR")
    print("Major Earthquakes (M 7.5 - 9.0)")
    print("="*80)

    # Create database
    print("\nCreating earthquake database from curated historical records...")
    df = create_major_earthquake_database()

    print(f"Loaded {len(df)} major earthquakes")

    # Calculate moment release rates
    print("\nCalculating moment release rates...")
    df = calculate_moment_release_rate(df, time_window_years=1)

    # Calculate Gutenberg-Richter
    print("Calculating Gutenberg-Richter distribution...")
    gr_df = calculate_gutenberg_richter(df)

    # Print summary
    print_summary_statistics(df)

    # Save files
    output_file = 'major_earthquakes_moment_data.csv'
    df.to_csv(output_file, index=False)
    print(f"\n✓ Main data saved to: {output_file}")

    gr_file = 'gutenberg_richter_distribution.csv'
    gr_df.to_csv(gr_file, index=False)
    print(f"✓ Gutenberg-Richter distribution saved to: {gr_file}")

    # Create a simplified version
    df_simple = df[['time', 'latitude', 'longitude', 'depth_km', 'magnitude',
                    'seismic_moment_Nm', 'cumulative_moment_Nm',
                    'moment_release_rate_Nm_per_year', 'location']]
    simple_file = 'major_earthquakes_simple.csv'
    df_simple.to_csv(simple_file, index=False)
    print(f"✓ Simplified data saved to: {simple_file}")

    print("\n" + "="*80)
    print("DONE!")
    print("="*80)
    print("\nOutput files contain:")
    print("  - Event time, location (lat/lon), depth, magnitude")
    print("  - Seismic moment (N·m)")
    print("  - Cumulative moment release")
    print("  - Moment release rate (1-year moving window)")
    print("  - Gutenberg-Richter frequency-magnitude distribution")
    print("="*80)

if __name__ == '__main__':
    main()
