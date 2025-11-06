#!/usr/bin/env python3
"""
Download seismic waveform data from GeoNet AWS S3 bucket

GeoNet provides open access to New Zealand seismic data via AWS S3.
This script downloads MiniSEED waveform files for specified stations and time periods.

Data source: s3://geonet-open-data/waveforms/miniseed/

Requirements:
    pip install obspy boto3
    # or use AWS CLI: apt-get install awscli (or brew install awscli on macOS)

Usage:
    python download_geonet_aws_data.py
    python download_geonet_aws_data.py --date 2025-11-06 --stations WEL,BFZ,PUZ
"""

import os
import sys
import argparse
import subprocess
from datetime import datetime, timedelta
from pathlib import Path


def check_aws_cli():
    """Check if AWS CLI is installed"""
    try:
        result = subprocess.run(['aws', '--version'],
                              capture_output=True,
                              text=True,
                              check=True)
        print(f"✓ AWS CLI found: {result.stdout.strip()}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("✗ AWS CLI not found")
        print("\nPlease install AWS CLI:")
        print("  Ubuntu/Debian: sudo apt-get install awscli")
        print("  macOS: brew install awscli")
        print("  pip: pip install awscli")
        return False


def list_geonet_s3_structure(date_str=None, max_depth=3):
    """
    List GeoNet S3 bucket structure

    Parameters:
    -----------
    date_str : str, optional
        Date in YYYY-MM-DD format to explore
    max_depth : int
        Maximum depth to explore
    """

    base_path = "s3://geonet-open-data/waveforms/miniseed/"

    if date_str:
        # Parse date and create path
        date_obj = datetime.strptime(date_str, "%Y-%m-%d")
        year = date_obj.strftime("%Y")
        doy = date_obj.strftime("%j")  # Day of year

        path = f"{base_path}{year}/{year}.{doy}/"
        print(f"\nExploring: {path}")
    else:
        path = base_path
        print(f"\nExploring: {path}")

    try:
        # List contents
        cmd = ['aws', 's3', 'ls', '--no-sign-request', path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        if result.stdout:
            print("Contents:")
            lines = result.stdout.strip().split('\n')
            for i, line in enumerate(lines[:50]):  # Show first 50 items
                print(f"  {line}")

            if len(lines) > 50:
                print(f"  ... and {len(lines) - 50} more items")

            return lines
        else:
            print("No contents found")
            return []

    except subprocess.CalledProcessError as e:
        print(f"Error listing S3 bucket: {e}")
        print(f"stderr: {e.stderr}")
        return []


def download_geonet_station_data(
    date_str,
    stations=None,
    network="NZ",
    location="10",
    channels=None,
    output_dir="./geonet_data",
    hours_before=0,
    hours_after=24
):
    """
    Download GeoNet waveform data from AWS S3

    GeoNet S3 structure:
    s3://geonet-open-data/waveforms/miniseed/YEAR/YEAR.DOY/NETWORK.STATION.LOCATION.CHANNEL.YEAR.DOY

    Example:
    s3://geonet-open-data/waveforms/miniseed/2025/2025.310/NZ.WEL.10.HHZ.D.2025.310

    Parameters:
    -----------
    date_str : str
        Date in YYYY-MM-DD format
    stations : list of str
        Station codes (e.g., ['WEL', 'BFZ', 'PUZ'])
    network : str
        Network code (default: 'NZ' for GeoNet)
    location : str
        Location code (default: '10')
    channels : list of str
        Channel codes (e.g., ['HHZ', 'HHN', 'HHE'])
        If None, downloads all available channels
    output_dir : str
        Output directory
    hours_before : int
        Hours before the specified date to download
    hours_after : int
        Hours after the specified date to download

    Returns:
    --------
    downloaded_files : list
        List of downloaded file paths
    """

    # Parse date
    date_obj = datetime.strptime(date_str, "%Y-%m-%d")

    # Calculate date range
    start_date = date_obj - timedelta(hours=hours_before)
    end_date = date_obj + timedelta(hours=hours_after)

    # Generate list of dates to download
    dates_to_download = []
    current = start_date
    while current <= end_date:
        dates_to_download.append(current)
        current += timedelta(days=1)

    print(f"Date range: {start_date.date()} to {end_date.date()}")
    print(f"Dates to process: {len(dates_to_download)}")

    # Default stations (major New Zealand seismic stations)
    if stations is None:
        stations = ['WEL', 'BFZ', 'PUZ', 'MXZ', 'THZ', 'KNZ', 'ODZ', 'DCZ']
        print(f"Using default stations: {', '.join(stations)}")
    else:
        print(f"Stations: {', '.join(stations)}")

    # Default channels (broadband)
    if channels is None:
        channels = ['HHZ', 'HHN', 'HHE', 'HH1', 'HH2']
        print(f"Channels: {', '.join(channels)}")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    downloaded_files = []
    total_attempts = 0
    successful_downloads = 0

    print("\n" + "="*70)
    print("DOWNLOADING DATA")
    print("="*70)

    # Download data for each date
    for date_to_download in dates_to_download:
        year = date_to_download.strftime("%Y")
        doy = date_to_download.strftime("%j")

        print(f"\nProcessing {date_to_download.date()} (DOY {doy})...")

        for station in stations:
            for channel in channels:
                # Construct S3 path
                # Format: NETWORK.STATION.LOCATION.CHANNEL.D.YEAR.DOY
                filename = f"{network}.{station}.{location}.{channel}.D.{year}.{doy}"
                s3_path = f"s3://geonet-open-data/waveforms/miniseed/{year}/{year}.{doy}/{filename}"

                # Local output path
                local_path = os.path.join(output_dir, f"{year}.{doy}")
                os.makedirs(local_path, exist_ok=True)
                local_file = os.path.join(local_path, filename)

                # Skip if already downloaded
                if os.path.exists(local_file):
                    print(f"  ⊙ {station}.{channel} - already exists")
                    downloaded_files.append(local_file)
                    continue

                # Download file
                total_attempts += 1
                cmd = [
                    'aws', 's3', 'cp',
                    '--no-sign-request',
                    s3_path,
                    local_file
                ]

                try:
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        timeout=60
                    )

                    if result.returncode == 0:
                        file_size = os.path.getsize(local_file)
                        print(f"  ✓ {station}.{channel} - {file_size/1024:.1f} KB")
                        downloaded_files.append(local_file)
                        successful_downloads += 1
                    else:
                        # File doesn't exist or other error
                        if "does not exist" in result.stderr or "404" in result.stderr:
                            print(f"  - {station}.{channel} - not available")
                        else:
                            print(f"  ✗ {station}.{channel} - error: {result.stderr.strip()[:50]}")

                except subprocess.TimeoutExpired:
                    print(f"  ✗ {station}.{channel} - timeout")
                except Exception as e:
                    print(f"  ✗ {station}.{channel} - {str(e)[:50]}")

    print("\n" + "="*70)
    print("DOWNLOAD SUMMARY")
    print("="*70)
    print(f"Total download attempts: {total_attempts}")
    print(f"Successful downloads: {successful_downloads}")
    print(f"Success rate: {successful_downloads/total_attempts*100:.1f}%")
    print(f"Total files downloaded: {len(downloaded_files)}")
    print(f"Output directory: {output_dir}")

    return downloaded_files


def merge_daily_files(output_dir, merge_output_dir="./geonet_merged"):
    """
    Merge downloaded daily files by station and channel using ObsPy

    Parameters:
    -----------
    output_dir : str
        Directory containing downloaded files
    merge_output_dir : str
        Directory to save merged files

    Returns:
    --------
    merged_files : list
        List of merged file paths
    """

    try:
        from obspy import read, Stream
    except ImportError:
        print("✗ ObsPy not installed. Cannot merge files.")
        print("  Install: pip install obspy")
        return []

    print("\n" + "="*70)
    print("MERGING FILES BY STATION")
    print("="*70)

    os.makedirs(merge_output_dir, exist_ok=True)

    # Find all downloaded files
    all_files = []
    for root, dirs, files in os.walk(output_dir):
        for f in files:
            if f.endswith('.D.2025.310') or f.endswith('.D.2025.309'):  # MiniSEED daily files
                all_files.append(os.path.join(root, f))

    if not all_files:
        print("No files to merge")
        return []

    print(f"Found {len(all_files)} files")

    # Group files by station and channel
    station_channel_files = {}
    for filepath in all_files:
        filename = os.path.basename(filepath)
        parts = filename.split('.')
        if len(parts) >= 4:
            network = parts[0]
            station = parts[1]
            location = parts[2]
            channel = parts[3]
            key = f"{network}.{station}.{location}.{channel}"

            if key not in station_channel_files:
                station_channel_files[key] = []
            station_channel_files[key].append(filepath)

    print(f"Found {len(station_channel_files)} unique station-channel combinations")

    merged_files = []

    for key, files in station_channel_files.items():
        print(f"\nMerging {key} ({len(files)} files)...")

        try:
            # Read all files
            st = Stream()
            for f in sorted(files):
                st += read(f)

            # Merge traces
            st.merge(method=1, fill_value='interpolate')

            # Save merged file
            output_file = os.path.join(merge_output_dir, f"{key}.mseed")
            st.write(output_file, format='MSEED')

            print(f"  ✓ Saved: {output_file}")
            print(f"    Traces: {len(st)}")
            print(f"    Duration: {st[0].stats.endtime - st[0].stats.starttime:.1f} seconds")

            merged_files.append(output_file)

        except Exception as e:
            print(f"  ✗ Error merging {key}: {e}")

    print(f"\n✓ Merged {len(merged_files)} station-channel files")
    print(f"Output directory: {merge_output_dir}")

    return merged_files


def main():
    """Main function"""

    parser = argparse.ArgumentParser(
        description='Download GeoNet seismic data from AWS S3',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List S3 bucket structure
  python download_geonet_aws_data.py --list

  # List data for specific date
  python download_geonet_aws_data.py --list --date 2025-11-06

  # Download data for default stations
  python download_geonet_aws_data.py --date 2025-11-06

  # Download data for specific stations
  python download_geonet_aws_data.py --date 2025-11-06 --stations WEL,BFZ,PUZ

  # Download with time range
  python download_geonet_aws_data.py --date 2025-11-06 --hours-before 2 --hours-after 12

  # Download and merge files
  python download_geonet_aws_data.py --date 2025-11-06 --merge
        """
    )

    parser.add_argument('--list', action='store_true',
                       help='List S3 bucket structure')

    parser.add_argument('--date', type=str,
                       default='2025-11-06',
                       help='Date in YYYY-MM-DD format (default: 2025-11-06)')

    parser.add_argument('--stations', type=str,
                       help='Comma-separated list of station codes (e.g., WEL,BFZ,PUZ)')

    parser.add_argument('--network', type=str, default='NZ',
                       help='Network code (default: NZ)')

    parser.add_argument('--location', type=str, default='10',
                       help='Location code (default: 10)')

    parser.add_argument('--channels', type=str,
                       help='Comma-separated channel codes (e.g., HHZ,HHN,HHE)')

    parser.add_argument('--hours-before', type=int, default=0,
                       help='Hours before date to download (default: 0)')

    parser.add_argument('--hours-after', type=int, default=24,
                       help='Hours after date to download (default: 24)')

    parser.add_argument('--output', type=str, default='./geonet_data',
                       help='Output directory (default: ./geonet_data)')

    parser.add_argument('--merge', action='store_true',
                       help='Merge downloaded files by station (requires ObsPy)')

    args = parser.parse_args()

    print("="*70)
    print("GEONET AWS S3 DATA DOWNLOADER")
    print("="*70)
    print("Data source: s3://geonet-open-data/waveforms/miniseed/")
    print("="*70)
    print()

    # Check AWS CLI
    if not check_aws_cli():
        sys.exit(1)

    # List mode
    if args.list:
        list_geonet_s3_structure(args.date if args.date else None)
        return

    # Parse stations
    stations = args.stations.split(',') if args.stations else None

    # Parse channels
    channels = args.channels.split(',') if args.channels else None

    # Download data
    downloaded_files = download_geonet_station_data(
        date_str=args.date,
        stations=stations,
        network=args.network,
        location=args.location,
        channels=channels,
        output_dir=args.output,
        hours_before=args.hours_before,
        hours_after=args.hours_after
    )

    # Merge files if requested
    if args.merge and downloaded_files:
        merge_daily_files(args.output)


if __name__ == "__main__":
    main()
