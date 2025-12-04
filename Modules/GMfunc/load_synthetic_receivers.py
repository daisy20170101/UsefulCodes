"""
Load SeisSol synthetic receiver data for comparison with observed records.

This module provides functions to load SeisSol receiver output files,
process them into ObsPy streams, and prepare them for comparison with
observed seismograms.
"""

import os
import numpy as np
from obspy import Stream, Trace
from obspy.core import UTCDateTime
import pandas as pd
from typing import Optional, Tuple, List


def load_seissol_receivers_for_comparison(
    folder: str,
    siteTable: pd.DataFrame,
    event_origin_time: UTCDateTime,
    num_receivers: int = 60,
    ista_start: int = 1,
    receiver_prefix: str = 'stf1-receiver',
    file_suffix: str = '-00000.dat',
    sample_rate: float = 200.0,
    velocity_cols: Tuple[int, int, int] = (7, 8, 9),
    time_col: int = 0,
    sta_col: str = 'sta',
    lon_col: str = 'lon',
    lat_col: str = 'lat',
    differentiate: bool = True,
    filter_freqs: Optional[Tuple[float, float]] = (0.01, 2.0),
    filter_corners: int = 4,
    zerophase: bool = True
) -> Stream:
    """
    Load SeisSol receiver data and prepare for comparison with observations.

    Reads velocity data from SeisSol receiver output files, converts to ObsPy
    stream format, optionally differentiates to acceleration, and applies
    bandpass filter.

    Parameters
    ----------
    folder : str
        Path to folder containing SeisSol receiver files
    siteTable : pd.DataFrame
        Site information table with columns for station codes and coordinates
        Required columns: sta (station code), lon, lat
    event_origin_time : obspy.UTCDateTime
        Event origin time for trace starttime
    num_receivers : int, optional
        Number of receiver files to load (default: 60)
        Files numbered from 1 to num_receivers
    receiver_prefix : str, optional
        Prefix for receiver files (default: 'stf1-receiver')
    file_suffix : str, optional
        Suffix for receiver files (default: '-00000.dat')
    sample_rate : float, optional
        Sampling rate in Hz (default: 200.0, corresponding to 0.005s dt)
    velocity_cols : tuple of int, optional
        Column indices for (vx, vy, vz) velocities (default: (7, 8, 9))
    time_col : int, optional
        Column index for time (default: 0)
    sta_col : str, optional
        Column name in siteTable for station codes (default: 'sta')
    lon_col : str, optional
        Column name in siteTable for longitude (default: 'lon')
    lat_col : str, optional
        Column name in siteTable for latitude (default: 'lat')
    differentiate : bool, optional
        Convert velocity to acceleration (default: True)
    filter_freqs : tuple of float, optional
        (freqmin, freqmax) for bandpass filter in Hz (default: (0.01, 2.0))
        Set to None to skip filtering
    filter_corners : int, optional
        Number of filter corners (default: 4)
    zerophase : bool, optional
        Use zero-phase filter (default: True)

    Returns
    -------
    obspy.Stream
        Processed stream with acceleration (if differentiate=True) or velocity
        Channel naming: HNE (East), HNN (North), HNZ (Up)
        Network: 'SeisSol'

    Examples
    --------
    >>> import pandas as pd
    >>> from obspy.core import UTCDateTime
    >>>
    >>> # Load site table
    >>> siteTable = pd.read_csv('site_info.csv')
    >>>
    >>> # Define event origin time
    >>> origin_time = UTCDateTime('2024-08-05T10:30:00')
    >>>
    >>> # Load SeisSol receivers
    >>> st_syn = load_seissol_receivers_for_comparison(
    ...     folder='/path/to/seissol/output/',
    ...     siteTable=siteTable,
    ...     event_origin_time=origin_time,
    ...     num_receivers=60,
    ...     differentiate=True,
    ...     filter_freqs=(0.01, 2.0)
    ... )
    >>>
    >>> print(f"Loaded {len(st_syn)} traces from {len(st_syn) // 3} stations")

    Notes
    -----
    File naming convention:
    - Files < 100: {prefix}-000{ista}{suffix}
    - Files >= 100: {prefix}-00{ista}{suffix}

    Processing steps:
    1. Load velocity data from .dat files
    2. Create ObsPy traces with proper metadata
    3. Differentiate to acceleration (if requested)
    4. Apply bandpass filter (if requested)

    The function assumes SeisSol output format with columns:
    - time (usually column 0)
    - velocities vx, vy, vz (usually columns 7, 8, 9)
    """
    # Channel names
    channels = ['HNE', 'HNN', 'HNZ']  # East, North, Up

    # Initialize empty stream
    st_syn = Stream()

    # Track loaded stations
    loaded_count = 0
    failed_files = []

    # Loop through receiver indices
    for ista in range(ista_start, ista_start + num_receivers ):
        # Construct filename based on receiver number
        if ista < 100:
            filename = os.path.join(
                folder,
                f'{receiver_prefix}-000{ista}{file_suffix}'
            )
        else:
            filename = os.path.join(
                folder,
                f'{receiver_prefix}-00{ista}{file_suffix}'
            )

        # Check if file exists
        if not os.path.isfile(filename):
            failed_files.append(filename)
            continue

        try:
            # Load data (skip header lines, extract time and velocities)
            data = np.loadtxt(
                filename,
                comments='#',
                skiprows=2,
                usecols=(time_col, velocity_cols[0], velocity_cols[1], velocity_cols[2]),
                unpack=True
            )

            t = data[0]
            vx = data[1]
            vy = data[2]
            vz = data[3]

            # Create stream with 3 components
            st = Stream([
                Trace(vx),
                Trace(vy),
                Trace(vz)
            ])

            # Get station info from siteTable
            # Receiver index starts at 1, but DataFrame index starts at 0
            site_idx = ista - ista_start  - 1

            if site_idx >= len(siteTable):
                print(f"Warning: Receiver {ista} exceeds siteTable length, skipping")
                continue

            # Set metadata for each component
            for tr_num in range(3):
                st[tr_num].stats['sampling_rate'] = sample_rate
                st[tr_num].stats['network'] = 'SeisSol'
                st[tr_num].stats['station'] = str(siteTable[sta_col].iloc[site_idx])
                st[tr_num].stats['channel'] = channels[tr_num]
                st[tr_num].stats['longitude'] = float(siteTable[lon_col].iloc[site_idx])
                st[tr_num].stats['latitude'] = float(siteTable[lat_col].iloc[site_idx])
                st[tr_num].stats['starttime'] = event_origin_time

            # Add to combined stream
            st_syn += st
            loaded_count += 1

        except Exception as e:
            print(f"Error loading {filename}: {str(e)}")
            failed_files.append(filename)
            continue

    print(f"Loaded {loaded_count} receivers successfully")
    if failed_files:
        print(f"Failed to load {len(failed_files)} files")

    # Process stream
    if len(st_syn) > 0:
        # Differentiate velocity to acceleration if requested
        if differentiate:
            st_syn.differentiate()

        # Apply bandpass filter if requested
        if filter_freqs is not None:
            freqmin, freqmax = filter_freqs
            st_syn.filter(
                'bandpass',
                freqmin=freqmin,
                freqmax=freqmax,
                corners=filter_corners,
                zerophase=zerophase
            )

    return st_syn


def load_synthetic_for_event(
    folder: str,
    siteTable_path: str,
    catalog_or_origin_time,
    modelname: str = 'm20240805',
    num_receivers: int = 60,
    filter_freqs: Tuple[float, float] = (0.01, 2.0),
    differentiate: bool = True
) -> Tuple[Stream, str]:
    """
    Simplified wrapper to load synthetic data for an event.

    Parameters
    ----------
    folder : str
        Path to SeisSol output folder
    siteTable_path : str
        Path to site table CSV file
    catalog_or_origin_time : obspy.Catalog or obspy.UTCDateTime
        Event catalog or origin time
    modelname : str, optional
        Model name for reference (default: 'm20240805')
    num_receivers : int, optional
        Number of receivers (default: 60)
    filter_freqs : tuple of float, optional
        (freqmin, freqmax) for filtering (default: (0.01, 2.0))
    differentiate : bool, optional
        Convert velocity to acceleration (default: True)

    Returns
    -------
    st_syn : obspy.Stream
        Processed synthetic stream
    modelname : str
        Model name (returned for convenience)

    Examples
    --------
    >>> # Using with catalog
    >>> from obspy import read_events
    >>> cat = read_events('event.xml')
    >>>
    >>> st_syn, model = load_synthetic_for_event(
    ...     folder='/path/to/output/',
    ...     siteTable_path='site_table.csv',
    ...     catalog_or_origin_time=cat,
    ...     modelname='m20240805'
    ... )

    >>> # Using with direct origin time
    >>> from obspy.core import UTCDateTime
    >>> origin_time = UTCDateTime('2024-08-05T10:30:00')
    >>>
    >>> st_syn, model = load_synthetic_for_event(
    ...     folder='/path/to/output/',
    ...     siteTable_path='site_table.csv',
    ...     catalog_or_origin_time=origin_time
    ... )
    """
    # Load site table
    siteTable = pd.read_csv(siteTable_path)

    # Extract origin time
    if hasattr(catalog_or_origin_time, 'events'):
        # It's a catalog
        origin_time = catalog_or_origin_time[0].origins[0].time
    else:
        # It's already a UTCDateTime
        origin_time = catalog_or_origin_time

    # Load receivers
    st_syn = load_seissol_receivers_for_comparison(
        folder=folder,
        siteTable=siteTable,
        event_origin_time=origin_time,
        num_receivers=num_receivers,
        differentiate=differentiate,
        filter_freqs=filter_freqs
    )

    return st_syn, modelname


# Example usage
if __name__ == "__main__":
    print("SeisSol Synthetic Receiver Loading Functions")
    print("=" * 70)
    print("\nAvailable functions:")
    print("  1. load_seissol_receivers_for_comparison() - Main function")
    print("  2. load_synthetic_for_event() - Simplified wrapper")
    print("\nImport with:")
    print("  from GMfunc.load_synthetic_receivers import load_seissol_receivers_for_comparison")
