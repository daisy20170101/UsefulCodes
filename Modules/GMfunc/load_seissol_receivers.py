"""
Load and process SeisSol receiver data for comparison with observed seismic records.

This module provides functions to load SeisSol synthetic receiver outputs,
convert them to ObsPy Stream objects, and prepare them for comparison with
real seismic data.
"""

import os
import numpy as np
from obspy import Stream, Trace
from typing import Optional, Tuple
import pandas as pd


def load_seissol_receivers(
    folder: str,
    siteTable: pd.DataFrame,
    event_origin_time,
    receiver_range: Tuple[int, int] = (181, 241),
    receiver_offset: int = 180,
    sample_rate: float = 200.0,
    channels: list = None,
    filter_params: Optional[dict] = None,
    differentiate: bool = True
) -> Stream:
    """
    Load SeisSol synthetic receiver data and convert to ObsPy Stream format.

    This function reads SeisSol receiver output files (.dat format), extracts velocity
    time series, converts them to ObsPy traces with proper metadata, optionally
    differentiates to acceleration, and applies bandpass filtering.

    Parameters
    ----------
    folder : str
        Directory containing SeisSol receiver output files
        Format expected: 'm7liu-receiver-XXXXX-00000.dat'
    siteTable : pd.DataFrame
        Site information table with columns: 'sta', 'lon', 'lat'
        Must have at least (receiver_range[1] - receiver_offset) rows
    event_origin_time : obspy.UTCDateTime
        Event origin time from earthquake catalog (e.g., cat[0].origins[0].time)
    receiver_range : tuple of int, optional
        (start, end) receiver numbers to load (default: (181, 241))
        End is exclusive (Python range convention)
    receiver_offset : int, optional
        Offset to apply when indexing siteTable (default: 180)
        siteTable_index = receiver_number - 1 - receiver_offset
    sample_rate : float, optional
        Data sample rate in Hz (default: 200.0)
        SeisSol default: 1/0.005 = 200 Hz
    channels : list of str, optional
        Channel codes for [E, N, Z] components (default: ['HNE','HNN','HNZ'])
    filter_params : dict, optional
        Bandpass filter parameters with keys:
        - 'freqmin': minimum frequency (Hz)
        - 'freqmax': maximum frequency (Hz)
        - 'corners': filter corners (default: 4)
        - 'zerophase': zero-phase filter (default: True)
        If None, no filtering is applied
    differentiate : bool, optional
        If True, differentiate velocity to acceleration (default: True)

    Returns
    -------
    obspy.Stream
        Stream object containing synthetic receiver data with metadata:
        - network: 'SeisSol'
        - station: from siteTable
        - channel: HNE/HNN/HNZ (or custom channels)
        - location: coordinates from siteTable
        - starttime: event origin time
        - sampling_rate: as specified

    Examples
    --------
    >>> import pandas as pd
    >>> from obspy import UTCDateTime
    >>>
    >>> # Setup
    >>> siteTable = pd.read_csv('sites.csv')  # Contains 'sta', 'lon', 'lat'
    >>> origin_time = UTCDateTime('2024-08-05T11:38:37.569936Z')
    >>> folder = '/path/to/seissol/output/'
    >>>
    >>> # Load with default settings (acceleration, 0.01-2 Hz)
    >>> st = load_seissol_receivers(
    ...     folder=folder,
    ...     siteTable=siteTable,
    ...     event_origin_time=origin_time,
    ...     filter_params={'freqmin': 0.01, 'freqmax': 2.0}
    ... )
    >>>
    >>> # Load velocity only, no filtering
    >>> st_vel = load_seissol_receivers(
    ...     folder=folder,
    ...     siteTable=siteTable,
    ...     event_origin_time=origin_time,
    ...     differentiate=False,
    ...     filter_params=None
    ... )

    Notes
    -----
    - SeisSol receiver files format:
      Column 0: time
      Columns 7,8,9: vx, vy, vz (velocity components in m/s)
    - File naming: 'm7liu-receiver-XXXXX-00000.dat' where XXXXX is zero-padded receiver number
    - Missing files are silently skipped
    - Velocity units: m/s
    - Acceleration units (after differentiate): m/s²
    """

    # Set default channels if not provided
    if channels is None:
        channels = ['HNE', 'HNN', 'HNZ']

    # Initialize empty stream
    st_syn = Stream()

    # Loop through receiver range
    for ista in receiver_range[0] + np.arange(receiver_range[1] - receiver_range[0]):

        # Construct filename with zero-padding
        if ista < 100:
            filename = os.path.join(folder, f'm7liu-receiver-000{ista}-00000.dat')
        else:
            filename = os.path.join(folder, f'm7liu-receiver-00{ista}-00000.dat')

        # Check if file exists
        if not os.path.isfile(filename):
            continue

        # Load SeisSol receiver data
        # Columns: 0=time, 7=vx, 8=vy, 9=vz
        t, vx, vy, vz = np.loadtxt(
            filename,
            comments='#',
            skiprows=2,
            usecols=(0, 7, 8, 9),
            unpack=True
        )

        # Create ObsPy traces for each component
        st = Stream([Trace(vx), Trace(vy), Trace(vz)])

        # Calculate siteTable index
        site_idx = ista - 1 - receiver_offset

        # Add metadata to each trace
        for tr_num in range(3):
            st[tr_num].stats['sampling_rate'] = sample_rate
            st[tr_num].stats['network'] = 'SeisSol'
            st[tr_num].stats['station'] = siteTable.sta.iloc[site_idx]
            st[tr_num].stats['channel'] = channels[tr_num]
            st[tr_num].stats['longitude'] = siteTable.lon.iloc[site_idx]
            st[tr_num].stats['latitude'] = siteTable.lat.iloc[site_idx]
            st[tr_num].stats['starttime'] = event_origin_time

        # Add to main stream
        st_syn += st

    # Process the stream
    st_processed = st_syn.copy()

    # Differentiate velocity to acceleration if requested
    if differentiate:
        st_processed.differentiate()

    # Apply bandpass filter if parameters provided
    if filter_params is not None:
        freqmin = filter_params.get('freqmin', 0.01)
        freqmax = filter_params.get('freqmax', 2.0)
        corners = filter_params.get('corners', 4)
        zerophase = filter_params.get('zerophase', True)

        st_processed.filter(
            'bandpass',
            freqmin=freqmin,
            freqmax=freqmax,
            corners=corners,
            zerophase=zerophase
        )

    return st_processed


def load_seissol_receivers_simple(
    folder: str,
    modelname: str,
    siteTable: pd.DataFrame,
    event_origin_time,
    freqmin: float = 0.01,
    freqmax: float = 2.0
) -> Stream:
    """
    Simplified interface to load SeisSol receivers with common defaults.

    This is a convenience wrapper around load_seissol_receivers() with
    commonly used parameters pre-set.

    Parameters
    ----------
    folder : str
        Directory containing SeisSol receiver files
    modelname : str
        Model name (for reference/logging only)
    siteTable : pd.DataFrame
        Site table with 'sta', 'lon', 'lat' columns
    event_origin_time : obspy.UTCDateTime
        Event origin time
    freqmin : float, optional
        Minimum frequency for bandpass filter in Hz (default: 0.01)
    freqmax : float, optional
        Maximum frequency for bandpass filter in Hz (default: 2.0)

    Returns
    -------
    obspy.Stream
        Processed stream with acceleration data (differentiated and filtered)

    Examples
    --------
    >>> st = load_seissol_receivers_simple(
    ...     folder='/path/to/data/',
    ...     modelname='m20240805',
    ...     siteTable=siteTable,
    ...     event_origin_time=cat[0].origins[0].time,
    ...     freqmin=0.01,
    ...     freqmax=2.0
    ... )
    """

    print(f"Loading SeisSol receivers for model: {modelname}")
    print(f"Folder: {folder}")
    print(f"Filter: {freqmin}-{freqmax} Hz")

    st = load_seissol_receivers(
        folder=folder,
        siteTable=siteTable,
        event_origin_time=event_origin_time,
        receiver_range=(181, 241),
        receiver_offset=180,
        sample_rate=200.0,
        filter_params={
            'freqmin': freqmin,
            'freqmax': freqmax,
            'corners': 4,
            'zerophase': True
        },
        differentiate=True
    )

    print(f"Loaded {len(st)} traces from {len(st) // 3} stations")

    return st


# Example usage
if __name__ == "__main__":
    import pandas as pd
    from obspy import UTCDateTime

    # Example setup (adjust paths as needed)
    folder = '/Users/DuoL/Documents/NSHM/Model_kinematic/mesh4/m20240805-att/'
    modelname = 'm20240805'

    # Load site table
    # siteTable = pd.read_csv('sites.csv')

    # Create dummy event time for testing
    event_time = UTCDateTime('2024-08-05T11:38:37.569936Z')

    # Example usage would be:
    # st = load_seissol_receivers_simple(
    #     folder=folder,
    #     modelname=modelname,
    #     siteTable=siteTable,
    #     event_origin_time=event_time
    # )

    print("Module loaded successfully")
    print("Use load_seissol_receivers() or load_seissol_receivers_simple()")
