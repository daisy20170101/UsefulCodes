"""
Compare observed and synthetic seismographs.

This module provides functions to compare observed seismic data with SeisSol
synthetic simulations, including waveform alignment, cross-correlation analysis,
and three-component comparison plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.geodetics import gps2dist_azimuth
import os
from typing import Optional, Tuple
import pandas as pd


def compare_station_waveforms(
    st_data,
    st_syn,
    station_code: str,
    event_lat: float,
    event_lon: float,
    site_lat: float,
    site_lon: float,
    max_lag: int = 1000,
    align_waveforms: bool = True,
    normalize: bool = True
) -> Tuple[float, float, float, float, float]:
    """
    Compare observed and synthetic waveforms for a single station.

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream (filtered and processed)
    st_syn : obspy.Stream
        Synthetic data stream (filtered and processed)
    station_code : str
        Station code to analyze
    event_lat : float
        Event latitude (degrees)
    event_lon : float
        Event longitude (degrees)
    site_lat : float
        Station latitude (degrees)
    site_lon : float
        Station longitude (degrees)
    max_lag : int, optional
        Maximum lag for cross-correlation (samples) (default: 1000)
    align_waveforms : bool, optional
        If True, align waveforms using cross-correlation (default: True)
    normalize : bool, optional
        If True, normalize waveforms by max horizontal amplitude (default: True)

    Returns
    -------
    distance_km : float
        Epicentral distance in km
    azimuth : float
        Azimuth from event to station (degrees)
    backazimuth : float
        Back-azimuth from station to event (degrees)
    shift : float
        Time shift in samples (from cross-correlation)
    cc_value : float
        Cross-correlation coefficient

    Raises
    ------
    ValueError
        If required traces are missing from streams
    """

    # Calculate distance and azimuth
    distance, azimuth, backazimuth = gps2dist_azimuth(
        event_lat, event_lon, site_lat, site_lon
    )
    distance_km = distance / 1e3

    # Get normalization factor from max horizontal amplitude
    if normalize:
        try:
            norm12 = np.max(np.abs(st_data.select(station=station_code, channel='??E')[0].data))
            norm13 = np.max(np.abs(st_data.select(station=station_code, channel='??N')[0].data))
            norm1 = np.max([norm12, norm13])

            norm22 = np.max(np.abs(st_syn.select(station=station_code, channel='??E')[0].data))
            norm23 = np.max(np.abs(st_syn.select(station=station_code, channel='??N')[0].data))
            norm2 = np.max([norm22, norm23])
        except (IndexError, AttributeError):
            raise ValueError(f"Missing horizontal components for station {station_code}")
    else:
        norm1 = 1.0
        norm2 = 1.0

    # Compute cross-correlation on E component
    if align_waveforms:
        try:
            data_e = st_data.select(station=station_code, channel='H?E')[0].data / norm1
            syn_e = st_syn.select(station=station_code, channel='?NE')[0].data / norm2

            cc = correlate(data_e, syn_e, max_lag)
            shift, cc_value = xcorr_max(cc, abs_max=False)
        except (IndexError, AttributeError):
            raise ValueError(f"Missing E component for station {station_code}")
    else:
        shift = 0
        cc_value = 0.0

    return distance_km, azimuth, backazimuth, shift, cc_value


def plot_station_comparison(
    st_data,
    st_syn,
    station_code: str,
    distance_km: float,
    backazimuth: float,
    shift: float,
    cc_value: float,
    output_path: Optional[str] = None,
    time_window: float = 30.0,
    normalize: bool = True,
    show_plot: bool = False
) -> plt.Figure:
    """
    Create three-component waveform comparison plot.

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream
    st_syn : obspy.Stream
        Synthetic data stream
    station_code : str
        Station code
    distance_km : float
        Epicentral distance in km
    backazimuth : float
        Back-azimuth in degrees
    shift : float
        Time shift in samples for alignment
    cc_value : float
        Cross-correlation coefficient
    output_path : str, optional
        Path to save figure. If None, figure not saved.
    time_window : float, optional
        Time window to display in seconds (default: 30.0)
    normalize : bool, optional
        Normalize by max horizontal amplitude (default: True)
    show_plot : bool, optional
        Show plot interactively (default: False)

    Returns
    -------
    matplotlib.figure.Figure
        Figure object

    Notes
    -----
    - Observed data plotted in black
    - Synthetic data plotted in red, offset by -2 for clarity
    - Synthetic data is shifted by cross-correlation lag
    """

    # Get normalization factors
    if normalize:
        norm12 = np.max(np.abs(st_data.select(station=station_code, channel='??E')[0].data))
        norm13 = np.max(np.abs(st_data.select(station=station_code, channel='??N')[0].data))
        norm1 = np.max([norm12, norm13])

        norm22 = np.max(np.abs(st_syn.select(station=station_code, channel='??E')[0].data))
        norm23 = np.max(np.abs(st_syn.select(station=station_code, channel='??N')[0].data))
        norm2 = np.max([norm22, norm23])
    else:
        norm1 = 1.0
        norm2 = 1.0

    # Create figure
    fig = plt.figure(figsize=(12, 1.5))

    # E component
    ax3 = plt.subplot(131)
    tr_data_e = st_data.select(station=station_code, channel='H?E')[0]
    tr_syn_e = st_syn.select(station=station_code, channel='?NE')[0]

    ax3.plot(tr_data_e.times(), tr_data_e.data / norm1, 'k', label='E')
    ax3.plot(tr_syn_e.times(),
             -2 + np.roll(tr_syn_e.data, int(shift)) / norm2, 'r')
    ax3.set_title(station_code)
    ax3.set_xlim([0, time_window])
    ax3.set_axis_off()
    ax3.set_yticklabels([])

    # N component
    ax2 = plt.subplot(132, sharey=ax3, sharex=ax3)
    tr_data_n = st_data.select(station=station_code, channel='H?N')[0]
    tr_syn_n = st_syn.select(station=station_code, channel='?NN')[0]

    ax2.plot(tr_data_n.times(), tr_data_n.data / norm1, 'k', label='N')
    ax2.plot(tr_syn_n.times(),
             -2 + np.roll(tr_syn_n.data, int(shift)) / norm2, 'r')
    ax2.set_yticklabels([])
    ax2.set_axis_off()
    ax2.set_title(f'D: {distance_km:.1f} km   BAZ: {backazimuth:.1f}°')

    # Z component
    ax1 = plt.subplot(133, sharey=ax3, sharex=ax3)
    tr_data_z = st_data.select(station=station_code, channel='H?Z')[0]
    tr_syn_z = st_syn.select(station=station_code, channel='?NZ')[0]

    ax1.plot(tr_data_z.times(), tr_data_z.data / norm1, 'k', label='Z')
    ax1.plot(tr_syn_z.times(),
             -2 + np.roll(tr_syn_z.data, int(shift)) / norm2, 'r')
    ax1.set_title(f'CC = {cc_value:.3f}')
    ax1.set_yticklabels([])
    ax1.set_axis_off()

    # Save or show
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=300, transparent=False, bbox_inches='tight')

    if show_plot:
        plt.show()
    else:
        plt.close()

    return fig


def compare_all_seismographs(
    st_data,
    st_syn,
    inventory,
    siteTable: pd.DataFrame,
    event,
    eventname: str,
    output_folder: str,
    time_window: float = 30.0,
    taper_fraction: float = 0.5,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Compare all seismographs between observed and synthetic data.

    This function processes all stations in the inventory, computes comparison
    metrics, and creates three-component waveform plots.

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream (already filtered)
    st_syn : obspy.Stream
        Synthetic data stream (already filtered)
    inventory : obspy.Inventory
        Station inventory with network and station information
    siteTable : pd.DataFrame
        Site information table with columns: 'sta', 'lat', 'lon'
    event : obspy.Event
        Event object containing origin information
    eventname : str
        Event identifier for file naming (e.g., '20240805')
    output_folder : str
        Directory to save comparison plots
    time_window : float, optional
        Time window for trimming and plotting in seconds (default: 30.0)
    taper_fraction : float, optional
        Taper fraction for data preprocessing (default: 0.5)
    verbose : bool, optional
        Print progress information (default: True)

    Returns
    -------
    pd.DataFrame
        Comparison results with columns:
        - station: station code
        - distance_km: epicentral distance
        - azimuth: azimuth angle
        - backazimuth: back-azimuth angle
        - shift_samples: time shift from cross-correlation
        - cc_value: cross-correlation coefficient
        - norm_data: normalization factor for observed data
        - norm_syn: normalization factor for synthetic data

    Examples
    --------
    >>> results = compare_all_seismographs(
    ...     st_data=st_data1,
    ...     st_syn=st_syn1,
    ...     inventory=inv2,
    ...     siteTable=siteTable,
    ...     event=cat[0],
    ...     eventname='20240805',
    ...     output_folder='/path/to/plots/m20240805/'
    ... )
    >>> print(f"Processed {len(results)} stations")
    >>> print(f"Mean CC: {results['cc_value'].mean():.3f}")
    """

    # Get event coordinates
    event_lat = event.origins[0].latitude
    event_lon = event.origins[0].longitude

    # Trim streams to same length
    st_data_proc = st_data.copy()
    st_syn_proc = st_syn.copy()

    t_init_data = st_data_proc[0].stats.starttime
    t_init_syn = st_syn_proc[0].stats.starttime

    st_data_proc.trim(t_init_data, t_init_data + time_window, pad=True, fill_value=0.0)
    st_syn_proc.trim(t_init_syn, t_init_syn + time_window, pad=True, fill_value=0.0)

    if verbose:
        print(f"Data start time: {t_init_data}")
        print(f"Synthetic start time: {t_init_syn}")
        print(f"Time window: {time_window} s")

    # Initialize results list
    results = []

    # Create output folder
    outfolder = os.path.join(output_folder, f'm{eventname}')
    os.makedirs(outfolder, exist_ok=True)

    # Loop through all networks and stations
    station_count = 0
    success_count = 0

    for network in inventory:
        for sta in network:
            station_count += 1
            station_code = sta.code

            if verbose:
                print(f"Processing: {station_code}")

            try:
                # Apply taper to data
                st_data_proc.select(station=station_code, channel='??E').taper(
                    max_percentage=taper_fraction
                )
                st_data_proc.select(station=station_code, channel='??N').taper(
                    max_percentage=taper_fraction
                )
                st_data_proc.select(station=station_code, channel='??Z').taper(
                    max_percentage=taper_fraction
                )

                # Get station coordinates from siteTable
                site_info = siteTable[siteTable['sta'] == station_code]
                if len(site_info) == 0:
                    if verbose:
                        print(f"  Station {station_code} not found in siteTable")
                    continue

                sta_lat = site_info.lat.values[0]
                sta_lon = site_info.lon.values[0]

                # Compute comparison metrics
                distance_km, azimuth, backazimuth, shift, cc_value = compare_station_waveforms(
                    st_data=st_data_proc,
                    st_syn=st_syn_proc,
                    station_code=station_code,
                    event_lat=event_lat,
                    event_lon=event_lon,
                    site_lat=sta_lat,
                    site_lon=sta_lon
                )

                if verbose:
                    print(f"  Distance: {distance_km:.1f} km, "
                          f"BAZ: {backazimuth:.1f}°, CC: {cc_value:.3f}")

                # Get normalization factors for storage
                norm_data = np.max([
                    np.max(np.abs(st_data_proc.select(station=station_code, channel='??E')[0].data)),
                    np.max(np.abs(st_data_proc.select(station=station_code, channel='??N')[0].data))
                ])
                norm_syn = np.max([
                    np.max(np.abs(st_syn_proc.select(station=station_code, channel='??E')[0].data)),
                    np.max(np.abs(st_syn_proc.select(station=station_code, channel='??N')[0].data))
                ])

                # Create comparison plot
                output_path = os.path.join(outfolder, f'ACC2-{station_code}.png')
                plot_station_comparison(
                    st_data=st_data_proc,
                    st_syn=st_syn_proc,
                    station_code=station_code,
                    distance_km=distance_km,
                    backazimuth=backazimuth,
                    shift=shift,
                    cc_value=cc_value,
                    output_path=output_path,
                    time_window=time_window
                )

                # Store results
                results.append({
                    'station': station_code,
                    'distance_km': distance_km,
                    'azimuth': azimuth,
                    'backazimuth': backazimuth,
                    'shift_samples': shift,
                    'cc_value': cc_value,
                    'norm_data': norm_data,
                    'norm_syn': norm_syn
                })

                success_count += 1

            except (IndexError, ValueError, KeyError) as e:
                if verbose:
                    print(f"  Skipped {station_code}: {str(e)}")
                continue

    # Create results DataFrame
    results_df = pd.DataFrame(results)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Comparison Summary:")
        print(f"{'='*60}")
        print(f"Total stations: {station_count}")
        print(f"Successful: {success_count}")
        print(f"Failed: {station_count - success_count}")
        if len(results_df) > 0:
            print(f"Mean CC: {results_df['cc_value'].mean():.3f}")
            print(f"Median CC: {results_df['cc_value'].median():.3f}")
            print(f"Mean shift: {results_df['shift_samples'].mean():.1f} samples")
        print(f"Output folder: {outfolder}")

    return results_df


# Example usage
if __name__ == "__main__":
    print("Module loaded successfully")
    print("Use compare_all_seismographs() to process all stations")
    print("Or use compare_station_waveforms() and plot_station_comparison() for individual stations")
