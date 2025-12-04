"""
Compare observed and synthetic waveforms with cross-correlation alignment.

This module provides functions to compare seismic waveforms by computing
cross-correlation to find optimal time shifts, then plotting aligned comparisons.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from obspy import Stream
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.geodetics import gps2dist_azimuth
from typing import Optional, Tuple, Dict


def compare_waveforms_all_stations(
    st_data: Stream,
    st_syn: Stream,
    inventory,
    siteTable: pd.DataFrame,
    event,
    eventname: str,
    output_folder: str,
    trim_duration: float = 30.0,
    taper_fraction: float = 0.5,
    xcorr_shift: int = 1000,
    vertical_offset: float = -2.0,
    figsize: Tuple[float, float] = (12, 1.5),
    dpi: int = 300,
    file_prefix: str = 'ACC2-',
    verbose: bool = True,
    sta_col: str = 'sta',
    lat_col: str = 'lat',
    lon_col: str = 'lon'
) -> pd.DataFrame:
    """
    Compare waveforms for all stations with cross-correlation alignment.

    Creates 3-panel plots (E, N, Z components) for each station showing:
    - Observed waveform (black)
    - Synthetic waveform shifted by cross-correlation lag (red, offset by -2)
    - Distance, back-azimuth, and cross-correlation coefficient

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream (already filtered)
    st_syn : obspy.Stream
        Synthetic data stream (already filtered)
    inventory : obspy.Inventory
        Station inventory with network and station information
    siteTable : pd.DataFrame
        Site table with station coordinates
        Required columns: sta (station code), lat, lon
    event : obspy.Event
        Event object with origin information
    eventname : str
        Event name for output file naming
    output_folder : str
        Folder to save plots (will be created if doesn't exist)
    trim_duration : float, optional
        Duration in seconds to trim waveforms (default: 30.0)
    taper_fraction : float, optional
        Taper fraction for cosine taper (default: 0.5)
    xcorr_shift : int, optional
        Maximum lag in samples for cross-correlation (default: 1000)
    vertical_offset : float, optional
        Vertical offset for synthetic trace (default: -2.0)
    figsize : tuple, optional
        Figure size (width, height) in inches (default: (12, 1.5))
    dpi : int, optional
        Resolution for saved figures (default: 300)
    file_prefix : str, optional
        Prefix for output filenames (default: 'ACC2-')
    verbose : bool, optional
        Print processing information (default: True)
    sta_col : str, optional
        Column name for station codes in siteTable (default: 'sta')
    lat_col : str, optional
        Column name for latitude in siteTable (default: 'lat')
    lon_col : str, optional
        Column name for longitude in siteTable (default: 'lon')

    Returns
    -------
    pd.DataFrame
        Results table with columns:
        - station: station code
        - distance_km: epicentral distance
        - azimuth: azimuth from event to station
        - backazimuth: back-azimuth from station to event
        - xcorr_coef: cross-correlation coefficient
        - xcorr_shift: time shift in samples
        - max_amp_obs: max horizontal amplitude in observed data
        - max_amp_syn: max horizontal amplitude in synthetic data
        - plot_file: path to saved plot

    Examples
    --------
    >>> # Load and filter data
    >>> st_data = read('observed.mseed')
    >>> st_syn = read('synthetic.mseed')
    >>> st_data.filter('bandpass', freqmin=0.01, freqmax=2.0)
    >>> st_syn.filter('bandpass', freqmin=0.01, freqmax=2.0)
    >>>
    >>> # Compare all stations
    >>> results = compare_waveforms_all_stations(
    ...     st_data=st_data,
    ...     st_syn=st_syn,
    ...     inventory=inv,
    ...     siteTable=siteTable,
    ...     event=cat[0],
    ...     eventname='20240805',
    ...     output_folder='./plots/'
    ... )
    >>>
    >>> # View results
    >>> print(results[['station', 'distance_km', 'xcorr_coef']])

    Notes
    -----
    - Streams are trimmed to trim_duration seconds
    - Cross-correlation is computed on the East component
    - The same time shift is applied to all three components
    - Normalization is by maximum horizontal amplitude (max of E and N)
    - Failed stations (missing data, errors) are skipped with warning
    """
    # Make copies to avoid modifying originals
    st_data1 = st_data.copy()
    st_syn1 = st_syn.copy()

    # Get event location
    evelat = event.origins[0].latitude
    evelon = event.origins[0].longitude

    # Trim streams to same duration
    t_init1 = st_data1[0].stats.starttime
    t_init2 = st_syn1[0].stats.starttime

    st_data1.trim(t_init1, t_init1 + trim_duration, pad=True, fill_value=0.0)
    st_syn1.trim(t_init2, t_init2 + trim_duration, pad=True, fill_value=0.0)

    if verbose:
        print(f"Trimmed streams to {trim_duration} seconds")
        print(f"Observed start: {t_init1}")
        print(f"Synthetic start: {t_init2}")

    # Create output folder
    os.makedirs(output_folder, exist_ok=True)

    # Results storage
    results = []

    # Process each station
    station_count = 0
    failed_count = 0

    for network in inventory:
        for sta in network:
            station_code = sta.code

            if verbose:
                print(f"\nProcessing station: {station_code}")

            try:
                # Apply taper to all components
                for comp in ['??E', '??N', '??Z']:
                    st_data1.select(station=station_code, channel=comp).taper(
                        max_percentage=taper_fraction
                    )

                # Get station coordinates from siteTable
                sta_info = siteTable[siteTable[sta_col] == station_code]
                if len(sta_info) == 0:
                    if verbose:
                        print(f"  Warning: {station_code} not found in siteTable, skipping")
                    failed_count += 1
                    continue

                sta_lat = sta_info[lat_col].to_numpy()[0]
                sta_lon = sta_info[lon_col].to_numpy()[0]

                # Compute distance and azimuth
                distance_m, azimuth, backazimuth = gps2dist_azimuth(
                    evelat, evelon, sta_lat, sta_lon
                )
                distance_km = distance_m / 1e3

                if verbose:
                    print(f"  Distance: {distance_km:.1f} km, BAZ: {backazimuth:.1f}°")

                # Get horizontal components for normalization
                tr_data_e = st_data1.select(station=station_code, channel='??E')
                tr_data_n = st_data1.select(station=station_code, channel='??N')
                tr_syn_e = st_syn1.select(station=station_code, channel='??E')
                tr_syn_n = st_syn1.select(station=station_code, channel='??N')

                if len(tr_data_e) == 0 or len(tr_data_n) == 0:
                    if verbose:
                        print(f"  Warning: Missing observed horizontal components, skipping")
                    failed_count += 1
                    continue

                if len(tr_syn_e) == 0 or len(tr_syn_n) == 0:
                    if verbose:
                        print(f"  Warning: Missing synthetic horizontal components, skipping")
                    failed_count += 1
                    continue

                # Normalize by max horizontal amplitude
                norm12 = np.max(np.abs(tr_data_e[0].data))
                norm13 = np.max(np.abs(tr_data_n[0].data))
                norm1 = np.max([norm12, norm13])  # Observed normalization

                norm22 = np.max(np.abs(tr_syn_e[0].data))
                norm23 = np.max(np.abs(tr_syn_n[0].data))
                norm2 = np.max([norm22, norm23])  # Synthetic normalization

                if verbose:
                    print(f"  Max amplitude - Obs: {norm1:.2e}, Syn: {norm2:.2e}")

                # Compute cross-correlation on E component
                # Note: Using different channel patterns for observed vs synthetic
                tr_obs_xcorr = st_data1.select(station=station_code, channel='H?E')
                tr_syn_xcorr = st_syn1.select(station=station_code, channel='?NE')

                if len(tr_obs_xcorr) == 0:
                    # Try alternative pattern
                    tr_obs_xcorr = st_data1.select(station=station_code, channel='??E')

                if len(tr_syn_xcorr) == 0:
                    # Try alternative pattern
                    tr_syn_xcorr = st_syn1.select(station=station_code, channel='??E')

                if len(tr_obs_xcorr) == 0 or len(tr_syn_xcorr) == 0:
                    if verbose:
                        print(f"  Warning: Cannot get E component for xcorr, skipping")
                    failed_count += 1
                    continue

                # Compute cross-correlation
                cc = correlate(
                    tr_obs_xcorr[0].data / norm1,
                    tr_syn_xcorr[0].data / norm2,
                    xcorr_shift
                )
                shift, value = xcorr_max(cc, abs_max=False)

                if verbose:
                    print(f"  Cross-correlation: shift={shift} samples, coef={value:.3f}")

                # Create figure
                fig = plt.figure(figsize=figsize)

                # Get all three components for plotting
                # East component
                tr_data_e_plot = st_data1.select(station=station_code, channel='H?E')
                tr_syn_e_plot = st_syn1.select(station=station_code, channel='?NE')
                if len(tr_data_e_plot) == 0:
                    tr_data_e_plot = st_data1.select(station=station_code, channel='??E')
                if len(tr_syn_e_plot) == 0:
                    tr_syn_e_plot = st_syn1.select(station=station_code, channel='??E')

                # North component
                tr_data_n_plot = st_data1.select(station=station_code, channel='H?N')
                tr_syn_n_plot = st_syn1.select(station=station_code, channel='?NN')
                if len(tr_data_n_plot) == 0:
                    tr_data_n_plot = st_data1.select(station=station_code, channel='??N')
                if len(tr_syn_n_plot) == 0:
                    tr_syn_n_plot = st_syn1.select(station=station_code, channel='??N')

                # Vertical component
                tr_data_z_plot = st_data1.select(station=station_code, channel='H?Z')
                tr_syn_z_plot = st_syn1.select(station=station_code, channel='?NZ')
                if len(tr_data_z_plot) == 0:
                    tr_data_z_plot = st_data1.select(station=station_code, channel='??Z')
                if len(tr_syn_z_plot) == 0:
                    tr_syn_z_plot = st_syn1.select(station=station_code, channel='??Z')

                # East component plot
                ax3 = plt.subplot(131)
                if len(tr_data_e_plot) > 0:
                    ax3.plot(
                        tr_data_e_plot[0].times(),
                        tr_data_e_plot[0].data / norm1,
                        'k',
                        label='E'
                    )
                if len(tr_syn_e_plot) > 0:
                    ax3.plot(
                        tr_syn_e_plot[0].times(),
                        vertical_offset + np.roll(tr_syn_e_plot[0].data, int(shift)) / norm2,
                        'r'
                    )
                ax3.set_title('sta:'+station_code)
                ax3.set_xlim([0, trim_duration])
                ax3.set_axis_off()

                # North component plot
                ax2 = plt.subplot(132, sharey=ax3, sharex=ax3)
                if len(tr_data_n_plot) > 0:
                    ax2.plot(
                        tr_data_n_plot[0].times(),
                        tr_data_n_plot[0].data / norm1,
                        'k',
                        label='N'
                    )
                if len(tr_syn_n_plot) > 0:
                    ax2.plot(
                        tr_syn_n_plot[0].times(),
                        vertical_offset + np.roll(tr_syn_n_plot[0].data, int(shift)) / norm2,
                        'r'
                    )
                ax2.set_title(f'D: {distance_km:.1f}   BAZ: {backazimuth:.1f}')
                ax2.set_axis_off()

                # Vertical component plot
                ax1 = plt.subplot(133, sharey=ax3, sharex=ax3)
                if len(tr_data_z_plot) > 0:
                    ax1.plot(
                        tr_data_z_plot[0].times(),
                        tr_data_z_plot[0].data / norm1,
                        'k',
                        label='Z'
                    )
                if len(tr_syn_z_plot) > 0:
                    ax1.plot(
                        tr_syn_z_plot[0].times(),
                        vertical_offset + np.roll(tr_syn_z_plot[0].data, int(shift)) / norm2,
                        'r'
                    )
                ax1.set_title(f'cc = {value:.3f}')
                ax1.set_axis_off()

                # Save figure
                plot_file = os.path.join(
                    output_folder,
                    f'{file_prefix}{station_code}.png'
                )
                plt.savefig(plot_file, dpi=dpi, transparent=False, bbox_inches='tight')
                plt.close()

                if verbose:
                    print(f"  Saved: {plot_file}")

                # Store results
                results.append({
                    'station': station_code,
                    'distance_km': distance_km,
                    'azimuth': azimuth,
                    'backazimuth': backazimuth,
                    'xcorr_coef': value,
                    'xcorr_shift': int(shift),
                    'max_amp_obs': norm1,
                    'max_amp_syn': norm2,
                    'plot_file': plot_file
                })

                station_count += 1

            except Exception as e:
                if verbose:
                    print(f"  Error processing {station_code}: {str(e)}")
                failed_count += 1
                continue

    if verbose:
        print(f"\n{'='*70}")
        print(f"Processing complete:")
        print(f"  Successfully processed: {station_count} stations")
        print(f"  Failed: {failed_count} stations")
        print(f"  Plots saved to: {output_folder}")

    return pd.DataFrame(results)


def plot_single_station_xcorr(
    st_data: Stream,
    st_syn: Stream,
    station_code: str,
    event_lat: float,
    event_lon: float,
    sta_lat: float,
    sta_lon: float,
    trim_duration: float = 30.0,
    taper_fraction: float = 0.5,
    xcorr_shift: int = 1000,
    vertical_offset: float = -2.0,
    figsize: Tuple[float, float] = (12, 1.5),
    output_path: Optional[str] = None,
    show_plot: bool = True
) -> Tuple[plt.Figure, Dict]:
    """
    Create cross-correlation aligned comparison for a single station.

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream
    st_syn : obspy.Stream
        Synthetic data stream
    station_code : str
        Station code
    event_lat : float
        Event latitude
    event_lon : float
        Event longitude
    sta_lat : float
        Station latitude
    sta_lon : float
        Station longitude
    trim_duration : float, optional
        Duration to trim (default: 30.0)
    taper_fraction : float, optional
        Taper fraction (default: 0.5)
    xcorr_shift : int, optional
        Max cross-correlation lag (default: 1000)
    vertical_offset : float, optional
        Vertical offset for synthetic (default: -2.0)
    figsize : tuple, optional
        Figure size (default: (12, 1.5))
    output_path : str, optional
        Path to save figure
    show_plot : bool, optional
        Display plot (default: True)

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    metrics : dict
        Dictionary with distance, azimuth, xcorr_coef, shift, amplitudes

    Examples
    --------
    >>> fig, metrics = plot_single_station_xcorr(
    ...     st_data, st_syn, 'WEL',
    ...     event_lat=-41.0, event_lon=174.0,
    ...     sta_lat=-41.3, sta_lon=174.8,
    ...     output_path='WEL_comparison.png'
    ... )
    >>> print(f"Cross-correlation: {metrics['xcorr_coef']:.3f}")
    """
    # Make copies
    st_data1 = st_data.copy()
    st_syn1 = st_syn.copy()

    # Trim
    t_init1 = st_data1[0].stats.starttime
    t_init2 = st_syn1[0].stats.starttime
    st_data1.trim(t_init1, t_init1 + trim_duration, pad=True, fill_value=0.0)
    st_syn1.trim(t_init2, t_init2 + trim_duration, pad=True, fill_value=0.0)

    # Apply taper
    for comp in ['??E', '??N', '??Z']:
        st_data1.select(station=station_code, channel=comp).taper(
            max_percentage=taper_fraction
        )

    # Compute distance and azimuth
    distance_m, azimuth, backazimuth = gps2dist_azimuth(
        event_lat, event_lon, sta_lat, sta_lon
    )
    distance_km = distance_m / 1e3

    # Normalize by max horizontal amplitude
    tr_data_e = st_data1.select(station=station_code, channel='??E')[0]
    tr_data_n = st_data1.select(station=station_code, channel='??N')[0]
    tr_syn_e = st_syn1.select(station=station_code, channel='??E')[0]
    tr_syn_n = st_syn1.select(station=station_code, channel='??N')[0]

    norm1 = np.max([np.max(np.abs(tr_data_e.data)), np.max(np.abs(tr_data_n.data))])
    norm2 = np.max([np.max(np.abs(tr_syn_e.data)), np.max(np.abs(tr_syn_n.data))])

    # Cross-correlation
    cc = correlate(tr_data_e.data / norm1, tr_syn_e.data / norm2, xcorr_shift)
    shift, value = xcorr_max(cc, abs_max=False)

    # Create plot
    fig = plt.figure(figsize=figsize)

    tr_data_z = st_data1.select(station=station_code, channel='??Z')[0]
    tr_syn_z = st_syn1.select(station=station_code, channel='??Z')[0]

    # East
    ax3 = plt.subplot(131)
    ax3.plot(tr_data_e.times(), tr_data_e.data / norm1, 'k')
    ax3.plot(tr_syn_e.times(), vertical_offset + np.roll(tr_syn_e.data, int(shift)) / norm2, 'r')
    ax3.set_title('sta:'+station_code)
    ax3.set_xlim([0, trim_duration])
    ax3.set_axis_off()

    # North
    ax2 = plt.subplot(132, sharey=ax3, sharex=ax3)
    ax2.plot(tr_data_n.times(), tr_data_n.data / norm1, 'k')
    ax2.plot(tr_syn_n.times(), vertical_offset + np.roll(tr_syn_n.data, int(shift)) / norm2, 'r')
    ax2.set_title(f'D: {distance_km:.1f}   BAZ: {backazimuth:.1f}')
    ax2.set_axis_off()

    # Vertical
    ax1 = plt.subplot(133, sharey=ax3, sharex=ax3)
    ax1.plot(tr_data_z.times(), tr_data_z.data / norm1, 'k')
    ax1.plot(tr_syn_z.times(), vertical_offset + np.roll(tr_syn_z.data, int(shift)) / norm2, 'r')
    ax1.set_title(f'cc = {value:.3f}')
    ax1.set_axis_off()

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        plt.savefig(output_path, dpi=300, transparent=False, bbox_inches='tight')

    if show_plot:
        plt.show()
    else:
        plt.close()

    metrics = {
        'distance_km': distance_km,
        'azimuth': azimuth,
        'backazimuth': backazimuth,
        'xcorr_coef': value,
        'xcorr_shift': int(shift),
        'max_amp_obs': norm1,
        'max_amp_syn': norm2
    }

    return fig, metrics


# Example usage
if __name__ == "__main__":
    print("Waveform Comparison with Cross-Correlation")
    print("=" * 70)
    print("\nAvailable functions:")
    print("  1. compare_waveforms_all_stations() - Process all stations")
    print("  2. plot_single_station_xcorr() - Single station comparison")
    print("\nImport with:")
    print("  from GMfunc.compare_waveforms_xcorr import compare_waveforms_all_stations")
