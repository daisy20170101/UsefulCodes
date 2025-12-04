"""
Goodness-of-Fit (GOF) analysis for Fourier Amplitude Spectra.

This module provides functions to compute and visualize GOF metrics comparing
observed and synthetic Fourier amplitude spectra across multiple events.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from obspy import read
from typing import Optional, Tuple, Callable, List


def compute_gof_multiple_events(
    event_list_path: str,
    station_code: str,
    data_folder: str,
    syn_folder: str,
    snzo_folder: str,
    fourierspec_cal_func: Callable,
    konno_ohmachi_func: Callable,
    rootfolder: str = '',
    data_prefix: str = 'p',
    data_suffix: str = '.mseed',
    syn_prefix: str = 'p',
    syn_suffix: str = '-2hzB.mseed',
    snzo_prefix: str = 'p',
    snzo_suffix: str = '-snzo.mseed',
    reference_station: str = 'snzo',
    reference_channel: str = '*E',
    bandwidth: int = 30,
    freq_points: Optional[np.ndarray] = None,
    norm_data: float = 1.0,
    verbose: bool = True
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Compute Goodness-of-Fit for Fourier Amplitude Spectra across multiple events.

    GOF is computed as: log10(FAS_data / FAS_synthetic) at discrete frequencies.
    Positive values indicate synthetic underprediction, negative values indicate
    overprediction.

    Parameters
    ----------
    event_list_path : str
        Path to CSV file with event list (columns: year, month, day)
    station_code : str
        Station code to analyze (e.g., 'tfss')
    data_folder : str
        Folder containing observed data files
    syn_folder : str
        Folder containing synthetic data files
    snzo_folder : str
        Folder containing SNZO reference data
    fourierspec_cal_func : callable
        Function to compute Fourier spectrum
        Signature: fourierspec_cal(data, delta) -> (freqlist, fas)
    konno_ohmachi_func : callable
        Konno-Ohmachi smoothing function
        Signature: konno_ohmachi_smoothing(fas, freqlist, bandwidth, normalize)
    rootfolder : str, optional
        Root folder path (default: '')
    data_prefix : str, optional
        Prefix for data files (default: 'p')
    data_suffix : str, optional
        Suffix for data files (default: '.mseed')
    syn_prefix : str, optional
        Prefix for synthetic files (default: 'p')
    syn_suffix : str, optional
        Suffix for synthetic files (default: '-2hzB.mseed')
    snzo_prefix : str, optional
        Prefix for SNZO files (default: 'p')
    snzo_suffix : str, optional
        Suffix for SNZO files (default: '-snzo.mseed')
    reference_station : str, optional
        Reference station for amplitude correction (default: 'snzo')
    reference_channel : str, optional
        Channel for amplitude correction (default: '*E')
    bandwidth : int, optional
        Konno-Ohmachi smoothing bandwidth (default: 30)
    freq_points : np.ndarray, optional
        Frequency points for GOF evaluation (default: logspace(-2, 1, 100))
    norm_data : float, optional
        Normalization factor for observed data (default: 1.0)
    verbose : bool, optional
        Print processing information (default: True)

    Returns
    -------
    gof : np.ndarray
        GOF matrix, shape (n_frequencies, n_events*3)
        Each event contributes 3 columns (Z, E, N components)
    freq_points : np.ndarray
        Frequency points where GOF was evaluated
    event_names : list of str
        List of event names (YYYYMMDD format)

    Examples
    --------
    >>> from SpecFunc.fourierspec_call import fourierspec_cal, konno_ohmachi_smoothing
    >>>
    >>> gof, freqs, events = compute_gof_multiple_events(
    ...     event_list_path='evelist.csv',
    ...     station_code='tfss',
    ...     data_folder='Data2/',
    ...     syn_folder='Syn/',
    ...     snzo_folder='Data2/SNZO/',
    ...     fourierspec_cal_func=fourierspec_cal,
    ...     konno_ohmachi_func=konno_ohmachi_smoothing,
    ...     rootfolder='/path/to/data/'
    ... )
    >>>
    >>> print(f"GOF computed for {len(events)} events")
    >>> print(f"Shape: {gof.shape}")

    Notes
    -----
    - GOF is computed for all 3 components (Z, E, N) for each event
    - Amplitude correction (Acorr) computed from reference station
    - GOF = log10(FAS_obs / FAS_syn) at discrete frequency points
    - Positive GOF: synthetic underpredicts observed
    - Negative GOF: synthetic overpredicts observed
    """
    # Load event list
    if verbose:
        print(f"Loading event list from: {event_list_path}")

    evelist = pd.read_csv(event_list_path)
    year = evelist.year.to_numpy()
    month = evelist.month.to_numpy()
    day = evelist.day.to_numpy()
    inum = len(year)

    if verbose:
        print(f"Found {inum} events")

    # Define frequency points for GOF evaluation
    if freq_points is None:
        freq_points = np.logspace(-2, 1, 100)  # 0.01 to 10 Hz, 100 points

    # Initialize GOF matrix (frequencies x events*3)
    gof = np.zeros((freq_points.shape[0], 3 * inum))

    event_names = []

    if verbose:
        print(f"Processing station: {station_code}")
        print(f"Frequency points: {len(freq_points)}")

    # Process each event
    for idd, iyr in enumerate(year):
        eventname = f"{iyr}{format(month[idd], '02d')}{format(day[idd], '02d')}"
        event_names.append(eventname)

        if verbose:
            print(f"\nProcessing event {idd+1}/{inum}: {eventname}")

        try:
            # Construct file paths
            data_file = os.path.join(rootfolder, data_folder, f"{data_prefix}{eventname}{data_suffix}")
            syn_file = os.path.join(rootfolder, syn_folder, f"{syn_prefix}{eventname}{syn_suffix}")
            snzo_file = os.path.join(rootfolder, snzo_folder, f"{snzo_prefix}{eventname}{snzo_suffix}")

            # Load data
            st_data1 = read(data_file, format='MSEED')
            st_syn1 = read(syn_file, format='MSEED')
            st_data0 = read(snzo_file, format='MSEED')

            # Compute amplitude correction from reference station
            max_syn = np.max(np.abs(st_syn1.select(station=reference_station, channel=reference_channel)[0].data))
            max_data = np.max(np.abs(st_data0.select(station=reference_station, channel=reference_channel)[0].data))
            Acorr = max_syn / max_data

            if verbose:
                print(f"  Amplitude correction (Acorr): {Acorr:.4f}")

            # Normalization factors
            norm1 = norm_data
            norm2 = Acorr

            # Compute FAS for Z component
            freqlist, fas_z_data = fourierspec_cal_func(
                st_data1.select(station=station_code, channel='H?Z')[0].data / norm1,
                st_data1.select(station=station_code)[0].stats['delta']
            )
            freqlist_syn, fas_z_syn = fourierspec_cal_func(
                st_syn1.select(station=station_code, channel='??Z')[0].data / norm2,
                st_syn1.select(station=station_code)[0].stats['delta']
            )

            # Compute FAS for E component
            freqlist, fas_e_data = fourierspec_cal_func(
                st_data1.select(station=station_code, channel='H?E')[0].data / norm1,
                st_data1.select(station=station_code)[0].stats['delta']
            )
            freqlist_syn, fas_e_syn = fourierspec_cal_func(
                st_syn1.select(station=station_code, channel='??E')[0].data / norm2,
                st_syn1.select(station=station_code)[0].stats['delta']
            )

            # Compute FAS for N component
            freqlist, fas_n_data = fourierspec_cal_func(
                st_data1.select(station=station_code, channel='H?N')[0].data / norm1,
                st_data1.select(station=station_code)[0].stats['delta']
            )
            freqlist_syn, fas_n_syn = fourierspec_cal_func(
                st_syn1.select(station=station_code, channel='??N')[0].data / norm2,
                st_syn1.select(station=station_code)[0].stats['delta']
            )

            # Apply Konno-Ohmachi smoothing
            fas_z_smooth_data = konno_ohmachi_func(fas_z_data, freqlist, bandwidth=bandwidth, normalize=True)
            fas_z_smooth_syn = konno_ohmachi_func(fas_z_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)

            fas_e_smooth_data = konno_ohmachi_func(fas_e_data, freqlist, bandwidth=bandwidth, normalize=True)
            fas_e_smooth_syn = konno_ohmachi_func(fas_e_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)

            fas_n_smooth_data = konno_ohmachi_func(fas_n_data, freqlist, bandwidth=bandwidth, normalize=True)
            fas_n_smooth_syn = konno_ohmachi_func(fas_n_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)

            # Compute GOF at discrete frequency points
            for iff, ff in enumerate(freq_points):
                # Find nearest frequency indices
                indx = np.abs(freqlist - ff).argmin()
                indy = np.abs(freqlist_syn - ff).argmin()

                # GOF = log10(data / synthetic)
                gof[iff, 3*idd] = np.log10(fas_z_smooth_data[indx] / fas_z_smooth_syn[indy])      # Z
                gof[iff, 3*idd + 1] = np.log10(fas_e_smooth_data[indx] / fas_e_smooth_syn[indy])  # E
                gof[iff, 3*idd + 2] = np.log10(fas_n_smooth_data[indx] / fas_n_smooth_syn[indy])  # N

            if verbose:
                print(f"  Successfully processed")

        except Exception as e:
            if verbose:
                print(f"  Error processing {eventname}: {str(e)}")
            # Fill with NaN for failed events
            gof[:, 3*idd:3*idd+3] = np.nan
            continue

    if verbose:
        print(f"\n{'='*70}")
        print(f"GOF computation complete")
        print(f"  Events processed: {inum}")
        print(f"  Station: {station_code}")
        print(f"  Components per event: 3 (Z, E, N)")
        print(f"  Total GOF curves: {inum * 3}")

    return gof, freq_points, event_names


def plot_gof(
    gof: np.ndarray,
    freq_points: np.ndarray,
    station_code: str,
    output_path: Optional[str] = None,
    freq_range: Tuple[float, float] = (0.0, 2.0),
    gof_range: Tuple[float, float] = (-5, 5),
    figsize: Tuple[float, float] = (4.5, 2.5),
    dpi: int = 300,
    show_plot: bool = False,
    invert_sign: bool = True
) -> plt.Figure:
    """
    Plot Goodness-of-Fit results with mean and standard deviation.

    Creates a figure showing individual GOF curves (transparent) and
    the mean ± 1 std deviation band.

    Parameters
    ----------
    gof : np.ndarray
        GOF matrix from compute_gof_multiple_events
    freq_points : np.ndarray
        Frequency points
    station_code : str
        Station code for title
    output_path : str, optional
        Path to save figure
    freq_range : tuple, optional
        (min_freq, max_freq) for x-axis (default: (0.0, 2.0))
    gof_range : tuple, optional
        (min_gof, max_gof) for y-axis (default: (-5, 5))
    figsize : tuple, optional
        Figure size (default: (4.5, 2.5))
    dpi : int, optional
        Plot resolution (default: 300)
    show_plot : bool, optional
        Display plot interactively (default: False)
    invert_sign : bool, optional
        Invert GOF sign for plotting (default: True)
        When True: positive = overprediction, negative = underprediction
        When False: positive = underprediction, negative = overprediction

    Returns
    -------
    matplotlib.figure.Figure
        Figure object

    Examples
    --------
    >>> fig = plot_gof(
    ...     gof, freq_points, 'tfss',
    ...     output_path='GOF-tfss.png',
    ...     show_plot=True
    ... )

    Notes
    -----
    - Individual curves shown in light blue (skyblue, alpha=0.3)
    - Mean curve shown in royal blue with ±1 std shaded band
    - Default sign convention (invert_sign=True):
      - Positive GOF: synthetic overpredicts observed
      - Negative GOF: synthetic underpredicts observed
    """
    # Invert sign if requested (for plotting convention)
    if invert_sign:
        gof_plot = -gof
    else:
        gof_plot = gof

    # Create figure
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)

    # Plot individual GOF curves (all events, all components)
    n_curves = gof_plot.shape[1]
    for ieve in range(n_curves):
        ax.plot(freq_points, gof_plot[:, ieve], 'skyblue', alpha=0.3)

    # Compute mean and std
    gof_mean = np.mean(gof_plot, axis=1)
    gof_std = np.std(gof_plot, axis=1)

    # Plot mean ± 1 std band
    ax.fill_between(
        freq_points,
        gof_mean + gof_std,
        gof_mean - gof_std,
        color='royalblue',
        alpha=0.5
    )

    # Plot mean line
    ax.plot(freq_points, gof_mean, 'royalblue', linewidth=2, label='mean')

    # Formatting
    ax.set_title(f'site: {station_code}')
    ax.set_ylim(gof_range)
    ax.set_xlabel('frequency (Hz)')
    ax.set_ylabel(r"log10($\frac{FAS^c_{mod}}{FAS_{data}}$)")
    ax.set_xlim(freq_range)

    # Set specific ticks if in standard range
    if freq_range == (0.0, 2.0):
        ax.set_xticks([0.5, 1, 1.5, 2])

    if gof_range == (-5, 5):
        ax.set_yticks([-6, -4, -2, 0, 2, 4, 6])

    # Add labels for overprediction/underprediction
    if invert_sign:
        ax.text(freq_range[0] + 0.1, gof_range[1] - 0.5, 'overpredicted')
        ax.text(freq_range[0] + 0.1, gof_range[0] + 0.5, 'underpredicted')
    else:
        ax.text(freq_range[0] + 0.1, gof_range[1] - 0.5, 'underpredicted')
        ax.text(freq_range[0] + 0.1, gof_range[0] + 0.5, 'overpredicted')

    ax.grid(linestyle=':', which='major')
    ax.legend()

    plt.tight_layout()

    # Save if requested
    if output_path:
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        plt.savefig(output_path, dpi=dpi, transparent=False, bbox_inches='tight')
        print(f"Saved GOF plot: {output_path}")

    # Show or close
    if show_plot:
        plt.show()
    else:
        plt.close()

    return fig


def compute_and_plot_gof(
    event_list_path: str,
    station_code: str,
    data_folder: str,
    syn_folder: str,
    snzo_folder: str,
    fourierspec_cal_func: Callable,
    konno_ohmachi_func: Callable,
    output_path: str,
    rootfolder: str = '',
    bandwidth: int = 30,
    freq_range: Tuple[float, float] = (0.0, 2.0),
    verbose: bool = True,
    **kwargs
) -> Tuple[np.ndarray, np.ndarray, List[str], plt.Figure]:
    """
    Compute GOF and create plot in one function call.

    Convenience function that combines compute_gof_multiple_events and plot_gof.

    Parameters
    ----------
    event_list_path : str
        Path to CSV file with event list
    station_code : str
        Station code to analyze
    data_folder : str
        Folder containing observed data
    syn_folder : str
        Folder containing synthetic data
    snzo_folder : str
        Folder containing SNZO reference data
    fourierspec_cal_func : callable
        Function to compute Fourier spectrum
    konno_ohmachi_func : callable
        Konno-Ohmachi smoothing function
    output_path : str
        Path to save plot
    rootfolder : str, optional
        Root folder path (default: '')
    bandwidth : int, optional
        Smoothing bandwidth (default: 30)
    freq_range : tuple, optional
        Frequency range for plotting (default: (0.0, 2.0))
    verbose : bool, optional
        Print information (default: True)
    **kwargs
        Additional arguments passed to compute_gof_multiple_events and plot_gof

    Returns
    -------
    gof : np.ndarray
        GOF matrix
    freq_points : np.ndarray
        Frequency points
    event_names : list of str
        Event names
    fig : matplotlib.figure.Figure
        Figure object

    Examples
    --------
    >>> from SpecFunc.fourierspec_call import fourierspec_cal, konno_ohmachi_smoothing
    >>>
    >>> gof, freqs, events, fig = compute_and_plot_gof(
    ...     event_list_path='evelist.csv',
    ...     station_code='tfss',
    ...     data_folder='Data2/',
    ...     syn_folder='Syn/',
    ...     snzo_folder='Data2/SNZO/',
    ...     fourierspec_cal_func=fourierspec_cal,
    ...     konno_ohmachi_func=konno_ohmachi_smoothing,
    ...     output_path='GOF-tfss.png',
    ...     rootfolder='/path/to/data/'
    ... )
    """
    # Compute GOF
    gof, freq_points, event_names = compute_gof_multiple_events(
        event_list_path=event_list_path,
        station_code=station_code,
        data_folder=data_folder,
        syn_folder=syn_folder,
        snzo_folder=snzo_folder,
        fourierspec_cal_func=fourierspec_cal_func,
        konno_ohmachi_func=konno_ohmachi_func,
        rootfolder=rootfolder,
        bandwidth=bandwidth,
        verbose=verbose,
        **{k: v for k, v in kwargs.items() if k in [
            'data_prefix', 'data_suffix', 'syn_prefix', 'syn_suffix',
            'snzo_prefix', 'snzo_suffix', 'reference_station', 'reference_channel',
            'freq_points', 'norm_data'
        ]}
    )

    # Create plot
    fig = plot_gof(
        gof=gof,
        freq_points=freq_points,
        station_code=station_code,
        output_path=output_path,
        freq_range=freq_range,
        **{k: v for k, v in kwargs.items() if k in [
            'gof_range', 'figsize', 'dpi', 'show_plot', 'invert_sign'
        ]}
    )

    return gof, freq_points, event_names, fig


# Example usage
if __name__ == "__main__":
    print("Goodness-of-Fit (GOF) Analysis for FAS")
    print("=" * 70)
    print("\nAvailable functions:")
    print("  1. compute_gof_multiple_events() - Compute GOF across multiple events")
    print("  2. plot_gof() - Plot GOF with mean and std")
    print("  3. compute_and_plot_gof() - Combined function (recommended)")
    print("\nImport with:")
    print("  from GMfunc.gof_analysis import compute_and_plot_gof")
