"""
Compare Fourier amplitude spectra between observed and synthetic seismograms.

This module provides functions to compute and compare Fourier amplitude spectra (FAS)
for observed and synthetic seismic data, with amplitude correction based on a reference station.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from typing import Optional, Tuple
import pandas as pd


def compute_amplitude_correction(
    st_data,
    st_syn,
    reference_station: str,
    channel: str = '*N',
    norm_factor: float = 1000.0
) -> float:
    """
    Compute amplitude correction factor from reference station.

    The correction factor is computed as:
    Acorr = max(|synthetic|) * norm_factor / max(|observed|)

    This accounts for differences in units and overall amplitude between
    observed and synthetic data.

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream
    st_syn : obspy.Stream
        Synthetic data stream
    reference_station : str
        Station code to use for amplitude calibration (e.g., 'snzo')
    channel : str, optional
        Channel pattern to use (default: '*N' for North component)
    norm_factor : float, optional
        Normalization factor for unit conversion (default: 1000.0)
        e.g., 1000.0 converts mm to m

    Returns
    -------
    Acorr : float
        Amplitude correction factor

    Examples
    --------
    >>> Acorr = compute_amplitude_correction(
    ...     st_data0,
    ...     st_syn1,
    ...     reference_station='snzo',
    ...     norm_factor=1000.0  # mm to m
    ... )
    >>> print(f"Amplitude correction: {Acorr:.2f}")
    """
    # Get traces for reference station
    tr_data = st_data.select(station=reference_station, channel=channel)[0]
    tr_syn = st_syn.select(station=reference_station, channel=channel)[0]

    # Compute max amplitudes
    max_syn = np.max(np.abs(tr_syn.data))
    max_data = np.max(np.abs(tr_data.data))

    # Compute correction factor
    Acorr = max_syn * norm_factor / max_data

    return Acorr


def compare_station_spectra(
    st_data,
    st_syn,
    station_code: str,
    fourierspec_cal_func,
    norm_data: float = 1.0,
    norm_syn: float = 1.0,
    freq_range: Tuple[float, float] = (0.0, 2.0),
    amp_range: Tuple[float, float] = (1e-10, 1e-3),
    output_path: Optional[str] = None,
    show_plot: bool = False
) -> Tuple[dict, plt.Figure]:
    """
    Compare Fourier amplitude spectra for a single station.

    Computes and plots FAS for all three components (E, N, Z) comparing
    observed data (solid lines) with synthetic data (dotted lines).

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream
    st_syn : obspy.Stream
        Synthetic data stream
    station_code : str
        Station code to analyze
    fourierspec_cal_func : callable
        Function to compute Fourier spectrum
        Should have signature: fourierspec_cal(data, delta) -> (freq, fas)
    norm_data : float, optional
        Normalization factor for observed data (default: 1.0)
    norm_syn : float, optional
        Normalization factor for synthetic data (default: 1.0)
    freq_range : tuple of float, optional
        (fmin, fmax) frequency range to display in Hz (default: (0.0, 2.0))
    amp_range : tuple of float, optional
        (amin, amax) amplitude range for y-axis (default: (1e-10, 1e-3))
    output_path : str, optional
        Path to save figure. If None, figure not saved.
    show_plot : bool, optional
        Display plot interactively (default: False)

    Returns
    -------
    spectra : dict
        Dictionary containing frequency and FAS arrays for each component:
        - 'freq': frequency array
        - 'fas_E_data': observed E component FAS
        - 'fas_E_syn': synthetic E component FAS
        - 'fas_N_data': observed N component FAS
        - 'fas_N_syn': synthetic N component FAS
        - 'fas_Z_data': observed Z component FAS
        - 'fas_Z_syn': synthetic Z component FAS
    fig : matplotlib.figure.Figure
        Figure object

    Examples
    --------
    >>> from SpecFunc.fourierspec_call import fourierspec_cal
    >>>
    >>> spectra, fig = compare_station_spectra(
    ...     st_data1,
    ...     st_syn1,
    ...     station_code='WEL',
    ...     fourierspec_cal_func=fourierspec_cal,
    ...     norm_data=1000.0,
    ...     norm_syn=Acorr,
    ...     output_path='spectra/WEL_FAS.png'
    ... )
    """
    # Get sampling interval
    delta = st_data.select(station=station_code)[0].stats['delta']

    # Compute FAS for E component
    tr_data_E = st_data.select(station=station_code, channel='H?E')[0]
    tr_syn_E = st_syn.select(station=station_code, channel='?NE')[0]
    freqlist, fas_E_data = fourierspec_cal_func(tr_data_E.data / norm_data, delta)
    freqlist_syn, fas_E_syn = fourierspec_cal_func(tr_syn_E.data / norm_syn, delta)

    # Compute FAS for N component
    tr_data_N = st_data.select(station=station_code, channel='H?N')[0]
    tr_syn_N = st_syn.select(station=station_code, channel='?NN')[0]
    _, fas_N_data = fourierspec_cal_func(tr_data_N.data / norm_data, delta)
    _, fas_N_syn = fourierspec_cal_func(tr_syn_N.data / norm_syn, delta)

    # Compute FAS for Z component
    tr_data_Z = st_data.select(station=station_code, channel='H?Z')[0]
    tr_syn_Z = st_syn.select(station=station_code, channel='?NZ')[0]
    _, fas_Z_data = fourierspec_cal_func(tr_data_Z.data / norm_data, delta)
    _, fas_Z_syn = fourierspec_cal_func(tr_syn_Z.data / norm_syn, delta)

    # Create figure
    fig = plt.figure(figsize=(5.0, 2.2))
    ax = plt.subplot(111)

    # Plot E component (East-West) - blue
    ax.plot(freqlist, fas_E_data, 'royalblue', label='EW')
    ax.plot(freqlist_syn, fas_E_syn, linestyle=':', c='royalblue')

    # Plot N component (North-South) - red
    ax.plot(freqlist, fas_N_data, 'tomato', label='NS')
    ax.plot(freqlist_syn, fas_N_syn, linestyle=':', c='tomato')

    # Plot Z component (Up-Down) - black
    ax.plot(freqlist, fas_Z_data, 'black', label='UD')
    ax.plot(freqlist_syn, fas_Z_syn, linestyle=':', c='black')

    # Formatting
    ax.legend(loc='lower right')
    ax.text(0.1, amp_range[0] * 10, station_code)
    ax.set_xlabel('freq (hz)')
    ax.set_ylabel('FAS (m/s/s)')
    ax.set_xlim(freq_range)
    ax.set_yscale('log', base=10)
    ax.grid(linestyle=':', which='both')
    ax.set_ylim(amp_range)
    ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])

    # Save or show
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=300, transparent=False, bbox_inches='tight')

    if show_plot:
        plt.show()
    else:
        plt.close()

    # Store spectra data
    spectra = {
        'freq': freqlist,
        'fas_E_data': fas_E_data,
        'fas_E_syn': fas_E_syn,
        'fas_N_data': fas_N_data,
        'fas_N_syn': fas_N_syn,
        'fas_Z_data': fas_Z_data,
        'fas_Z_syn': fas_Z_syn
    }

    return spectra, fig


def compare_all_spectra(
    st_data,
    st_syn,
    inventory,
    fourierspec_cal_func,
    reference_station: str,
    eventname: str,
    output_folder: str,
    norm_data: float = 1000.0,
    freq_range: Tuple[float, float] = (0.0, 2.0),
    amp_range: Tuple[float, float] = (1e-10, 1e-3),
    verbose: bool = True
) -> pd.DataFrame:
    """
    Compare Fourier spectra for all stations in inventory.

    Automatically computes amplitude correction from reference station,
    then processes all stations to create FAS comparison plots.

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream
    st_syn : obspy.Stream
        Synthetic data stream
    inventory : obspy.Inventory
        Station inventory
    fourierspec_cal_func : callable
        Function to compute Fourier spectrum
        Import from: SpecFunc.fourierspec_call import fourierspec_cal
    reference_station : str
        Station code for amplitude calibration (e.g., 'snzo')
    eventname : str
        Event identifier for output naming
    output_folder : str
        Base directory for output plots
    norm_data : float, optional
        Normalization for observed data, typically for unit conversion
        (default: 1000.0 for mm to m)
    freq_range : tuple of float, optional
        Frequency range to display (default: (0.0, 2.0) Hz)
    amp_range : tuple of float, optional
        Amplitude range for plot (default: (1e-10, 1e-3))
    verbose : bool, optional
        Print progress information (default: True)

    Returns
    -------
    pd.DataFrame
        Summary statistics with columns:
        - station: station code
        - Acorr: amplitude correction factor used
        - success: whether processing succeeded

    Examples
    --------
    >>> from SpecFunc.fourierspec_call import fourierspec_cal
    >>>
    >>> results = compare_all_spectra(
    ...     st_data=st_data1,
    ...     st_syn=st_syn1,
    ...     inventory=inv2,
    ...     fourierspec_cal_func=fourierspec_cal,
    ...     reference_station='snzo',
    ...     eventname='20240805',
    ...     output_folder='/path/to/plots5/',
    ...     norm_data=1000.0
    ... )
    >>> print(f"Processed {results['success'].sum()} stations successfully")

    Notes
    -----
    - Amplitude correction is computed once from reference station
    - Plots saved as: {output_folder}/m{eventname}-att/FASc-{station}.png
    - Solid lines: observed data
    - Dotted lines: synthetic data
    - Colors: blue (EW), red (NS), black (UD)
    """
    # Compute amplitude correction from reference station
    if verbose:
        print(f"Computing amplitude correction from reference station: {reference_station}")

    Acorr = compute_amplitude_correction(
        st_data,
        st_syn,
        reference_station=reference_station,
        channel='*N',
        norm_factor=norm_data
    )

    if verbose:
        print(f"Amplitude correction factor (Acorr): {Acorr:.4f}")
        print(f"Event: {eventname}")
        print(f"Data normalization: {norm_data}")
        print(f"Synthetic normalization: {Acorr}")

    # Create output folder
    outfolder = os.path.join(output_folder, f'm{eventname}-att')
    os.makedirs(outfolder, exist_ok=True)

    # Initialize results
    results = []
    success_count = 0
    fail_count = 0

    # Loop through stations
    for network in inventory:
        for sta in network:
            station_code = sta.code

            if verbose:
                print(f"Processing: {station_code}")

            try:
                # Compute and plot spectra
                output_path = os.path.join(outfolder, f'FASc-{station_code}.png')

                spectra, fig = compare_station_spectra(
                    st_data=st_data,
                    st_syn=st_syn,
                    station_code=station_code,
                    fourierspec_cal_func=fourierspec_cal_func,
                    norm_data=norm_data,
                    norm_syn=Acorr,
                    freq_range=freq_range,
                    amp_range=amp_range,
                    output_path=output_path,
                    show_plot=False
                )

                results.append({
                    'station': station_code,
                    'Acorr': Acorr,
                    'success': True
                })

                success_count += 1

            except Exception as e:
                if verbose:
                    print(f"  Failed: {str(e)}")

                results.append({
                    'station': station_code,
                    'Acorr': Acorr,
                    'success': False
                })

                fail_count += 1

    # Create results DataFrame
    results_df = pd.DataFrame(results)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Fourier Spectrum Comparison Summary:")
        print(f"{'='*60}")
        print(f"Reference station: {reference_station}")
        print(f"Amplitude correction: {Acorr:.4f}")
        print(f"Total stations: {len(results)}")
        print(f"Successful: {success_count}")
        print(f"Failed: {fail_count}")
        print(f"Output folder: {outfolder}")

    return results_df


# Example usage
if __name__ == "__main__":
    print("Module loaded successfully")
    print("Available functions:")
    print("  - compute_amplitude_correction(): Compute Acorr from reference station")
    print("  - compare_station_spectra(): Compare FAS for single station")
    print("  - compare_all_spectra(): Process all stations with automatic Acorr")
