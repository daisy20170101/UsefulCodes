"""
Compute spectral ratios for site characterization and validation.

This module provides functions to compute:
1. Site-to-reference spectral ratios (comparing each site to a reference station)
2. Horizontal-to-vertical spectral ratios (HVSR)

Both use Konno-Ohmachi smoothing for stable spectral estimates.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from obspy import Stream
from typing import Optional, Tuple, Callable, List

from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing


def compute_site_to_reference_ratio(
    st_data: Stream,
    st_syn: Stream,
    inventory,
    reference_station: str,
    fourierspec_cal_func: Callable,
    konno_ohmachi_func: Callable,
    eventname: str,
    output_folder: str,
    bandwidth: int = 40,
    freq_range: Tuple[float, float] = (0.05, 2.0),
    ratio_range: Tuple[float, float] = (0, 10),
    figsize: Tuple[float, float] = (4.5, 2.2),
    dpi: int = 300,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Compute site-to-reference spectral ratios for all stations.

    This function computes the ratio of horizontal Fourier amplitude spectra
    at each station relative to a reference station. This helps identify
    site amplification patterns.

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream
    st_syn : obspy.Stream
        Synthetic data stream
    inventory : obspy.Inventory
        Station inventory
    reference_station : str
        Reference station code (e.g., 'ptos')
    fourierspec_cal_func : callable
        Function to compute Fourier spectrum
        Should have signature: fourierspec_cal(data, delta) -> (freqlist, fas)
    konno_ohmachi_func : callable
        Konno-Ohmachi smoothing function
        Should have signature: konno_ohmachi_smoothing(fas, freqlist, bandwidth, normalize)
    eventname : str
        Event name for output files
    output_folder : str
        Folder to save plots
    bandwidth : int, optional
        Konno-Ohmachi smoothing bandwidth (default: 40)
    freq_range : tuple, optional
        (min_freq, max_freq) for plotting (default: (0.05, 2.0))
    ratio_range : tuple, optional
        (min_ratio, max_ratio) for y-axis (default: (0, 10))
    figsize : tuple, optional
        Figure size (default: (4.5, 2.2))
    dpi : int, optional
        Plot resolution (default: 300)
    verbose : bool, optional
        Print processing information (default: True)

    Returns
    -------
    pd.DataFrame
        Results with columns:
        - station: station code
        - mean_ratio_data: mean spectral ratio for observed data
        - mean_ratio_syn: mean spectral ratio for synthetic data
        - max_ratio_data: maximum spectral ratio for observed data
        - max_ratio_syn: maximum spectral ratio for synthetic data
        - freq_at_max_data: frequency of maximum ratio (observed)
        - freq_at_max_syn: frequency of maximum ratio (synthetic)
        - plot_file: path to saved plot

    Examples
    --------
    >>> from SpecFunc.fourierspec_call import fourierspec_cal, konno_ohmachi_smoothing
    >>>
    >>> results = compute_site_to_reference_ratio(
    ...     st_data=st_data1,
    ...     st_syn=st_syn1,
    ...     inventory=inv2,
    ...     reference_station='ptos',
    ...     fourierspec_cal_func=fourierspec_cal,
    ...     konno_ohmachi_func=konno_ohmachi_smoothing,
    ...     eventname='20240805',
    ...     output_folder='./spectral_ratios/'
    ... )
    >>>
    >>> # View amplification factors
    >>> print(results[['station', 'max_ratio_data', 'freq_at_max_data']])

    Notes
    -----
    - Computes horizontal component as geometric mean: H = sqrt(0.5*(E^2 + N^2))
    - Reference station spectra computed once, then used for all comparisons
    - Konno-Ohmachi smoothing applied to stabilize spectral estimates
    - Ratio = H_site / H_reference
    """
    # Create output folder
    os.makedirs(output_folder, exist_ok=True)

    if verbose:
        print(f"Computing site-to-reference ratios")
        print(f"Reference station: {reference_station}")
        print(f"Event: {eventname}")

    # Compute reference station spectra (E and N components)
    if verbose:
        print(f"\nProcessing reference station: {reference_station}")

    try:
        # Observed data - reference station
        freqlist, ref1_data = fourierspec_cal_func(
            st_data.select(station=reference_station, channel='H?E')[0].data,
            st_data.select(station=reference_station)[0].stats['delta']
        )
        freqlist, ref2_data = fourierspec_cal_func(
            st_data.select(station=reference_station, channel='H?N')[0].data,
            st_data.select(station=reference_station)[0].stats['delta']
        )

        # Synthetic data - reference station
        freqlist_syn, ref1_syn = fourierspec_cal_func(
            st_syn.select(station=reference_station, channel='?NE')[0].data,
            st_syn.select(station=reference_station)[0].stats['delta']
        )
        freqlist_syn, ref2_syn = fourierspec_cal_func(
            st_syn.select(station=reference_station, channel='?NN')[0].data,
            st_syn.select(station=reference_station)[0].stats['delta']
        )

        # Apply Konno-Ohmachi smoothing
        fas_SA1_data = konno_ohmachi_func(ref1_data, freqlist, bandwidth=bandwidth, normalize=True)
        fas_SA2_data = konno_ohmachi_func(ref2_data, freqlist, bandwidth=bandwidth, normalize=True)

        fas_SA1_syn = konno_ohmachi_func(ref1_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)
        fas_SA2_syn = konno_ohmachi_func(ref2_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)

        # Compute horizontal component (geometric mean)
        href_data = np.sqrt(0.5 * (fas_SA2_data**2 + fas_SA1_data**2))
        href_syn = np.sqrt(0.5 * (fas_SA2_syn**2 + fas_SA1_syn**2))

        # Store reference frequency lists for interpolation
        freqlist_ref_data = freqlist.copy()
        freqlist_ref_syn = freqlist_syn.copy()

        if verbose:
            print(f"  Reference spectra computed successfully")
            print(f"  Reference freq array size: data={len(freqlist_ref_data)}, syn={len(freqlist_ref_syn)}")

    except Exception as e:
        print(f"Error computing reference station spectra: {str(e)}")
        return pd.DataFrame()

    # Process all stations
    results = []
    station_count = 0
    failed_count = 0

    for network in inventory:
        for sta in network:
            station_code = sta.code

            if verbose:
                print(f"\nProcessing station: {station_code}")

            try:
                # Observed data - station
                freqlist, fasZ_data = fourierspec_cal_func(
                    st_data.select(station=station_code, channel='H?Z')[0].data,
                    st_data.select(station=station_code)[0].stats['delta']
                )
                freqlist, fas1_data = fourierspec_cal_func(
                    st_data.select(station=station_code, channel='H?E')[0].data,
                    st_data.select(station=station_code)[0].stats['delta']
                )
                freqlist, fas2_data = fourierspec_cal_func(
                    st_data.select(station=station_code, channel='H?N')[0].data,
                    st_data.select(station=station_code)[0].stats['delta']
                )

                # Synthetic data - station
                freqlist_syn, fasZ_syn = fourierspec_cal_func(
                    st_syn.select(station=station_code, channel='?NZ')[0].data,
                    st_syn.select(station=station_code)[0].stats['delta']
                )
                freqlist_syn, fas1_syn = fourierspec_cal_func(
                    st_syn.select(station=station_code, channel='?NE')[0].data,
                    st_syn.select(station=station_code)[0].stats['delta']
                )
                freqlist_syn, fas2_syn = fourierspec_cal_func(
                    st_syn.select(station=station_code, channel='?NN')[0].data,
                    st_syn.select(station=station_code)[0].stats['delta']
                )

                # Apply Konno-Ohmachi smoothing
                fasZ_f_data = konno_ohmachi_func(fasZ_data, freqlist, bandwidth=bandwidth, normalize=True)
                fas_SA1_data = konno_ohmachi_func(fas1_data, freqlist, bandwidth=bandwidth, normalize=True)
                fas_SA2_data = konno_ohmachi_func(fas2_data, freqlist, bandwidth=bandwidth, normalize=True)

                fasZ_f_syn = konno_ohmachi_func(fasZ_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)
                fas_SA1_syn = konno_ohmachi_func(fas1_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)
                fas_SA2_syn = konno_ohmachi_func(fas2_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)

                # Compute horizontal component (geometric mean)
                hf_data = np.sqrt(0.5 * (fas_SA1_data**2 + fas_SA2_data**2))
                hf_syn = np.sqrt(0.5 * (fas_SA1_syn**2 + fas_SA2_syn**2))

                # Interpolate reference spectrum to match station frequency grid (if needed)
                # This handles cases where traces have different lengths
                if len(freqlist) != len(freqlist_ref_data):
                    if verbose:
                        print(f"  Interpolating reference (data): {len(freqlist_ref_data)} -> {len(freqlist)} points")
                    href_data_interp = np.interp(freqlist, freqlist_ref_data, href_data)
                else:
                    href_data_interp = href_data

                if len(freqlist_syn) != len(freqlist_ref_syn):
                    if verbose:
                        print(f"  Interpolating reference (syn): {len(freqlist_ref_syn)} -> {len(freqlist_syn)} points")
                    href_syn_interp = np.interp(freqlist_syn, freqlist_ref_syn, href_syn)
                else:
                    href_syn_interp = href_syn

                # Compute spectral ratios (site / reference)
                ratio_data = hf_data / href_data_interp
                ratio_syn = hf_syn / href_syn_interp

                # Find frequency range for statistics
                freq_mask = (freqlist >= freq_range[0]) & (freqlist <= freq_range[1])
                freq_mask_syn = (freqlist_syn >= freq_range[0]) & (freqlist_syn <= freq_range[1])

                # Compute statistics
                mean_ratio_data = np.mean(ratio_data[freq_mask])
                mean_ratio_syn = np.mean(ratio_syn[freq_mask_syn])

                max_idx_data = np.argmax(ratio_data[freq_mask])
                max_idx_syn = np.argmax(ratio_syn[freq_mask_syn])

                max_ratio_data = ratio_data[freq_mask][max_idx_data]
                max_ratio_syn = ratio_syn[freq_mask_syn][max_idx_syn]

                freq_at_max_data = freqlist[freq_mask][max_idx_data]
                freq_at_max_syn = freqlist_syn[freq_mask_syn][max_idx_syn]

                if verbose:
                    print(f"  Max ratio - Data: {max_ratio_data:.2f} at {freq_at_max_data:.2f} Hz")
                    print(f"  Max ratio - Syn: {max_ratio_syn:.2f} at {freq_at_max_syn:.2f} Hz")

                # Determine y-axis limits based on max amplitude in frequency range
                max_amp_in_range = max(max_ratio_data, max_ratio_syn)
                ylim_max = max_amp_in_range * 1.15  # Add 15% margin above max

                # Create plot
                fig = plt.figure(figsize=figsize)
                ax = plt.subplot(111)

                ax.plot(freqlist, ratio_data, 'black', label='Data')
                ax.plot(freqlist_syn, ratio_syn, 'tomato', label='Model')

                ax.legend(loc='upper right')
                ax.set_title( 'sta:'+station_code)

                ax.set_xlabel('frequency (Hz)')
                ax.set_ylabel('amplitude ratio')
                ax.grid(linestyle=':', which='both')
                ax.set_ylim(0, ylim_max)  # Auto-scale based on data
                ax.set_xlim(freq_range)

                # Save plot
                plot_file = os.path.join(output_folder, f'StoR-Amp-{station_code}.png')
                plt.savefig(plot_file, dpi=dpi, transparent=False, bbox_inches='tight')
                plt.close()

                if verbose:
                    print(f"  Saved: {plot_file}")

                # Store results
                results.append({
                    'station': station_code,
                    'mean_ratio_data': mean_ratio_data,
                    'mean_ratio_syn': mean_ratio_syn,
                    'max_ratio_data': max_ratio_data,
                    'max_ratio_syn': max_ratio_syn,
                    'freq_at_max_data': freq_at_max_data,
                    'freq_at_max_syn': freq_at_max_syn,
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
        print(f"  Reference station: {reference_station}")

    return pd.DataFrame(results)


def compute_hvsr(
    st_data: Stream,
    st_syn: Stream,
    inventory,
    fourierspec_cal_func: Callable,
    konno_ohmachi_func: Callable,
    eventname: str,
    output_folder: str,
    bandwidth: int = 30,
    freq_range: Tuple[float, float] = (0.05, 2.0),
    ratio_range: Tuple[float, float] = (0, 10),
    figsize: Tuple[float, float] = (4.5, 2.2),
    dpi: int = 300,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Compute Horizontal-to-Vertical Spectral Ratios (HVSR) for all stations.

    HVSR is commonly used for site characterization and identifying
    fundamental resonance frequencies.

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
    konno_ohmachi_func : callable
        Konno-Ohmachi smoothing function
    eventname : str
        Event name for output files
    output_folder : str
        Folder to save plots
    bandwidth : int, optional
        Konno-Ohmachi smoothing bandwidth (default: 30)
    freq_range : tuple, optional
        (min_freq, max_freq) for plotting (default: (0.05, 2.0))
    ratio_range : tuple, optional
        (min_ratio, max_ratio) for y-axis (default: (0, 10))
    figsize : tuple, optional
        Figure size (default: (4.5, 2.2))
    dpi : int, optional
        Plot resolution (default: 300)
    verbose : bool, optional
        Print processing information (default: True)

    Returns
    -------
    pd.DataFrame
        Results with columns:
        - station: station code
        - peak_freq_data: frequency of HVSR peak (observed)
        - peak_freq_syn: frequency of HVSR peak (synthetic)
        - peak_amplitude_data: HVSR peak amplitude (observed)
        - peak_amplitude_syn: HVSR peak amplitude (synthetic)
        - mean_hvsr_data: mean HVSR (observed)
        - mean_hvsr_syn: mean HVSR (synthetic)
        - plot_file: path to saved plot

    Examples
    --------
    >>> from SpecFunc.fourierspec_call import fourierspec_cal, konno_ohmachi_smoothing
    >>>
    >>> results = compute_hvsr(
    ...     st_data=st_data1,
    ...     st_syn=st_syn1,
    ...     inventory=inv2,
    ...     fourierspec_cal_func=fourierspec_cal,
    ...     konno_ohmachi_func=konno_ohmachi_smoothing,
    ...     eventname='20240805',
    ...     output_folder='./hvsr/'
    ... )
    >>>
    >>> # View fundamental frequencies
    >>> print(results[['station', 'peak_freq_data', 'peak_amplitude_data']])

    Notes
    -----
    - Horizontal component: H = sqrt(0.5*(E^2 + N^2))
    - Vertical component: Z
    - HVSR = H / Z
    - Peak frequency often indicates site fundamental resonance
    """
    # Create output folder
    os.makedirs(output_folder, exist_ok=True)

    if verbose:
        print(f"Computing Horizontal-to-Vertical Spectral Ratios (HVSR)")
        print(f"Event: {eventname}")

    results = []
    station_count = 0
    failed_count = 0

    for network in inventory:
        for sta in network:
            station_code = sta.code

            if verbose:
                print(f"\nProcessing station: {station_code}")

            try:
                # Observed data
                freqlist, fasZ_data = fourierspec_cal_func(
                    st_data.select(station=station_code, channel='H?Z')[0].data,
                    st_data.select(station=station_code)[0].stats['delta']
                )
                freqlist, fas1_data = fourierspec_cal_func(
                    st_data.select(station=station_code, channel='H?E')[0].data,
                    st_data.select(station=station_code)[0].stats['delta']
                )
                freqlist, fas2_data = fourierspec_cal_func(
                    st_data.select(station=station_code, channel='H?N')[0].data,
                    st_data.select(station=station_code)[0].stats['delta']
                )

                # Synthetic data
                freqlist_syn, fasZ_syn = fourierspec_cal_func(
                    st_syn.select(station=station_code, channel='?NZ')[0].data,
                    st_syn.select(station=station_code)[0].stats['delta']
                )
                freqlist_syn, fas1_syn = fourierspec_cal_func(
                    st_syn.select(station=station_code, channel='?NE')[0].data,
                    st_syn.select(station=station_code)[0].stats['delta']
                )
                freqlist_syn, fas2_syn = fourierspec_cal_func(
                    st_syn.select(station=station_code, channel='?NN')[0].data,
                    st_syn.select(station=station_code)[0].stats['delta']
                )

                # Apply Konno-Ohmachi smoothing
                fasZ_f_data = konno_ohmachi_func(fasZ_data, freqlist, bandwidth=bandwidth, normalize=True)
                fas_SA1_data = konno_ohmachi_func(fas1_data, freqlist, bandwidth=bandwidth, normalize=True)
                fas_SA2_data = konno_ohmachi_func(fas2_data, freqlist, bandwidth=bandwidth, normalize=True)

                fasZ_f_syn = konno_ohmachi_func(fasZ_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)
                fas_SA1_syn = konno_ohmachi_func(fas1_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)
                fas_SA2_syn = konno_ohmachi_func(fas2_syn, freqlist_syn, bandwidth=bandwidth, normalize=True)

                # Compute horizontal component (geometric mean)
                hf_data = np.sqrt(0.5 * (fas_SA1_data**2 + fas_SA2_data**2))
                hf_syn = np.sqrt(0.5 * (fas_SA1_syn**2 + fas_SA2_syn**2))

                # Compute HVSR
                hvsr_data = hf_data / fasZ_f_data
                hvsr_syn = hf_syn / fasZ_f_syn

                # Find frequency range for statistics
                freq_mask = (freqlist >= freq_range[0]) & (freqlist <= freq_range[1])
                freq_mask_syn = (freqlist_syn >= freq_range[0]) & (freqlist_syn <= freq_range[1])

                # Find peak
                peak_idx_data = np.argmax(hvsr_data[freq_mask])
                peak_idx_syn = np.argmax(hvsr_syn[freq_mask_syn])

                peak_freq_data = freqlist[freq_mask][peak_idx_data]
                peak_freq_syn = freqlist_syn[freq_mask_syn][peak_idx_syn]

                peak_amplitude_data = hvsr_data[freq_mask][peak_idx_data]
                peak_amplitude_syn = hvsr_syn[freq_mask_syn][peak_idx_syn]

                mean_hvsr_data = np.mean(hvsr_data[freq_mask])
                mean_hvsr_syn = np.mean(hvsr_syn[freq_mask_syn])

                if verbose:
                    print(f"  Peak - Data: {peak_amplitude_data:.2f} at {peak_freq_data:.2f} Hz")
                    print(f"  Peak - Syn: {peak_amplitude_syn:.2f} at {peak_freq_syn:.2f} Hz")

                # Determine y-axis limits based on max amplitude in frequency range
                max_amp_in_range = max(peak_amplitude_data, peak_amplitude_syn)
                ylim_max = max_amp_in_range * 1.15  # Add 15% margin above max

                # Create plot
                fig = plt.figure(figsize=figsize)
                ax = plt.subplot(111)

                ax.plot(freqlist, hvsr_data, 'black', label='Data')
                ax.plot(freqlist_syn, hvsr_syn, 'tomato', label='Model')

                ax.legend(loc='upper right')
                ax.set_title('sta:' + station_code)

                ax.set_xlabel('frequency (Hz)')
                ax.set_ylabel('HV spectral amplitude ratio')
                ax.grid(linestyle=':', which='both')
                ax.set_ylim(0, ylim_max)  # Auto-scale based on data
                ax.set_xlim(freq_range)

                # Save plot
                plot_file = os.path.join(output_folder, f'HVSR-{station_code}.png')
                plt.savefig(plot_file, dpi=dpi, transparent=False, bbox_inches='tight')
                plt.close()

                if verbose:
                    print(f"  Saved: {plot_file}")

                # Store results
                results.append({
                    'station': station_code,
                    'peak_freq_data': peak_freq_data,
                    'peak_freq_syn': peak_freq_syn,
                    'peak_amplitude_data': peak_amplitude_data,
                    'peak_amplitude_syn': peak_amplitude_syn,
                    'mean_hvsr_data': mean_hvsr_data,
                    'mean_hvsr_syn': mean_hvsr_syn,
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

    return pd.DataFrame(results)


def compute_hvsr_two_models(
    st_data: Stream,
    st_syn: Stream,
    st_syn2: Stream,
    inventory,
    fourierspec_cal_func: Callable,
    konno_ohmachi_func: Callable,
    eventname: str,
    output_folder: str,
    bandwidth: int = 30,
    freq_range: Tuple[float, float] = (0.05, 2.0),
    ratio_range: Tuple[float, float] = (0, 10),
    model1_label: str = 'Model 1',
    model2_label: str = 'Model 2',
    model1_color: str = 'tomato',
    model2_color: str = 'dodgerblue',
    figsize: Tuple[float, float] = (4.5, 2.2),
    dpi: int = 300,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Compute HVSR for observed data and two synthetic models for comparison.

    This function extends compute_hvsr to compare observed HVSR with two
    different synthetic models (e.g., different velocity structures, attenuation
    models, or source parameters).

    Parameters
    ----------
    st_data : obspy.Stream
        Observed data stream
    st_syn : obspy.Stream
        First synthetic model stream
    st_syn2 : obspy.Stream
        Second synthetic model stream
    inventory : obspy.Inventory
        Station inventory
    fourierspec_cal_func : callable
        Function to compute Fourier spectrum
    konno_ohmachi_func : callable
        Konno-Ohmachi smoothing function
    eventname : str
        Event name for output files
    output_folder : str
        Folder to save plots
    bandwidth : int, optional
        Konno-Ohmachi smoothing bandwidth (default: 30)
    freq_range : tuple, optional
        (min_freq, max_freq) for plotting (default: (0.05, 2.0))
    ratio_range : tuple, optional
        (min_ratio, max_ratio) for y-axis (default: (0, 10))
    model1_label : str, optional
        Legend label for first model (default: 'Model 1')
    model2_label : str, optional
        Legend label for second model (default: 'Model 2')
    model1_color : str, optional
        Color for first model (default: 'tomato')
    model2_color : str, optional
        Color for second model (default: 'dodgerblue')
    figsize : tuple, optional
        Figure size (default: (4.5, 2.2))
    dpi : int, optional
        Plot resolution (default: 300)
    verbose : bool, optional
        Print processing information (default: True)

    Returns
    -------
    pd.DataFrame
        Results with columns:
        - station: station code
        - peak_freq_data: frequency of HVSR peak (observed)
        - peak_freq_syn1: frequency of HVSR peak (model 1)
        - peak_freq_syn2: frequency of HVSR peak (model 2)
        - peak_amplitude_data: HVSR peak amplitude (observed)
        - peak_amplitude_syn1: HVSR peak amplitude (model 1)
        - peak_amplitude_syn2: HVSR peak amplitude (model 2)
        - mean_hvsr_data: mean HVSR (observed)
        - mean_hvsr_syn1: mean HVSR (model 1)
        - mean_hvsr_syn2: mean HVSR (model 2)
        - plot_file: path to saved plot

    Examples
    --------
    >>> from SpecFunc.fourierspec_call import fourierspec_cal, konno_ohmachi_smoothing
    >>>
    >>> results = compute_hvsr_two_models(
    ...     st_data=st_data1,
    ...     st_syn=st_syn1,
    ...     st_syn2=st_syn2,
    ...     inventory=inv2,
    ...     fourierspec_cal_func=fourierspec_cal,
    ...     konno_ohmachi_func=konno_ohmachi_smoothing,
    ...     eventname='20240805',
    ...     output_folder='./hvsr_comparison/',
    ...     model1_label='1D Model',
    ...     model2_label='3D Model'
    ... )

    Notes
    -----
    - Horizontal component: H = sqrt(0.5*(E^2 + N^2))
    - Vertical component: Z
    - HVSR = H / Z
    - Useful for comparing different velocity models or attenuation parameters
    """
    # Create output folder
    os.makedirs(output_folder, exist_ok=True)

    if verbose:
        print(f"Computing HVSR for observed data and two synthetic models")
        print(f"Event: {eventname}")
        print(f"Model 1: {model1_label}")
        print(f"Model 2: {model2_label}")

    results = []
    station_count = 0
    failed_count = 0

    for network in inventory:
        for sta in network:
            station_code = sta.code

            if verbose:
                print(f"\nProcessing station: {station_code}")

            try:
                # Observed data
                freqlist, fasZ_data = fourierspec_cal_func(
                    st_data.select(station=station_code, channel='H?Z')[0].data,
                    st_data.select(station=station_code)[0].stats['delta']
                )
                freqlist, fas1_data = fourierspec_cal_func(
                    st_data.select(station=station_code, channel='H?E')[0].data,
                    st_data.select(station=station_code)[0].stats['delta']
                )
                freqlist, fas2_data = fourierspec_cal_func(
                    st_data.select(station=station_code, channel='H?N')[0].data,
                    st_data.select(station=station_code)[0].stats['delta']
                )

                # Synthetic model 1
                freqlist_syn1, fasZ_syn1 = fourierspec_cal_func(
                    st_syn.select(station=station_code, channel='?NZ')[0].data,
                    st_syn.select(station=station_code)[0].stats['delta']
                )
                freqlist_syn1, fas1_syn1 = fourierspec_cal_func(
                    st_syn.select(station=station_code, channel='?NE')[0].data,
                    st_syn.select(station=station_code)[0].stats['delta']
                )
                freqlist_syn1, fas2_syn1 = fourierspec_cal_func(
                    st_syn.select(station=station_code, channel='?NN')[0].data,
                    st_syn.select(station=station_code)[0].stats['delta']
                )

                # Synthetic model 2
                freqlist_syn2, fasZ_syn2 = fourierspec_cal_func(
                    st_syn2.select(station=station_code, channel='?NZ')[0].data,
                    st_syn2.select(station=station_code)[0].stats['delta']
                )
                freqlist_syn2, fas1_syn2 = fourierspec_cal_func(
                    st_syn2.select(station=station_code, channel='?NE')[0].data,
                    st_syn2.select(station=station_code)[0].stats['delta']
                )
                freqlist_syn2, fas2_syn2 = fourierspec_cal_func(
                    st_syn2.select(station=station_code, channel='?NN')[0].data,
                    st_syn2.select(station=station_code)[0].stats['delta']
                )

                # Apply Konno-Ohmachi smoothing - Observed
                fasZ_f_data = konno_ohmachi_func(fasZ_data, freqlist, bandwidth=bandwidth, normalize=True)
                fas_SA1_data = konno_ohmachi_func(fas1_data, freqlist, bandwidth=bandwidth, normalize=True)
                fas_SA2_data = konno_ohmachi_func(fas2_data, freqlist, bandwidth=bandwidth, normalize=True)

                # Apply Konno-Ohmachi smoothing - Model 1
                fasZ_f_syn1 = konno_ohmachi_func(fasZ_syn1, freqlist_syn1, bandwidth=bandwidth, normalize=True)
                fas_SA1_syn1 = konno_ohmachi_func(fas1_syn1, freqlist_syn1, bandwidth=bandwidth, normalize=True)
                fas_SA2_syn1 = konno_ohmachi_func(fas2_syn1, freqlist_syn1, bandwidth=bandwidth, normalize=True)

                # Apply Konno-Ohmachi smoothing - Model 2
                fasZ_f_syn2 = konno_ohmachi_func(fasZ_syn2, freqlist_syn2, bandwidth=bandwidth, normalize=True)
                fas_SA1_syn2 = konno_ohmachi_func(fas1_syn2, freqlist_syn2, bandwidth=bandwidth, normalize=True)
                fas_SA2_syn2 = konno_ohmachi_func(fas2_syn2, freqlist_syn2, bandwidth=bandwidth, normalize=True)

                # Compute horizontal component (geometric mean)
                hf_data = np.sqrt(0.5 * (fas_SA1_data**2 + fas_SA2_data**2))
                hf_syn1 = np.sqrt(0.5 * (fas_SA1_syn1**2 + fas_SA2_syn1**2))
                hf_syn2 = np.sqrt(0.5 * (fas_SA1_syn2**2 + fas_SA2_syn2**2))

                # Compute HVSR
                hvsr_data = hf_data / fasZ_f_data
                hvsr_syn1 = hf_syn1 / fasZ_f_syn1
                hvsr_syn2 = hf_syn2 / fasZ_f_syn2

                # Find frequency range for statistics
                freq_mask = (freqlist >= freq_range[0]) & (freqlist <= freq_range[1])
                freq_mask_syn1 = (freqlist_syn1 >= freq_range[0]) & (freqlist_syn1 <= freq_range[1])
                freq_mask_syn2 = (freqlist_syn2 >= freq_range[0]) & (freqlist_syn2 <= freq_range[1])

                # Find peak - Observed
                peak_idx_data = np.argmax(hvsr_data[freq_mask])
                peak_freq_data = freqlist[freq_mask][peak_idx_data]
                peak_amplitude_data = hvsr_data[freq_mask][peak_idx_data]
                mean_hvsr_data = np.mean(hvsr_data[freq_mask])

                # Find peak - Model 1
                peak_idx_syn1 = np.argmax(hvsr_syn1[freq_mask_syn1])
                peak_freq_syn1 = freqlist_syn1[freq_mask_syn1][peak_idx_syn1]
                peak_amplitude_syn1 = hvsr_syn1[freq_mask_syn1][peak_idx_syn1]
                mean_hvsr_syn1 = np.mean(hvsr_syn1[freq_mask_syn1])

                # Find peak - Model 2
                peak_idx_syn2 = np.argmax(hvsr_syn2[freq_mask_syn2])
                peak_freq_syn2 = freqlist_syn2[freq_mask_syn2][peak_idx_syn2]
                peak_amplitude_syn2 = hvsr_syn2[freq_mask_syn2][peak_idx_syn2]
                mean_hvsr_syn2 = np.mean(hvsr_syn2[freq_mask_syn2])

                if verbose:
                    print(f"  Peak - Data: {peak_amplitude_data:.2f} at {peak_freq_data:.2f} Hz")
                    print(f"  Peak - {model1_label}: {peak_amplitude_syn1:.2f} at {peak_freq_syn1:.2f} Hz")
                    print(f"  Peak - {model2_label}: {peak_amplitude_syn2:.2f} at {peak_freq_syn2:.2f} Hz")

                # Determine y-axis limits based on max amplitude in frequency range (all three curves)
                max_amp_in_range = max(peak_amplitude_data, peak_amplitude_syn1, peak_amplitude_syn2)
                ylim_max = max_amp_in_range * 1.15  # Add 15% margin above max

                # Create plot with three curves
                fig = plt.figure(figsize=figsize)
                ax = plt.subplot(111)

                ax.plot(freqlist, hvsr_data, 'black', linewidth=1.5, label='Data')
                ax.plot(freqlist_syn1, hvsr_syn1, model1_color, linewidth=1.5, label=model1_label)
                ax.plot(freqlist_syn2, hvsr_syn2, model2_color, linewidth=1.5, label=model2_label)

                ax.legend(loc='upper right', fontsize=8)
                ax.set_title('sta:' + station_code)

                ax.set_xlabel('frequency (Hz)')
                ax.set_ylabel('HV spectral amplitude ratio')
                ax.grid(linestyle=':', which='both')
                ax.set_ylim(0, ylim_max)  # Auto-scale based on data
                ax.set_xlim(freq_range)

                # Save plot
                plot_file = os.path.join(output_folder, f'HVSR-compare-{station_code}.png')
                plt.savefig(plot_file, dpi=dpi, transparent=False, bbox_inches='tight')
                plt.close()

                if verbose:
                    print(f"  Saved: {plot_file}")

                # Store results
                results.append({
                    'station': station_code,
                    'peak_freq_data': peak_freq_data,
                    'peak_freq_syn1': peak_freq_syn1,
                    'peak_freq_syn2': peak_freq_syn2,
                    'peak_amplitude_data': peak_amplitude_data,
                    'peak_amplitude_syn1': peak_amplitude_syn1,
                    'peak_amplitude_syn2': peak_amplitude_syn2,
                    'mean_hvsr_data': mean_hvsr_data,
                    'mean_hvsr_syn1': mean_hvsr_syn1,
                    'mean_hvsr_syn2': mean_hvsr_syn2,
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

    return pd.DataFrame(results)


# Example usage
if __name__ == "__main__":
    print("Spectral Ratio Analysis Functions")
    print("=" * 70)
    print("\nAvailable functions:")
    print("  1. compute_site_to_reference_ratio() - Site-to-reference spectral ratios")
    print("  2. compute_hvsr() - Horizontal-to-vertical spectral ratios")
    print("  3. compute_hvsr_two_models() - HVSR comparison with two synthetic models")
    print("\nImport with:")
    print("  from GMfunc.spectral_ratios import compute_site_to_reference_ratio, compute_hvsr, compute_hvsr_two_models")
