"""
Compute ground motion intensity measures from waveforms.

This module provides functions to compute:
- Peak Ground Velocity (PGV) from acceleration time histories
- Spectral Acceleration (SA) at specified periods
Both for observed and synthetic data, with amplitude correction.

Author: Refactored from SeisDataProcess3.ipynb
Date: 2025-11-24
"""

import numpy as np
import pandas as pd
from obspy.signal.differentiate_and_integrate import integrate_cumtrapz


def compute_intensity_measures(
    st_data,
    st_syn,
    siteTable,
    periods,
    compute_psa_func,
    Acorr=1.0,
    sampling_rate=200.0,
    damping=0.05,
    taper_fraction=0.1,
    verbose=True
):
    """
    Compute PGV and SA at specified periods for all stations.

    This function processes acceleration time histories to compute:
    1. Peak Ground Velocity (PGV) for all three components (Z, E, N)
       - Computed by integrating acceleration using trapezoidal rule
       - PGV = max(|integral(acceleration)|)
    2. Spectral Acceleration (SA) at user-specified periods
       - Computed using SDOF oscillator response
       - Default damping = 5%

    Parameters
    ----------
    st_data : obspy.Stream
        Stream containing observed acceleration data.
        Expected channels: H?E (E), H?N (N), H?Z (Z)
    st_syn : obspy.Stream
        Stream containing synthetic acceleration data.
        Expected channels: ?NE (E), ?NN (N), ?NZ (Z)
    siteTable : pd.DataFrame
        DataFrame with station information. Must have 'sta' column.
    periods : array-like
        Periods (in seconds) at which to compute SA.
        Example: [1.0, 3.0, 0.3]
    compute_psa_func : callable
        Function to compute pseudo-spectral acceleration.
        Signature: compute_psa(data, dt, period, damping)
        Should return SA value in same units as input data.
    Acorr : float, optional
        Amplitude correction factor to apply to synthetic data.
        Default: 1.0 (no correction)
    sampling_rate : float, optional
        Sampling rate in Hz. Default: 200.0
    damping : float, optional
        Damping ratio for SDOF oscillator. Default: 0.05 (5%)
    taper_fraction : float, optional
        Fraction of trace to taper at edges. Default: 0.1
    verbose : bool, optional
        Print progress messages. Default: True

    Returns
    -------
    results : pd.DataFrame
        DataFrame with columns:
        - station: Station code
        - pgv_z_obs, pgv_e_obs, pgv_n_obs: Observed PGV (m/s) for Z, E, N
        - pgv_z_syn, pgv_e_syn, pgv_n_syn: Synthetic PGV (m/s) for Z, E, N
        - sa_{period}_z_obs, sa_{period}_e_obs, sa_{period}_n_obs: Observed SA
        - sa_{period}_z_syn, sa_{period}_e_syn, sa_{period}_n_syn: Synthetic SA

        Example column names for periods=[1.0, 3.0, 0.3]:
        'sa_1.0_z_obs', 'sa_3.0_e_syn', etc.

    Notes
    -----
    - Data streams are tapered and detrended before processing
    - Synthetic data is corrected by Acorr (divided) for amplitude normalization
    - PGV units: same as input acceleration integrated once (typically m/s if input is m/s²)
    - SA units: same as input acceleration (typically m/s²)
    - Stations that fail processing are skipped (logged if verbose=True)

    Examples
    --------
    >>> from SpecFunc.compute_psa import compute_psa
    >>> periods = np.array([1.0, 3.0, 0.3])
    >>> results = compute_intensity_measures(
    ...     st_data=st_data,
    ...     st_syn=st_syn,
    ...     siteTable=siteTable,
    ...     periods=periods,
    ...     compute_psa_func=compute_psa,
    ...     Acorr=1.5,
    ...     damping=0.05
    ... )
    >>> print(results[['station', 'pgv_e_obs', 'sa_1.0_e_obs']])
    """
    dt = 1.0 / sampling_rate
    results = []

    for id, sta in enumerate(siteTable['sta']):
        if verbose:
            print(f"Processing station {id}: {sta}")

        try:
            # Select and prepare traces
            st_data_sta = st_data.select(station=sta)
            st_syn_sta = st_syn.select(station=sta)

            # Taper all components
            for channel in ['??E', '??N', '??Z']:
                if st_data_sta.select(channel=channel):
                    st_data_sta.select(channel=channel).taper(max_percentage=taper_fraction)

            # Detrend
            st_data_sta.detrend('simple')

            # Get observed traces
            tr_data_z = st_data_sta.select(channel='?NZ')[0] if st_data_sta.select(channel='?NZ') else st_data_sta.select(channel='H?Z')[0]
            tr_data_e = st_data_sta.select(channel='?NE')[0] if st_data_sta.select(channel='?NE') else st_data_sta.select(channel='H?E')[0]
            tr_data_n = st_data_sta.select(channel='?NN')[0] if st_data_sta.select(channel='?NN') else st_data_sta.select(channel='H?N')[0]

            # Get synthetic traces
            tr_syn_z = st_syn_sta.select(channel='?NZ')[0]
            tr_syn_e = st_syn_sta.select(channel='?NE')[0]
            tr_syn_n = st_syn_sta.select(channel='?NN')[0]

            # Compute PGV (integrate acceleration to get velocity, then take max abs)
            pgv_z_obs = np.max(np.abs(integrate_cumtrapz(tr_data_z.data, dt)))
            pgv_e_obs = np.max(np.abs(integrate_cumtrapz(tr_data_e.data, dt)))
            pgv_n_obs = np.max(np.abs(integrate_cumtrapz(tr_data_n.data, dt)))

            pgv_z_syn = np.max(np.abs(integrate_cumtrapz(tr_syn_z.data, dt))) / Acorr
            pgv_e_syn = np.max(np.abs(integrate_cumtrapz(tr_syn_e.data, dt))) / Acorr
            pgv_n_syn = np.max(np.abs(integrate_cumtrapz(tr_syn_n.data, dt))) / Acorr

            # Initialize result dict for this station
            result = {
                'station': sta,
                'pgv_z_obs': pgv_z_obs,
                'pgv_e_obs': pgv_e_obs,
                'pgv_n_obs': pgv_n_obs,
                'pgv_z_syn': pgv_z_syn,
                'pgv_e_syn': pgv_e_syn,
                'pgv_n_syn': pgv_n_syn
            }

            # Compute SA at each period
            for period in periods:
                # Observed SA
                sa_z_obs = compute_psa_func(
                    tr_data_z.data, tr_data_z.stats['delta'],
                    period=period, damping=damping
                )
                sa_e_obs = compute_psa_func(
                    tr_data_e.data, tr_data_e.stats['delta'],
                    period=period, damping=damping
                )
                sa_n_obs = compute_psa_func(
                    tr_data_n.data, tr_data_n.stats['delta'],
                    period=period, damping=damping
                )

                # Synthetic SA (with amplitude correction)
                sa_z_syn = compute_psa_func(
                    tr_syn_z.data / Acorr, tr_syn_z.stats['delta'],
                    period=period, damping=damping
                )
                sa_e_syn = compute_psa_func(
                    tr_syn_e.data / Acorr, tr_syn_e.stats['delta'],
                    period=period, damping=damping
                )
                sa_n_syn = compute_psa_func(
                    tr_syn_n.data / Acorr, tr_syn_n.stats['delta'],
                    period=period, damping=damping
                )

                # Store with period in column name
                result[f'sa_{period}_z_obs'] = sa_z_obs
                result[f'sa_{period}_e_obs'] = sa_e_obs
                result[f'sa_{period}_n_obs'] = sa_n_obs
                result[f'sa_{period}_z_syn'] = sa_z_syn
                result[f'sa_{period}_e_syn'] = sa_e_syn
                result[f'sa_{period}_n_syn'] = sa_n_syn

            results.append(result)

        except Exception as e:
            if verbose:
                print(f"  WARNING: Failed to process station {sta}: {e}")
            continue

    return pd.DataFrame(results)


def save_intensity_measures(
    results_df,
    output_folder,
    eventname,
    periods,
    obs_file_suffix='Data',
    syn_file_suffix='',
    verbose=True
):
    """
    Save intensity measures to CSV files in the format expected by legacy code.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame returned by compute_intensity_measures()
    output_folder : str
        Directory to save CSV files
    eventname : str
        Event identifier (e.g., '20240805')
    periods : array-like
        Same periods used in compute_intensity_measures()
        Order matters: [period1, period2, period3] maps to columns in output
    obs_file_suffix : str, optional
        Suffix for observed data file. Default: 'Data'
    syn_file_suffix : str, optional
        Suffix for synthetic data file. Default: ''
    verbose : bool, optional
        Print file paths. Default: True

    Returns
    -------
    obs_file : str
        Path to observed data CSV file
    syn_file : str
        Path to synthetic data CSV file

    Notes
    -----
    Output CSV format (13 columns):
    site, PGV vertical, PGV horizontal 1, PGV horizontal 2,
    SA {period1} vertical, SA {period1} horizontal 1, SA {period1} horizontal 2,
    SA {period2} vertical, SA {period2} horizontal 1, SA {period2} horizontal 2,
    SA {period3} vertical, SA {period3} horizontal 1, SA {period3} horizontal 2

    Legacy mapping (for periods=[1.0, 3.0, 0.3]):
    - Column 1: PGV Z
    - Column 2: PGV E
    - Column 3: PGV N
    - Columns 4-6: SA 1.0s (Z, E, N)
    - Columns 7-9: SA 3.0s (Z, E, N)
    - Columns 10-12: SA 0.3s (Z, E, N)

    Examples
    --------
    >>> obs_file, syn_file = save_intensity_measures(
    ...     results_df=results,
    ...     output_folder='/path/to/output/',
    ...     eventname='20240805',
    ...     periods=[1.0, 3.0, 0.3]
    ... )
    """
    import os
    os.makedirs(output_folder, exist_ok=True)

    # Construct file paths
    obs_file = os.path.join(output_folder, f'intensity{obs_file_suffix}-m{eventname}.csv')
    syn_file = os.path.join(output_folder, f'intensity{syn_file_suffix}-m{eventname}.csv')

    # Verify required periods
    if len(periods) != 3:
        raise ValueError(f"Expected 3 periods, got {len(periods)}. Legacy format requires exactly 3 periods.")

    # Observed data
    with open(obs_file, 'w') as f:
        # Header
        f.write("site,PGV vertical,PGV horizontal 1,PGV horizontal 2,")
        f.write(f"SA {periods[0]} vertical,SA {periods[0]} horizontal 1,SA {periods[0]} horizontal 2,")
        f.write(f"SA {periods[1]} vertical,SA {periods[1]} horizontal 1,SA {periods[1]} horizontal 2,")
        f.write(f"SA {periods[2]} vertical,SA {periods[2]} horizontal 1,SA {periods[2]} horizontal 2\n")

        # Data rows
        for _, row in results_df.iterrows():
            f.write(f"{row['station']},")
            f.write(f"{row['pgv_z_obs']},{row['pgv_e_obs']},{row['pgv_n_obs']},")
            f.write(f"{row[f'sa_{periods[0]}_z_obs']},{row[f'sa_{periods[0]}_e_obs']},{row[f'sa_{periods[0]}_n_obs']},")
            f.write(f"{row[f'sa_{periods[1]}_z_obs']},{row[f'sa_{periods[1]}_e_obs']},{row[f'sa_{periods[1]}_n_obs']},")
            f.write(f"{row[f'sa_{periods[2]}_z_obs']},{row[f'sa_{periods[2]}_e_obs']},{row[f'sa_{periods[2]}_n_obs']}\n")

    # Synthetic data
    with open(syn_file, 'w') as f:
        # Header
        f.write("site,PGV vertical,PGV horizontal 1,PGV horizontal 2,")
        f.write(f"SA {periods[0]} vertical,SA {periods[0]} horizontal 1,SA {periods[0]} horizontal 2,")
        f.write(f"SA {periods[1]} vertical,SA {periods[1]} horizontal 1,SA {periods[1]} horizontal 2,")
        f.write(f"SA {periods[2]} vertical,SA {periods[2]} horizontal 1,SA {periods[2]} horizontal 2\n")

        # Data rows
        for _, row in results_df.iterrows():
            f.write(f"{row['station']},")
            f.write(f"{row['pgv_z_syn']},{row['pgv_e_syn']},{row['pgv_n_syn']},")
            f.write(f"{row[f'sa_{periods[0]}_z_syn']},{row[f'sa_{periods[0]}_e_syn']},{row[f'sa_{periods[0]}_n_syn']},")
            f.write(f"{row[f'sa_{periods[1]}_z_syn']},{row[f'sa_{periods[1]}_e_syn']},{row[f'sa_{periods[1]}_n_syn']},")
            f.write(f"{row[f'sa_{periods[2]}_z_syn']},{row[f'sa_{periods[2]}_e_syn']},{row[f'sa_{periods[2]}_n_syn']}\n")

    if verbose:
        print(f"Saved observed intensity measures to: {obs_file}")
        print(f"Saved synthetic intensity measures to: {syn_file}")

    return obs_file, syn_file


def compute_and_save_intensity_measures(
    st_data,
    st_syn,
    siteTable,
    periods,
    compute_psa_func,
    output_folder,
    eventname,
    Acorr=1.0,
    sampling_rate=200.0,
    damping=0.05,
    verbose=True
):
    """
    Convenience function to compute and save intensity measures in one step.

    This combines compute_intensity_measures() and save_intensity_measures()
    for the most common use case.

    Parameters
    ----------
    st_data : obspy.Stream
        Stream containing observed acceleration data
    st_syn : obspy.Stream
        Stream containing synthetic acceleration data
    siteTable : pd.DataFrame
        DataFrame with station information
    periods : array-like
        Periods at which to compute SA (must be length 3 for legacy format)
    compute_psa_func : callable
        Function to compute pseudo-spectral acceleration
    output_folder : str
        Directory to save CSV files
    eventname : str
        Event identifier
    Acorr : float, optional
        Amplitude correction for synthetic data. Default: 1.0
    sampling_rate : float, optional
        Sampling rate in Hz. Default: 200.0
    damping : float, optional
        Damping ratio. Default: 0.05
    verbose : bool, optional
        Print progress. Default: True

    Returns
    -------
    results_df : pd.DataFrame
        DataFrame with all computed intensity measures
    obs_file : str
        Path to observed data CSV
    syn_file : str
        Path to synthetic data CSV

    Examples
    --------
    >>> from SpecFunc.compute_psa import compute_psa
    >>> periods = [1.0, 3.0, 0.3]
    >>> results_df, obs_file, syn_file = compute_and_save_intensity_measures(
    ...     st_data=st_data,
    ...     st_syn=st_syn,
    ...     siteTable=siteTable,
    ...     periods=periods,
    ...     compute_psa_func=compute_psa,
    ...     output_folder='/path/to/output/',
    ...     eventname='20240805',
    ...     Acorr=1.5
    ... )
    """
    # Compute intensity measures
    results_df = compute_intensity_measures(
        st_data=st_data,
        st_syn=st_syn,
        siteTable=siteTable,
        periods=periods,
        compute_psa_func=compute_psa_func,
        Acorr=Acorr,
        sampling_rate=sampling_rate,
        damping=damping,
        verbose=verbose
    )

    # Save to CSV files
    obs_file, syn_file = save_intensity_measures(
        results_df=results_df,
        output_folder=output_folder,
        eventname=eventname,
        periods=periods,
        verbose=verbose
    )

    return results_df, obs_file, syn_file
