#!/usr/bin/env python3
"""
Power Spectral Decay Analysis
Calculates power spectral density and fits decay exponent for time series data from CSV
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit
from scipy.stats import linregress
import argparse
import os


def power_law(f, a, beta):
    """Power law function: S(f) = a * f^(-beta)"""
    return a * f ** (-beta)


def log_power_law(log_f, log_a, beta):
    """Linear form: log(S) = log(a) - beta * log(f)"""
    return log_a - beta * log_f


def calculate_psd(data, sampling_rate, method="welch", nperseg=None):
    """
    Calculate Power Spectral Density

    Parameters:
    -----------
    data : array_like
        Time series data
    sampling_rate : float
        Sampling frequency in Hz
    method : str
        Method for PSD calculation ('welch', 'periodogram', 'multitaper')
    nperseg : int
        Length of each segment for Welch's method

    Returns:
    --------
    frequencies : array
        Frequency bins
    psd : array
        Power spectral density values
    """

    if nperseg is None:
        nperseg = min(len(data) // 8, 1024)

    if method == "welch":
        frequencies, psd = signal.welch(
            data, fs=sampling_rate, nperseg=nperseg, window="hanning", overlap=0.5
        )
    elif method == "periodogram":
        frequencies, psd = signal.periodogram(data, fs=sampling_rate, window="hanning")
    elif method == "multitaper":
        frequencies, psd = signal.multitaper.multitaper(data, fs=sampling_rate)
    else:
        raise ValueError("Method must be 'welch', 'periodogram', or 'multitaper'")

    # Remove DC component
    frequencies = frequencies[1:]
    psd = psd[1:]

    return frequencies, psd


def fit_power_law_decay(frequencies, psd, fit_range=None, method="linear"):
    """
    Fit power law decay to PSD

    Parameters:
    -----------
    frequencies : array
        Frequency values
    psd : array
        Power spectral density values
    fit_range : tuple
        (f_min, f_max) frequency range for fitting
    method : str
        Fitting method ('linear' for log-log regression, 'nonlinear' for curve_fit)

    Returns:
    --------
    beta : float
        Decay exponent
    amplitude : float
        Amplitude coefficient
    r_squared : float
        Coefficient of determination
    fit_freq : array
        Frequencies used for fitting
    fit_psd : array
        Fitted PSD values
    """

    # Apply frequency range filter
    if fit_range is not None:
        mask = (frequencies >= fit_range[0]) & (frequencies <= fit_range[1])
        fit_freq = frequencies[mask]
        fit_psd_data = psd[mask]
    else:
        fit_freq = frequencies
        fit_psd_data = psd

    if method == "linear":
        # Linear regression in log-log space
        log_freq = np.log10(fit_freq)
        log_psd = np.log10(fit_psd_data)

        slope, intercept, r_value, p_value, std_err = linregress(log_freq, log_psd)

        beta = -slope  # Negative because we want positive decay exponent
        amplitude = 10**intercept
        r_squared = r_value**2

        # Generate fitted curve
        fit_psd = amplitude * fit_freq ** (-beta)

    elif method == "nonlinear":
        # Nonlinear curve fitting
        try:
            popt, pcov = curve_fit(
                power_law, fit_freq, fit_psd_data, p0=[fit_psd_data[0], 2.0]
            )
            amplitude, beta = popt

            # Calculate R-squared
            fit_psd = power_law(fit_freq, amplitude, beta)
            ss_res = np.sum((fit_psd_data - fit_psd) ** 2)
            ss_tot = np.sum((fit_psd_data - np.mean(fit_psd_data)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)

        except RuntimeError:
            print("Nonlinear fitting failed, using linear method")
            return fit_power_law_decay(frequencies, psd, fit_range, "linear")

    else:
        raise ValueError("Method must be 'linear' or 'nonlinear'")

    return beta, amplitude, r_squared, fit_freq, fit_psd


def analyze_spectral_decay(
    csv_file,
    column_name=None,
    sampling_rate=1.0,
    psd_method="welch",
    fit_range=None,
    fit_method="linear",
    plot=True,
    save_results=True,
):
    """
    Complete spectral decay analysis

    Parameters:
    -----------
    csv_file : str
        Path to CSV file
    column_name : str
        Name of column to analyze (if None, uses first numeric column)
    sampling_rate : float
        Sampling frequency in Hz
    psd_method : str
        Method for PSD calculation
    fit_range : tuple
        Frequency range for power law fitting
    fit_method : str
        Method for power law fitting
    plot : bool
        Whether to create plots
    save_results : bool
        Whether to save results to file

    Returns:
    --------
    results : dict
        Analysis results
    """

    # Load data
    print(f"Loading data from {csv_file}")
    df = pd.read_csv(csv_file)

    if column_name is None:
        # Find first numeric column
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) == 0:
            raise ValueError("No numeric columns found in CSV file")
        column_name = numeric_cols[0]
        print(f"Using column: {column_name}")

    data = df[column_name].dropna().values
    print(f"Data length: {len(data)} points")

    # Calculate PSD
    print("Calculating power spectral density...")
    frequencies, psd = calculate_psd(data, sampling_rate, method=psd_method)

    # Fit power law decay
    print("Fitting power law decay...")
    beta, amplitude, r_squared, fit_freq, fit_psd = fit_power_law_decay(
        frequencies, psd, fit_range, fit_method
    )

    # Prepare results
    results = {
        "decay_exponent": beta,
        "amplitude": amplitude,
        "r_squared": r_squared,
        "frequencies": frequencies,
        "psd": psd,
        "fit_frequencies": fit_freq,
        "fit_psd": fit_psd,
        "data_length": len(data),
        "sampling_rate": sampling_rate,
        "psd_method": psd_method,
        "fit_method": fit_method,
        "fit_range": fit_range,
    }

    # Print results
    print("\n" + "=" * 50)
    print("POWER SPECTRAL DECAY ANALYSIS RESULTS")
    print("=" * 50)
    print(f"Decay exponent (β): {beta:.3f}")
    print(f"Amplitude: {amplitude:.2e}")
    print(f"R-squared: {r_squared:.4f}")
    print(f"Power law: S(f) = {amplitude:.2e} × f^(-{beta:.3f})")

    if fit_range:
        print(f"Fitted frequency range: {fit_range[0]:.3f} - {fit_range[1]:.3f} Hz")
    else:
        print(
            f"Fitted frequency range: {frequencies[0]:.3f} - {frequencies[-1]:.3f} Hz"
        )

    # Create plots
    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Time series plot
        time = np.arange(len(data)) / sampling_rate
        ax1.plot(time, data, "b-", alpha=0.7)
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Amplitude")
        ax1.set_title(f"Time Series: {column_name}")
        ax1.grid(True, alpha=0.3)

        # Power spectral density plot
        ax2.loglog(frequencies, psd, "b-", alpha=0.7, label="Data PSD")
        ax2.loglog(
            fit_freq,
            fit_psd,
            "r-",
            linewidth=2,
            label=f"Power Law Fit: f^(-{beta:.2f})",
        )
        ax2.set_xlabel("Frequency (Hz)")
        ax2.set_ylabel("Power Spectral Density")
        ax2.set_title(f"Power Spectral Decay Analysis (R² = {r_squared:.4f})")
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plt.tight_layout()

        if save_results:
            plot_file = csv_file.replace(".csv", "_spectral_analysis.png")
            plt.savefig(plot_file, dpi=300, bbox_inches="tight")
            print(f"Plot saved as: {plot_file}")

        plt.show()

    # Save numerical results
    if save_results:
        results_file = csv_file.replace(".csv", "_spectral_results.txt")
        with open(results_file, "w") as f:
            f.write("POWER SPECTRAL DECAY ANALYSIS RESULTS\n")
            f.write("=" * 50 + "\n")
            f.write(f"Input file: {csv_file}\n")
            f.write(f"Column analyzed: {column_name}\n")
            f.write(f"Data length: {len(data)} points\n")
            f.write(f"Sampling rate: {sampling_rate} Hz\n")
            f.write(f"PSD method: {psd_method}\n")
            f.write(f"Fit method: {fit_method}\n")
            if fit_range:
                f.write(f"Fit range: {fit_range[0]:.3f} - {fit_range[1]:.3f} Hz\n")
            f.write(f"\nDecay exponent (β): {beta:.6f}\n")
            f.write(f"Amplitude: {amplitude:.6e}\n")
            f.write(f"R-squared: {r_squared:.6f}\n")
            f.write(f"Power law: S(f) = {amplitude:.6e} × f^(-{beta:.6f})\n")

        print(f"Results saved as: {results_file}")

    return results


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(description="Power Spectral Decay Analysis")
    parser.add_argument("csv_file", help="Path to CSV file")
    parser.add_argument("-c", "--column", help="Column name to analyze")
    parser.add_argument(
        "-sr",
        "--sampling_rate",
        type=float,
        default=1.0,
        help="Sampling rate in Hz (default: 1.0)",
    )
    parser.add_argument(
        "-m",
        "--method",
        choices=["welch", "periodogram", "multitaper"],
        default="welch",
        help="PSD calculation method",
    )
    parser.add_argument(
        "-fr",
        "--fit_range",
        nargs=2,
        type=float,
        help="Frequency range for fitting (f_min f_max)",
    )
    parser.add_argument(
        "-fm",
        "--fit_method",
        choices=["linear", "nonlinear"],
        default="linear",
        help="Power law fitting method",
    )
    parser.add_argument("--no_plot", action="store_true", help="Disable plotting")
    parser.add_argument("--no_save", action="store_true", help="Disable saving results")

    args = parser.parse_args()

    # Check if file exists
    if not os.path.exists(args.csv_file):
        print(f"Error: File {args.csv_file} not found")
        return

    # Run analysis
    try:
        results = analyze_spectral_decay(
            csv_file=args.csv_file,
            column_name=args.column,
            sampling_rate=args.sampling_rate,
            psd_method=args.method,
            fit_range=tuple(args.fit_range) if args.fit_range else None,
            fit_method=args.fit_method,
            plot=not args.no_plot,
            save_results=not args.no_save,
        )

    except Exception as e:
        print(f"Error during analysis: {e}")


if __name__ == "__main__":
    main()

# Example usage:
"""
# Basic usage
python power_spectral_decay.py data.csv

# Advanced usage
python power_spectral_decay.py data.csv -c "acceleration" -sr 100 -fr 1 50 -m welch -fm linear

# In Python script:
results = analyze_spectral_decay('data.csv', 
                               column_name='signal', 
                               sampling_rate=100,
                               fit_range=(1, 50))
print(f"Decay exponent: {results['decay_exponent']}")
"""
