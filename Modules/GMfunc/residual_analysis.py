"""
Residual analysis functions for ground motion prediction equations (GMPE).

This module provides functions to compute and visualize residuals between
observed and predicted ground motions, including binned statistics and
residual plots.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple, Optional, List


def bin_mean_std(nbins: int, sa1: np.ndarray, rrup: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute binned mean and standard deviation in log10 space.

    Divides distance range into bins and computes mean and standard deviation
    of log10(ground motion) in each bin.

    Parameters
    ----------
    nbins : int
        Number of distance bins
    sa1 : np.ndarray
        Ground motion intensity values (linear space, not log)
    rrup : np.ndarray
        Rupture distances (km)

    Returns
    -------
    mean_bins : np.ndarray
        Mean of log10(sa1) in each bin (length: nbins-1)
    std_bins : np.ndarray
        Standard deviation of log10(sa1) in each bin (length: nbins-1)
    rjb_bins : np.ndarray
        Distance bin edges (left edge of each bin, length: nbins-1)

    Examples
    --------
    >>> # Compute binned statistics
    >>> mean_bin, std_bins, rhypo_bins = bin_mean_std(
    ...     nbins=20,
    ...     sa1=10**residuals,  # Convert from log to linear
    ...     rrup=distances
    ... )
    >>>
    >>> # Plot mean ± 1 sigma
    >>> plt.plot(rhypo_bins, mean_bin, '-', label='Mean')
    >>> plt.fill_between(
    ...     rhypo_bins,
    ...     mean_bin - std_bins,
    ...     mean_bin + std_bins,
    ...     alpha=0.25,
    ...     label='±1σ'
    ... )

    Notes
    -----
    - Empty bins are filled with NaN
    - Input sa1 should be in linear space (will be converted to log10 internally)
    - Returns mean and std in log10 space
    """
    nbins = int(nbins)
    rjb_bins = np.linspace(rrup.min(), rrup.max(), nbins)
    mean_bins = np.zeros(nbins - 1)
    std_bins = np.zeros(nbins - 1)

    for ik, rjb_bin in enumerate(rjb_bins[0:-1]):
        # Find data in current bin
        mask = (rrup > rjb_bin) & (rrup < rjb_bins[ik + 1])
        bin_data = sa1[mask]

        # Handle empty bins gracefully
        if len(bin_data) > 0:
            log_data = np.log10(bin_data)
            mean_bins[ik] = np.mean(log_data)
            std_bins[ik] = np.std(log_data) if len(bin_data) > 1 else 0.0
        else:
            # Empty bin - set to NaN so it won't be plotted
            mean_bins[ik] = np.nan
            std_bins[ik] = np.nan

    return mean_bins, std_bins, rjb_bins[0:-1]


def bin_mean_err(nbins: int, sa1: np.ndarray, rrup: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute binned mean and standard error in log10 space.

    Similar to bin_mean_std but returns standard error (std/sqrt(n)) instead
    of standard deviation.

    Parameters
    ----------
    nbins : int
        Number of distance bins
    sa1 : np.ndarray
        Ground motion intensity values (linear space)
    rrup : np.ndarray
        Rupture distances (km)

    Returns
    -------
    mean_bins : np.ndarray
        Mean of log10(sa1) in each bin
    err_bins : np.ndarray
        Standard error of log10(sa1) in each bin (std/n)
    rjb_bins : np.ndarray
        All distance bin edges (length: nbins)

    Examples
    --------
    >>> mean_bin, err_bins, rhypo_bins = bin_mean_err(20, data, distances)
    >>> plt.errorbar(rhypo_bins[:-1], mean_bin, yerr=err_bins, fmt='o')
    """
    nbins = int(nbins)
    rjb_bins = np.linspace(rrup.min(), rrup.max(), nbins)
    mean_bins = np.zeros(nbins - 1)
    std_bins = np.zeros(nbins - 1)

    for ik, rjb_bin in enumerate(rjb_bins[0:-1]):
        # Find data in current bin
        mask = (rrup > rjb_bin) & (rrup < rjb_bins[ik + 1])
        bin_data = sa1[mask]
        n_data = len(bin_data)

        # Handle empty bins gracefully
        if n_data > 0:
            log_data = np.log10(bin_data)
            mean_bins[ik] = np.mean(log_data)
            if n_data > 1:
                std_bins[ik] = np.std(log_data) / n_data  # Standard error
            else:
                std_bins[ik] = 0.0
        else:
            # Empty bin - set to NaN
            mean_bins[ik] = np.nan
            std_bins[ik] = np.nan

    return mean_bins, std_bins, rjb_bins


def plot_residuals_vs_distance(
    df: pd.DataFrame,
    residual_cols: List[str],
    distance_col: str = 'r_rup',
    nbins: int = 20,
    labels: Optional[List[str]] = None,
    colors: Optional[List[str]] = None,
    figsize: Tuple[float, float] = (6, 8),
    output_path: Optional[str] = None,
    show_plot: bool = True
) -> plt.Figure:
    """
    Plot multiple residual columns vs distance with binned statistics.

    Creates a subplot for each residual column, showing scatter points,
    binned mean, and ±1 sigma bands.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing residuals and distances
    residual_cols : list of str
        Column names containing residuals (in log10 space)
    distance_col : str, optional
        Column name for distance (default: 'r_rup')
    nbins : int, optional
        Number of distance bins (default: 20)
    labels : list of str, optional
        Labels for each subplot (default: use column names)
    colors : list of str, optional
        Colors for scatter and mean line (default: ['lightblue', 'royalblue'])
    figsize : tuple, optional
        Figure size (width, height) in inches (default: (6, 8))
    output_path : str, optional
        Path to save figure
    show_plot : bool, optional
        Display plot interactively (default: True)

    Returns
    -------
    matplotlib.figure.Figure
        Figure object

    Examples
    --------
    >>> # Plot SA at 1.0s and 3.0s
    >>> fig = plot_residuals_vs_distance(
    ...     df_all_resid,
    ...     residual_cols=['resid_pSA_1.0', 'resid_pSA_3.0'],
    ...     distance_col='r_rup',
    ...     labels=['SA(1.0s)', 'SA(3.0s)'],
    ...     output_path='residuals.png'
    ... )
    """
    n_plots = len(residual_cols)

    # Default labels and colors
    if labels is None:
        labels = residual_cols
    if colors is None:
        colors = ['lightblue', 'royalblue'] * n_plots

    # Create subplots
    fig, axes = plt.subplots(n_plots, 1, figsize=figsize, sharex=True)

    # Handle single subplot case
    if n_plots == 1:
        axes = [axes]

    for idx, (resid_col, label) in enumerate(zip(residual_cols, labels)):
        ax = axes[idx]

        # Convert residuals from log to linear space for bin_mean_std
        residuals_linear = 10**df[resid_col].values
        distances = df[distance_col].values

        # Compute binned statistics
        mean_bin, std_bins, rhypo_bins = bin_mean_std(nbins, residuals_linear, distances)

        # Scatter plot of raw residuals
        ax.plot(
            df[distance_col],
            df[resid_col],
            '.',
            alpha=0.5,
            color=colors[0],
            markersize=4
        )

        # Binned mean
        ax.plot(
            rhypo_bins,
            mean_bin,
            '-',
            color=colors[1],
            linewidth=2.5,
            label='Mean',
            zorder=5
        )

        # ±1 sigma band
        ax.fill_between(
            rhypo_bins,
            mean_bin - std_bins,
            mean_bin + std_bins,
            color=colors[1],
            alpha=0.25,
            label=r'$\pm 1\sigma$'
        )

        # Zero reference line
        ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)

        # Formatting
        ax.set_ylabel(f'Residual {label}\nlog10(Obs/Pred)')
        ax.legend(loc='best', framealpha=0.9)
        ax.grid(which='both', linestyle=':', alpha=0.5)

        # Only show x-label on bottom plot
        if idx == n_plots - 1:
            ax.set_xlabel(f'{distance_col} ')

    plt.tight_layout()

    # Save if requested
    if output_path:
        import os
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    # Show or close
    if show_plot:
        plt.show()
    else:
        plt.close()

    return fig


def compute_residual_statistics(
    df: pd.DataFrame,
    residual_col: str,
    distance_col: str = 'r_rup',
    distance_bins: Optional[np.ndarray] = None
) -> pd.DataFrame:
    """
    Compute residual statistics in distance bins.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with residuals and distances
    residual_col : str
        Column name for residuals (log10 space)
    distance_col : str, optional
        Column name for distance (default: 'r_rup')
    distance_bins : np.ndarray, optional
        Distance bin edges. If None, uses 20 equal bins.

    Returns
    -------
    pd.DataFrame
        Statistics table with columns:
        - distance_min, distance_max: bin edges
        - distance_center: bin center
        - count: number of observations
        - mean: mean residual
        - std: standard deviation
        - median: median residual
        - q25, q75: 25th and 75th percentiles

    Examples
    --------
    >>> stats = compute_residual_statistics(
    ...     df_all_resid,
    ...     'resid_pSA_1.0',
    ...     distance_col='r_rup'
    ... )
    >>> print(stats)
    """
    if distance_bins is None:
        distance_bins = np.linspace(
            df[distance_col].min(),
            df[distance_col].max(),
            21  # 20 bins
        )

    results = []

    for i in range(len(distance_bins) - 1):
        d_min = distance_bins[i]
        d_max = distance_bins[i + 1]

        # Filter data in bin
        mask = (df[distance_col] >= d_min) & (df[distance_col] < d_max)
        bin_data = df.loc[mask, residual_col].values

        if len(bin_data) > 0:
            results.append({
                'distance_min': d_min,
                'distance_max': d_max,
                'distance_center': (d_min + d_max) / 2,
                'count': len(bin_data),
                'mean': np.mean(bin_data),
                'std': np.std(bin_data) if len(bin_data) > 1 else 0,
                'median': np.median(bin_data),
                'q25': np.percentile(bin_data, 25),
                'q75': np.percentile(bin_data, 75)
            })

    return pd.DataFrame(results)


def plot_residual_comparison(
    df_list: List[pd.DataFrame],
    residual_col: str,
    distance_col: str = 'r_rup',
    labels: Optional[List[str]] = None,
    nbins: int = 20,
    figsize: Tuple[float, float] = (8, 5),
    output_path: Optional[str] = None,
    show_plot: bool = True
) -> plt.Figure:
    """
    Compare residuals from multiple models/datasets.

    Parameters
    ----------
    df_list : list of pd.DataFrame
        List of DataFrames to compare
    residual_col : str
        Column name for residuals (must exist in all DataFrames)
    distance_col : str, optional
        Distance column name (default: 'r_rup')
    labels : list of str, optional
        Labels for each dataset
    nbins : int, optional
        Number of bins for statistics (default: 20)
    figsize : tuple, optional
        Figure size (default: (8, 5))
    output_path : str, optional
        Path to save figure
    show_plot : bool, optional
        Display plot (default: True)

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> # Compare multiple events
    >>> fig = plot_residual_comparison(
    ...     [df_event1, df_event2, df_event3],
    ...     residual_col='resid_pSA_1.0',
    ...     labels=['Event 1', 'Event 2', 'Event 3']
    ... )
    """
    if labels is None:
        labels = [f'Dataset {i+1}' for i in range(len(df_list))]

    fig, ax = plt.subplots(figsize=figsize)

    colors = plt.cm.tab10(np.linspace(0, 1, len(df_list)))

    for idx, (df, label) in enumerate(zip(df_list, labels)):
        # Compute binned statistics
        residuals_linear = 10**df[residual_col].values
        distances = df[distance_col].values

        mean_bin, std_bins, rhypo_bins = bin_mean_std(nbins, residuals_linear, distances)

        # Plot mean
        ax.plot(rhypo_bins, mean_bin, '-', label=label, color=colors[idx], linewidth=2)

        # ±1 sigma band
        ax.fill_between(
            rhypo_bins,
            mean_bin - std_bins,
            mean_bin + std_bins,
            color=colors[idx],
            alpha=0.15
        )

    # Zero reference
    ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)

    # Formatting
    ax.set_xlabel(f'{distance_col} (km)')
    ax.set_ylabel(f'Residual log10(Obs/Pred)')
    ax.legend(loc='best')
    ax.grid(which='both', linestyle=':', alpha=0.5)

    plt.tight_layout()

    if output_path:
        import os
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    if show_plot:
        plt.show()
    else:
        plt.close()

    return fig


# Example usage
if __name__ == "__main__":
    print("Residual Analysis Functions")
    print("=" * 70)
    print("\nAvailable functions:")
    print("  1. bin_mean_std() - Compute binned mean and std dev")
    print("  2. bin_mean_err() - Compute binned mean and std error")
    print("  3. plot_residuals_vs_distance() - Plot multiple residual columns")
    print("  4. compute_residual_statistics() - Detailed statistics by distance bin")
    print("  5. plot_residual_comparison() - Compare multiple datasets")
    print("\nImport with:")
    print("  from GMfunc.residual_analysis import plot_residuals_vs_distance")
