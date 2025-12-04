"""
Source time function and seismic moment plotting utilities.

This module provides functions for plotting source time functions (STF),
moment rates, and seismic moment evolution from earthquake simulations.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Tuple, List


def calculate_moment_rate(
    seismic_moment: np.ndarray,
    time: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate moment rate from seismic moment time series.

    Parameters:
    -----------
    seismic_moment : np.ndarray
        Seismic moment values (Nm)
    time : np.ndarray
        Time values (s)

    Returns:
    --------
    tuple
        (time_diff, moment_rate) - Time points and moment rate values
    """
    moment_rate = np.diff(seismic_moment) / np.diff(time)
    time_diff = time[1:]

    return time_diff, moment_rate


def calculate_magnitude(seismic_moment: float) -> float:
    """
    Calculate moment magnitude (Mw) from seismic moment.

    Uses Hanks & Kanamori (1979) formula:
    Mw = (2/3) * log10(M0) - 6.07
    where M0 is in Nm

    Parameters:
    -----------
    seismic_moment : float
        Final seismic moment in Nm

    Returns:
    --------
    float
        Moment magnitude (Mw)
    """
    return (2.0 / 3.0) * np.log10(seismic_moment) - 6.07


def load_energy_file(
    filepath: str,
    variable: str = 'seismic_moment'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load seismic moment or energy data from CSV file.

    Parameters:
    -----------
    filepath : str
        Path to the energy CSV file
    variable : str, optional
        Variable name to extract. Default: 'seismic_moment'

    Returns:
    --------
    tuple
        (time, values) - Time array and values for the requested variable
    """
    data = pd.read_csv(filepath)

    # Filter for the requested variable
    data_filtered = data[data['variable'] == variable]

    time = data_filtered['time'].to_numpy()
    values = data_filtered['measurement'].to_numpy()

    return time, values


def plot_stf_single(
    time: np.ndarray,
    moment_rate: np.ndarray,
    magnitude: float,
    label: str = None,
    ax: plt.Axes = None,
    figsize: Tuple[float, float] = (10, 6),
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    save_path: str = None
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot source time function (moment rate) for a single event.

    Parameters:
    -----------
    time : np.ndarray
        Time values (s)
    moment_rate : np.ndarray
        Moment rate values (Nm/s)
    magnitude : float
        Moment magnitude (Mw)
    label : str, optional
        Label for the plot
    ax : plt.Axes, optional
        Matplotlib axes to plot on. If None, creates new figure
    figsize : tuple, optional
        Figure size. Default: (10, 6)
    xlim : tuple, optional
        X-axis limits (xmin, xmax)
    ylim : tuple, optional
        Y-axis limits (ymin, ymax)
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    tuple
        (fig, ax) - matplotlib figure and axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Plot moment rate (normalized to 10^20 Nm/s)
    plot_label = label if label else f'Mw {magnitude:.2f}'
    ax.plot(time, moment_rate / 1e20, '-', linewidth=2, label=plot_label)

    ax.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Moment Rate (×10²⁰ Nm/s)', fontsize=12, fontweight='bold')
    ax.set_title('Source Time Function', fontsize=13, fontweight='bold')
    ax.legend(loc='best', frameon=True, shadow=True, fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=10)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig, ax


def plot_stf_loop(
    eventlist_file: str,
    energy_folder: str = 'Joint4/energy/',
    energy_suffix: str = '-energy.csv',
    figsize: Tuple[float, float] = (10, 7),
    xlim: Optional[Tuple[float, float]] = (0, 80),
    ylim: Optional[Tuple[float, float]] = (0, 20),
    m0scale: float=1e19,
    save_path: str = './slab-mag.png',
    event_column: str = 'event'
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot source time functions for multiple events from a list.

    Parameters:
    -----------
    eventlist_file : str
        Path to CSV file containing event list
    energy_folder : str, optional
        Folder containing energy files. Default: 'Joint4/energy/'
    energy_suffix : str, optional
        Suffix for energy files. Default: '-energy.csv'
    figsize : tuple, optional
        Figure size. Default: (10, 7)
    xlim : tuple, optional
        X-axis limits. Default: (0, 80)
    ylim : tuple, optional
        Y-axis limits. Default: (0, 20)
    save_path : str, optional
        Path to save figure. Default: './slab-mag.png'
    event_column : str, optional
        Column name for event IDs. Default: 'event'

    Returns:
    --------
    tuple
        (fig, ax) - matplotlib figure and axes objects

    Example:
    --------
    >>> fig, ax = plot_stf_loop('eventlist.csv', save_path='stf_comparison.png')
    """
    # Read event list
    eventfile = pd.read_csv(eventlist_file)
    print(f"Event file columns: {eventfile.keys()}")

    eventlist = eventfile[event_column]
    print(f"Processing {len(eventlist)} events")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Loop through events
    for ieve, eve in enumerate(eventlist):
        print(f"Processing event {ieve + 1}/{len(eventlist)}: {eve}")

        # Construct energy file path
        energy_file = f'{energy_folder}{eve}{energy_suffix}'

        try:
            # Load data
            time, seismic_moment = load_energy_file(energy_file, 'seismic_moment')

            # Calculate moment rate
            time_rate, moment_rate = calculate_moment_rate(seismic_moment, time)

            # Calculate magnitude
            magnitude = calculate_magnitude(seismic_moment[-1])

            print(f"  Event: {eve}, Mw = {magnitude:.2f}")

            # Plot
            ax.plot(time_rate, moment_rate / m0scale, '-',
                   label=f'{eve}, Mw {magnitude:.2f}',
                   linewidth=1.5)

        except FileNotFoundError:
            print(f"  Warning: Energy file not found: {energy_file}")
            continue
        except Exception as e:
            print(f"  Error processing event {eve}: {e}")
            continue

    # Format plot
    ax.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Moment Rate (×{m0scale} Nm/s)', fontsize=12, fontweight='bold')
    ax.set_title('Source Time Functions Comparison', fontsize=13, fontweight='bold')
    ax.legend(loc='best', frameon=True, shadow=True, fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=10)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    plt.tight_layout()

    # Save figure
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig, ax


def plot_seismic_moment_evolution(
    time: np.ndarray,
    seismic_moment: np.ndarray,
    magnitude: float,
    label: str = None,
    figsize: Tuple[float, float] = (10, 6),
    log_scale: bool = False,
    save_path: str = None
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot seismic moment evolution over time.

    Parameters:
    -----------
    time : np.ndarray
        Time values (s)
    seismic_moment : np.ndarray
        Seismic moment values (Nm)
    magnitude : float
        Final moment magnitude (Mw)
    label : str, optional
        Label for the plot
    figsize : tuple, optional
        Figure size. Default: (10, 6)
    log_scale : bool, optional
        Use logarithmic scale for y-axis. Default: False
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    tuple
        (fig, ax) - matplotlib figure and axes objects
    """
    fig, ax = plt.subplots(figsize=figsize)

    plot_label = label if label else f'Mw {magnitude:.2f}'
    ax.plot(time, seismic_moment / 1e18, '-', linewidth=2, label=plot_label)

    ax.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Seismic Moment (×10¹⁸ Nm)', fontsize=12, fontweight='bold')
    ax.set_title('Seismic Moment Evolution', fontsize=13, fontweight='bold')
    ax.legend(loc='best', frameon=True, shadow=True, fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=10)

    if log_scale:
        ax.set_yscale('log')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig, ax


def plot_stf_and_moment(
    time: np.ndarray,
    seismic_moment: np.ndarray,
    magnitude: float,
    label: str = None,
    figsize: Tuple[float, float] = (12, 5),
    save_path: str = None
) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Plot both source time function and seismic moment evolution side by side.

    Parameters:
    -----------
    time : np.ndarray
        Time values (s)
    seismic_moment : np.ndarray
        Seismic moment values (Nm)
    magnitude : float
        Moment magnitude (Mw)
    label : str, optional
        Label for the plots
    figsize : tuple, optional
        Figure size. Default: (12, 5)
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    tuple
        (fig, axes) - matplotlib figure and list of axes objects
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    plot_label = label if label else f'Mw {magnitude:.2f}'

    # Calculate moment rate
    time_rate, moment_rate = calculate_moment_rate(seismic_moment, time)

    # Plot moment rate (STF)
    axes[0].plot(time_rate, moment_rate / 1e20, '-', linewidth=2,
                 label=plot_label, color='royalblue')
    axes[0].set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Moment Rate (×10²⁰ Nm/s)', fontsize=12, fontweight='bold')
    axes[0].set_title('Source Time Function', fontsize=13, fontweight='bold')
    axes[0].legend(loc='best', frameon=True, shadow=True)
    axes[0].grid(True, alpha=0.3, linestyle='--')

    # Plot seismic moment
    axes[1].plot(time, seismic_moment / 1e18, '-', linewidth=2,
                 label=plot_label, color='tomato')
    axes[1].set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Seismic Moment (×10¹⁸ Nm)', fontsize=12, fontweight='bold')
    axes[1].set_title('Seismic Moment Evolution', fontsize=13, fontweight='bold')
    axes[1].legend(loc='best', frameon=True, shadow=True)
    axes[1].grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig, axes


