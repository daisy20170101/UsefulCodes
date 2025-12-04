import numpy as np
import seissolxdmf
import matplotlib.pyplot as plt
import pandas as pd


def extract_sr_vr_t0_from_table(model_prefix: str,
                                point_indices: np.ndarray,
                                ndt: int,
                                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract strike-slip rate (SRs) and dip-slip rate (SRd) at each time step from XDMF file.

    Parameters:
    -----------
    xdmfFilename : str
        Path to the XDMF file containing the seismic data

    Returns:
    --------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        - times: Array of time values for each time step
        - SRs_data: Array of strike-slip rates at each time step (shape: n_timesteps x n_elements)
        - SRd_data: Array of dip-slip rates at each time step (shape: n_timesteps x n_elements)

    Example:
    --------
    >>> times, SRs, SRd = extract_slip_rates('output-fault.xdmf')
    >>> print(f"Number of time steps: {len(times)}")
    >>> print(f"SRs shape: {SRs.shape}")
    """
    
    
    times = np.arange(0,0.5*ndt,0.5)

    # Initialize arrays to store data
    t0_data = np.zeros((ndt, len(point_indices)))
    sr_data = np.zeros((ndt, len(point_indices)))
    vr_data = np.zeros((ndt, len(point_indices)))
    rt_data = np.zeros((ndt, len(point_indices)))

    # Extract data at each time step
    for istep in range(ndt):

        df_series = pd.read_csv(model_prefix + str(istep) +'.csv')

        df_series_select = df_series.iloc[point_indices]

        t0_data[istep,:] = df_series_select['T0']
        sr_data[istep,:] = df_series_select['sr']
        vr_data[istep,:] = df_series_select['Vr']
        rt_data[istep,:] = df_series_select['RT']


    return times, sr_data, t0_data, vr_data,rt_data



def extract_slip_rates(xdmfFilename: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract strike-slip rate (SRs) and dip-slip rate (SRd) at each time step from XDMF file.

    Parameters:
    -----------
    xdmfFilename : str
        Path to the XDMF file containing the seismic data

    Returns:
    --------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        - times: Array of time values for each time step
        - SRs_data: Array of strike-slip rates at each time step (shape: n_timesteps x n_elements)
        - SRd_data: Array of dip-slip rates at each time step (shape: n_timesteps x n_elements)

    Example:
    --------
    >>> times, SRs, SRd = extract_slip_rates('output-fault.xdmf')
    >>> print(f"Number of time steps: {len(times)}")
    >>> print(f"SRs shape: {SRs.shape}")
    """

    # Initialize seissolxdmf reader
    sx = seissolxdmf.seissolxdmf(xdmfFilename)

    # Get number of time steps
    ndt = sx.ReadNdt()
    print(f"Number of time steps: {ndt}")

    # Get geometry information
    surfxyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    print(f"Geometry shape: {surfxyz.shape}, Connectivity shape: {connect.shape}")

    # Get time values for all time steps
    times = np.array(sx.ReadTimes())
    print(f"Time range: {times[0]:.3f} to {times[-1]:.3f} seconds")

    # Initialize arrays to store data
    SRs_data = []
    SRd_data = []

    # Extract data at each time step
    for i in range(ndt):
        # Read SRs (strike-slip rate) at this time step
        SRs = sx.ReadData('SRs', idt=i)
        SRs_data.append(SRs)

        # Read SRd (dip-slip rate) at this time step
        SRd = sx.ReadData('SRd', idt=i)
        SRd_data.append(SRd)

        if (i + 1) % 100 == 0:
            print(f"Processed {i + 1}/{ndt} time steps")

    # Convert lists to numpy arrays
    SRs_data = np.array(SRs_data)
    SRd_data = np.array(SRd_data)

    print(f"Extraction complete!")
    print(f"Times shape: {times.shape}")
    print(f"SRs data shape: {SRs_data.shape}")
    print(f"SRd data shape: {SRd_data.shape}")

    return times, SRs_data, SRd_data


def extract_slip_rates_at_points(xdmfFilename: str, point_indices: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract strike-slip rate (SRs) and dip-slip rate (SRd) at specific points for all time steps.

    Parameters:
    -----------
    xdmfFilename : str
        Path to the XDMF file containing the seismic data
    point_indices : np.ndarray
        Array of element indices to extract data for

    Returns:
    --------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        - times: Array of time values for each time step
        - SRs_data: Array of strike-slip rates at selected points (shape: n_timesteps x n_points)
        - SRd_data: Array of dip-slip rates at selected points (shape: n_timesteps x n_points)

    Example:
    --------
    >>> indices = np.array([0, 100, 500, 1000])
    >>> times, SRs, SRd = extract_slip_rates_at_points('output-fault.xdmf', indices)
    """

    # Initialize seissolxdmf reader
    sx = seissolxdmf.seissolxdmf(xdmfFilename)

    # Get number of time steps
    ndt = sx.ReadNdt()
    print(f"Number of time steps: {ndt}")
    print(f"Number of points to extract: {len(point_indices)}")

    # Get time values for all time steps
    times = np.array(sx.ReadTimes())
    print(f"Time range: {times[0]:.3f} to {times[-1]:.3f} seconds")

    # Initialize arrays to store data
    SRs_data = np.zeros((ndt, len(point_indices)))
    SRd_data = np.zeros((ndt, len(point_indices)))

    # Extract data at each time step
    for i in range(ndt):
        # Read SRs and SRd at this time step
        SRs_full = sx.ReadData('SRs', idt=i)
        SRd_full = sx.ReadData('SRd', idt=i)

        # Extract only the selected points
        SRs_data[i, :] = SRs_full[point_indices]
        SRd_data[i, :] = SRd_full[point_indices]

        if (i + 1) % 100 == 0:
            print(f"Processed {i + 1}/{ndt} time steps")

    print(f"Extraction complete!")
    print(f"Times shape: {times.shape}")
    print(f"SRs data shape: {SRs_data.shape}")
    print(f"SRd data shape: {SRd_data.shape}")

    return times, SRs_data, SRd_data



def extract_traction_three_at_points(xdmfFilename: str, point_indices: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract strike-tranctiona normal at specific points for all time steps.

    Parameters:
    -----------
    xdmfFilename : str
        Path to the XDMF file containing the seismic data
    point_indices : np.ndarray
        Array of element indices to extract data for

    Returns:
    --------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        - times: Array of time values for each time step
        - SRs_data: Array of strike-slip rates at selected points (shape: n_timesteps x n_points)
        - SRd_data: Array of dip-slip rates at selected points (shape: n_timesteps x n_points)

    Example:
    --------
    >>> indices = np.array([0, 100, 500, 1000])
    >>> times, SRs, SRd = extract_slip_rates_at_points('output-fault.xdmf', indices)
    """

    # Initialize seissolxdmf reader
    sx = seissolxdmf.seissolxdmf(xdmfFilename)

    # Get number of time steps
    ndt = sx.ReadNdt()
    print(f"Number of time steps: {ndt}")
    print(f"Number of points to extract: {len(point_indices)}")

    # Get time values for all time steps
    times = np.array(sx.ReadTimes())
    print(f"Time range: {times[0]:.3f} to {times[-1]:.3f} seconds")

    # Initialize arrays to store data
    SRs_data = np.zeros((ndt, len(point_indices)))
    SRd_data = np.zeros((ndt, len(point_indices)))
    Pn0_data = np.zeros((ndt, len(point_indices)))

    # Extract data at each time step
    for i in range(ndt):
        # Read SRs and SRd at this time step
        SRs_full = sx.ReadData('Ts0', idt=i)
        SRd_full = sx.ReadData('Td0', idt=i)
        Pn0_full = sx.ReadData('Pn0', idt=i)

        # Extract only the selected points
        SRs_data[i, :] = SRs_full[point_indices]
        SRd_data[i, :] = SRd_full[point_indices]
        Pn0_data[i, :] = Pn0_full[point_indices]

        if (i + 1) % 100 == 0:
            print(f"Processed {i + 1}/{ndt} time steps")

    print(f"Extraction complete!")
    print(f"Times shape: {times.shape}")
    print(f"SRs data shape: {SRs_data.shape}")
    print(f"SRd data shape: {SRd_data.shape}")

    return times, SRs_data,SRd_data,Pn0_data



def plot_slip_rate_timeseries(
    times: np.ndarray,
    SRs_data: np.ndarray = None,
    SRd_data: np.ndarray = None,
    point_indices: list = None,
    labels: list = None,
    figsize: tuple = (6, 5),
    save_path: str = None,
    title: str = None
) -> tuple:
    """
    Plot time series of strike-slip rate (SRs) and/or dip-slip rate (SRd).

    Parameters:
    -----------
    times : np.ndarray
        Array of time values
    SRs_data : np.ndarray, optional
        Strike-slip rate data (shape: n_timesteps x n_points)
    SRd_data : np.ndarray, optional
        Dip-slip rate data (shape: n_timesteps x n_points)
    point_indices : list, optional
        List of point indices to plot. If None, plots all points (or first 5 if >5)
    labels : list, optional
        List of labels for each point. If None, uses "Point 0", "Point 1", etc.
    figsize : tuple, optional
        Figure size (width, height). Default: (10, 6)
    save_path : str, optional
        Path to save the figure. If None, figure is not saved
    title : str, optional
        Plot title. If None, auto-generates based on data type

    Returns:
    --------
    tuple
        (fig, axes) - matplotlib figure and axes objects

    Example:
    --------
    >>> times, SRs, SRd = extract_slip_rates_at_points('output.xdmf', np.array([0, 100, 500]))
    >>> fig, axes = plot_slip_rate_timeseries(times, SRs, SRd,
    ...                                        labels=['Point A', 'Point B', 'Point C'])
    """

    # Validate inputs
    if SRs_data is None and SRd_data is None:
        raise ValueError("At least one of SRs_data or SRd_data must be provided")

    # Determine which data to plot
    plot_SRs = SRs_data is not None
    plot_SRd = SRd_data is not None

    # Set up number of subplots
    n_plots = int(plot_SRs) + int(plot_SRd)

    # Create figure
    fig, axes = plt.subplots(n_plots, 1, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    # Determine points to plot
    if plot_SRs:
        n_points = SRs_data.shape[1]
    else:
        n_points = SRd_data.shape[1]

    if point_indices is None:
        # Plot all points if 5 or fewer, otherwise first 5
        if n_points <= 5:
            point_indices = list(range(n_points))
        else:
            point_indices = list(range(5))
            print(f"Plotting first 5 of {n_points} points")

    # Generate labels if not provided
    if labels is None:
        labels = [f"Point {i}" for i in point_indices]

    # Plot counter
    plot_idx = 0

    # Plot SRs if provided
    if plot_SRs:
        ax = axes[plot_idx]
        for i, point_idx in enumerate(point_indices):
            ax.plot(times, SRs_data[:, point_idx], label=labels[i], linewidth=1.5)

        ax.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Strike-Slip Rate (m/s)', fontsize=12, fontweight='bold')
        ax.set_title('Strike-Slip Rate (SRs) Time Series', fontsize=13, fontweight='bold')
        ax.legend(loc='best', frameon=True, shadow=True)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.tick_params(axis='both', which='major', labelsize=10)
        plot_idx += 1

    # Plot SRd if provided
    if plot_SRd:
        ax = axes[plot_idx]
        for i, point_idx in enumerate(point_indices):
            ax.plot(times, SRd_data[:, point_idx], label=labels[i], linewidth=1.5)

        ax.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Dip-Slip Rate (m/s)', fontsize=12, fontweight='bold')
        ax.set_title('Dip-Slip Rate (SRd) Time Series', fontsize=13, fontweight='bold')
        ax.legend(loc='best', frameon=True, shadow=True)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.tick_params(axis='both', which='major', labelsize=10)

    # Add overall title if provided
    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold', y=0.995)

    plt.tight_layout()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig, axes


def plot_slip_rate_comparison(
    times: np.ndarray,
    SRs_data: np.ndarray,
    SRd_data: np.ndarray,
    point_idx: int = 0,
    label: str = None,
    figsize: tuple = (6, 5),
    save_path: str = None
) -> tuple:
    """
    Plot SRs and SRd for a single point on the same axes for comparison.

    Parameters:
    -----------
    times : np.ndarray
        Array of time values
    SRs_data : np.ndarray
        Strike-slip rate data (shape: n_timesteps x n_points)
    SRd_data : np.ndarray
        Dip-slip rate data (shape: n_timesteps x n_points)
    point_idx : int, optional
        Index of the point to plot. Default: 0
    label : str, optional
        Label for the point. If None, uses "Point {point_idx}"
    figsize : tuple, optional
        Figure size (width, height). Default: (10, 4)
    save_path : str, optional
        Path to save the figure. If None, figure is not saved

    Returns:
    --------
    tuple
        (fig, ax) - matplotlib figure and axis objects

    Example:
    --------
    >>> times, SRs, SRd = extract_slip_rates_at_points('output.xdmf', np.array([100]))
    >>> fig, ax = plot_slip_rate_comparison(times, SRs, SRd, point_idx=0, label='Hypocenter')
    """

    if label is None:
        label = f"Point {point_idx}"

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Plot both SRs and SRd
    ax.plot(times, SRs_data[:, point_idx], label=f'SRs ({label})',
            linewidth=2, color='royalblue')
    ax.plot(times, SRd_data[:, point_idx], label=f'SRd ({label})',
            linewidth=2, color='tomato', linestyle='--')

    ax.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Slip Rate (m/s)', fontsize=12, fontweight='bold')
    ax.set_title(f'Slip Rate Comparison - {label}', fontsize=13, fontweight='bold')
    ax.legend(loc='best', frameon=True, shadow=True, fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=10)

    plt.tight_layout()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig, ax


def plot_slip_rate_comparison_all(
    times: np.ndarray,
    SRs_data: np.ndarray,
    SRd_data: np.ndarray,
    max_elements: int = None,
    label: str = 'slip rate(m/s)',
    figsize: tuple = (6, 5),
    save_path: str = None,
    alpha: float = 0.7
) -> tuple:
    """
    Plot SRs and SRd for all elements on the same plot for comparison.

    Parameters:
    -----------
    times : np.ndarray
        Array of time values
    SRs_data : np.ndarray
        Strike-slip rate data (shape: n_timesteps x n_elements)
    SRd_data : np.ndarray
        Dip-slip rate data (shape: n_timesteps x n_elements)
    max_elements : int, optional
        Maximum number of elements to plot. If None, plots all elements
    figsize : tuple, optional
        Figure size (width, height). Default: (12, 5)
    save_path : str, optional
        Path to save the figure. If None, figure is not saved
    alpha : float, optional
        Transparency for individual element lines. Default: 0.7

    Returns:
    --------
    tuple
        (fig, axes) - matplotlib figure and axes objects

    Example:
    --------
    >>> times, SRs, SRd = extract_slip_rates('output.xdmf')
    >>> fig, axes = plot_slip_rate_comparison_all(times, SRs, SRd, max_elements=100)
    """


    import matplotlib as mpl
    import numpy as np
    import matplotlib.colors as mcolors

    
    n_elements = SRs_data.shape[1]

    cmap1 =  mpl.colormaps['plasma'].resampled(n_elements)

    # Limit number of elements if specified
    if max_elements is not None and n_elements > max_elements:
        print(f"Limiting to first {max_elements} of {n_elements} elements")
        n_elements = max_elements
        SRs_data = SRs_data[:, :max_elements]
        SRd_data = SRd_data[:, :max_elements]

    fig, axes = plt.subplots(2, 1, figsize=figsize)

    # Plot SRs for all elements
    for elem_idx in range(n_elements):
        axes[0].plot(times, SRs_data[:, elem_idx], linewidth=0.8, c=cmap1(elem_idx/n_elements))

    axes[0].set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    axes[0].set_ylabel(f'{label}', fontsize=12, fontweight='bold')
    # axes[0].set_title(f'SRs - All Elements (n={n_elements})', fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3, linestyle='--')
    axes[0].tick_params(axis='both', which='major', labelsize=10)

    # Plot SRd for all elements
    for elem_idx in range(n_elements):
        axes[1].plot(times, SRd_data[:, elem_idx], linewidth=0.8, c=cmap1(elem_idx/n_elements))

    axes[1].set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    axes[1].set_ylabel(f'{label}', fontsize=12, fontweight='bold')
    # axes[1].set_title(f'SRd - All Elements (n={n_elements})', fontsize=13, fontweight='bold')
    axes[1].grid(True, alpha=0.3, linestyle='--')
    axes[1].tick_params(axis='both', which='major', labelsize=10)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig, axes


def plot_slip_rate_comparison_overlay(
    times: np.ndarray,
    SRs_data: np.ndarray,
    SRd_data: np.ndarray,
    max_elements: int = None,
    figsize: tuple = (6, 5),
    save_path: str = None,
    alpha: float = 0.5
) -> tuple:
    """
    Plot SRs and SRd for all elements overlaid on a single plot.

    Parameters:
    -----------
    times : np.ndarray
        Array of time values
    SRs_data : np.ndarray
        Strike-slip rate data (shape: n_timesteps x n_elements)
    SRd_data : np.ndarray
        Dip-slip rate data (shape: n_timesteps x n_elements)
    max_elements : int, optional
        Maximum number of elements to plot. If None, plots all elements
    figsize : tuple, optional
        Figure size (width, height). Default: (10, 6)
    save_path : str, optional
        Path to save the figure. If None, figure is not saved
    alpha : float, optional
        Transparency for individual element lines. Default: 0.5

    Returns:
    --------
    tuple
        (fig, ax) - matplotlib figure and axis objects

    Example:
    --------
    >>> times, SRs, SRd = extract_slip_rates('output.xdmf')
    >>> fig, ax = plot_slip_rate_comparison_overlay(times, SRs, SRd, max_elements=100)
    """

    n_elements = SRs_data.shape[1]
    
    # Limit number of elements if specified
    if max_elements is not None and n_elements > max_elements:
        print(f"Limiting to first {max_elements} of {n_elements} elements")
        n_elements = max_elements
        SRs_data = SRs_data[:, :max_elements]
        SRd_data = SRd_data[:, :max_elements]

    fig, ax = plt.subplots(figsize=figsize)

    # Plot SRs for all elements in blue
    for elem_idx in range(n_elements):
        if elem_idx == 0:
            ax.plot(times, SRs_data[:, elem_idx], color='royalblue', 
                   linewidth=0.8, alpha=alpha, label='SRs')
        else:
            ax.plot(times, SRs_data[:, elem_idx], color='royalblue', 
                   linewidth=0.8, alpha=alpha)

    # Plot SRd for all elements in red
    for elem_idx in range(n_elements):
        if elem_idx == 0:
            ax.plot(times, SRd_data[:, elem_idx], color='tomato', 
                   linewidth=0.8, alpha=alpha, label='SRd', linestyle='--')
        else:
            ax.plot(times, SRd_data[:, elem_idx], color='tomato', 
                   linewidth=0.8, alpha=alpha, linestyle='--')

    ax.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Slip Rate (m/s)', fontsize=12, fontweight='bold')
    ax.set_title(f'Slip Rate Comparison - All Elements (n={n_elements})', 
                fontsize=13, fontweight='bold')
    ax.legend(loc='best', frameon=True, shadow=True, fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=10)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    return fig, ax
