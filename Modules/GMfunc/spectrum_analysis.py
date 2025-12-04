"""
Seismic signal spectrum analysis utilities.

This module provides functions to compute and visualize 2D spectrograms
of seismic waveforms using Short-Time Fourier Transform (STFT).
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from obspy import Stream
from typing import Optional, Tuple
import os


def compute_spectrogram(
    data: np.ndarray,
    sampling_rate: float,
    nperseg: int = 256,
    noverlap: int = 200,
    window: str = 'hann'
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute 2D spectrogram using Short-Time Fourier Transform (STFT).

    Parameters
    ----------
    data : np.ndarray
        Time series data (e.g., acceleration, velocity)
    sampling_rate : float
        Sampling rate in Hz (samples per second)
    nperseg : int, optional
        Length of each segment for STFT (default: 256)
        Larger values give better frequency resolution but worse time resolution
    noverlap : int, optional
        Number of points to overlap between segments (default: 200)
        Should be less than nperseg
    window : str, optional
        Window function to use (default: 'hann')
        Options: 'hann', 'hamming', 'blackman', 'boxcar', etc.

    Returns
    -------
    f : np.ndarray
        Array of frequency values (Hz)
    t : np.ndarray
        Array of time values (s)
    Zxx : np.ndarray
        Complex STFT matrix (frequency × time)
        Use np.abs(Zxx) for magnitude spectrogram

    Examples
    --------
    >>> import numpy as np
    >>> # Generate test signal
    >>> fs = 200.0
    >>> t = np.linspace(0, 10, int(10*fs))
    >>> data = np.sin(2*np.pi*5*t) + 0.5*np.sin(2*np.pi*10*t)
    >>>
    >>> # Compute spectrogram
    >>> f, t, Zxx = compute_spectrogram(data, fs)
    >>> print(f"Frequency range: {f[0]:.2f} - {f[-1]:.2f} Hz")
    >>> print(f"Time range: {t[0]:.2f} - {t[-1]:.2f} s")

    Notes
    -----
    The STFT represents a signal in the time-frequency domain by computing
    discrete Fourier transforms (DFT) over short overlapping windows.

    Time resolution: Better with smaller nperseg
    Frequency resolution: Better with larger nperseg
    """
    f, t, Zxx = signal.stft(
        data,
        fs=sampling_rate,
        nperseg=nperseg,
        noverlap=noverlap,
        window=window
    )

    return f, t, Zxx


def plot_spectrogram(
    f: np.ndarray,
    t: np.ndarray,
    Zxx: np.ndarray,
    title: str = "2D Spectrogram",
    freq_range: Optional[Tuple[float, float]] = None,
    time_range: Optional[Tuple[float, float]] = None,
    cmap: str = 'plasma',
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    log_scale: bool = False,
    output_path: Optional[str] = None,
    show_plot: bool = True,
    figsize: Tuple[float, float] = (10, 4),
    dpi: int = 300
) -> plt.Figure:
    """
    Plot 2D spectrogram.

    Parameters
    ----------
    f : np.ndarray
        Frequency array from compute_spectrogram()
    t : np.ndarray
        Time array from compute_spectrogram()
    Zxx : np.ndarray
        Complex STFT matrix from compute_spectrogram()
    title : str, optional
        Plot title (default: "2D Spectrogram")
    freq_range : tuple of float, optional
        (fmin, fmax) frequency range to display in Hz
        If None, shows all frequencies
    time_range : tuple of float, optional
        (tmin, tmax) time range to display in seconds
        If None, shows all times
    cmap : str, optional
        Matplotlib colormap name (default: 'plasma')
        Options: 'viridis', 'plasma', 'inferno', 'magma', 'hot', 'jet', etc.
    vmin : float, optional
        Minimum value for color scale
    vmax : float, optional
        Maximum value for color scale
    log_scale : bool, optional
        Use logarithmic color scale (default: False)
    output_path : str, optional
        Path to save figure. If None, figure not saved.
    show_plot : bool, optional
        Display plot interactively (default: True)
    figsize : tuple of float, optional
        Figure size (width, height) in inches (default: (10, 4))
    dpi : int, optional
        Resolution for saved figure (default: 300)

    Returns
    -------
    matplotlib.figure.Figure
        Figure object

    Examples
    --------
    >>> # After computing spectrogram
    >>> fig = plot_spectrogram(
    ...     f, t, Zxx,
    ...     title="Acceleration Spectrogram",
    ...     freq_range=(0, 10),
    ...     cmap='viridis',
    ...     output_path='spectrogram.png'
    ... )
    """
    # Compute magnitude
    magnitude = np.abs(Zxx)

    # Apply log scale if requested
    if log_scale:
        magnitude = 10 * np.log10(magnitude + 1e-10)  # Add small value to avoid log(0)

    # Create figure
    fig = plt.figure(figsize=figsize)

    # Plot spectrogram
    mesh = plt.pcolormesh(
        t, f, magnitude,
        shading='auto',
        cmap=cmap,
        vmin=vmin,
        vmax=vmax
    )

    # Set labels and title
    plt.xlabel("Time (s)", fontsize=12)
    plt.ylabel("Frequency (Hz)", fontsize=12)
    plt.title(title, fontsize=14)

    # Add colorbar
    cbar = plt.colorbar(mesh, label="Amplitude")

    # Set axis limits
    if freq_range is not None:
        plt.ylim(freq_range)

    if time_range is not None:
        plt.xlim(time_range)

    # Tight layout
    plt.tight_layout()

    # Save if path provided
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')

    # Show or close
    if show_plot:
        plt.show()
    else:
        plt.close()

    return fig


def compute_and_plot_spectrogram(
    trace_data: np.ndarray,
    sampling_rate: float,
    title: str = "2D Spectrogram of acceleration",
    freq_range: Tuple[float, float] = (0, 2.0),
    nperseg: int = 256,
    noverlap: int = 200,
    cmap: str = 'plasma',
    output_path: Optional[str] = None,
    show_plot: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, plt.Figure]:
    """
    Convenience function to compute and plot spectrogram in one call.

    Parameters
    ----------
    trace_data : np.ndarray
        Seismic trace data
    sampling_rate : float
        Sampling rate in Hz
    title : str, optional
        Plot title (default: "2D Spectrogram of acceleration")
    freq_range : tuple of float, optional
        (fmin, fmax) frequency range to display (default: (0, 2.0))
    nperseg : int, optional
        STFT segment length (default: 256)
    noverlap : int, optional
        STFT overlap (default: 200)
    cmap : str, optional
        Colormap (default: 'plasma')
    output_path : str, optional
        Path to save figure
    show_plot : bool, optional
        Display plot (default: True)

    Returns
    -------
    f : np.ndarray
        Frequency array
    t : np.ndarray
        Time array
    Zxx : np.ndarray
        STFT matrix
    fig : matplotlib.figure.Figure
        Figure object

    Examples
    --------
    >>> # From ObsPy stream
    >>> trace = st.select(station='VUWS', channel='HNN')[0]
    >>> f, t, Zxx, fig = compute_and_plot_spectrogram(
    ...     trace.data,
    ...     trace.stats.sampling_rate,
    ...     title=f"Spectrogram - {trace.stats.station}",
    ...     freq_range=(0, 10)
    ... )
    """
    # Compute STFT
    f, t, Zxx = compute_spectrogram(
        trace_data,
        sampling_rate,
        nperseg=nperseg,
        noverlap=noverlap
    )

    # Plot
    fig = plot_spectrogram(
        f, t, Zxx,
        title=title,
        freq_range=freq_range,
        cmap=cmap,
        output_path=output_path,
        show_plot=show_plot
    )

    return f, t, Zxx, fig


def analyze_stream_spectrogram(
    st: Stream,
    station: str,
    channel: str = '??N',
    title: Optional[str] = None,
    freq_range: Tuple[float, float] = (0, 2.0),
    nperseg: int = 256,
    noverlap: int = 200,
    output_path: Optional[str] = None,
    show_plot: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, plt.Figure]:
    """
    Compute and plot spectrogram directly from ObsPy Stream.

    Parameters
    ----------
    st : obspy.Stream
        ObsPy stream object
    station : str
        Station code to analyze
    channel : str, optional
        Channel pattern (default: '??N' for North component)
        Use wildcards: ?? matches any 2 chars, * matches any length
    title : str, optional
        Plot title. If None, auto-generated from station/channel
    freq_range : tuple of float, optional
        Frequency range to display (default: (0, 2.0) Hz)
    nperseg : int, optional
        STFT segment length (default: 256)
    noverlap : int, optional
        STFT overlap (default: 200)
    output_path : str, optional
        Path to save figure
    show_plot : bool, optional
        Display plot (default: True)

    Returns
    -------
    f : np.ndarray
        Frequency array
    t : np.ndarray
        Time array
    Zxx : np.ndarray
        STFT matrix
    fig : matplotlib.figure.Figure
        Figure object

    Examples
    --------
    >>> # Basic usage
    >>> f, t, Zxx, fig = analyze_stream_spectrogram(
    ...     st_syn1,
    ...     station='VUWS',
    ...     channel='??N',
    ...     freq_range=(0, 2.0)
    ... )

    >>> # Save multiple stations
    >>> for sta in ['VUWS', 'SNZO', 'WEL']:
    ...     analyze_stream_spectrogram(
    ...         st_syn1,
    ...         station=sta,
    ...         output_path=f'spectrograms/{sta}_spectrogram.png',
    ...         show_plot=False
    ...     )
    """
    # Get trace
    try:
        tr = st.select(station=station, channel=channel)[0]
    except IndexError:
        raise ValueError(f"No trace found for station={station}, channel={channel}")

    # Get sampling rate
    sampling_rate = tr.stats.sampling_rate

    # Auto-generate title if not provided
    if title is None:
        title = f"Spectrogram - {tr.stats.station} {tr.stats.channel}"

    # Compute and plot
    f, t, Zxx, fig = compute_and_plot_spectrogram(
        tr.data,
        sampling_rate,
        title=title,
        freq_range=freq_range,
        nperseg=nperseg,
        noverlap=noverlap,
        output_path=output_path,
        show_plot=show_plot
    )

    return f, t, Zxx, fig


# Example usage
if __name__ == "__main__":
    print("Module loaded successfully")
    print("Available functions:")
    print("  - compute_spectrogram(): Compute STFT")
    print("  - plot_spectrogram(): Plot STFT results")
    print("  - compute_and_plot_spectrogram(): Combined convenience function")
    print("  - analyze_stream_spectrogram(): Direct analysis from ObsPy Stream")
