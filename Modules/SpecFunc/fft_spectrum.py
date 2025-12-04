import numpy as np
import pandas as pd

def _infer_fs(series):
    """Infer sampling frequency (Hz) from Series index."""
    idx = series.index
    if isinstance(idx, pd.DatetimeIndex):
        # seconds between samples → Hz
        dt = np.median(np.diff(idx.view("int64"))) / 1e9  # ns → s
    else:
        # assume numeric axis is in seconds
        dt = float(np.median(np.diff(idx.values)))
    if dt <= 0 or not np.isfinite(dt):
        raise ValueError("Cannot infer positive sampling interval from index.")
    return 1.0 / dt

def fft_spectrum(series, fs=None, detrend="mean", window="hann"):
    """
    Compute one-sided amplitude and power spectra via FFT.
    
    Parameters
    ----------
    series : pd.Series
        Time series with regular sampling.
    fs : float or None
        Sampling frequency (Hz). If None, inferred from index spacing.
    detrend : {'mean','linear',None}
        Remove mean or best-fit line before FFT.
    window : {'hann', None}
        Optional taper to reduce leakage.
        
    Returns
    -------
    f : ndarray
        Frequencies (Hz).
    A : ndarray
        One-sided amplitude spectrum (same units as input per √Hz is NOT applied).
    P : ndarray
        One-sided power spectrum (|FFT|^2) scaled to preserve Parseval energy.
    """
    x = series.dropna().values.astype(float)
    if fs is None:
        fs = _infer_fs(series)

    N = x.size
    t = np.arange(N) / fs

    # detrend
    if detrend == "mean":
        x = x - x.mean()
    elif detrend == "linear":
        p = np.polyfit(t, x, 1)
        x = x - np.polyval(p, t)

    # window
    if window == "hann":
        w = np.hanning(N)
    else:
        w = np.ones(N)

    xw = x * w
    # coherent gain for amplitude scaling
    cg = w.mean()

    # FFT (one-sided)
    X = np.fft.rfft(xw)
    f = np.fft.rfftfreq(N, d=1/fs)

    # Amplitude spectrum (one-sided, correct for window + factor 2 except DC/Nyquist)
    A = np.abs(X) / (cg * N)
    if N % 2 == 0:
        A[1:-1] *= 2.0
    else:
        A[1:] *= 2.0

    # Power spectrum (one-sided), scale so that sum(P)*df ~= signal power after window
    # Window power correction:
    U = (w**2).mean()
    P = (np.abs(X)**2) / (fs * N**2 * U)
    if N % 2 == 0:
        P[1:-1] *= 2.0
    else:
        P[1:] *= 2.0

    return f, A, P

# --- Optional: Welch PSD (requires SciPy) ---
def welch_psd(series, fs=None, nperseg=None, noverlap=None, detrend='constant'):
    """
    Power Spectral Density via Welch’s method (if SciPy is available).
    Returns (f, PSD).
    """
    try:
        from scipy.signal import welch
    except Exception as e:
        raise ImportError("SciPy is required for welch_psd().") from e

    if fs is None:
        fs = _infer_fs(series)

    f, Pxx = welch(series.dropna().values.astype(float), fs=fs,
                   nperseg=nperseg, noverlap=noverlap, detrend=detrend,
                   window='hann', return_onesided=True, scaling='density')
    return f, Pxx
