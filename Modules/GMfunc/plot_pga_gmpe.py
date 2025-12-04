"""
Ground Motion Prediction Equation (GMPE) plotting utilities.

This module provides functions to compare GMPE predictions with observed data.
For residual analysis and binned statistics, see residual_analysis.py module.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_pga_gmpe_vs_observed(
    pred_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    *,
    distance_col: str = "distance_km",  # column name in BOTH tables
    mean_ln_col: str = "mean_ln",  # ln(PGA in g) from GMPE
    sigma_col: str = "sigma_tot",  # optional; ln-stddev
    obs_im_col: str = "pga",  # observed PGA column
    obs_unit: str = "g",  # 'g' | 'm/s^2' | 'cm/s^2'
    filter_mw: float | None = None,  # if pred_df has 'Mw' column, select one
    distance_type_label: str = "Rrup",  # axis label only (e.g., Rrup, Rhypo)
    label_pred: str = "GMPE mean",
    label_obs: str = "Observed",
    show_residuals: bool = True,
    figsize=(7, 6),
):
    """
    Plot GMPE predictions (mean ±1σ if available) against observed PGA.

    pred_df required columns:
        distance_col, mean_ln_col; optionally sigma_col and 'Mw'.
        mean_ln must be natural log of PGA in g.

    obs_df required columns:
        distance_col, obs_im_col (PGA in obs_unit).

    Returns
    -------
    (fig, (ax_main, ax_resid_or_None))
    """
    dfp = pred_df.copy()
    dfo = obs_df.copy()

    # If multiple magnitudes are present, filter to one
    if filter_mw is not None and "Mw" in dfp.columns:
        dfp = dfp[dfp["Mw"] == filter_mw].copy()
        if dfp.empty:
            raise ValueError(f"No prediction rows found for Mw={filter_mw}")

    # Rename to standard working names
    if distance_col not in dfp or distance_col not in dfo:
        raise KeyError(f"'{distance_col}' must exist in both pred_df and obs_df")

    # Convert observed PGA to g
    unit = obs_unit.strip().lower()
    if unit in ("g",):
        dfo["_pga_g"] = dfo[obs_im_col].astype(float)
    elif unit in ("m/s^2", "mps2", "m s^-2", "m s^-2"):
        dfo["_pga_g"] = dfo[obs_im_col].astype(float) / 9.80665
    elif unit in ("cm/s^2", "cmps2", "gal"):
        dfo["_pga_g"] = dfo[obs_im_col].astype(float) / 981.0
    else:
        raise ValueError("obs_unit must be 'g', 'm/s^2', or 'cm/s^2'")

    # Compute GMPE mean and bands in g
    if mean_ln_col not in dfp:
        raise KeyError(f"'{mean_ln_col}' (ln(PGA in g)) is required in pred_df")
    dfp["_mean_g"] = np.exp(dfp[mean_ln_col].astype(float))

    has_sigma = sigma_col in dfp.columns and np.isfinite(dfp[sigma_col]).any()
    if has_sigma:
        s = dfp[sigma_col].astype(float)
        dfp["_lo_g"] = np.exp(dfp[mean_ln_col] - s)
        dfp["_hi_g"] = np.exp(dfp[mean_ln_col] + s)

    # Clean/validate for log plotting
    def _posfinite(a):
        return np.isfinite(a) & (a > 0)

    m_pred = _posfinite(dfp[distance_col].values) & _posfinite(dfp["_mean_g"].values)
    m_obs = _posfinite(dfo[distance_col].values) & _posfinite(dfo["_pga_g"].values)

    dfp = dfp.loc[m_pred].sort_values(distance_col)
    dfo = dfo.loc[m_obs].sort_values(distance_col)

    if dfp.empty or dfo.empty:
        raise ValueError(
            "After filtering nonpositive/nonfinite values, there is no data to plot."
        )

    # Figure
    if show_residuals:
        fig, (ax, axr) = plt.subplots(
            2, 1, figsize=figsize, sharex=True, gridspec_kw={"height_ratios": [3, 1]}
        )
    else:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        axr = None

    # Observations
    ax.scatter(
        dfo[distance_col], dfo["_pga_g"], s=28, alpha=0.8, c="gray", label=label_obs
    )

    # GMPE mean and ±1σ
    ax.plot(dfp[distance_col], dfp["_mean_g"], lw=2, label=label_pred)
    if has_sigma:
        ax.fill_between(
            dfp[distance_col], dfp["_lo_g"], dfp["_hi_g"], alpha=0.25, label="±1σ"
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"{distance_type_label} (km)")
    ax.set_ylabel("PGA (g)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    # Residuals: ln(obs) − ln(pred) (interpolate mean_ln onto obs distances in log-space)
    if show_residuals:
        # Interpolate predicted ln(PGA) vs ln(distance)
        from scipy.interpolate import interp1d

        f_ln_mean = interp1d(
            np.log(dfp[distance_col].values),
            np.log(dfp["_mean_g"].values),
            kind="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        ln_res = np.log(dfo["_pga_g"].values) - f_ln_mean(
            np.log(dfo[distance_col].values)
        )
        ok = np.isfinite(ln_res)
        axr.axhline(0.0, color="k", lw=1)
        axr.scatter(dfo[distance_col].values[ok], ln_res[ok], s=20, alpha=0.9)
        axr.set_xscale("log")
        axr.set_ylabel("ln(obs) − ln(pred)")
        axr.set_xlabel(f"{distance_type_label} (km)")
        axr.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    return fig, (ax, axr)


def to_g(series, unit="g"):
    u = unit.lower()
    if u == "g":
        return series.astype(float)
    if u in ("m/s^2", "mps2", "m s^-2"):
        return series.astype(float) / 9.80665
    if u in ("cm/s^2", "cmps2", "gal"):
        return series.astype(float) / 981.0
    raise ValueError("obs_unit must be 'g', 'm/s^2', or 'cm/s^2'.")


# --- 3) Plot GMPE vs observed (log–log), with ±1σ if available ---
def plot_gmpe_vs_obs(
    pred_df,
    obs_df,
    distance_col="distance_km",
    mean_ln_col="mean_ln",
    sigma_col="sigma_tot",
    obs_im_col="pga",
    obs_unit="g",
    label_pred="GMPE mean",
    label_obs="Observed",
    distance_label="Distance (km)",
):
    # Prepare curves
    p = pred_df[
        [distance_col, mean_ln_col]
        + ([sigma_col] if sigma_col in pred_df.columns else [])
    ].dropna()
    p = p.sort_values(distance_col)
    mean_g = np.exp(p[mean_ln_col].values)
    d_pred = p[distance_col].values
    # Observations to g
    o = obs_df[[distance_col, obs_im_col]].copy()
    o["_g"] = to_g(o[obs_im_col], obs_unit)
    o = o.replace([np.inf, -np.inf], np.nan).dropna()
    o = o[(o[distance_col] > 0) & (o["_g"] > 0)]
    d_obs = o[distance_col].values
    y_obs = o["_g"].values

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(d_obs, y_obs, s=12, alpha=0.7, label=label_obs, c="gray")
    ax.plot(d_pred, mean_g, lw=2, label=label_pred, c="royalblue")

    if sigma_col in p.columns:
        s = p[sigma_col].values.astype(float)
        ax.fill_between(
            d_pred,
            np.exp(np.log(mean_g) - s),
            np.exp(np.log(mean_g) + s),
            alpha=0.25,
            label="±1σ",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(distance_label)
    ax.set_ylabel("PGA (g)")
    ax.grid(True, which="both", alpha=0.3, linestyle=":")
    ax.legend()
    plt.tight_layout()
    return fig, ax


def bin_mean_std(nbins, sa1, rrup):
    """
    Compute binned mean and standard deviation in log10 space.

    NOTE: This function is also available in GMfunc.residual_analysis with
    enhanced documentation and additional plotting functions.

    Parameters
    ----------
    nbins : int
        Number of distance bins
    sa1 : array-like
        Ground motion data (linear space, will be converted to log10)
    rrup : array-like
        Distance to hypocenter or rupture (km)

    Returns
    -------
    mean_bins : np.ndarray
        Mean of log10(sa1) in each bin
    std_bins : np.ndarray
        Standard deviation of log10(sa1) in each bin
    rjb_bins : np.ndarray
        Distance bin edges (left edges)
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


def bin_mean_err(nbins, sa1, rrup):
    """
    Compute binned mean and standard error in log10 space.

    NOTE: This function is also available in GMfunc.residual_analysis with
    enhanced documentation and additional plotting functions.

    Parameters
    ----------
    nbins : int
        Number of distance bins
    sa1 : array-like
        Ground motion data (linear space, will be converted to log10)
    rrup : array-like
        Distance to hypocenter or rupture (km)

    Returns
    -------
    mean_bins : np.ndarray
        Mean of log10(sa1) in each bin
    err_bins : np.ndarray
        Standard error of log10(sa1) in each bin (std/n)
    rjb_bins : np.ndarray
        All distance bin edges
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
                std_bins[ik] = np.std(log_data) / n_data
            else:
                std_bins[ik] = 0.0
        else:
            # Empty bin - set to NaN so it won't be plotted
            mean_bins[ik] = np.nan
            std_bins[ik] = np.nan

    return mean_bins, std_bins, rjb_bins
