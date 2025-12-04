import numpy as np
import pandas as pd


# ----------------------------
# 1) Units and residuals
# ----------------------------
def _to_g(x: pd.Series, unit: str) -> pd.Series:
    unit = unit.strip().lower()
    if unit in ("g",):
        return x.astype(float)
    if unit in ("m/s^2", "mps2", "m s^-2"):
        return x.astype(float) / 9.80665
    if unit in ("cm/s^2", "cmps2", "gal"):
        return x.astype(float) / 981.0
    raise ValueError("obs_unit must be 'g', 'm/s^2', or 'cm/s^2'.")


def build_residual_table(
    pred_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    *,
    on=("event_id", "site_id"),  # join keys (edit as needed)
    pred_mean_ln_col="mean_ln",  # ln(IM_pred in g)
    pred_sigma_col="sigma_tot",  # optional
    obs_im_col="pga",  # observed IM
    obs_unit="g",  # 'g' | 'm/s^2' | 'cm/s^2'
) -> pd.DataFrame:
    """
    Merge prediction and observation tables and compute residuals:
        resid = ln(IM_obs_g) - mean_ln_pred
    Returns a clean DataFrame with columns:
        ['event_id','site_id','ln_obs_g','mean_ln','resid','sigma_tot'(if any)] + any carry-over cols.
    """
    df = pd.merge(
        obs_df.copy(),
        pred_df.copy(),
        on=list(on),
        how="inner",
        suffixes=("_obs", "_pred"),
    )
    if pred_mean_ln_col not in df:
        raise KeyError(f"Prediction column '{pred_mean_ln_col}' not found.")
    if obs_im_col not in df:
        raise KeyError(f"Observation column '{obs_im_col}' not found.")

    df["ln_obs_g"] = np.log(_to_g(df[obs_im_col], obs_unit))
    df["mean_ln"] = df[pred_mean_ln_col].astype(float)
    df["resid"] = df["ln_obs_g"] - df["mean_ln"]

    # keep sigma if present
    if pred_sigma_col in df.columns:
        df["sigma_tot"] = df[pred_sigma_col].astype(float)
    # drop non-finite rows
    df = df[np.isfinite(df["resid"])].reset_index(drop=True)
    return df


# --------------------------------------------------
# 2) ANOVA-style decomposition (event & site terms)
# --------------------------------------------------
def decompose_residuals_simple(
    df: pd.DataFrame, *, event_col="event_id", site_col="site_id", resid_col="resid"
):
    """
    Decompose residuals into:
      - event terms e_i (between-event),
      - site terms s_k (site-to-site, after removing event terms),
      - within-event single-station residuals ε_ijk.
    Estimates:
      tau      = std(e_i)
      phi_S2S  = std(s_k)
      phi0     = std(ε_ijk)
    Returns dict with arrays and scalars.
    """
    if not {event_col, site_col, resid_col}.issubset(df.columns):
        raise KeyError("Required columns missing.")

    r = df[resid_col].to_numpy()

    # event terms (means per event)
    e = df.groupby(event_col)[resid_col].mean()
    df = df.join(e.rename("event_term"), on=event_col)

    # remove event term, then site terms (means per site of event-detrended residual)
    df["r_no_event"] = df[resid_col] - df["event_term"]
    s = df.groupby(site_col)["r_no_event"].mean()
    df = df.join(s.rename("site_term"), on=site_col)

    # within-event single-station residuals
    df["eps_single"] = df[resid_col] - df["event_term"] - df["site_term"]

    # component standard deviations (unbiased)
    tau = float(np.std(e.values, ddof=1)) if e.size > 1 else 0.0
    phi_S2S = float(np.std(s.values, ddof=1)) if s.size > 1 else 0.0
    phi0 = float(np.std(df["eps_single"].values, ddof=1)) if df.shape[0] > 1 else 0.0

    return {
        "event_terms": e,  # pd.Series, index=event_id
        "site_terms": s,  # pd.Series, index=site_id
        "eps_single": df["eps_single"],  # pd.Series aligned with df rows
        "tau": tau,
        "phi_S2S": phi_S2S,
        "phi0": phi0,
        "phi_total": float(np.sqrt(tau**2 + phi_S2S**2 + phi0**2)),
    }


# --------------------------------------------------------
# 3) Mixed-effects estimate (random event & site effects)
# --------------------------------------------------------
def decompose_residuals_mixedlm(
    df: pd.DataFrame,
    *,
    event_col="event_id",
    site_col="site_id",
    resid_col="resid",
    reml=True,
):
    """
    Linear mixed model:
        resid = μ + (1|event) + (1|site) + ε
    Uses statsmodels MixedLM with variance components for site.
    Returns tau, phi_S2S, phi0 and fitted random effects.
    """
    import statsmodels.api as sm

    if not {event_col, site_col, resid_col}.issubset(df.columns):
        raise KeyError("Required columns missing.")
    d = df[[event_col, site_col, resid_col]].dropna().copy()

    endog = d[resid_col].to_numpy()
    exog = np.ones((len(d), 1))  # intercept
    groups = d[event_col].astype("category")  # group by event
    # site variance component
    vc = {"site": f"0 + C({site_col})"}

    model = sm.MixedLM(
        endog, exog, groups=groups, exog_re=None, vc_formula=vc, re_formula="1"
    )
    res = model.fit(method="lbfgs", reml=reml, disp=False)

    # Components:
    # - Between-event variance (random intercept): cov_re[0,0]
    tau = float(np.sqrt(res.cov_re.iloc[0, 0]))
    # - Residual variance (within-event single-station):
    phi0 = float(np.sqrt(res.scale))
    # - Site-to-site variance: extract from variance components 'site'
    #   res.vcomp is a Series with one entry per variance component + group var.
    #   Safer: compute as max(total var minus others, or parse index)
    phi_S2S = 0.0
    if hasattr(res, "vcomp") and res.vcomp is not None:
        # try to find the 'site' component
        try:
            vc_ser = res.vcomp
            # statsmodels >= 0.13: MultiIndex with names ('Component','Group')
            if hasattr(vc_ser, "index"):
                idx = [
                    i
                    for i, name in enumerate(vc_ser.index)
                    if "site" in str(name).lower()
                ]
                if idx:
                    phi_S2S = float(np.sqrt(vc_ser.iloc[idx[0]]))
        except Exception:
            pass

    # Random effects (BLUPs)
    event_re = pd.Series({k: float(v.values[0]) for k, v in res.random_effects.items()})
    # site effects not returned explicitly; recompute as best linear unbiased predictions:
    # site term = group-specific VC BLUP: using design matrix for vc; for brevity, provide zeros if unavailable.
    site_re = pd.Series(dtype=float)
    if "site" in vc:
        try:
            # approximate site effects as mean residual per site after removing event BLUPs and fixed μ
            mu = float(res.fe_params[0])
            tmp = d.copy()
            tmp["event_re"] = event_re.reindex(tmp[event_col]).to_numpy()
            tmp["site_res"] = tmp[resid_col] - mu - tmp["event_re"]
            site_re = tmp.groupby(site_col)["site_res"].mean()
        except Exception:
            pass

    return {
        "tau": tau,
        "phi_S2S": phi_S2S,
        "phi0": phi0,
        "phi_total": float(np.sqrt(tau**2 + phi_S2S**2 + phi0**2)),
        "event_random_effects": event_re,  # pd.Series
        "site_random_effects": site_re,  # pd.Series (approx.)
    }


# --------------------------------------------
# 4) Standardization & bin statistics
# --------------------------------------------
def standardize_residuals(
    df: pd.DataFrame, resid_col="resid", tau=0.0, phi_S2S=0.0, phi0=0.0
):
    """
    Add standardized components:
      - z_total   = resid / sqrt(tau^2 + phi_S2S^2 + phi0^2)
      - z_event   = event_term / tau                (if available)
      - z_site    = site_term  / phi_S2S            (if available)
      - z_single  = eps_single / phi0               (if available)
    """
    out = df.copy()
    denom = (
        np.sqrt(tau**2 + phi_S2S**2 + phi0**2) if (tau + phi_S2S + phi0) > 0 else np.nan
    )
    out["z_total"] = out[resid_col] / denom
    for name, comp, sd in [
        ("event_term", "z_event", tau),
        ("site_term", "z_site", phi_S2S),
        ("eps_single", "z_single", phi0),
    ]:
        if name in out.columns and sd > 0:
            out[comp] = out[name] / sd
    return out


def binned_stats(x, y, bins):
    """
    Compute count, mean, and std of y in bins of x.
    Returns DataFrame with ['bin_left','bin_right','count','mean','std'].
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    idx = np.digitize(x, bins) - 1
    out = []
    for i in range(len(bins) - 1):
        sel = idx == i
        if np.any(sel):
            out.append(
                {
                    "bin_left": bins[i],
                    "bin_right": bins[i + 1],
                    "count": int(sel.sum()),
                    "mean": float(np.mean(y[sel])),
                    "std": float(np.std(y[sel], ddof=1)) if sel.sum() > 1 else np.nan,
                }
            )
    return pd.DataFrame(out)


import numpy as np
import pandas as pd
from typing import Tuple, Iterable, Optional
from scipy.stats import t as student_t


def binned_residual_stats(
    distance: Iterable[float],
    residual: Iterable[float],
    *,
    n_bins: int = 10,
    bin_edges: Optional[Iterable[float]] = None,
    log_bins: bool = True,
    min_count: int = 1,
    ci_level: float = 0.95,
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Compute summary statistics of residuals in distance bins.

    Parameters
    ----------
    distance : array-like
        Distances (km). Must be positive if log_bins=True.
    residual : array-like
        Residuals r = ln(obs) − ln(pred).
    n_bins : int, default 10
        Number of bins (ignored if bin_edges is provided).
    bin_edges : array-like, optional
        Explicit bin edges. If given, overrides n_bins/log_bins.
    log_bins : bool, default True
        If True, construct log-spaced bins from min..max(distance);
        otherwise use linear spacing.
    min_count : int, default 1
        Minimum number of samples required to report a bin.
    ci_level : float, default 0.95
        Confidence level for the half-width of the mean (uses Student-t).

    Returns
    -------
    bin_df : pd.DataFrame
        Columns:
        - bin_left, bin_right : float
        - center_km           : geometric mean (log_bins) or arithmetic mean (linear)
        - count               : int
        - mean_resid          : float
        - std_resid           : float (sample std, ddof=1; NaN if count<2)
        - se_resid            : float (std/sqrt(n); NaN if count<2)
        - ci_halfwidth        : float (Student-t * SE; NaN if count<2)
    bins : np.ndarray
        The bin edges actually used (length = number_of_bins + 1).
    """
    d = np.asarray(distance, dtype=float)
    r = np.asarray(residual, dtype=float)

    # Basic cleaning
    m = np.isfinite(d) & np.isfinite(r)
    if log_bins:
        m &= d > 0
    d, r = d[m], r[m]
    if d.size == 0:
        return pd.DataFrame(
            columns=[
                "bin_left",
                "bin_right",
                "center_km",
                "count",
                "mean_resid",
                "std_resid",
                "se_resid",
                "ci_halfwidth",
            ]
        ), np.array([])

    # Bins
    if bin_edges is not None:
        bins = np.asarray(bin_edges, dtype=float)
        if np.any(~np.isfinite(bins)) or bins.ndim != 1 or bins.size < 2:
            raise ValueError("bin_edges must be a 1D finite array with length ≥ 2.")
    else:
        lo = d.min()
        hi = d.max()
        if log_bins:
            bins = np.logspace(np.log10(lo), np.log10(hi), n_bins + 1)
        else:
            bins = np.linspace(lo, hi, n_bins + 1)

    # Assign to bins
    idx = np.digitize(d, bins) - 1  # 0..n_bins-1
    rows = []
    for i in range(bins.size - 1):
        sel = idx == i
        if not np.any(sel) or sel.sum() < min_count:
            continue
        x = d[sel]
        y = r[sel]
        n = int(sel.sum())

        center = float(np.exp(np.mean(np.log(x)))) if log_bins else float(np.mean(x))
        mean = float(np.mean(y))
        std = float(np.std(y, ddof=1)) if n > 1 else np.nan
        se = float(std / np.sqrt(n)) if n > 1 else np.nan

        # Student-t half-width for the mean
        if n > 1 and np.isfinite(se):
            alpha = 1.0 - ci_level
            t_mult = student_t.ppf(1.0 - alpha / 2.0, df=n - 1)
            ci_half = float(t_mult * se)
        else:
            ci_half = np.nan

        rows.append(
            {
                "bin_left": bins[i],
                "bin_right": bins[i + 1],
                "center_km": center,
                "count": n,
                "mean_resid": mean,
                "std_resid": std,
                "se_resid": se,
                "ci_halfwidth": ci_half,
            }
        )

    return pd.DataFrame(rows), bins
