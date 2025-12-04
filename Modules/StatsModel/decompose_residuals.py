import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


def mixed_decompose(
    df,
    event_col="event_id",
    site_col="site_id",
    resid_col="resid",
    reml=True,
    min_events_per_site=5,
    recenter_site_terms=True,
):
    # one value per (event, site)
    es = (
        df[[event_col, site_col, resid_col]]
        .dropna()
        .groupby([event_col, site_col], as_index=False)[resid_col]
        .mean()
        .rename(columns={resid_col: "r_es"})
    )
    es[event_col] = es[event_col].astype("category")
    es[site_col] = es[site_col].astype("category")

    md = smf.mixedlm(
        "r_es ~ 1",
        es,
        groups=es[event_col],  # random intercept: event
        re_formula="1",
        vc_formula={"site": f"0 + C({site_col})"},  # variance component: site
    )

    # res = md.fit(method="bfgs", reml=reml, disp=True)
    res_ev = md.fit(method="lbfgs", reml=reml, disp=False)

    # Fixed effect (Intercept)
    mu = res_ev.fe_params["Intercept"]

    # Random-effect SD for events (tau): from cov_re (1x1 for intercept)
    tau = float(np.sqrt(res_ev.cov_re.iloc[0, 0]))

    # Station VC SD (phi_s2s): take the (only) variance component
    if getattr(res_ev, "vcomp", None) is None or len(res_ev.vcomp) == 0:
        raise RuntimeError(
            "statsmodels did not return vcomp; please upgrade to >=0.13."
        )
    if len(res_ev.vcomp) == 1:
        phi_S2S = float(np.sqrt(res_ev.vcomp[0]))
    else:
        # Fallback for multiple VC terms: take the first; adjust if you later add more VC terms.
        phi_S2S = float(np.sqrt(res_ev.vcomp[0]))

    # The key in vcomp_map is typically the name you used in vc_formula (sta_col)

    # Residual SD (phi): statsmodels uses 'scale' for residual variance
    phi0_model = float(np.sqrt(res_ev.scale))

    # mu = float(res.fe_params["Intercept"])
    # tau = float(np.sqrt(res.cov_re.iloc[0, 0]))

    # phi0_model = float(np.sqrt(res.scale))
    # phi_S2S = extract_phi_S2S(res)  # from helper above

    # event BLUPs
    event_eff = pd.Series(
        {k: float(v.values[0]) for k, v in res_ev.random_effects.items()},
        name="event_term",
    ).sort_index()

    # site EBLUPs via shrinkage
    es = es.join(event_eff, on=event_col)
    es["r_corr"] = es["r_es"] - mu - es["event_term"]
    m_site = es.groupby(site_col, observed=True)["r_corr"].mean()
    n_site = es.groupby(site_col, observed=True).size().astype(float)

    if phi_S2S > 0 and phi0_model > 0:
        shrink = (phi_S2S**2) / (phi_S2S**2 + (phi0_model**2) / n_site)
        site_eff = shrink * m_site
    else:
        site_eff = 0.0 * m_site

    site_eff = site_eff[n_site >= float(min_events_per_site)]
    if recenter_site_terms and not site_eff.empty:
        site_eff = site_eff - site_eff.mean()
    site_eff.name = "site_term"

    out = (
        df[[event_col, site_col, resid_col]]
        .dropna()
        .join(event_eff, on=event_col)
        .join(site_eff, on=site_col)
    )
    out["site_term"] = out["site_term"].fillna(0.0)
    out["eps_single"] = out[resid_col] - mu - out["event_term"] - out["site_term"]
    phi0_empirical = (
        float(np.std(out["eps_single"].dropna().values, ddof=1))
        if len(out) > 1
        else 0.0
    )

    return {
        "mu": mu,
        "tau": tau,
        "phi_S2S": phi_S2S,
        "phi0_model": phi0_model,
        "phi0_empirical": phi0_empirical,
        "event_effects": event_eff.sort_index(),
        "site_effects": site_eff.sort_index(),
        "df_components": out[
            [event_col, site_col, resid_col, "event_term", "site_term", "eps_single"]
        ],
    }


def extract_phi_S2S(res, label_substr="site"):
    import numpy as np

    v = getattr(res, "vcomp", None)
    if v is None:
        return 0.0
    # pandas Series case
    if hasattr(v, "index"):
        for name, val in v.items():
            if label_substr in str(name).lower():
                return float(np.sqrt(val))
        if len(v) == 1:  # unnamed single component
            return float(np.sqrt(v.iloc[0]))
        return 0.0
    # ndarray case
    arr = np.asarray(v).ravel()
    names = getattr(res, "vcomp_names", None) or getattr(
        getattr(res, "model", None), "vcomp_names", None
    )
    if names and len(names) == len(arr):
        for nm, val in zip(names, arr):
            if label_substr in str(nm).lower():
                return float(np.sqrt(val))
    return float(np.sqrt(arr[0])) if arr.size == 1 else 0.0


# ==============================================================================
# APPROACH 4: Sequential Fitting
# ==============================================================================


def method_4_sequential_fitting(
    df,
    event_col="event_id",
    site_col="site_id",
    resid_col="resid",
    reml=True,
    min_events_per_site=5,
    recenter_site_terms=True,
):
    """
    Fit effects sequentially and combine
    """

    print("=" * 50)
    print("METHOD 4: Sequential Fitting")
    print("=" * 50)

    # one value per (event, site)
    es = (
        df[[event_col, site_col, resid_col]]
        .dropna()
        .groupby([event_col, site_col], as_index=False)[resid_col]
        .mean()
        .rename(columns={resid_col: "r_es"})
    )
    es[event_col] = es[event_col].astype("category")
    es[site_col] = es[site_col].astype("category")

    # Step 1: Fit first random effect
    model1 = smf.mixedlm(
        "r_es ~ 1",
        es,
        groups=es[event_col],  # random intercept: event
    )
    result1 = model1.fit()

    # Get residuals from first model
    df_temp = es.copy()
    df_temp["residuals1"] = df_temp["r_es"] - result1.fittedvalues

    # Step 2: Fit second random effect on residuals
    model2 = smf.mixedlm("residuals1 ~ 1", df_temp, groups=df_temp["site_id"])
    result2 = model2.fit()

    # Combine results
    fixed_effect = result1.params["Intercept"]
    var_evid = result1.cov_re.iloc[0, 0]
    var_sta = result2.cov_re.iloc[0, 0]
    var_residual = result2.scale

    print(f"Fixed intercept: {fixed_effect:.6f}")
    print(f"Random effect variance (evid): {var_evid:.6f}")
    print(f"Random effect variance (sta): {var_sta:.6f}")
    print(f"Residual variance: {var_residual:.6f}")

    # return {
    #     'fixed_intercept': fixed_effect,
    #     'var_evid': var_evid,
    #     'var_sta': var_sta,
    #     'var_residual': var_residual,
    #     'model1': result1,
    #     'model2': result2
    # }

    return {
        "mu": fixed_effect,
        "tau": var_evid,
        "phi_S2S": var_sta,
        "phi0_model": var_residual,
        "phi0_empirical": var_residual,
    }


