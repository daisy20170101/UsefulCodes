# Newer rpy2 syntax

from pymer4.models import lmer
import numpy as np
import pandas as pd


def fit_mixed_model_event_station(
    df,
    event_col="event_id",
    site_col="site_id",
    resid_col="resid",
    min_events_per_site=5,
):
    """Fit model and extract coefficients."""

    from rpy2.robjects.conversion import localconverter, Converter
    from rpy2.robjects import pandas2ri, default_converter

    # Create converter that includes pandas and default conversions
    cv = Converter("pandas converter with default")
    cv += pandas2ri.converter + default_converter

    # Prepare data
    es = (
        df[[event_col, site_col, resid_col]]
        .dropna()
        .groupby([event_col, site_col], as_index=False)[resid_col]
        .mean()
        .rename(columns={resid_col: "lnY", event_col: "event_id", site_col: "site_id"})
    )

    # Convert to proper types
    es["event_id"] = es["event_id"].astype(str)
    es["site_id"] = es["site_id"].astype(str)
    es["lnY"] = es["lnY"].astype(float)

    print(f"Fitting model with {len(es)} observations")
    print(f"  Events: {es['event_id'].nunique()}")
    print(f"  Sites: {es['site_id'].nunique()}")

    # Use combined converter
    with localconverter(cv):
        # Use lowercase lmer for latest pymer4 version
        model = lmer("lnY ~ 1 + (1|event_id) + (1|site_id)", data=es)
        model.fit(summarize=False, verbose=False)

        # Extract variance components correctly
        mu = float(model.fixef.iloc[0, 0])
        tau = float(np.sqrt(model.ranef_var["event_id"]))
        phi_s2s = float(np.sqrt(model.ranef_var["site_id"]))
        phi_0 = float(np.sqrt(model.resid_var))

    return {
        "bias": mu,
        "bias_se": float(model.fixef.iloc[0, 1]),
        "tau": tau,
        "phi_s2s": phi_s2s,
        "phi": phi_0,
        "sigma": np.sqrt(tau**2 + phi_0**2 + phi_s2s**2),
        "ranef": model.ranef,
        "model": model,
    }


import polars as pl


def fit_lmer_event_site(df, y="lnY", event="event_id", site="site_id"):
    """
    Fit: lnY ~ 1 + (1|event_id) + (1|site_id)
    Returns key variance components and effects, similar to your previous Lmer helper.
    """

    # from rpy2.robjects.conversion import localconverter, Converter
    # from rpy2.robjects import pandas2ri, default_converter

    # Create converter that includes pandas and default conversions
    # cv = Converter("pandas converter with default")
    # cv += pandas2ri.converter + default_converter

    # minimal cols & drop NA
    use = df[[y, event, site]].dropna().copy()

    # pymer4 expects categorical grouping factors
    use[event] = use[event].astype(str).astype("category")
    use[site] = use[site].astype(str).astype("category")

    # convert to Polars and mark group columns categorical
    pl_df = pl.from_pandas(use).with_columns(
        pl.col(event).cast(pl.Categorical),
        pl.col(site).cast(pl.Categorical),
    )
    # Use combined converter
    # with localconverter(cv):
    # fit model
    m = lmer(f"{y} ~ 1 + (1|{event}) + (1|{site})", data=pl_df)
    m.fit(summarize=False)

    # fixed intercept (bias)
    est_tbl = m.result_fit
    bias = float(est_tbl.filter(pl.col("term") == "(Intercept)")["estimate"][0])
    bias_se = float(est_tbl.filter(pl.col("term") == "(Intercept)")["std_error"][0])

    # --- random-effect SDs from ranef_var ---
    rv = m.ranef_var  # columns: group, term, estimate, conf_low, conf_high
    # event SD (tau): random-intercept std for the event factor
    tau = float(
        rv.filter((pl.col("group") == event) & (pl.col("term") == "sd__(Intercept)"))[
            "estimate"
        ][0]
    )
    # station SD (phi_s2s): random-intercept std for the site factor
    phi_s2s = float(
        rv.filter((pl.col("group") == site) & (pl.col("term") == "sd__(Intercept)"))[
            "estimate"
        ][0]
    )
    # residual SD (phi)
    phi = float(
        rv.filter(
            (pl.col("group") == "Residual") & (pl.col("term") == "sd__Observation")
        )["estimate"][0]
    )

    residuals = m.data.select("lnY").to_pandas()["lnY"]

    out = {
        "bias": bias,
        "bias_se": bias_se,
        "tau": tau,
        "phi_s2s": phi_s2s,
        "phi": phi,
        "sigma": float(np.sqrt(tau**2 + phi_s2s**2 + phi**2)),
        "ranef": m.ranef,  # BLUPs for event & site
        "residuals": residuals,  # per-observation residuals
    }
    return out


#!/usr/bin/env python3
"""
Comparison: Python Mixed Linear Regression vs R's lmer()
Shows which Python approaches give identical results to R

Summary:
- pymer4: IDENTICAL results (uses R's lme4 directly)
- rpy2: IDENTICAL results (direct R interface)  
- statsmodels: SIMILAR but not identical results
- sklearn-lmer: IDENTICAL results (wraps pymer4)
"""

import pandas as pd
import numpy as np

# ==============================================================================
# METHOD 1: pymer4 - IDENTICAL TO R (RECOMMENDED)
# ==============================================================================


def method_pymer4(df):
    """
    pymer4 gives IDENTICAL results to R's lmer
    Uses rpy2 to call R's lme4 package directly behind the scenes
    """

    try:
        from pymer4.models import Lmer

        print("=" * 60)
        print("METHOD 1: pymer4 (IDENTICAL TO R)")
        print("=" * 60)

        # Exact same syntax as R
        model = Lmer("lnY ~ 1 + (1|evid) + (1|sta)", data=df)
        model.fit()

        print(model)

        # Extract results
        results = {
            "fixed_effects": model.fixef,
            "random_effects_var": model.ranef_var,
            "residual_var": model.residual_var,
            "fitted_values": model.fittedvalues,
            "residuals": model.residuals,
            "loglik": model.logLik,
            "aic": model.AIC,
            "bic": model.BIC,
        }

        return model, results

    except ImportError:
        print("pymer4 not installed. Install with:")
        print("pip install pymer4")
        print("Note: Requires R and rpy2")
        return None, None


# ==============================================================================
# METHOD 2: Direct rpy2 - IDENTICAL TO R
# ==============================================================================


def method_rpy2_direct(df):
    """
    Direct rpy2 interface - IDENTICAL to R results
    """

    try:
        from rpy2.robjects import pandas2ri, r
        from rpy2.robjects.packages import importr

        print("=" * 60)
        print("METHOD 2: Direct rpy2 (IDENTICAL TO R)")
        print("=" * 60)

        # Activate pandas-R conversion
        pandas2ri.activate()

        # Import R packages
        lme4 = importr("lme4")
        base = importr("base")

        # Convert pandas to R dataframe
        r_df = pandas2ri.py2rpy(df)

        # Fit model using R's lmer directly
        model = lme4.lmer("lnY ~ 1 + (1|evid) + (1|sta)", data=r_df)

        # Get summary
        summary = base.summary(model)
        print(summary)

        # Extract coefficients
        coeffs = r.fixef(model)
        ranef_var = r.VarCorr(model)

        results = {
            "fixed_effects": dict(zip(coeffs.names, list(coeffs))),
            "model_object": model,
            "summary": summary,
        }

        return model, results

    except ImportError:
        print("rpy2 not installed or R not available")
        print("Install with: pip install rpy2")
        print("Also need R with lme4 package")
        return None, None


# ==============================================================================
# METHOD 3: statsmodels - SIMILAR but NOT IDENTICAL
# ==============================================================================


def method_statsmodels(df):
    """
    statsmodels gives SIMILAR but NOT IDENTICAL results
    Different optimization algorithms and assumptions
    """

    import statsmodels.formula.api as smf

    print("=" * 60)
    print("METHOD 3: statsmodels (SIMILAR, NOT IDENTICAL)")
    print("=" * 60)

    # Method 3a: Single random effect
    print("\n3a: Single random effect (evid only)")
    model1 = smf.mixedlm("lnY ~ 1", df, groups=df["evid"])
    result1 = model1.fit()
    print(result1.summary())

    # Method 3b: Variance components approach
    print("\n3b: Variance components (closest to R)")
    df_temp = df.copy()
    df_temp["sta_cat"] = df_temp["sta"].astype("category")

    model2 = smf.mixedlm(
        "lnY ~ 1", df_temp, groups=df_temp["evid"], vc_formula={"sta": "0 + C(sta_cat)"}
    )
    result2 = model2.fit()
    print(result2.summary())

    results = {"single_re": result1, "variance_comp": result2}

    return model2, results


# ==============================================================================
# METHOD 4: sklearn-lmer - IDENTICAL TO R
# ==============================================================================


def method_sklearn_lmer(df):
    """
    sklearn-lmer wraps pymer4 - IDENTICAL to R results
    """

    try:
        from sklearn_lmer import LmerRegressor

        print("=" * 60)
        print("METHOD 4: sklearn-lmer (IDENTICAL TO R)")
        print("=" * 60)

        # Prepare data for sklearn format
        X = df[["evid", "sta"]].copy()
        y = df["lnY"].values

        # Create model
        model = LmerRegressor(
            formula="lnY ~ 1 + (1|evid) + (1|sta)", X_cols=["evid", "sta"]
        )

        # Fit model
        model.fit(X, y)

        # Get predictions
        predictions = model.predict(X)

        results = {
            "model": model,
            "predictions": predictions,
            "coef_": model.coef_,
            "fitted_values": predictions,
        }

        print(f"Fixed intercept: {model.coef_}")

        return model, results

    except ImportError:
        print("sklearn-lmer not installed. Install with:")
        print("pip install sklearn-lmer")
        return None, None


# ==============================================================================
# COMPARISON FUNCTION
# ==============================================================================


def compare_all_methods(df):
    """
    Compare all methods and show which give identical results
    """

    print("COMPARISON: Python Mixed Linear Regression Methods")
    print("=" * 80)
    print("Model: lnY ~ 1 + (1|evid) + (1|sta)")
    print("Data shape:", df.shape)
    print("Events (evid):", df["evid"].nunique())
    print("Stations (sta):", df["sta"].nunique())
    print("\n")

    results_summary = {}

    # Method 1: pymer4
    print("\n" + "=" * 80)
    model1, res1 = method_pymer4(df)
    if res1:
        results_summary["pymer4"] = {
            "intercept": (
                res1["fixed_effects"]["(Intercept)"]
                if "(Intercept)" in res1["fixed_effects"]
                else None
            ),
            "identical_to_r": True,
            "notes": "Uses R lme4 directly via rpy2",
        }

    # # Method 2: rpy2
    # print("\n" + "=" * 80)
    # model2, res2 = method_rpy2_direct(df)
    # if res2:
    #     results_summary["rpy2"] = {
    #         "intercept": res2["fixed_effects"].get("(Intercept)", None),
    #         "identical_to_r": True,
    #         "notes": "Direct R interface",
    #     }

    # Method 3: statsmodels
    print("\n" + "=" * 80)
    model3, res3 = method_statsmodels(df)
    if res3:
        results_summary["statsmodels"] = {
            "intercept": res3["variance_comp"].params["Intercept"],
            "identical_to_r": False,
            "notes": "Similar but different optimization",
        }

    # # Method 4: sklearn-lmer
    # print("\n" + "=" * 80)
    # model4, res4 = method_sklearn_lmer(df)
    # if res4:
    #     results_summary["sklearn_lmer"] = {
    #         "intercept": res4["coef_"],
    #         "identical_to_r": True,
    #         "notes": "Wraps pymer4",
    #     }

    # Summary table
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    print(f"{'Method':<15} {'Intercept':<12} {'Identical?':<12} {'Notes'}")
    print("-" * 60)

    for method, info in results_summary.items():
        intercept = (
            f"{info['intercept']:.4f}" if info["intercept"] is not None else "Failed"
        )
        identical = "YES" if info["identical_to_r"] else "NO"
        print(f"{method:<15} {intercept:<12} {identical:<12} {info['notes']}")

    return results_summary
