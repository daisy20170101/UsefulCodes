import numpy as np
import pandas as pd


def create_resid_table(
    df="example.csv",
    gmpefile="gmpe_example.csv",
    output_dir="atk22_resd",
    backarc_flag=True,
    magnitude=5.0,
):
    """
    Prepare residual data by calculating between-event and within-event residuals

    Parameters:
    df: DataFrame with columns ['event_id', 'site_id', 'magnitude', 'r_rup','PGA']
    gmpefile2: GMPE model from OpenQuake: columns ['rrup_km','mean_ln','sigma_tot']
    output_dir: dir to output residual table
    backarc_flag: True or False
    magnitude: e.g. 5.0

    """

    mag = np.round(magnitude, 1)

    # backarc_flag = False

    if backarc_flag:
        arcflag = "backarc"
    else:
        arcflag = "forearc"

    print(arcflag)

    # An example of CSV files
    # gmpefile = '/Users/DuoL/Documents/NSHM/Attenuation/atk22_'+ arcflag +'/atk22_slab_pga_Mw' + str(mag)+'.csv'
    db4 = pd.read_csv(gmpefile)

    print(db4.keys())

    # db3_s5a is forearc while db3_s5 is backarc; output, db3_s4 for selected magnitude
    db3_s6 = df[df["mag"] < mag + 0.25]
    db3_s4 = db3_s6[db3_s6["mag"] > mag - 0.25]

    db4_dist = []
    db4_mean = []
    db4_evid = []
    db4_staid = []
    db4_sigma = []

    for ik, rrup in enumerate(db3_s4["r_rup"]):

        id_db4 = np.abs(db4["rrup_km"] - rrup).idxmin()
        db4_dist.append(rrup)
        db4_mean.append(db4["mean_ln"][id_db4])
        db4_evid.append(db3_s4["evid"].values[ik])
        db4_staid.append(db3_s4["sta"].values[ik])
        db4_sigma.append(db4["sigma_tot"][id_db4])

    obs = pd.DataFrame(
        {
            "distance_km": db3_s4[
                "r_rup"
            ],  # if atk22, r_rup is used, while bc hydro uses r_hypo
            "pga": db3_s4["PGA"],
            "event_id": db3_s4["evid"],
            "magnitude": db3_s4["mag"],
            "site_id": db3_s4["sta"],
            "Vs30": db3_s4["Vs30"],
        }
    )

    pred = pd.DataFrame(
        {
            "distance_km": db4_dist,
            "mean_ln": db4_mean,  # ln(PGA in g)
            "sigma_tot": db4_sigma,
        }
    )

    obs.to_csv(output_dir + "/obs_df_m" + str(mag) + "_" + arcflag + ".csv")
    pred.to_csv(output_dir + "/pred_df_m" + str(mag) + "_" + arcflag + ".csv")

    print("mag:", mag)

    # --- 2) Merge and compute residuals: r = ln(obs_g) - mean_ln ---
    def to_g(series, unit="g"):
        u = unit.lower()
        if u == "g":
            return series.astype(float)
        if u in ("m/s^2", "mps2", "m s^-2"):
            return series.astype(float) / 9.80665
        if u in ("cm/s^2", "cmps2", "gal"):
            return series.astype(float) / 981.0
        raise ValueError("obs_unit must be 'g', 'm/s^2', or 'cm/s^2'.")

    df = pd.merge(
        obs, pred, on=["distance_km"], how="inner", suffixes=("_obs", "_pred")
    )
    df["ln_obs_g"] = np.log(to_g(df["pga"], "g"))

    df["resid"] = df["ln_obs_g"] - df["mean_ln"]  # natural-log residual

    print(df["resid"].mean())

    # Harmonize distance column once:
    if "distance_km_pred" in df.columns:
        df["distance_km"] = df["distance_km_pred"]
    elif "distance_km_obs" in df.columns:
        df["distance_km"] = df["distance_km_obs"]
    elif "distance_km" in df.columns:
        pass
    else:
        # last resort: look for common alternatives
        for cand in [
            "rrup",
            "rjb",
            "rx",
            "rhypo",
            "rrup_km",
            "rjb_km",
            "rx_km",
            "rhypo_km",
        ]:
            if cand in df.columns:
                df["distance_km"] = df[cand].astype(float)
                break
        else:
            raise KeyError("No distance column found after merge.")

    df = df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["distance_km", "ln_obs_g", "mean_ln", "resid"]
    )

    # Export residuals
    df[
        [
            "event_id",
            "site_id",
            "magnitude",
            "distance_km",
            "ln_obs_g",
            "mean_ln",
            "resid",
            "Vs30",
        ]
        + (["sigma_tot"] if "sigma_tot" in df.columns else [])
    ].to_csv(output_dir + "/resid_m" + str(mag) + "-" + arcflag + ".csv", index=False)


def create_resid_table_band(
    df_input="example.csv",
    gmpefile_dir="/Users/DuoL/Documents/NSHM/Attenuation/atk22_",
    modelname="ag20",
    model_tag="slab",
    backarc_flag=True,
    sa_list=["sa01", "sa02"],
    psa_list=["pSA_0.1", "pSA_0.2"],
    magnitude=5.0,
):
    """
    Prepare residual data by calculating between-event and within-event residuals

    Parameters:
    df: DataFrame with columns ['event_id', 'site_id', 'magnitude', 'r_rup','PGA']
    gmpefile2: GMPE model from OpenQuake: columns ['rrup_km','mean_ln','sigma_tot']
    output_dir: dir to output residual table
    backarc_flag: True or False
    magnitude: e.g. 5.0

    """

    mag1 = np.round(magnitude, 1)

    if backarc_flag:
        arcflag = "backarc"
    else:
        arcflag = "forearc"

    print(arcflag)

    output_dir = modelname + "_resd"

    mask1 = df_input["mag"].between(mag1 - 0.25, mag1 + 0.25, inclusive="both")
    db3_s4 = df_input[mask1].copy()

    for isa, sa in enumerate(sa_list):

        print(sa)

        # An example of CSV files
        gmpefile = (
            gmpefile_dir
            + arcflag
            + "/"
            + modelname
            + "_"
            + model_tag
            + "_"
            + sa
            + "_Mw"
            + str(mag1)
            + ".csv"
        )

        print(gmpefile)

        db4 = pd.read_csv(gmpefile)

        # Handle both 'mag' and 'magnitude' column names

        db4_dist = []
        db4_mean = []
        db4_evid = []
        db4_staid = []
        db4_sigma = []

        period, db3_s4_sub = get_period_and_filter(psa_list[isa], db3_s4)

        for ik, rrup in enumerate(db3_s4_sub["r_rup"]):

            id_db4 = np.abs(db4["rrup_km"] - rrup).idxmin()
            db4_dist.append(rrup)
            db4_mean.append(db4["mean_ln"][id_db4])
            db4_evid.append(db3_s4_sub["evid"].values[ik])
            db4_staid.append(db3_s4_sub["sta"].values[ik])
            db4_sigma.append(db4["sigma_tot"][id_db4])

        obs = pd.DataFrame(
            {
                "distance_km": db3_s4_sub[
                    "r_rup"
                ],  # if atk22, r_rup is used, while bc hydro uses r_hypo
                f"{psa_list[isa]}": db3_s4_sub[psa_list[isa]],
                "event_id": db3_s4_sub["evid"],
                "magnitude": db3_s4_sub["mag"],
                "site_id": db3_s4_sub["sta"],
                "Vs30": db3_s4_sub["Vs30"],
            }
        )

        pred = pd.DataFrame(
            {
                "distance_km": db4_dist,
                "mean_ln": db4_mean,  # ln(PGA in g)
                # "sigma_tot": db4_sigma,
            }
        )

        # obs.to_csv(output_dir + "/obs_df_m" + str(mag) + "_" + arcflag + ".csv")
        # pred.to_csv(output_dir + "/pred_df_m" + str(mag) + "_" + arcflag + ".csv")

        print(obs.keys(), pred.keys())

        # --- 2) Merge and compute residuals: r = ln(obs_g) - mean_ln ---
        def to_g(series, unit="g"):
            u = unit.lower()
            if u == "g":
                return series.astype(float)
            if u in ("m/s^2", "mps2", "m s^-2"):
                return series.astype(float) / 9.80665
            if u in ("cm/s^2", "cmps2", "gal"):
                return series.astype(float) / 981.0
            raise ValueError("obs_unit must be 'g', 'm/s^2', or 'cm/s^2'.")

        df = pd.merge(
            obs, pred, on=["distance_km"], how="inner", suffixes=("_obs", "_pred")
        )
        df["ln_obs_g"] = np.log(to_g(df[psa_list[isa]], "g"))

        df[psa_list[isa]] = (
            df["ln_obs_g"] - df["mean_ln"]
        )  # natural-log residual # replace data with residual

        print("pSA:", sa, df[psa_list[isa]].mean())

        df = df.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["distance_km", "ln_obs_g", "mean_ln", f"{psa_list[isa]}"]
        )

        # Harmonize distance column once:
        if "distance_km_pred" in df.columns:
            df["distance_km"] = df["distance_km_pred"]
        elif "distance_km_obs" in df.columns:
            df["distance_km"] = df["distance_km_obs"]
        elif "distance_km" in df.columns:
            pass
        else:
            # last resort: look for common alternatives
            for cand in [
                "rrup",
                "rjb",
                "rx",
                "rhypo",
                "rrup_km",
                "rjb_km",
                "rx_km",
                "rhypo_km",
            ]:
                if cand in df.columns:
                    df["distance_km"] = df[cand].astype(float)
                    break
            else:
                raise KeyError("No distance column found after merge.")

        # Export residuals
        df[
            [
                "event_id",
                "site_id",
                "magnitude",
                "distance_km",
                f"{psa_list[isa]}",
            ]
            + (["sigma_tot"] if "sigma_tot" in df.columns else [])
        ].to_csv(
            output_dir
            + "/resid_"
            + model_tag
            + "_"
            + sa
            + "_m"
            + str(mag1)
            + "-"
            + arcflag
            + ".csv",
            index=False,
        )


def get_period_and_filter(col_name, df):
    """Extract period and apply frequency filter based on column name."""
    if col_name == "PGA":
        return 0, df
    elif col_name == "PGV":
        return -1, df
    else:
        period = float(col_name.split("_")[1])
        fosc = 1 / period
        return period, df[df["fmin"] <= fosc]


def get_period_and_filter_v2(col_name, df):
    """Extract period and apply frequency filter based on column name."""
    if col_name == "resid_PGA":
        return 0, df
    elif col_name == "resid_PGV":
        return -1, df
    else:
        period = float(col_name.split("_")[2])
        fosc = 1 / period
        return period, df[df["fmin"] <= fosc]
