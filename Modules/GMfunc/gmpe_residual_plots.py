import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# import seaborn as sns
from scipy import stats
import warnings

warnings.filterwarnings("ignore")

# Set plotting style
plt.style.use("default")
# sns.set_palette("husl")


def prepare_residual_data(df):
    """
    Prepare residual data by calculating between-event and within-event residuals

    Parameters:
    df: DataFrame with columns ['event_id', 'site_id', 'residual', 'magnitude', 'distance']
    """
    # Calculate between-event residuals (δB)
    event_means = df.groupby("event_id")["residual"].mean()
    df["delta_B"] = df["event_id"].map(event_means)

    # Calculate within-event residuals (δW) - before removing site effects
    df["delta_W"] = df["residual"] - df["delta_B"]

    # Calculate site-to-site residuals (δS2S)
    # This is the average within-event residual for each site across all events
    site_means = df.groupby("site_id")["delta_W"].mean()
    # Re-center to zero mean
    site_means_centered = site_means - site_means.mean()
    df["delta_S2S"] = df["site_id"].map(site_means_centered)

    # Calculate within-event, single-station residual (δWS)
    # This is the remaining residual after removing both event and site effects
    df["delta_WS"] = df["delta_W"] - df["delta_S2S"]

    return df


def plot_residual_analysis(df, figsize=(15, 9)):
    """
    Create comprehensive residual analysis plots for GMPE studies
    """
    # Prepare data
    df = prepare_residual_data(df)

    fig, axes = plt.subplots(2, 3, figsize=figsize)
    fig.suptitle("GMPE Residual Analysis", fontsize=16, fontweight="bold")

    # 1. Between-event residuals vs Magnitude
    ax1 = axes[0, 0]
    event_data = (
        df.groupby("event_id")
        .agg({"delta_B": "first", "magnitude": "first"})
        .reset_index()
    )

    ax1.scatter(
        event_data["magnitude"],
        event_data["delta_B"],
        alpha=0.7,
        s=50,
        color="red",
        edgecolor="darkred",
        linewidth=0.5,
    )

    # Add trend line
    z = np.polyfit(event_data["magnitude"], event_data["delta_B"], 1)
    p = np.poly1d(z)
    ax1.plot(
        event_data["magnitude"],
        p(event_data["magnitude"]),
        "r--",
        alpha=0.8,
        linewidth=2,
        label=f"Trend: y={z[0]:.3f}x+{z[1]:.3f}",
    )

    ax1.axhline(y=0, color="black", linestyle="-", alpha=0.3)
    ax1.set_xlabel("Magnitude")
    ax1.set_ylabel("Between-event Residuals (δB)")
    ax1.set_title("δB vs Magnitude\n(Magnitude Scaling Check)")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2. Within-event residuals vs Distance
    ax2 = axes[0, 1]
    # Sample data for better visualization if too many points
    if len(df) > 5000:
        df_sample = df.sample(n=5000, random_state=42)
    else:
        df_sample = df

    scatter = ax2.scatter(
        df_sample["distance"],
        df_sample["delta_W"],
        alpha=0.5,
        s=20,
        c=df_sample["magnitude"],
        cmap="viridis",
        edgecolor="none",
    )

    # Add trend line
    z = np.polyfit(df["distance"], df["delta_W"], 1)
    p = np.poly1d(z)
    ax2.plot(
        df["distance"],
        p(df["distance"]),
        "r--",
        alpha=0.8,
        linewidth=2,
        label=f"Trend: y={z[0]:.4f}x+{z[1]:.3f}",
    )

    ax2.axhline(y=0, color="black", linestyle="-", alpha=0.3)
    ax2.set_xlabel("Distance (km)")
    ax2.set_ylabel("Within-event Residuals (δW)")
    ax2.set_title("δW vs Distance\n(Distance Attenuation Check)")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax2)
    cbar.set_label("Magnitude", rotation=270, labelpad=15)

    # 3. Within-event residuals vs Magnitude
    ax3 = axes[0, 2]
    scatter2 = ax3.scatter(
        df_sample["magnitude"],
        df_sample["delta_W"],
        alpha=0.5,
        s=20,
        c=df_sample["distance"],
        cmap="plasma",
        edgecolor="none",
    )

    # Add trend line
    z = np.polyfit(df["magnitude"], df["delta_W"], 1)
    p = np.poly1d(z)
    ax3.plot(
        df["magnitude"],
        p(df["magnitude"]),
        "r--",
        alpha=0.8,
        linewidth=2,
        label=f"Trend: y={z[0]:.3f}x+{z[1]:.3f}",
    )

    ax3.axhline(y=0, color="black", linestyle="-", alpha=0.3)
    ax3.set_xlabel("Magnitude")
    ax3.set_ylabel("Within-event Residuals (δW)")
    ax3.set_title("δW vs Magnitude\n(Mag-Distance Interaction)")
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # Add colorbar
    cbar2 = plt.colorbar(scatter2, ax=ax3)
    cbar2.set_label("Distance (km)", rotation=270, labelpad=15)

    # 4. Binned residuals vs Magnitude
    ax4 = axes[1, 0]
    mag_bins = np.arange(df["magnitude"].min(), df["magnitude"].max() + 0.5, 0.5)
    binned_stats = []
    bin_centers = []

    for i in range(len(mag_bins) - 1):
        mask = (df["magnitude"] >= mag_bins[i]) & (df["magnitude"] < mag_bins[i + 1])
        if mask.sum() > 0:
            residuals_in_bin = df.loc[mask, "residual"]
            bin_centers.append((mag_bins[i] + mag_bins[i + 1]) / 2)
            binned_stats.append(
                {
                    "mean": residuals_in_bin.mean(),
                    "std": residuals_in_bin.std(),
                    "count": len(residuals_in_bin),
                }
            )

    bin_centers = np.array(bin_centers)
    means = [stat["mean"] for stat in binned_stats]
    stds = [stat["std"] for stat in binned_stats]
    counts = [stat["count"] for stat in binned_stats]

    ax4.errorbar(
        bin_centers,
        means,
        yerr=stds,
        fmt="o-",
        capsize=5,
        capthick=2,
        linewidth=2,
        markersize=6,
    )
    ax4.axhline(y=0, color="black", linestyle="-", alpha=0.3)
    ax4.set_xlabel("Magnitude")
    ax4.set_ylabel("Mean Residual ± 1σ")
    ax4.set_title("Binned Residuals vs Magnitude")
    ax4.grid(True, alpha=0.3)

    # Add count annotations
    for i, (x, y, count) in enumerate(zip(bin_centers, means, counts)):
        ax4.annotate(
            f"n={count}",
            (x, y),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
            alpha=0.7,
        )

    # 5. Binned residuals vs Distance
    ax5 = axes[1, 1]
    # Use log-spaced bins for distance
    dist_bins = np.logspace(
        np.log10(max(1, df["distance"].min())), np.log10(df["distance"].max()), 15
    )
    binned_stats_dist = []
    bin_centers_dist = []

    for i in range(len(dist_bins) - 1):
        mask = (df["distance"] >= dist_bins[i]) & (df["distance"] < dist_bins[i + 1])
        if mask.sum() > 0:
            residuals_in_bin = df.loc[mask, "residual"]
            bin_centers_dist.append(np.sqrt(dist_bins[i] * dist_bins[i + 1]))
            binned_stats_dist.append(
                {
                    "mean": residuals_in_bin.mean(),
                    "std": residuals_in_bin.std(),
                    "count": len(residuals_in_bin),
                }
            )

    bin_centers_dist = np.array(bin_centers_dist)
    means_dist = [stat["mean"] for stat in binned_stats_dist]
    stds_dist = [stat["std"] for stat in binned_stats_dist]
    counts_dist = [stat["count"] for stat in binned_stats_dist]

    ax5.errorbar(
        bin_centers_dist,
        means_dist,
        yerr=stds_dist,
        fmt="o-",
        capsize=5,
        capthick=2,
        linewidth=2,
        markersize=6,
    )
    ax5.axhline(y=0, color="black", linestyle="-", alpha=0.3)
    ax5.set_xlabel("Distance (km)")
    ax5.set_ylabel("Mean Residual ± 1σ")
    ax5.set_title("Binned Residuals vs Distance")
    ax5.set_xscale("log")
    ax5.grid(True, alpha=0.3)

    # Add count annotations
    for i, (x, y, count) in enumerate(zip(bin_centers_dist, means_dist, counts_dist)):
        ax5.annotate(
            f"n={count}",
            (x, y),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
            alpha=0.7,
        )

    # 6. Residual distribution
    ax6 = axes[1, 2]
    ax6.hist(
        df["residual"],
        bins=30,
        alpha=0.7,
        color="skyblue",
        edgecolor="black",
        density=True,
    )

    # Overlay normal distribution
    mu, sigma = df["residual"].mean(), df["residual"].std()
    x = np.linspace(df["residual"].min(), df["residual"].max(), 100)
    ax6.plot(
        x,
        stats.norm.pdf(x, mu, sigma),
        "r-",
        linewidth=2,
        label=f"Normal(μ={mu:.3f}, σ={sigma:.3f})",
    )

    ax6.axvline(x=0, color="black", linestyle="-", alpha=0.3)
    ax6.set_xlabel("Residuals")
    ax6.set_ylabel("Density")
    ax6.set_title("Residual Distribution")
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig, axes


def plot_resid_s2s(df, figsize=(9, 4)):
    """
    Create comprehensive residual analysis plots for GMPE studies
    """
    # Prepare data
    df = prepare_residual_data(df)

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # 1. Between-site residuals vs Magnitude
    ax1 = axes[0]
    event_data = (
        df.groupby("site_id").agg({"delta_S2S": "first", "Vs30": "first"}).reset_index()
    )

    ax1.scatter(
        event_data["Vs30"],
        event_data["delta_S2S"],
        alpha=0.7,
        s=50,
        color="red",
        edgecolor="darkred",
        linewidth=0.5,
    )

    # Add trend line
    z = np.polyfit(event_data["Vs30"], event_data["delta_S2S"], 1)
    p = np.poly1d(z)
    ax1.plot(
        event_data["Vs30"],
        p(event_data["Vs30"]),
        "r--",
        alpha=0.8,
        linewidth=2,
        label=f"Trend: y={z[0]:.3f}x+{z[1]:.3f}",
    )

    ax1.axhline(y=0, color="black", linestyle="-", alpha=0.3)
    ax1.set_xlabel("Vs30 (m/s)")
    ax1.set_ylabel("Site-to-site Residuals (δS2S)")
    ax1.set_title("δS2S vs Vs30)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc=3)

    # 2. Within-event residuals vs Distance
    ax2 = axes[1]
    # Sample data for better visualization if too many points
    if len(df) > 5000:
        df_sample = df.sample(n=5000, random_state=42)
    else:
        df_sample = df

    scatter = ax2.scatter(
        df_sample["distance"],
        df_sample["delta_S2S"],
        alpha=0.5,
        s=20,
        c=df_sample["magnitude"],
        cmap="viridis",
        edgecolor="none",
    )

    # Add trend line
    z = np.polyfit(df["distance"], df["delta_S2S"], 1)
    p = np.poly1d(z)
    ax2.plot(
        df["distance"],
        p(df["distance"]),
        "r--",
        alpha=0.8,
        linewidth=2,
        label=f"Trend: y={z[0]:.4f}x+{z[1]:.3f}",
    )

    ax2.axhline(y=0, color="black", linestyle="-", alpha=0.3)
    ax2.set_xlabel("Distance (km)")
    ax2.set_ylabel("Within-site Residuals (δS2S)")
    ax2.set_title("δS2S vs Distance\n(Distance Attenuation Check)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc=3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax2)
    cbar.set_label("Magnitude", rotation=270, labelpad=15)

    return fig, axes


def prepare_df_dict_for_period_plot(df, period_columns, pga_column='PGA'):
    """
    Prepare a dictionary of DataFrames from a single DataFrame with multiple period columns

    Parameters:
    -----------
    df : DataFrame
        DataFrame containing residuals for multiple periods
        Must have columns: 'event_id', 'site_id', and period columns
    period_columns : list
        List of column names for different periods (e.g., ['SA_0.05', 'SA_0.1', 'SA_0.2', ...])
        Or list of tuples [(column_name, period_string), ...] for custom naming
    pga_column : str, optional
        Column name for PGA. Default is 'PGA'

    Returns:
    --------
    df_dict : dict
        Dictionary with period keys and DataFrames as values
        Each DataFrame has columns: 'event_id', 'site_id', 'residual', and any other columns from original df

    Examples:
    ---------
    # Example 1: Simple column names
    period_cols = ['SA_0.05', 'SA_0.1', 'SA_0.2', 'SA_0.5', 'SA_1.0', 'SA_2.0', 'SA_5.0', 'SA_10.0']
    df_dict = prepare_df_dict_for_period_plot(df, period_cols, 'PGA')

    # Example 2: Custom period naming
    period_cols = [('res_005', '0.05'), ('res_01', '0.1'), ('res_02', '0.2')]
    df_dict = prepare_df_dict_for_period_plot(df, period_cols, 'res_pga')
    """
    df_dict = {}

    # Get all non-period columns to preserve
    if isinstance(period_columns[0], tuple):
        # Custom naming: extract column names
        col_names = [col[0] for col in period_columns]
    else:
        col_names = period_columns

    # Add PGA column to the list
    all_period_cols = [pga_column] + col_names

    # Get columns to preserve (everything except the period residual columns)
    preserve_cols = [col for col in df.columns if col not in all_period_cols]

    # Process PGA
    if pga_column in df.columns:
        df_pga = df[preserve_cols].copy()
        df_pga['residual'] = df[pga_column]
        df_dict['PGA'] = df_pga

    # Process each period
    for period_col in period_columns:
        if isinstance(period_col, tuple):
            col_name, period_key = period_col
        else:
            col_name = period_col
            # Try to extract period from column name
            # Assumes format like 'SA_0.05', 'SA_0.1', 'res_0.05', etc.
            period_key = col_name.split('_')[-1]

        if col_name in df.columns:
            df_period = df[preserve_cols].copy()
            df_period['residual'] = df[col_name]
            df_dict[period_key] = df_period

    return df_dict


def plot_resid_band(df_dict, periods=None, figsize=(12, 8)):
    """
    Plot residual statistics (mean and std) as a function of period

    Parameters:
    -----------
    df_dict : dict
        Dictionary with period keys (e.g., 'PGA', '0.05', '0.1', etc.)
        Each value is a DataFrame with 'residual', 'event_id', 'site_id' columns
        Can be created using prepare_df_dict_for_period_plot()
    periods : list, optional
        List of period values to plot. If None, uses all periods from df_dict
    figsize : tuple
        Figure size (width, height)

    Returns:
    --------
    fig, axes : matplotlib figure and axes objects

    Examples:
    ---------
    # Method 1: From wide format DataFrame
    df_dict = prepare_df_dict_for_period_plot(df, period_columns, 'PGA')
    fig, axes = plot_resid_band(df_dict)

    # Method 2: Manual dictionary creation
    df_dict = {
        'PGA': df_pga,
        '0.05': df_005,
        '0.1': df_01,
        '0.2': df_02,
        '0.5': df_05,
        '1.0': df_10,
        '2.0': df_20,
        '5.0': df_50,
        '10.0': df_100
    }
    fig, axes = plot_resid_band(df_dict)
    """
    if periods is None:
        periods = sorted([k for k in df_dict.keys() if k != 'PGA'],
                        key=lambda x: float(x))
        if 'PGA' in df_dict:
            periods = ['PGA'] + periods

    # Prepare data - convert periods to numeric values for plotting
    period_values = []
    means = []
    stds = []
    taus = []
    phis = []

    for period in periods:
        if period not in df_dict:
            continue

        df = df_dict[period].copy()
        df = prepare_residual_data(df)

        # Get period value (PGA = 0.01 for plotting purposes)
        if period == 'PGA':
            period_val = 0.01
        else:
            period_val = float(period)

        period_values.append(period_val)
        means.append(df['residual'].mean())
        stds.append(df['residual'].std())

        # Calculate tau and phi
        event_residuals = df.groupby("event_id")["delta_B"].first()
        taus.append(event_residuals.std())
        phis.append(df['delta_W'].std())

    period_values = np.array(period_values)
    means = np.array(means)
    stds = np.array(stds)
    taus = np.array(taus)
    phis = np.array(phis)

    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    fig.suptitle("Residual Statistics vs Period", fontsize=16, fontweight="bold")

    # Plot 1: Mean residual vs period
    ax1 = axes[0, 0]
    ax1.semilogx(period_values, means, 'o-', linewidth=2, markersize=8, color='blue')
    ax1.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    ax1.set_xlabel("Period (s)")
    ax1.set_ylabel("Mean Residual")
    ax1.set_title("Mean Residual vs Period")
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0.008, 12)

    # Add PGA label
    if periods[0] == 'PGA':
        ax1.text(period_values[0], means[0], 'PGA',
                ha='right', va='bottom', fontsize=9)

    # Plot 2: Total std vs period
    ax2 = axes[0, 1]
    ax2.semilogx(period_values, stds, 'o-', linewidth=2, markersize=8, color='red')
    ax2.set_xlabel("Period (s)")
    ax2.set_ylabel("Total Standard Deviation (σ)")
    ax2.set_title("Total Std vs Period")
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0.008, 12)

    # Add PGA label
    if periods[0] == 'PGA':
        ax2.text(period_values[0], stds[0], 'PGA',
                ha='right', va='bottom', fontsize=9)

    # Plot 3: Tau and Phi vs period
    ax3 = axes[1, 0]
    ax3.semilogx(period_values, taus, 'o-', linewidth=2, markersize=8,
                label='τ (between-event)', color='red')
    ax3.semilogx(period_values, phis, 's-', linewidth=2, markersize=8,
                label='φ (within-event)', color='blue')
    ax3.set_xlabel("Period (s)")
    ax3.set_ylabel("Standard Deviation")
    ax3.set_title("τ and φ vs Period")
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0.008, 12)

    # Add PGA labels
    if periods[0] == 'PGA':
        ax3.text(period_values[0], taus[0], 'PGA',
                ha='right', va='bottom', fontsize=9)

    # Plot 4: Ratio phi/tau vs period
    ax4 = axes[1, 1]
    ratio = phis / taus
    ax4.semilogx(period_values, ratio, 'o-', linewidth=2, markersize=8, color='green')
    ax4.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='φ/τ = 1')
    ax4.set_xlabel("Period (s)")
    ax4.set_ylabel("φ/τ Ratio")
    ax4.set_title("Within-event to Between-event Ratio")
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(0.008, 12)

    # Add PGA label
    if periods[0] == 'PGA':
        ax4.text(period_values[0], ratio[0], 'PGA',
                ha='right', va='bottom', fontsize=9)

    plt.tight_layout()
    return fig, axes


def print_residual_statistics(df):
    """Print key statistics for residual analysis"""
    df = prepare_residual_data(df)

    print("=== RESIDUAL ANALYSIS STATISTICS ===")
    print(f"Total records: {len(df)}")
    print(f"Number of events: {df['event_id'].nunique()}")
    print(f"Number of sites: {df['site_id'].nunique()}")
    print()

    print("Overall Residuals:")
    print(f"  Mean: {df['residual'].mean():.4f}")
    print(f"  Std:  {df['residual'].std():.4f}")
    print()

    print("Between-event Residuals (δB):")
    event_residuals = df.groupby("event_id")["delta_B"].first()
    print(f"  Mean: {event_residuals.mean():.4f}")
    print(f"  Std (τ - between-event):  {event_residuals.std():.4f}")
    print()

    print("Within-event Residuals (δW):")
    print(f"  Mean: {df['delta_W'].mean():.4f}")
    print(f"  Std (φ - within-event):  {df['delta_W'].std():.4f}")
    print()

    print("Between-site Residuals (δS2S):")
    site_residuals = df.groupby("site_id")["delta_S2S"].first()
    print(f"  Mean: {site_residuals.mean():.4f}")
    print(f"  Std (τ(S2S) - between-site):  {site_residuals.std():.4f}")
    print()

    # Calculate total standard deviation from components (non-ergodic)
    tau = event_residuals.std()
    phi = df["delta_W"].std()

    sigma_total = np.sqrt(tau**2 + phi**2)

    print("Standard Deviation Components:")
    print(f"  τ (between-event std): {tau:.4f}")
    print(f"  φ (within-event std):  {phi:.4f}")
    print(f"  σ_total = √(τ² + φ²): {sigma_total:.4f}")

    print(f"  (Observed total std:   {df['residual'].std():.4f})")

    # Calculate total standard deviation from components (fully ergodic)
    tau = event_residuals.std()
    phi = df["delta_WS"].std()
    taus2s = site_residuals.std()

    sigma_total = np.sqrt(tau**2 + phi**2 + taus2s**2)

    print("Standard Deviation Components:")
    print(f"  τ (between-event std): {tau:.4f}")
    print(f"  φ (remaining within-event std):  {phi:.4f}")
    print(f"  τ(S2S) (between-site std): {taus2s:.4f}")
    print(f"  σ_total = √(τ² + φ²+s2s^2): {sigma_total:.4f}")
    print()

    # Correlation analysis
    print("Correlations with δB (between-event):")
    event_data = df.groupby("event_id").agg({"delta_B": "first", "magnitude": "first"})
    mag_corr = event_data["delta_B"].corr(event_data["magnitude"])
    print(f"  δB vs Magnitude: {mag_corr:.4f}")
    print()

    print("Correlations with δW (within-event):")
    dist_corr = df["delta_W"].corr(df["distance"])
    mag_corr_w = df["delta_W"].corr(df["magnitude"])
    print(f"  δW vs Distance:  {dist_corr:.4f}")
    print(f"  δW vs Magnitude: {mag_corr_w:.4f}")


# Example usage:
if __name__ == "__main__":
    # Create example data (replace this with your actual data loading)
    np.random.seed(42)
    n_records = 2000
    n_events = 50
    n_sites = 200

    # Generate synthetic data
    event_ids = np.random.choice(range(1, n_events + 1), n_records)
    site_ids = np.random.choice(range(1, n_sites + 1), n_records)
    magnitudes = np.random.uniform(4.0, 8.0, n_records)
    distances = np.random.lognormal(2.0, 1.5, n_records)

    # Generate residuals with some realistic bias patterns
    residuals = (
        0.1 * (magnitudes - 6.0)  # magnitude bias
        + 0.001 * (distances - 50)  # distance bias
        + np.random.normal(0, 0.3, n_records)
    )  # random scatter

    # Create DataFrame
    data = pd.DataFrame(
        {
            "event_id": event_ids,
            "site_id": site_ids,
            "magnitude": magnitudes,
            "distance": distances,
            "residual": residuals,
        }
    )

    # To use with your data, replace the above with:
    # data = pd.read_csv('your_residual_data.csv')
    # Make sure your DataFrame has columns: ['event_id', 'site_id', 'magnitude', 'distance', 'residual']

    # Generate plots
    fig, axes = plot_residual_analysis(data)
    plt.show()

    # Print statistics
    print_residual_statistics(data)
