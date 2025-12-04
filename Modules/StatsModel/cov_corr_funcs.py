import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def cov_to_corr(
    cov: pd.DataFrame | np.ndarray, labels: list[str] | None = None
) -> pd.DataFrame:
    """
    Convert a covariance matrix to a correlation matrix:
    R = D^{-1/2} Σ D^{-1/2}, where D = diag(Σ).
    Handles near-zero variances robustly.
    """
    if isinstance(cov, pd.DataFrame):
        arr = cov.to_numpy()
        labels = list(cov.index) if labels is None else labels
    else:
        cov = np.asarray(cov)
        arr = cov
        if labels is None:
            labels = [f"v{i}" for i in range(arr.shape[0])]

    if arr.shape[0] != arr.shape[1]:
        raise ValueError("Covariance matrix must be square.")

    # Diagonal (variances) and safe inverse sqrt
    d = np.diag(arr).astype(float)
    # Avoid division by ~0
    eps = np.finfo(float).eps
    inv_sqrt_d = 1.0 / np.sqrt(np.maximum(d, eps))

    # R = D^{-1/2} Σ D^{-1/2}
    R = (arr * inv_sqrt_d).T * inv_sqrt_d  # broadcasting
    # Numerical cleanup to keep diag exactly 1
    # np.fill_diagonal(R, 1.0)

    return pd.DataFrame(R, index=labels, columns=labels)


def plot_corr_heatmap(
    R: pd.DataFrame | np.ndarray,
    labels: list[str] | None = None,
    annotate: bool = True,
    fmt: str = ".2f",
    figsize=(6, 5),
    title: str = "Correlation matrix",
):
    """
    Plot a correlation matrix as a heat map using matplotlib (single axis).
    Shows only the upper triangle of the matrix.
    Keeps color scale fixed to [-1, 1] for comparability across runs.
    """
    if isinstance(R, pd.DataFrame):
        data = R.to_numpy().copy()
        labels = list(R.index) if labels is None else labels
    else:
        R = np.asarray(R)
        data = R.copy()
        if labels is None:
            labels = [f"v{i}" for i in range(data.shape[0])]

    n = data.shape[0]
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        data, vmin=-1, vmax=1, cmap="RdBu_r"
    )  # default colormap; no explicit colors
    # ax.set_title(title, pad=10)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)
    ax.set_xlabel("Variables")
    ax.set_ylabel("Variables")
    ax.yaxis.tick_right()
    # ax.xaxis.tick_bottom()

    # Add a colorbar
    # cbar = fig.colorbar(im, ax=ax,shrink=0.6)
    # cbar.set_label("Correlation")

    # Optional numeric annotations (only for upper triangle)
    if annotate:
        for i in range(n):
            for j in range(n):  # Only upper triangle (i <= j)
                ax.text(j, i, format(data[i, j], fmt), ha="center", va="center")

    # Tight layout to prevent label cutoff
    fig.tight_layout()
    return fig, ax


def plot_corr_heatmap_half(
    R: pd.DataFrame | np.ndarray,
    labels: list[str] | None = None,
    annotate: bool = True,
    fmt: str = ".2f",
    figsize=(6, 5),
    title: str = "Correlation matrix",
):
    """
    Plot a correlation matrix as a heat map using matplotlib (single axis).
    Shows only the upper triangle of the matrix.
    Keeps color scale fixed to [-1, 1] for comparability across runs.
    """
    if isinstance(R, pd.DataFrame):
        data = R.to_numpy().copy()
        labels = list(R.index) if labels is None else labels
    else:
        R = np.asarray(R)
        data = R.copy()
        if labels is None:
            labels = [f"v{i}" for i in range(data.shape[0])]

    # Mask the lower triangle (set to NaN for transparency)
    n = data.shape[0]
    mask = np.tril(np.ones((n, n), dtype=bool), k=-1)
    data_masked = data.copy()
    data_masked[mask] = np.nan

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        data_masked, vmin=-1, vmax=1, cmap="RdBu_r"
    )  # default colormap; no explicit colors
    # ax.set_title(title, pad=10)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels,rotation=45,ha="right")
    ax.set_xlabel("Variables")
    ax.set_ylabel("Variables")
    ax.yaxis.tick_right()
    ax.xaxis.tick_top()

    # Add a colorbar
    # cbar = fig.colorbar(im, ax=ax)
    # cbar.set_label("Correlation")

    # Optional numeric annotations (only for upper triangle)
    if annotate:
        for i in range(n):
            for j in range(i, n):  # Only upper triangle (i <= j)
                ax.text(j, i, format(data[i, j], fmt), ha="center", va="center")

    # Tight layout to prevent label cutoff
    fig.tight_layout()
    plt.show()

    return fig, ax


def plot_density_2d(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    bins: int | tuple[int, int] = 50,
    figsize=(8, 6),
    cmap: str = "viridis",
    title: str | None = None,
    show_scatter: bool = True,
    scatter_alpha: float = 0.3,
    scatter_size: float = 10,
):
    """
    Plot the 2D density distribution between two variables using hexbin or histogram.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the data
    x_col : str
        Name of the column for x-axis variable
    y_col : str
        Name of the column for y-axis variable
    bins : int or tuple of (int, int)
        Number of bins for 2D histogram. If int, same bins for both dimensions.
        Default is 50.
    figsize : tuple
        Figure size (width, height)
    cmap : str
        Colormap name for density visualization
    title : str or None
        Plot title. If None, generates title from variable names.
    show_scatter : bool
        Whether to overlay scatter points on top of density
    scatter_alpha : float
        Transparency of scatter points (0-1)
    scatter_size : float
        Size of scatter points

    Returns
    -------
    fig, ax : matplotlib figure and axis objects
    """
    # Extract data and remove NaN
    data = df[[x_col, y_col]].dropna()
    x = data[x_col].to_numpy()
    y = data[y_col].to_numpy()

    if len(x) == 0:
        raise ValueError(f"No valid data found for columns '{x_col}' and '{y_col}'")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create 2D histogram for density
    if isinstance(bins, int):
        bins = (bins, bins)

    # Overlay scatter plot if requested
    if show_scatter:
        ax.scatter(
            x,
            y,
            c="gray",
            s=scatter_size,
            alpha=scatter_alpha,
            edgecolors="none",
            rasterized=True,
        )
    # Plot 2D histogram with transparency
    h = ax.hist2d(x, y, bins=bins, cmap=cmap, alpha=0.5)

    # Add colorbar
    cbar = fig.colorbar(h[3], ax=ax)
    cbar.set_label("Count", rotation=270, labelpad=20)

    # Labels and title
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    if title is None:
        title = f"Density: {x_col} vs {y_col}"
    ax.set_title(title, pad=10)

    # Add correlation coefficient to plot
    corr = np.corrcoef(x, y)[0, 1]
    ax.text(
        0.05,
        0.95,
        f"r = {corr:.3f}",
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    fig.tight_layout()
    # plt.show()

    return fig, ax


def plot_density_hexbin(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    gridsize: int = 30,
    figsize=(10, 8),
    cmap: str = "viridis",
    title: str | None = None,
    reduce_C_function=np.mean,
    show_marginals: bool = True,
    hist_bins: int = 30,
    hist_color: str = "steelblue",
    hist_alpha: float = 0.7,
):
    """
    Plot the 2D density distribution using hexagonal binning (hexbin)
    with optional marginal histograms showing x and y distributions.

    Hexbin is often better for large datasets and creates a smoother appearance.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the data
    x_col : str
        Name of the column for x-axis variable
    y_col : str
        Name of the column for y-axis variable
    gridsize : int
        Number of hexagons in the x-direction (default 30)
    figsize : tuple
        Figure size (width, height). Default (10, 8) for marginal plots.
    cmap : str
        Colormap name for density visualization
    title : str or None
        Plot title. If None, generates title from variable names.
    reduce_C_function : callable
        Function to aggregate values in each hexagon (default: np.mean)
    show_marginals : bool
        If True, show marginal histograms for x and y distributions (default True)
    hist_bins : int
        Number of bins for marginal histograms (default 30)
    hist_color : str
        Color for marginal histograms (default "steelblue")
    hist_alpha : float
        Transparency for marginal histograms (default 0.7)

    Returns
    -------
    fig, axes : matplotlib figure and dictionary of axis objects
        If show_marginals=True, axes = {'main': ax_main, 'top': ax_top, 'right': ax_right}
        If show_marginals=False, axes = {'main': ax_main}
    """
    # Extract data and remove NaN
    data = df[[x_col, y_col]].dropna()
    x = data[x_col].to_numpy()
    y = data[y_col].to_numpy()

    if len(x) == 0:
        raise ValueError(f"No valid data found for columns '{x_col}' and '{y_col}'")

    if show_marginals:
        # Create figure with GridSpec for marginal plots
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(
            2, 2,
            width_ratios=(4, 1),
            height_ratios=(1, 4),
            wspace=0.05,
            hspace=0.05
        )

        # Create axes
        ax_main = fig.add_subplot(gs[1, 0])      # Main hexbin plot (bottom-left)
        ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)   # Top histogram (x distribution)
        ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main) # Right histogram (y distribution)

        # Hide tick labels for marginal plots
        ax_top.tick_params(labelbottom=False)
        ax_right.tick_params(labelleft=False)

        # Create hexbin plot on main axis
        hexbin = ax_main.hexbin(
            x,
            y,
            gridsize=gridsize,
            cmap=cmap,
            reduce_C_function=reduce_C_function,
            mincnt=1,
        )

        # Add colorbar
        cbar = fig.colorbar(hexbin, ax=ax_right, pad=0.1)
        cbar.set_label("Count", rotation=270, labelpad=20)

        # Labels for main plot
        ax_main.set_xlabel(x_col, fontsize=11)
        ax_main.set_ylabel(y_col, fontsize=11)

        # Add correlation coefficient to main plot
        corr = np.corrcoef(x, y)[0, 1]
        ax_main.text(
            0.05,
            0.95,
            f"r = {corr:.3f}",
            transform=ax_main.transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
            fontsize=10
        )

        # Top histogram (x distribution)
        ax_top.hist(
            x,
            bins=hist_bins,
            color=hist_color,
            alpha=hist_alpha,
            edgecolor='black',
            linewidth=0.5
        )
        ax_top.set_ylabel("Count", fontsize=9)
        ax_top.tick_params(labelsize=8)

        # Add mean line to top histogram
        ax_top.axvline(np.mean(x), color='red', linestyle='--',
                      linewidth=1.5, alpha=0.7, label=f'μ={np.mean(x):.2f}')
        ax_top.legend(fontsize=8, loc='upper right')

        # Right histogram (y distribution)
        ax_right.hist(
            y,
            bins=hist_bins,
            orientation='horizontal',
            color=hist_color,
            alpha=hist_alpha,
            edgecolor='black',
            linewidth=0.5
        )
        ax_right.set_xlabel("Count", fontsize=9)
        ax_right.tick_params(labelsize=8)

        # Add mean line to right histogram
        ax_right.axhline(np.mean(y), color='red', linestyle='--',
                        linewidth=1.5, alpha=0.7, label=f'μ={np.mean(y):.2f}')
        ax_right.legend(fontsize=8, loc='upper right')

        # Overall title
        if title is None:
            title = f"Density (hexbin): {x_col} vs {y_col}"
        fig.suptitle(title, fontsize=13, y=0.98)

        axes_dict = {'main': ax_main, 'top': ax_top, 'right': ax_right}

    else:
        # Create simple figure without marginals (original behavior)
        fig, ax_main = plt.subplots(figsize=figsize)

        # Create hexbin plot
        hexbin = ax_main.hexbin(
            x,
            y,
            gridsize=gridsize,
            cmap=cmap,
            reduce_C_function=reduce_C_function,
            mincnt=1,
        )

        # Add colorbar
        cbar = fig.colorbar(hexbin, ax=ax_main)
        cbar.set_label("Count", rotation=270, labelpad=20)

        # Labels and title
        ax_main.set_xlabel(x_col)
        ax_main.set_ylabel(y_col)
        if title is None:
            title = f"Density (hexbin): {x_col} vs {y_col}"
        ax_main.set_title(title, pad=10)

        # Add correlation coefficient to plot
        corr = np.corrcoef(x, y)[0, 1]
        ax_main.text(
            0.05,
            0.95,
            f"r = {corr:.3f}",
            transform=ax_main.transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        fig.tight_layout()

        axes_dict = {'main': ax_main}

    return fig, axes_dict


def plot_histogram_1d(
    df: pd.DataFrame,
    col: str,
    bins: int = 30,
    figsize=(6, 5),
    color: str = "steelblue",
    edgecolor: str = "black",
    alpha: float = 0.7,
    title: str | None = None,
    xlabel: str | None = None,
    density: bool = False,
    show_stats: bool = True,
    show_kde: bool = True,
    kde_color: str = "darkred",
    kde_linewidth: float = 2.5,
):
    """
    Plot a 1D histogram for a single variable with optional KDE overlay.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the data
    col : str
        Name of the column to plot
    bins : int
        Number of bins for the histogram (default 30)
    figsize : tuple
        Figure size (width, height)
    color : str
        Bar color (default "steelblue")
    edgecolor : str
        Edge color for bars (default "black")
    alpha : float
        Transparency of bars (0-1, default 0.7)
    title : str or None
        Plot title. If None, uses column name.
    xlabel : str or None
        X-axis label. If None, uses column name.
    density : bool
        If True, plot probability density instead of counts (default False)
    show_stats : bool
        If True, display mean and std on plot (default True)
    show_kde : bool
        If True, overlay a kernel density estimate curve (default True)
    kde_color : str
        Color of the KDE curve (default "darkred")
    kde_linewidth : float
        Width of the KDE curve line (default 2.5)

    Returns
    -------
    fig, ax : matplotlib figure and axis objects
    """
    # Extract data and remove NaN
    data = df[col].dropna().to_numpy()

    if len(data) == 0:
        raise ValueError(f"No valid data found for column '{col}'")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram
    n, bins_edges, patches = ax.hist(
        data,
        bins=bins,
        color=color,
        edgecolor=edgecolor,
        alpha=alpha,
        density=density,
    )

    # Add KDE (kernel density estimate) curve
    if show_kde:
        from scipy.stats import gaussian_kde

        # Compute KDE
        kde = gaussian_kde(data)
        x_range = np.linspace(data.min(), data.max(), 300)
        kde_values = kde(x_range)

        # If not density mode, scale KDE to match histogram counts
        if not density:
            # Scale KDE to match histogram by multiplying by n_samples * bin_width
            bin_width = bins_edges[1] - bins_edges[0]
            kde_values = kde_values * len(data) * bin_width

        # Plot KDE curve
        ax.plot(
            x_range,
            kde_values,
            color=kde_color,
            linewidth=kde_linewidth,
            label="Density",
            zorder=10,
        )

    # ax.legend(loc="upper right")

    # Labels and title
    if xlabel is None:
        xlabel = col
    ax.set_xlabel(xlabel)

    if density:
        ax.set_ylabel("Probability Density")
    else:
        ax.set_ylabel("Count")

    if title is None:
        title = f"Distribution of {col}"
    ax.set_title(title, pad=10)

    # Add statistics text box
    if show_stats:
        mean_val = np.mean(data)
        std_val = np.std(data)
        median_val = np.median(data)
        p20 = np.percentile(data, 20)
        p80 = np.percentile(data, 80)
        stats_text = (
            f"n = {len(data)}\n"
            f"Mean = {mean_val:.2f}\n"
            f"Std = {std_val:.2f}\n"
            f"Median = {median_val:.2f}\n"
            f"P20 = {p20:.2f}\n"
            f"P80 = {p80:.2f}"
        )

        ax.text(
            0.97,
            0.97,
            stats_text,
            transform=ax.transAxes,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
            fontsize=9,
        )

    # Add grid for better readability
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)

    fig.tight_layout()
    # plt.show()

    return fig, ax


def plot_histogram_circular(
    df: pd.DataFrame,
    col: str,
    bins: int = 36,
    figsize=(8, 8),
    color: str = "steelblue",
    alpha: float = 0.7,
    title: str | None = None,
    degrees: bool = True,
):
    """
    Plot a circular histogram (rose diagram) for directional data like rake angles.

    This is more appropriate than a linear histogram for directional/circular data
    as it properly represents the periodic nature (e.g., 0° = 360°).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the data
    col : str
        Name of the column with directional data
    bins : int
        Number of angular bins (default 36, i.e., 10° bins)
    figsize : tuple
        Figure size (width, height)
    color : str
        Bar color (default "steelblue")
    alpha : float
        Transparency of bars (0-1, default 0.7)
    title : str or None
        Plot title. If None, uses column name.
    degrees : bool
        If True, input data is in degrees. If False, in radians (default True)

    Returns
    -------
    fig, ax : matplotlib figure and axis objects
    """
    # Extract data and remove NaN
    data = df[col].dropna().to_numpy()

    if len(data) == 0:
        raise ValueError(f"No valid data found for column '{col}'")

    # Convert to radians if needed
    if degrees:
        theta = np.deg2rad(data)
    else:
        theta = data

    # Wrap to [0, 2π)
    theta = np.mod(theta, 2 * np.pi)

    # Create figure with polar projection
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection="polar"))

    # Create histogram
    counts, bin_edges = np.histogram(theta, bins=bins, range=(0, 2 * np.pi))

    # Bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    width = 2 * np.pi / bins

    # Plot bars
    bars = ax.bar(
        bin_centers,
        counts,
        width=width,
        color=color,
        edgecolor="black",
        alpha=alpha,
        linewidth=0.5,
    )

    # Configure polar plot
    ax.set_theta_zero_location("N")  # 0° at top
    ax.set_theta_direction(-1)  # Clockwise

    # Set degree labels
    ax.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
    ax.set_xticklabels(["0°", "45°", "90°", "135°", "180°", "225°", "270°", "315°"])

    # Title
    if title is None:
        title = f"Circular Distribution of {col}"
    ax.set_title(title, pad=20)

    # Add statistics text
    mean_angle = np.arctan2(np.sin(theta).mean(), np.cos(theta).mean())
    if mean_angle < 0:
        mean_angle += 2 * np.pi
    mean_angle_deg = np.rad2deg(mean_angle)

    # Circular variance (R statistic)
    R = np.sqrt(np.sin(theta).mean() ** 2 + np.cos(theta).mean() ** 2)
    circular_var = 1 - R

    stats_text = f"n = {len(data)}\nMean direction = {mean_angle_deg:.1f}°\nR = {R:.3f}\nCirc. var = {circular_var:.3f}"

    # Add text box outside the polar plot
    fig.text(
        0.02,
        0.02,
        stats_text,
        verticalalignment="bottom",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        fontsize=9,
    )

    fig.tight_layout()
    plt.show()

    return fig, ax


# --- Example usage ---
# Suppose you already computed a covariance matrix, e.g.:
# cov_aug = covariance_augmented(df, "SHmax_deg", "s2ratio", "xnuc")
# Convert to correlation and plot:
# R = cov_to_corr(cov_aug)
# plot_corr_heatmap(R, title="Correlation (augmented variables)")


import numpy as np
import pandas as pd


# ---------- Circular helpers (axial angle = orientation, period = 180°) ----------
def _axial_radians(shmax_deg: pd.Series) -> np.ndarray:
    """
    Convert SHmax in degrees to axial radians (mod π). SHmax is an orientation
    with 180° ambiguity, so 0° ≡ 180°.
    """
    theta = np.deg2rad(pd.Series(shmax_deg).to_numpy())
    return np.mod(theta, np.pi)


def axial_components(shmax_deg: pd.Series) -> pd.DataFrame:
    """
    Map axial angle θ to the unit circle using the double-angle embedding:
    u = cos(2θ), v = sin(2θ). This preserves the 180° symmetry.
    Returns a DataFrame with columns ['SHmax_c2','SHmax_s2'].
    """
    theta = _axial_radians(shmax_deg)
    return pd.DataFrame(
        {
            "SHmax_c2": np.cos(2.0 * theta),
            "SHmax_s2": np.sin(2.0 * theta),
        }
    )


def _to_radians_directional(angles_deg: pd.Series) -> np.ndarray:
    """
    Wrap directional angles given in degrees to [0, 2π) and return radians.
    Suitable for rake in [0, 360).
    """
    theta = np.deg2rad(pd.Series(angles_deg).to_numpy())
    return np.mod(theta, 2 * np.pi)


def covariance_with_cos_rake(
    df: pd.DataFrame,
    rake_col: str = "rake_deg",
    extra_cols: list[str] = ("s2ratio", "xnuc"),
    harmonic: int = 1,  # k = 1 uses cos(rake); try k=2 for cos(2*rake) if you expect 180° periodicity
    ddof: int = 1,
) -> pd.DataFrame:
    """
    Build covariance of [ cos(k·rake), <extra parameters...> ] where rake is directional (0–360°).
    Parameters
    ----------
    df : DataFrame containing rake and extra parameter columns.
    rake_col : name of the rake column in degrees (0–360).
    extra_cols : list of additional (linear) parameter columns to include.
    harmonic : integer k ≥ 1. Use 1 for cos(rake). Use 2 to test 180° periodic patterns, etc.
    ddof : degrees of freedom for covariance (1 gives sample covariance).
    """
    if harmonic < 1:
        raise ValueError("harmonic must be an integer ≥ 1")

    cols_needed = [rake_col] + list(extra_cols)
    data = df[cols_needed].dropna().copy()

    # either use cos(angle) or angle
    # theta = _to_radians_directional(data[rake_col])
    # cos_term = np.cos(harmonic * theta)
    cos_term = data[rake_col]

    X = np.column_stack([cos_term] + [data[c].to_numpy() for c in extra_cols])
    cov = np.cov(X, rowvar=False, ddof=ddof)

    colnames = [f"cos({harmonic}*{rake_col})"] + list(extra_cols)
    return pd.DataFrame(cov, index=colnames, columns=colnames)


def circular_strength_axial(shmax_deg: pd.Series) -> dict:
    """
    Resultant length R (on the doubled angles) and circular variance (1-R)
    for axial data. R≈1 means tightly clustered orientations; R≈0 diffuse.
    """
    comp = axial_components(shmax_deg)
    c_bar, s_bar = comp["SHmax_c2"].mean(), comp["SHmax_s2"].mean()
    R = float(np.hypot(c_bar, s_bar))
    return {"R_double_angle": R, "circular_variance_axial": 1.0 - R}


# ---------- Covariance matrices ----------
def covariance_augmented(
    df: pd.DataFrame, sh_col="SHmax_deg", s2_col="s2ratio", x_col="xnuc", ddof=1
) -> pd.DataFrame:
    """
    Circular-aware covariance of the augmented vector:
    [cos(2*SHmax), sin(2*SHmax), s2ratio, xnuc].
    """
    data = df[[sh_col, s2_col, x_col]].dropna()
    comp = axial_components(data[sh_col])
    X = np.column_stack(
        [
            comp["SHmax_c2"],
            comp["SHmax_s2"],
            data[s2_col].to_numpy(),
            data[x_col].to_numpy(),
        ]
    )
    cov = np.cov(X, rowvar=False, ddof=ddof)
    cols = ["SHmax_c2", "SHmax_s2", s2_col, x_col]
    return pd.DataFrame(cov, index=cols, columns=cols)


def covariance_naive_linear(
    df: pd.DataFrame, sh_col="SHmax_deg", s2_col="s2ratio", x_col="xnuc", ddof=1
) -> pd.DataFrame:
    """
    Naïve (reference-only) covariance treating SHmax as linear after unwrapping.
    Not recommended for inference, but useful to see the effect of circularity.
    """
    data = df[[sh_col, s2_col, x_col]].dropna()
    theta = _axial_radians(data[sh_col])
    theta_unwrap = np.unwrap(theta)  # unwrap with period π
    sh_unwrap_deg = np.rad2deg(theta_unwrap)
    X = np.column_stack(
        [sh_unwrap_deg, data[s2_col].to_numpy(), data[x_col].to_numpy()]
    )
    cov = np.cov(X, rowvar=False, ddof=ddof)
    cols = [f"{sh_col}_unwrapped_deg", s2_col, x_col]
    return pd.DataFrame(cov, index=cols, columns=cols)


# ---------- Circular–linear association (Mardia, 1976) ----------
def circular_linear_association(shmax_deg: pd.Series, x: pd.Series) -> dict:
    """
    Mardia's circular–linear association:
    ρ_cl^2 = (r_xc^2 + r_xs^2 - 2 r_xc r_xs r_cs) / (1 - r_cs^2),
    where c=cos(θ), s=sin(θ) with θ axial (mod π).
    Returns components and ρ_cl (≥0).
    """
    theta = _axial_radians(shmax_deg)
    c, s = np.cos(theta), np.sin(theta)
    x = pd.Series(x).to_numpy()

    def _corr(a, b):
        A = a - a.mean()
        B = b - b.mean()
        denom = np.sqrt((A**2).sum() * (B**2).sum())
        return float((A @ B) / denom) if denom > 0 else np.nan

    r_xc = _corr(x, c)
    r_xs = _corr(x, s)
    r_cs = _corr(c, s)
    denom = max(1.0 - r_cs**2, 1e-12)
    rho2 = max((r_xc**2 + r_xs**2 - 2.0 * r_xc * r_xs * r_cs) / denom, 0.0)
    return {"r_xc": r_xc, "r_xs": r_xs, "r_cs": r_cs, "rho_cl": float(np.sqrt(rho2))}


def summarize_interdependency(
    df: pd.DataFrame, sh_col="SHmax_deg", s2_col="s2ratio", x_col="xnuc"
) -> pd.DataFrame:
    """
    Small summary table combining:
    - circular–linear association of SHmax with s2ratio and xnuc;
    - Pearson correlation between s2ratio and xnuc.
    """
    a1 = circular_linear_association(df[sh_col], df[s2_col])
    a2 = circular_linear_association(df[sh_col], df[x_col])
    r_lin = df[[s2_col, x_col]].dropna().corr().iloc[0, 1]
    return pd.DataFrame(
        [
            {"pair": f"{sh_col} ↔ {s2_col}", "rho_cl": a1["rho_cl"]},
            {"pair": f"{sh_col} ↔ {x_col}", "rho_cl": a2["rho_cl"]},
            {"pair": f"{s2_col} ↔ {x_col}", "pearson_r": r_lin},
        ]
    )


from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def plot_pca_analysis(data_dict, variables=['stress', 'slip', 'vr'],
                      labels=['stress', 'slip', 'vr']):
    """
    Perform PCA and plot results.

    Parameters
    ----------
    data_dict : dict
        Dictionary with variable names as keys and numpy arrays as values
    variables : list of str
        List of variable names to include in PCA

    Returns
    -------
    fig, axes : matplotlib objects
    pca : sklearn PCA object

    Interpretation of PC1 vs PC2 plot:
    - Clustering: Indicates distinct physical regimes or failure modes
    - Spread along PC1: Main mode of variation (highest variance direction)
    - Spread along PC2: Second mode of variation (orthogonal to PC1)
    - Patterns reveal relationships between variables that are hidden in high dimensions
    """
    # Prepare data matrix
    data_list = [data_dict[var].flatten() for var in variables]
    X = np.column_stack(data_list)

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # PCA
    pca = PCA()
    X_pca = pca.fit_transform(X_scaled)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # 1. Explained variance
    axes[0].bar(range(1, len(variables)+1), pca.explained_variance_ratio_,color='tomato')
    axes[0].set_xlabel('Principal Component', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Explained Variance Ratio', fontsize=12, fontweight='bold')
    axes[0].set_title('Scree Plot', fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3, axis='y')

    # Add cumulative variance line
    cumvar = np.cumsum(pca.explained_variance_ratio_)
    ax_twin = axes[0].twinx()
    ax_twin.plot(range(1, len(variables)+1), cumvar, 'ro-', linewidth=2, markersize=6)
    ax_twin.set_ylabel('Cumulative Variance', fontsize=11, color='tomato')
    ax_twin.tick_params(axis='y', labelcolor='tomato')

    # 2. Component loadings
    loadings = pca.components_.T
    im = axes[1].imshow(loadings, cmap='BrBG_r', aspect='auto', vmin=-1, vmax=1)
    axes[1].set_xticks(range(len(variables)))
    axes[1].set_xticklabels([f'PC{i+1}' for i in range(len(variables))])
    axes[1].set_yticks(range(len(variables)))
    axes[1].set_yticklabels(labels,rotation=45,ha='right')

    axes[1].set_title('Component Loadings', fontsize=13, fontweight='bold')
    plt.colorbar(im, ax=axes[1], label='Loading')

    # Add text annotations for loadings
    for i in range(len(variables)):
        for j in range(min(3, len(variables))):  # Only first 3 PCs
            text = axes[1].text(j, i, f'{loadings[i, j]:.2f}',
                              ha="center", va="center", color="black", fontsize=9)

    # 3. PC1 vs PC2
    axes[2].scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.5, s=30, c='forestgreen', edgecolors='none')
    axes[2].axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)
    axes[2].axvline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)
    axes[2].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})',
                      fontsize=12, fontweight='bold')
    axes[2].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})',
                      fontsize=12, fontweight='bold')
    axes[2].set_title('First Two Principal Components', fontsize=13, fontweight='bold')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    return fig, axes, pca


def interpret_pca_loadings(pca, variables, n_components=2):
    """
    Print interpretation of PCA loadings to understand what each PC represents.

    Parameters
    ----------
    pca : sklearn PCA object
        Fitted PCA model
    variables : list of str
        Variable names
    n_components : int
        Number of components to interpret (default: 2)

    Example
    -------
    >>> fig, axes, pca = plot_pca_analysis(data_dict, variables)
    >>> interpret_pca_loadings(pca, variables)
    """
    print("="*70)
    print("PCA COMPONENT INTERPRETATION")
    print("="*70)

    loadings = pca.components_.T

    for i in range(min(n_components, len(variables))):
        print(f"\nPC{i+1} (Explains {pca.explained_variance_ratio_[i]:.1%} of variance):")
        print("-" * 70)

        # Sort variables by absolute loading
        loading_order = np.argsort(np.abs(loadings[:, i]))[::-1]

        for j in loading_order:
            loading_val = loadings[j, i]
            var_name = variables[j]

            # Interpretation
            if abs(loading_val) > 0.5:
                strength = "STRONG"
            elif abs(loading_val) > 0.3:
                strength = "moderate"
            else:
                strength = "weak"

            direction = "positive" if loading_val > 0 else "negative"

            print(f"  {var_name:15s}: {loading_val:+.3f}  ({strength:8s} {direction:8s})")

        # Suggest interpretation
        print(f"\n  → Interpretation:")
        strong_vars = [variables[j] for j in range(len(variables))
                      if abs(loadings[j, i]) > 0.5]
        if strong_vars:
            print(f"    PC{i+1} primarily represents variation in: {', '.join(strong_vars)}")

    # Overall summary
    total_var = np.sum(pca.explained_variance_ratio_[:n_components])
    print(f"\n{'='*70}")
    print(f"First {n_components} PCs explain {total_var:.1%} of total variance")
    print(f"{'='*70}\n")


# nonlinear pca decomposition

from sklearn.decomposition import KernelPCA

def plot_kernel_pca_analysis(data_dict, variables=['stress', 'slip', 'vr'],
                              kernel='rbf', gamma=None, n_components=None):
    """
    Perform Kernel PCA (non-linear dimensionality reduction) and plot results.

    Kernel PCA can capture non-linear relationships that linear PCA might miss.

    Parameters
    ----------
    data_dict : dict
        Dictionary with variable names as keys and numpy arrays as values
    variables : list of str
        List of variable names to include in Kernel PCA
    kernel : str
        Kernel type: 'rbf' (default), 'poly', 'sigmoid', 'cosine', 'linear'
        - 'rbf': Radial Basis Function, good for general non-linear patterns
        - 'poly': Polynomial, good for polynomial relationships
        - 'sigmoid': Sigmoid, similar to neural network activation
        - 'linear': Equivalent to standard PCA
    gamma : float, optional
        Kernel coefficient for rbf, poly and sigmoid kernels.
        If None, defaults to 1/n_features
    n_components : int, optional
        Number of components to compute. If None, uses len(variables)

    Returns
    -------
    fig : matplotlib figure
        The figure containing all plots
    axes : array of matplotlib axes
        Array of subplot axes for customization
    kpca : sklearn.decomposition.KernelPCA object
        The FITTED Kernel PCA model. This is the "transformer" object.
        Contains:
        - kpca.lambdas_ : eigenvalues (importance of each component)
        - kpca.alphas_ : eigenvectors in kernel space
        Use this to: transform NEW data with same kernel transformation
    X_kpca : numpy array, shape (n_samples, n_components)
        The TRANSFORMED data in kernel principal component space.
        This is your original data projected onto the kernel PCs.
        - Each row = one data point (e.g., one earthquake fault element)
        - Each column = one kernel principal component (KPC1, KPC2, ...)
        - X_kpca[:, 0] = all KPC1 values
        - X_kpca[:, 1] = all KPC2 values
        Use this for: visualization, clustering, regression, classification

    Notes
    -----
    Kernel PCA doesn't have explicit loadings like linear PCA because the
    transformation is performed in a high-dimensional feature space. However,
    this function approximates "pseudo-loadings" by computing the correlation
    between each kernel PC and the original variables. These correlations help
    interpret which original variables contribute most to each kernel component.

    The three plots show:
    1. Scree plot: Variance explained by each kernel PC
    2. Component loadings: Correlation heatmap (how original variables relate to KPCs)
    3. PC scatter: KPC1 vs KPC2 visualization

    What are kpca and X_kpca?
    -------------------------
    Think of it like this:

    kpca (the model):
        - Like a "recipe" or "transformer"
        - Stores how to do the kernel transformation
        - Use: kpca.transform(new_data) to transform new data
        - Example: Train on 2020 earthquakes, apply to 2021 earthquakes

    X_kpca (the transformed data):
        - Your actual data in the new kernel PC space
        - Shape: (n_samples, n_components)
        - Ready to use for analysis, plotting, modeling
        - Example: Use for clustering similar fault behaviors

    Example
    -------
    >>> # Basic usage
    >>> fig, axes, kpca, X_kpca = plot_kernel_pca_analysis(
    ...     data_dict, variables=['rake', 'ASl', 'Vr', 'PSR', 'T0', 'Pf'],
    ...     kernel='rbf'
    ... )
    >>>
    >>> # What is X_kpca?
    >>> print(X_kpca.shape)  # (100000, 6) = 100k samples, 6 kernel PCs
    >>> print(X_kpca[:, 0])  # KPC1 values for all samples
    >>> print(X_kpca[:, 1])  # KPC2 values for all samples
    >>>
    >>> # Use X_kpca for clustering
    >>> from sklearn.cluster import KMeans
    >>> kmeans = KMeans(n_clusters=3)
    >>> clusters = kmeans.fit_predict(X_kpca[:, :2])  # Cluster using KPC1, KPC2
    >>>
    >>> # Use kpca to transform NEW data
    >>> new_data = np.column_stack([new_rake, new_ASl, new_Vr, ...])
    >>> from sklearn.preprocessing import StandardScaler
    >>> scaler = StandardScaler()
    >>> scaler.fit(X)  # Fit on original data
    >>> new_data_scaled = scaler.transform(new_data)
    >>> new_data_kpca = kpca.transform(new_data_scaled)  # Apply same transformation
    """
    # Prepare data matrix
    data_list = [data_dict[var].flatten() for var in variables]
    X = np.column_stack(data_list)

    # Standardize (important for kernel methods)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Set default number of components
    if n_components is None:
        n_components = min(len(variables), X_scaled.shape[0] - 1)

    # Set default gamma for rbf kernel
    if gamma is None and kernel in ['rbf', 'poly', 'sigmoid']:
        gamma = 1.0 / len(variables)

    # Kernel PCA with error handling
    try:
        kpca = KernelPCA(n_components=n_components, kernel=kernel, gamma=gamma,
                         fit_inverse_transform=True, eigen_solver='dense')
        X_kpca = kpca.fit_transform(X_scaled)
    except np.linalg.LinAlgError:
        # Try without inverse transform if singular
        print(f"Warning: Singular matrix detected, disabling inverse_transform")
        kpca = KernelPCA(n_components=n_components, kernel=kernel, gamma=gamma,
                         fit_inverse_transform=False, eigen_solver='arpack')
        X_kpca = kpca.fit_transform(X_scaled)

    # Try to compute explained variance (not always available for all kernels)
    try:
        # Reconstruct to estimate variance
        X_reconstructed = kpca.inverse_transform(X_kpca)
        reconstruction_error = np.mean((X_scaled - X_reconstructed) ** 2, axis=0)

        # Compute variance explained by each component
        total_var = np.var(X_scaled, axis=0).sum()
        var_per_component = np.var(X_kpca, axis=0)
        explained_var_ratio = var_per_component / total_var

        has_variance_info = True
    except:
        explained_var_ratio = None
        has_variance_info = False
        print(f"Note: Explained variance not available for kernel='{kernel}'")

    # Create plots
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # 1. Explained variance (if available)
    if has_variance_info:
        axes[0].bar(range(1, n_components+1), explained_var_ratio, color='steelblue')
        axes[0].set_xlabel('Kernel PC', fontsize=12, fontweight='bold')
        axes[0].set_ylabel('Explained Variance Ratio', fontsize=12, fontweight='bold')
        axes[0].set_title(f'Kernel PCA Scree Plot ({kernel} kernel)', fontsize=13, fontweight='bold')
        axes[0].grid(True, alpha=0.3, axis='y')

        # Add cumulative variance line
        cumvar = np.cumsum(explained_var_ratio)
        ax_twin = axes[0].twinx()
        ax_twin.plot(range(1, n_components+1), cumvar, 'ro-', linewidth=2, markersize=6)
        ax_twin.set_ylabel('Cumulative Variance', fontsize=11, color='red')
        ax_twin.tick_params(axis='y', labelcolor='red')
    else:
        axes[0].text(0.5, 0.5, f'Variance information\nnot available\nfor {kernel} kernel',
                    ha='center', va='center', transform=axes[0].transAxes,
                    fontsize=12, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        axes[0].set_title(f'Kernel PCA ({kernel} kernel)', fontsize=13, fontweight='bold')
        axes[0].axis('off')

    # 2. Component loadings (correlation between KPCs and original variables)
    # For Kernel PCA, we approximate "loadings" by computing correlations
    # between kernel PCs and original (scaled) variables
    from scipy.stats import pearsonr

    # Compute correlations (pseudo-loadings)
    n_components_to_show = min(n_components, len(variables))
    loadings_approx = np.zeros((len(variables), n_components_to_show))

    for i in range(len(variables)):
        for j in range(n_components_to_show):
            # Correlation between original variable i and kernel PC j
            corr, _ = pearsonr(X_scaled[:, i], X_kpca[:, j])
            loadings_approx[i, j] = corr

    # Plot loadings heatmap
    im = axes[1].imshow(loadings_approx, cmap='BrBG_r', aspect='auto', vmin=-1, vmax=1)
    axes[1].set_xticks(range(n_components_to_show))
    axes[1].set_xticklabels([f'KPC{i+1}' for i in range(n_components_to_show)], fontsize=10)
    axes[1].set_yticks(range(len(variables)))
    axes[1].set_yticklabels(variables, fontsize=10)
    axes[1].set_title('Component Loadings\n(Correlation with Original Variables)',
                     fontsize=13, fontweight='bold')

    # Add colorbar
    cbar = plt.colorbar(im, ax=axes[1])
    cbar.set_label('Correlation', fontsize=10)

    # Add correlation values as text annotations
    for i in range(len(variables)):
        for j in range(n_components_to_show):
            text_color = 'white' if abs(loadings_approx[i, j]) > 0.5 else 'black'
            axes[1].text(j, i, f'{loadings_approx[i, j]:.2f}',
                        ha='center', va='center', color=text_color, fontsize=9)

    # 3. KPC1 vs KPC2 with density
    if has_variance_info:
        title_suffix = f'({explained_var_ratio[0]:.1%}, {explained_var_ratio[1]:.1%})'
    else:
        title_suffix = ''

    axes[2].scatter(X_kpca[:, 0], X_kpca[:, 1], alpha=0.5, s=30,
                   c='forestgreen', edgecolors='none')
    axes[2].axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)
    axes[2].axvline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)
    axes[2].set_xlabel(f'KPC1', fontsize=12, fontweight='bold')
    axes[2].set_ylabel(f'KPC2', fontsize=12, fontweight='bold')
    axes[2].set_title(f'First Two Kernel PCs {title_suffix}', fontsize=13, fontweight='bold')
    axes[2].grid(True, alpha=0.3)

    # Add kernel info to figure
    kernel_info = f"Kernel: {kernel}"
    if gamma is not None:
        kernel_info += f", γ={gamma:.4f}"
    fig.suptitle(kernel_info, fontsize=11, y=0.98)

    plt.tight_layout()
    return fig, axes, kpca, X_kpca


def compare_linear_vs_kernel_pca(data_dict, variables=['stress', 'slip', 'vr'],
                                  kernel='rbf', gamma=None):
    """
    Compare linear PCA vs Kernel PCA side by side.

    This helps identify whether non-linear relationships exist in the data.

    Parameters
    ----------
    data_dict : dict
        Dictionary with variable names as keys and numpy arrays as values
    variables : list of str
        List of variable names to include
    kernel : str
        Kernel type for Kernel PCA (default: 'rbf')
    gamma : float, optional
        Kernel coefficient

    Returns
    -------
    fig : matplotlib figure
    (pca, kpca) : tuple of fitted models
    (X_pca, X_kpca) : tuple of transformed data

    Example
    -------
    >>> fig, (pca, kpca), (X_pca, X_kpca) = compare_linear_vs_kernel_pca(
    ...     data_dict, variables=['rake', 'ASl', 'Vr', 'PSR', 'T0', 'Pf'])
    """
    # Prepare data
    data_list = [data_dict[var].flatten() for var in variables]
    X = np.column_stack(data_list)

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Set number of components
    n_comp = min(len(variables), X_scaled.shape[0] - 1)

    # Linear PCA
    pca = PCA(n_components=n_comp)
    X_pca = pca.fit_transform(X_scaled)

    # Set default gamma for rbf kernel
    if gamma is None and kernel in ['rbf', 'poly', 'sigmoid']:
        gamma = 1.0 / len(variables)

    # Kernel PCA with error handling
    try:
        kpca = KernelPCA(n_components=n_comp, kernel=kernel, gamma=gamma,
                         fit_inverse_transform=True, eigen_solver='dense')
        X_kpca = kpca.fit_transform(X_scaled)
    except np.linalg.LinAlgError:
        # Try without inverse transform if singular
        print(f"Warning: Singular matrix detected, disabling inverse_transform")
        kpca = KernelPCA(n_components=n_comp, kernel=kernel, gamma=gamma,
                         fit_inverse_transform=False, eigen_solver='arpack')
        X_kpca = kpca.fit_transform(X_scaled)

    # Create comparison plot
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # Row 1: Linear PCA
    # Scree plot
    var_ratio_pca = pca.explained_variance_ratio_
    axes[0, 0].bar(range(1, len(variables)+1), var_ratio_pca, color='tomato', alpha=0.7)
    axes[0, 0].set_xlabel('PC', fontsize=11, fontweight='bold')
    axes[0, 0].set_ylabel('Variance Ratio', fontsize=11, fontweight='bold')
    axes[0, 0].set_title('Linear PCA - Scree Plot', fontsize=12, fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3, axis='y')

    # Loadings
    loadings = pca.components_.T
    im1 = axes[0, 1].imshow(loadings, cmap='BrBG_r', aspect='auto', vmin=-1, vmax=1)
    axes[0, 1].set_xticks(range(len(variables)))
    axes[0, 1].set_xticklabels([f'PC{i+1}' for i in range(len(variables))], fontsize=9)
    axes[0, 1].set_yticks(range(len(variables)))
    axes[0, 1].set_yticklabels(variables, fontsize=9)
    axes[0, 1].set_title('Linear PCA - Loadings', fontsize=12, fontweight='bold')
    plt.colorbar(im1, ax=axes[0, 1], label='Loading')

    # PC1 vs PC2
    axes[0, 2].scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.5, s=5, c='royalblue', edgecolors='none')
    axes[0, 2].axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)
    axes[0, 2].axvline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)
    axes[0, 2].set_xlabel(f'PC1 ({var_ratio_pca[0]:.1%})', fontsize=11, fontweight='bold')
    axes[0, 2].set_ylabel(f'PC2 ({var_ratio_pca[1]:.1%})', fontsize=11, fontweight='bold')
    axes[0, 2].set_title('Linear PCA - PC1 vs PC2', fontsize=12, fontweight='bold')
    axes[0, 2].grid(True, alpha=0.3)

    # Row 2: Kernel PCA
    try:
        # Estimate variance
        X_reconstructed = kpca.inverse_transform(X_kpca)
        total_var = np.var(X_scaled, axis=0).sum()
        var_per_component = np.var(X_kpca, axis=0)
        var_ratio_kpca = var_per_component / total_var

        axes[1, 0].bar(range(1, len(variables)+1), var_ratio_kpca, color='steelblue', alpha=0.7)
        axes[1, 0].set_xlabel('KPC', fontsize=11, fontweight='bold')
        axes[1, 0].set_ylabel('Variance Ratio', fontsize=11, fontweight='bold')
        axes[1, 0].set_title(f'Kernel PCA ({kernel}) - Scree Plot', fontsize=12, fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3, axis='y')
    except:
        axes[1, 0].text(0.5, 0.5, 'Variance not available', ha='center', va='center',
                       transform=axes[1, 0].transAxes, fontsize=11)
        axes[1, 0].set_title(f'Kernel PCA ({kernel})', fontsize=12, fontweight='bold')
        var_ratio_kpca = None

    # No loadings for kernel PCA (non-linear transformation)
    axes[1, 1].text(0.5, 0.5, 'Kernel PCA uses\nnon-linear transformation\n(no explicit loadings)',
                   ha='center', va='center', transform=axes[1, 1].transAxes,
                   fontsize=11, bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    axes[1, 1].set_title('Kernel PCA - Feature Space', fontsize=12, fontweight='bold')
    axes[1, 1].axis('off')

    # KPC1 vs KPC2
    axes[1, 2].scatter(X_kpca[:, 0], X_kpca[:, 1], alpha=0.5, s=5, c='forestgreen', edgecolors='none')
    axes[1, 2].axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)
    axes[1, 2].axvline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)

    if var_ratio_kpca is not None:
        axes[1, 2].set_xlabel(f'KPC1 ({var_ratio_kpca[0]:.1%})', fontsize=11, fontweight='bold')
        axes[1, 2].set_ylabel(f'KPC2 ({var_ratio_kpca[1]:.1%})', fontsize=11, fontweight='bold')
    else:
        axes[1, 2].set_xlabel('KPC1', fontsize=11, fontweight='bold')
        axes[1, 2].set_ylabel('KPC2', fontsize=11, fontweight='bold')

    axes[1, 2].set_title('Kernel PCA - KPC1 vs KPC2', fontsize=12, fontweight='bold')
    axes[1, 2].grid(True, alpha=0.3)

    plt.tight_layout()

    # Print comparison summary
    print("="*70)
    print("LINEAR PCA vs KERNEL PCA COMPARISON")
    print("="*70)
    print(f"\nLinear PCA:")
    print(f"  PC1 + PC2 variance: {var_ratio_pca[0] + var_ratio_pca[1]:.1%}")

    if var_ratio_kpca is not None:
        print(f"\nKernel PCA ({kernel}):")
        print(f"  KPC1 + KPC2 variance: {var_ratio_kpca[0] + var_ratio_kpca[1]:.1%}")

        if var_ratio_kpca[0] + var_ratio_kpca[1] > var_ratio_pca[0] + var_ratio_pca[1] + 0.05:
            print(f"\n  → Kernel PCA captures MORE variance")
            print(f"    This suggests significant NON-LINEAR relationships in data")
        else:
            print(f"\n  → Similar variance captured by both methods")
            print(f"    Data relationships are mostly LINEAR")

    print("="*70 + "\n")

    return fig, (pca, kpca), (X_pca, X_kpca)


import time

def benchmark_kernel_pca(data_dict, variables=['stress', 'slip', 'vr'],
                         kernels=['linear', 'poly', 'rbf', 'sigmoid'],
                         n_samples=None):
    """
    Benchmark different Kernel PCA kernels to compare speed and variance captured.

    Parameters
    ----------
    data_dict : dict
        Dictionary with variable names as keys and numpy arrays as values
    variables : list of str
        List of variable names to include
    kernels : list of str
        List of kernels to test (default: ['linear', 'poly', 'rbf', 'sigmoid'])
    n_samples : int, optional
        Number of samples to use for benchmarking. If None, uses all data.
        Use smaller values (e.g., 1000-10000) for faster benchmarking.

    Returns
    -------
    results_df : pd.DataFrame
        Benchmark results with columns: kernel, time_seconds, variance_pc1, variance_pc2

    Notes
    -----
    Typical speed ranking (fastest to slowest):
    1. 'linear' - Equivalent to standard PCA, very fast
    2. 'poly' (degree=2) - Fast, polynomial kernel
    3. 'rbf' - Moderate speed, most commonly used non-linear kernel
    4. 'sigmoid' - Similar speed to rbf

    Memory usage increases with n_samples^2 for kernel matrix computation.

    Example
    -------
    >>> results = benchmark_kernel_pca(data_dict, variables=['rake', 'ASl', 'Vr'],
    ...                                 n_samples=5000)
    >>> print(results)
    """
    import time

    # Prepare data
    data_list = [data_dict[var].flatten() for var in variables]
    X = np.column_stack(data_list)

    # Subsample if requested
    if n_samples is not None and n_samples < X.shape[0]:
        indices = np.random.choice(X.shape[0], n_samples, replace=False)
        X = X[indices]
        print(f"Using {n_samples} samples for benchmarking (out of {len(data_list[0].flatten())})")

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    n_comp = min(len(variables), X_scaled.shape[0] - 1)

    results = []

    print("="*70)
    print("KERNEL PCA BENCHMARK")
    print("="*70)
    print(f"Data shape: {X_scaled.shape}")
    print(f"Number of components: {n_comp}")
    print(f"\nTesting kernels: {kernels}")
    print("="*70)

    for kernel in kernels:
        print(f"\nTesting kernel: '{kernel}'...")

        # Set gamma for non-linear kernels
        gamma = 1.0 / len(variables) if kernel in ['rbf', 'poly', 'sigmoid'] else None

        try:
            # Time the fitting
            start_time = time.time()

            kpca = KernelPCA(n_components=n_comp, kernel=kernel, gamma=gamma,
                            fit_inverse_transform=False, eigen_solver='auto')
            X_kpca = kpca.fit_transform(X_scaled)

            elapsed_time = time.time() - start_time

            # Try to compute variance
            try:
                # For kernel PCA, compute variance in KERNEL space
                total_var_kernel = np.var(X_kpca, axis=0).sum()
                var_per_component = np.var(X_kpca, axis=0)

                # Variance relative to kernel space (always ≤100%)
                var_ratio_kernel = var_per_component / total_var_kernel

                # Also show variance relative to original space (can be >100% for poly)
                total_var_original = np.var(X_scaled, axis=0).sum()
                var_ratio_original = var_per_component / total_var_original

                var_pc1 = var_ratio_original[0]
                var_pc2 = var_ratio_original[1] if len(var_ratio_original) > 1 else 0.0
                var_total = var_pc1 + var_pc2

                # Add warning if >100%
                if var_total > 1.0 and kernel == 'poly':
                    print(f"      Note: Variance >100% is EXPECTED for polynomial kernel")
                    print(f"      This indicates feature space expansion (not an error)")
            except:
                var_pc1 = np.nan
                var_pc2 = np.nan
                var_total = np.nan

            results.append({
                'kernel': kernel,
                'time_seconds': elapsed_time,
                'variance_pc1': var_pc1,
                'variance_pc2': var_pc2,
                'variance_pc1+pc2': var_total,
                'status': 'success'
            })

            print(f"  ✓ Time: {elapsed_time:.4f}s | PC1+PC2 variance: {var_total:.1%}")

        except Exception as e:
            print(f"  ✗ Error: {e}")
            results.append({
                'kernel': kernel,
                'time_seconds': np.nan,
                'variance_pc1': np.nan,
                'variance_pc2': np.nan,
                'variance_pc1+pc2': np.nan,
                'status': f'error: {str(e)[:50]}'
            })

    # Create results DataFrame
    results_df = pd.DataFrame(results)

    # Print summary
    print("\n" + "="*70)
    print("BENCHMARK RESULTS SUMMARY")
    print("="*70)
    print(results_df.to_string(index=False))

    # Find fastest
    if results_df['time_seconds'].notna().any():
        fastest = results_df.loc[results_df['time_seconds'].idxmin()]
        print(f"\nFastest kernel: '{fastest['kernel']}' ({fastest['time_seconds']:.4f}s)")

        # Speed comparison
        print("\nSpeed comparison (relative to fastest):")
        for _, row in results_df.iterrows():
            if pd.notna(row['time_seconds']):
                speedup = row['time_seconds'] / fastest['time_seconds']
                print(f"  {row['kernel']:10s}: {speedup:.2f}x slower")

    print("="*70 + "\n")

    return results_df


def plot_kernel_comparison(data_dict, variables=['stress', 'slip', 'vr'],
                           kernels=['linear', 'poly', 'rbf', 'sigmoid'],
                           n_samples=None):
    """
    Visualize speed and variance comparison for different kernels.

    Parameters
    ----------
    data_dict : dict
        Dictionary with variable names as keys and numpy arrays as values
    variables : list of str
        List of variable names to include
    kernels : list of str
        List of kernels to test
    n_samples : int, optional
        Number of samples for benchmarking

    Returns
    -------
    fig : matplotlib figure
    results_df : pd.DataFrame

    Example
    -------
    >>> fig, results = plot_kernel_comparison(data_dict,
    ...                                        variables=['rake', 'ASl', 'Vr'],
    ...                                        n_samples=5000)
    """
    # Run benchmark
    results_df = benchmark_kernel_pca(data_dict, variables, kernels, n_samples)

    # Filter successful results
    results_valid = results_df[results_df['status'] == 'success'].copy()

    if len(results_valid) == 0:
        print("No valid results to plot")
        return None, results_df

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # 1. Speed comparison
    colors = ['steelblue' if k == 'linear' else 'tomato' if k == 'rbf'
              else 'orange' if k == 'poly' else 'purple'
              for k in results_valid['kernel']]

    axes[0].barh(results_valid['kernel'], results_valid['time_seconds'], color=colors, alpha=0.7)
    axes[0].set_xlabel('Time (seconds)', fontsize=12, fontweight='bold')
    axes[0].set_title('Computation Speed', fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3, axis='x')

    # Add time labels
    for i, (idx, row) in enumerate(results_valid.iterrows()):
        axes[0].text(row['time_seconds'], i, f"  {row['time_seconds']:.3f}s",
                    va='center', fontsize=10)

    # 2. Variance explained (PC1+PC2)
    axes[1].barh(results_valid['kernel'], results_valid['variance_pc1+pc2']*100,
                color=colors, alpha=0.7)
    axes[1].set_xlabel('Variance Explained (%)', fontsize=12, fontweight='bold')
    axes[1].set_title('Total Variance (PC1+PC2)', fontsize=13, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='x')

    # Add percentage labels
    for i, (idx, row) in enumerate(results_valid.iterrows()):
        if pd.notna(row['variance_pc1+pc2']):
            axes[1].text(row['variance_pc1+pc2']*100, i,
                        f"  {row['variance_pc1+pc2']*100:.1f}%",
                        va='center', fontsize=10)

    # 3. PC1 vs PC2 variance breakdown
    x_pos = np.arange(len(results_valid))
    width = 0.35

    axes[2].bar(x_pos - width/2, results_valid['variance_pc1']*100, width,
               label='PC1', color='royalblue', alpha=0.7)
    axes[2].bar(x_pos + width/2, results_valid['variance_pc2']*100, width,
               label='PC2', color='coral', alpha=0.7)

    axes[2].set_xlabel('Kernel', fontsize=12, fontweight='bold')
    axes[2].set_ylabel('Variance Explained (%)', fontsize=12, fontweight='bold')
    axes[2].set_title('PC1 vs PC2 Variance', fontsize=13, fontweight='bold')
    axes[2].set_xticks(x_pos)
    axes[2].set_xticklabels(results_valid['kernel'])
    axes[2].legend()
    axes[2].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    return fig, results_df


def explain_kernel_pca_variance():
    """
    Explain why Kernel PCA can show >100% variance explained.

    This is a common confusion when using non-linear kernels.

    Returns
    -------
    str
        Detailed explanation

    Notes
    -----
    Why Kernel PCA variance can exceed 100%:

    1. **Different Feature Spaces**:
       - Linear PCA: variance in ORIGINAL space (n features)
       - Kernel PCA: variance in KERNEL space (potentially infinite dimensions!)

    2. **Polynomial Kernel Example**:
       - Original space: 2D with variance = 1.0
       - Poly kernel (degree=2): Maps to ~6D space
       - New space can have total variance > 1.0 (e.g., 2.5)
       - First PC might capture 1.5 units → 1.5/1.0 = 150% !

    3. **The Math**:
       For polynomial kernel with degree d:
       φ(x) creates cross-products: x₁, x₂, x₁², x₂², x₁x₂, ...
       These new dimensions ADD variance!

    4. **What This Means**:
       - Variance >100% is NORMAL for non-linear kernels
       - It means the kernel space has MORE total variance than original space
       - The kernel is "expanding" the feature representation
       - NOT an error - just a different reference frame

    5. **Correct Interpretation**:
       - For Kernel PCA, compare components to EACH OTHER
       - Don't compare to original space variance
       - Focus on: "PC1 captures 60% of KERNEL space variance"
       - NOT: "PC1 captures 60% of original space variance"

    6. **Which Kernels Show This**:
       - 'poly': YES - creates polynomial features (expands variance)
       - 'rbf': NO - typically stays ≤100% (maps to infinite but bounded space)
       - 'sigmoid': SOMETIMES - depends on parameters
       - 'linear': NO - equivalent to regular PCA (always ≤100%)

    Recommendation
    --------------
    For interpretability, use the RELATIVE variance between components,
    not the absolute percentage compared to original space.

    Example
    -------
    >>> # Polynomial kernel might show:
    >>> PC1: 85% of kernel space (but 150% of original!)
    >>> PC2: 15% of kernel space (but 30% of original!)
    >>> Total: 100% of kernel space (but 180% of original!)
    >>>
    >>> # This is CORRECT behavior!
    """
    explanation = """
    ╔════════════════════════════════════════════════════════════════════╗
    ║           WHY KERNEL PCA VARIANCE CAN EXCEED 100%                  ║
    ╚════════════════════════════════════════════════════════════════════╝

    The issue: You're seeing variance >100% with polynomial kernel.

    Why this happens:
    ─────────────────
    1. Kernel PCA operates in a TRANSFORMED feature space
    2. Polynomial kernels EXPAND the dimensionality
    3. This expansion can INCREASE total variance

    Example with degree=2 polynomial:
    ──────────────────────────────────
    Original space (2D):
      Features: [x₁, x₂]
      Variance: 1.0 units

    Kernel space (~6D):
      Features: [x₁, x₂, x₁², x₂², x₁x₂, √2·x₁, ...]
      Variance: 2.5 units (EXPANDED!)

    Result: PC1 might capture 1.8 units of variance
            1.8 / 1.0 = 180% ← Appears >100%!

    Is this wrong? NO!
    ──────────────────
    ✓ It's mathematically correct
    ✓ It means the kernel is "expanding" your features
    ✓ The variance is measured in DIFFERENT spaces

    How to interpret:
    ─────────────────
    ✗ DON'T compare to original space (meaningless)
    ✓ DO compare kernel PCs to each other
    ✓ DO use for dimensionality reduction
    ✓ DO check if kernel PCA captures MORE patterns than linear PCA

    Summary:
    ────────
    Variance >100% is EXPECTED for polynomial kernels.
    It indicates feature space expansion, not an error.
    """
    return explanation


def explain_kpca_outputs():
    """
    Print a simple explanation of what kpca and X_kpca are.

    This is a helper function to understand the outputs of Kernel PCA functions.

    Returns
    -------
    None
        Prints explanation to console

    Example
    -------
    >>> from cov_corr_funcs import explain_kpca_outputs
    >>> explain_kpca_outputs()
    """
    explanation = """
    ╔════════════════════════════════════════════════════════════════════╗
    ║         UNDERSTANDING KERNEL PCA OUTPUTS: kpca vs X_kpca           ║
    ╚════════════════════════════════════════════════════════════════════╝

    When you run:
        fig, axes, kpca, X_kpca = plot_kernel_pca_analysis(data_dict, variables)

    You get TWO outputs:

    1. kpca - The Fitted MODEL (the "transformer")
    ══════════════════════════════════════════════
    What it is:
      • A trained Kernel PCA model object
      • Stores the transformation "recipe"
      • Like a saved function that can be applied to new data

    What it contains:
      • kpca.lambdas_ = eigenvalues (component importance)
      • kpca.alphas_ = eigenvectors (in kernel space)
      • The kernel type and parameters used

    How to use it:
      • Transform NEW data: new_transformed = kpca.transform(new_data_scaled)
      • Apply same transformation to test set or future data
      • Save it with pickle to use later

    Example:
      # Train on 2020 earthquake data
      fig, axes, kpca, X_kpca = plot_kernel_pca_analysis(data_2020, vars)

      # Apply to 2021 data (must scale first!)
      data_2021_scaled = scaler.transform(data_2021)
      data_2021_kpca = kpca.transform(data_2021_scaled)


    2. X_kpca - The TRANSFORMED DATA (the result)
    ═══════════════════════════════════════════════
    What it is:
      • Your original data in the new kernel PC space
      • A numpy array with shape (n_samples, n_components)
      • Each row = one sample (e.g., one fault element)
      • Each column = one kernel principal component

    Structure:
      X_kpca = [
        [KPC1_sample1, KPC2_sample1, KPC3_sample1, ...],  ← Sample 1
        [KPC1_sample2, KPC2_sample2, KPC3_sample2, ...],  ← Sample 2
        [KPC1_sample3, KPC2_sample3, KPC3_sample3, ...],  ← Sample 3
        ...
      ]

    How to access:
      • X_kpca[:, 0]  →  All KPC1 values (first component)
      • X_kpca[:, 1]  →  All KPC2 values (second component)
      • X_kpca[0, :]  →  All KPCs for first sample

    How to use it:
      • Clustering: kmeans.fit(X_kpca[:, :3])  # Use first 3 KPCs
      • Visualization: plt.scatter(X_kpca[:, 0], X_kpca[:, 1])
      • Regression: model.fit(X_kpca, target_variable)
      • Classification: classifier.fit(X_kpca, labels)


    QUICK REFERENCE
    ═══════════════

    ┌──────────┬────────────────────────┬─────────────────────────────┐
    │ Object   │ What it is             │ When to use                 │
    ├──────────┼────────────────────────┼─────────────────────────────┤
    │ kpca     │ The fitted model       │ Transform NEW data          │
    │          │ (transformer)          │ Save for later use          │
    │          │                        │ Apply to test sets          │
    ├──────────┼────────────────────────┼─────────────────────────────┤
    │ X_kpca   │ Transformed data       │ Immediate analysis          │
    │          │ (results)              │ Clustering/visualization    │
    │          │                        │ Downstream modeling         │
    └──────────┴────────────────────────┴─────────────────────────────┘


    ANALOGY
    ═══════

    Think of PCA like translating a book:

    • kpca = The translator (the person/tool)
      - Can translate MORE books later
      - Stores the translation rules

    • X_kpca = The translated book (the result)
      - Ready to read and use
      - The actual content you work with

    ════════════════════════════════════════════════════════════════════
    """
    print(explanation)


def demonstrate_poly_kernel_variance():
    """
    Create a simple demonstration of why polynomial kernel gives >100% variance.

    Returns
    -------
    None
        Prints explanation and shows numerical example

    Example
    -------
    >>> demonstrate_poly_kernel_variance()
    """
    print("\n" + "="*70)
    print("DEMONSTRATION: Why Polynomial Kernel Shows >100% Variance")
    print("="*70 + "\n")

    # Create simple 2D example
    np.random.seed(42)
    X = np.random.randn(100, 2)  # 2D data

    # Standardize
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Original space variance
    var_original = np.var(X_scaled, axis=0).sum()
    print(f"1. Original Space (2D)")
    print(f"   Variables: x₁, x₂")
    print(f"   Total variance: {var_original:.3f}")

    # Apply polynomial features explicitly (degree=2)
    from sklearn.preprocessing import PolynomialFeatures
    poly = PolynomialFeatures(degree=2, include_bias=False)
    X_poly = poly.fit_transform(X_scaled)

    var_poly = np.var(X_poly, axis=0).sum()
    print(f"\n2. Polynomial Space (degree=2)")
    print(f"   Features: {poly.get_feature_names_out(['x1', 'x2'])}")
    print(f"   Dimensionality: {X_poly.shape[1]}D")
    print(f"   Total variance: {var_poly:.3f}")
    print(f"   Variance expansion: {var_poly/var_original:.2f}x")

    # Apply Kernel PCA
    kpca = KernelPCA(n_components=2, kernel='poly', degree=2, gamma=1.0)
    X_kpca = kpca.fit_transform(X_scaled)

    var_kpca = np.var(X_kpca, axis=0)
    var_kpca_total = var_kpca.sum()

    print(f"\n3. Kernel PCA (polynomial, degree=2)")
    print(f"   KPC1 variance: {var_kpca[0]:.3f}")
    print(f"   KPC2 variance: {var_kpca[1]:.3f}")
    print(f"   Total: {var_kpca_total:.3f}")

    print(f"\n4. Variance Ratios")
    print(f"   Compared to ORIGINAL space:")
    print(f"     KPC1: {var_kpca[0]/var_original*100:.1f}%")
    print(f"     KPC2: {var_kpca[1]/var_original*100:.1f}%")
    print(f"     Total: {var_kpca_total/var_original*100:.1f}% ← Can be >100%!")

    print(f"\n   Compared to KERNEL space:")
    print(f"     KPC1: {var_kpca[0]/var_kpca_total*100:.1f}%")
    print(f"     KPC2: {var_kpca[1]/var_kpca_total*100:.1f}%")
    print(f"     Total: {100.0:.1f}% ← Always 100%!")

    print("\n" + "="*70)
    print("CONCLUSION:")
    print("="*70)
    print("✓ Polynomial kernel expands feature space (2D → 5D)")
    print("✓ Expanded space has MORE total variance than original")
    print("✓ Variance >100% relative to original is MATHEMATICALLY CORRECT")
    print("✓ Use RELATIVE variance between components for interpretation")
    print("="*70 + "\n")


import networkx as nx

def plot_correlation_network(corr_matrix, variables=None, threshold=0.3):
    """
    Plot correlation matrix as a network graph.

    Parameters
    ----------
    corr_matrix : pd.DataFrame or np.ndarray
        Correlation matrix (n x n). If DataFrame, variable names are extracted from index.
    variables : list of str, optional
        Variable names. If None and corr_matrix is DataFrame, uses index/columns.
    threshold : float
        Minimum absolute correlation to display as edge (default 0.3)

    Returns
    -------
    fig, ax : matplotlib figure and axis objects

    Example
    -------
    >>> # Direct usage with cov_to_corr output
    >>> R = cov_to_corr(cov_matrix)
    >>> fig, ax = plot_correlation_network(R, threshold=0.3)
    >>>
    >>> # Or convert manually
    >>> fig, ax = plot_correlation_network(R.to_numpy(), variables=R.columns.tolist())
    """
    # Handle DataFrame input
    if isinstance(corr_matrix, pd.DataFrame):
        if variables is None:
            variables = corr_matrix.columns.tolist()
        corr_matrix = corr_matrix.to_numpy()
    elif variables is None:
        # If numpy array without variable names, create generic labels
        n = corr_matrix.shape[0]
        variables = [f'v{i}' for i in range(n)]

    fig, ax = plt.subplots(figsize=(6, 6))

    # Create graph
    G = nx.Graph()
    G.add_nodes_from(variables)
    
    # Add edges for correlations above threshold
    for i in range(len(variables)):
        for j in range(i+1, len(variables)):
            corr = corr_matrix[i, j]
            if abs(corr) > threshold:
                G.add_edge(variables[i], variables[j], 
                          weight=abs(corr), correlation=corr)
    
    # Layout
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=3000, node_color='lightblue',
                          alpha=0.9, ax=ax)
    
    # Draw edges with color based on correlation sign
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    corrs = [G[u][v]['correlation'] for u, v in edges]
    
    colors = ['tomato' if c < 0 else 'blue' for c in corrs]
    
    nx.draw_networkx_edges(G, pos, width=[w*5 for w in weights],
                          edge_color=colors, alpha=0.6, ax=ax)
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold', ax=ax)
    
    # Add edge labels
    edge_labels = {(u, v): f'{G[u][v]["correlation"]:.2f}' 
                  for u, v in edges}
    nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=9, ax=ax)
    
    ax.set_title('Correlation Network (threshold={:.2f})'.format(threshold),
                fontsize=14, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()
    return fig, ax


def add_coordinates_from_xdmf(
    fault_data: pd.DataFrame,
    xdmf_file: str,
    method: str = 'centroid'
) -> pd.DataFrame:
    """
    Add x, y, z coordinates to fault data by extracting from SeisSol XDMF file.

    Uses seissolxdmf to read geometry and connectivity, then computes element
    positions (centroids or first vertex).

    Parameters
    ----------
    fault_data : pd.DataFrame
        Fault data DataFrame (e.g., from stress CSV files).
        Must have same number of rows as XDMF elements.
    xdmf_file : str
        Path to SeisSol XDMF file (e.g., 'model-fault.xdmf')
    method : str
        Method to compute element position (default: 'centroid')
        - 'centroid': Average of all vertices (recommended)
        - 'first': Use first vertex of each element
        - 'mean': Same as centroid (alias)

    Returns
    -------
    fault_data_with_coords : pd.DataFrame
        Input DataFrame with added columns: 'x', 'y', 'z'

    Examples
    --------
    >>> # Basic usage
    >>> ssTable1 = pd.read_csv('stress_jp3z_final.csv')
    >>> ssTable1 = add_coordinates_from_xdmf(
    ...     ssTable1,
    ...     '/path/to/model-fault.xdmf'
    ... )
    >>> print(ssTable1[['x', 'y', 'z']].head())

    >>> # Then calculate seismic moment
    >>> df_moment = calculate_seismic_moment_by_fault(
    ...     ssTable1,
    ...     'mat3d_fault.csv'
    ... )

    Notes
    -----
    - Requires seissolxdmf module to be installed
    - XDMF file must correspond to the same simulation as fault_data
    - Element order must match between CSV and XDMF
    - For triangular elements: 3 vertices, centroid = (v0+v1+v2)/3
    - Coordinates are extracted as-is from XDMF (check units!)
    """
    try:
        import seissolxdmf
    except ImportError:
        raise ImportError(
            "seissolxdmf module is required. Install from:\n"
            "  https://github.com/SeisSol/SeisSol/tree/master/postprocessing/science/seissolxdmf"
        )

    print("="*70)
    print("ADDING COORDINATES FROM XDMF")
    print("="*70)
    print(f"\nLoading XDMF file: {xdmf_file}")

    # Read XDMF
    try:
        sx = seissolxdmf.seissolxdmf(xdmf_file)
    except Exception as e:
        raise RuntimeError(f"Failed to read XDMF file: {e}")

    # Get geometry (vertices) and connectivity (elements)
    print("  Reading geometry (vertices)...")
    vertices = sx.ReadGeometry()  # Shape: (n_vertices, 3)

    print("  Reading connectivity (elements)...")
    connect = sx.ReadConnect()    # Shape: (n_elements, 3 or 4)

    n_vertices = len(vertices)
    n_elements = len(connect)

    print(f"\nXDMF data:")
    print(f"  Vertices: {n_vertices}")
    print(f"  Elements: {n_elements}")
    print(f"  Vertices per element: {connect.shape[1]}")

    # Check dimensions match
    if len(fault_data) != n_elements:
        print(f"\n⚠ WARNING: Row count mismatch!")
        print(f"  fault_data rows: {len(fault_data)}")
        print(f"  XDMF elements: {n_elements}")
        print(f"  Proceeding anyway (may cause issues)...")

    # Compute element positions
    print(f"\nComputing element positions using method='{method}'...")

    if method in ['centroid', 'mean']:
        # Compute centroid of each element
        coords = np.zeros((n_elements, 3))

        for i, elem_vertices in enumerate(connect):
            # Get coordinates of all vertices in this element
            elem_coords = vertices[elem_vertices]
            # Centroid = average of vertex positions
            coords[i] = elem_coords.mean(axis=0)

            if (i + 1) % 50000 == 0:
                print(f"  Processed {i+1}/{n_elements} elements...")

    elif method == 'first':
        # Use first vertex of each element
        coords = vertices[connect[:, 0]]

    else:
        raise ValueError(f"Unknown method: {method}. Use 'centroid' or 'first'")

    # Add coordinates to DataFrame
    fault_data_with_coords = fault_data.copy()
    fault_data_with_coords['x'] = coords[:, 0]
    fault_data_with_coords['y'] = coords[:, 1]
    fault_data_with_coords['z'] = coords[:, 2]

    print("\nCoordinate ranges:")
    print(f"  x: [{coords[:, 0].min():.2f}, {coords[:, 0].max():.2f}]")
    print(f"  y: [{coords[:, 1].min():.2f}, {coords[:, 1].max():.2f}]")
    print(f"  z: [{coords[:, 2].min():.2f}, {coords[:, 2].max():.2f}]")

    print(f"\n✓ Added coordinates to DataFrame")
    print(f"  New columns: ['x', 'y', 'z']")
    print("="*70 + "\n")

    return fault_data_with_coords


def calculate_seismic_moment_by_fault(
    fault_data: pd.DataFrame,
    material_file: str,
    x_col: str = 'x',
    y_col: str = 'y',
    z_col: str = 'z',
    asl_col: str = 'ASl',
    area_col: str = 'Area',
    fault_tag_col: str = 'fault-tag',
    mu_col: str = 'mu',
    use_kdtree: bool = True,
    max_distance: float = None
) -> pd.DataFrame:
    """
    Calculate accumulated seismic moment for each fault segment.

    Seismic moment is calculated as: M0 = mu * ASl * Area
    where mu (shear modulus) is found from the closest point in material properties file.

    Parameters
    ----------
    fault_data : pd.DataFrame
        Fault data with columns for coordinates, slip, area, and fault-tag.
        Required columns: x, y, z, ASl, Area, fault-tag (or custom names)
    material_file : str
        Path to CSV file with material properties (x, y, z, mu, lambda, rho)
    x_col, y_col, z_col : str
        Column names for coordinates in both DataFrames (default: 'x', 'y', 'z')
    asl_col : str
        Column name for accumulated slip (default: 'ASl')
    area_col : str
        Column name for element area (default: 'Area')
    fault_tag_col : str
        Column name for fault segment identifier (default: 'fault-tag')
    mu_col : str
        Column name for shear modulus in material file (default: 'mu')
    use_kdtree : bool
        Use KDTree for fast nearest neighbor search (default: True)
    max_distance : float, optional
        Maximum distance for nearest neighbor search. If closest point is
        farther than this, use a default mu value (default: None, no limit)

    Returns
    -------
    df_moment : pd.DataFrame
        Summary DataFrame with columns:
        - fault-tag: Fault segment identifier
        - count: Number of elements in segment
        - total_M0: Total seismic moment (N⋅m)
        - mean_mu: Mean shear modulus (Pa)
        - total_slip: Total slip (sum of ASl, m)
        - mean_slip: Mean slip (m)
        - total_area: Total fault area (m²)
        - Mw: Moment magnitude (calculated from M0)
        - mean_M0_per_element: Average moment per element (N⋅m)

    Examples
    --------
    >>> # Basic usage
    >>> df_moment = calculate_seismic_moment_by_fault(
    ...     ssTable1,
    ...     'mat3d_fault.csv'
    ... )

    >>> # With custom column names
    >>> df_moment = calculate_seismic_moment_by_fault(
    ...     fault_df,
    ...     'materials.csv',
    ...     x_col='X', y_col='Y', z_col='Z',
    ...     asl_col='slip', area_col='area'
    ... )

    >>> # Display results
    >>> print(df_moment)
    >>> print(f"Total moment: {df_moment['total_M0'].sum():.3e} N⋅m")
    >>> print(f"Equivalent Mw: {df_moment['Mw'].max():.2f}")

    Notes
    -----
    - Seismic moment: M0 = μ × ΔA × d, where:
      * μ (mu): shear modulus (Pa)
      * ΔA (Area): fault area (m²)
      * d (ASl): average slip (m)

    - Moment magnitude: Mw = (2/3) × log10(M0) - 6.07
      where M0 is in N⋅m (Hanks & Kanamori, 1979)

    - For each fault element, the closest mu value from the material
      properties file is used based on 3D Euclidean distance

    - If coordinates have different units (e.g., x,y in km, z in m),
      normalize them before calling this function
    """
    from scipy.spatial import KDTree

    print("="*70)
    print("SEISMIC MOMENT CALCULATION")
    print("="*70)

    # Load material properties
    print(f"\nLoading material properties from: {material_file}")
    try:
        df_mat = pd.read_csv(material_file)
    except FileNotFoundError:
        raise FileNotFoundError(f"Material file not found: {material_file}")

    print(f"  Loaded {len(df_mat)} material property points")
    print(f"  Available columns: {df_mat.columns.tolist()}")

    # Check required columns
    required_fault = [x_col, y_col, z_col, asl_col, area_col, fault_tag_col]
    required_mat = [x_col, y_col, z_col, mu_col]

    missing_fault = [col for col in required_fault if col not in fault_data.columns]
    missing_mat = [col for col in required_mat if col not in df_mat.columns]

    if missing_fault:
        raise ValueError(f"Missing columns in fault_data: {missing_fault}")
    if missing_mat:
        raise ValueError(f"Missing columns in material file: {missing_mat}")

    # Extract coordinates and properties
    fault_coords = fault_data[[x_col, y_col, z_col]].values
    mat_coords = df_mat[[x_col, y_col, z_col]].values
    mat_mu = df_mat[mu_col].values

    print(f"\nFault data:")
    print(f"  Points: {len(fault_coords)}")
    print(f"  x range: [{fault_coords[:, 0].min():.1f}, {fault_coords[:, 0].max():.1f}]")
    print(f"  y range: [{fault_coords[:, 1].min():.1f}, {fault_coords[:, 1].max():.1f}]")
    print(f"  z range: [{fault_coords[:, 2].min():.1f}, {fault_coords[:, 2].max():.1f}]")

    print(f"\nMaterial properties:")
    print(f"  Points: {len(mat_coords)}")
    print(f"  mu range: [{mat_mu.min():.3e}, {mat_mu.max():.3e}] Pa")

    # Find closest mu for each fault point
    print("\nFinding closest material properties for each fault point...")

    if use_kdtree:
        # Build KDTree for fast nearest neighbor search
        tree = KDTree(mat_coords)
        distances, indices = tree.query(fault_coords)
        closest_mu = mat_mu[indices]

        print(f"  Using KDTree (fast method)")
        print(f"  Mean distance to nearest material point: {distances.mean():.1f}")
        print(f"  Max distance to nearest material point: {distances.max():.1f}")

        # Handle max_distance constraint
        if max_distance is not None:
            far_points = distances > max_distance
            if far_points.any():
                default_mu = mat_mu.mean()
                print(f"  Warning: {far_points.sum()} points farther than {max_distance:.1f}")
                print(f"           Using default mu = {default_mu:.3e} Pa for these points")
                closest_mu[far_points] = default_mu
    else:
        # Brute force search (slower but simple)
        print(f"  Using brute force search (may be slow for large datasets)")
        closest_mu = np.zeros(len(fault_coords))

        for i, fault_pt in enumerate(fault_coords):
            distances = np.sqrt(np.sum((mat_coords - fault_pt)**2, axis=1))
            closest_idx = np.argmin(distances)
            closest_mu[i] = mat_mu[closest_idx]

            if (i + 1) % 10000 == 0:
                print(f"    Processed {i+1}/{len(fault_coords)} points...")

    # Add mu to fault data
    fault_data_with_mu = fault_data.copy()
    fault_data_with_mu['mu_closest'] = closest_mu

    # Calculate seismic moment for each element: M0 = mu * ASl * Area
    fault_data_with_mu['M0_element'] = (
        fault_data_with_mu['mu_closest'] *
        fault_data_with_mu[asl_col] *
        fault_data_with_mu[area_col]
    )

    print("\nCalculating seismic moment by fault segment...")

    # Group by fault-tag and calculate statistics
    grouped = fault_data_with_mu.groupby(fault_tag_col)

    df_moment = pd.DataFrame({
        'fault-tag': grouped[fault_tag_col].first(),
        'count': grouped.size(),
        'total_M0': grouped['M0_element'].sum(),
        'mean_mu': grouped['mu_closest'].mean(),
        'std_mu': grouped['mu_closest'].std(),
        'total_slip': grouped[asl_col].sum(),
        'mean_slip': grouped[asl_col].mean(),
        'total_area': grouped[area_col].sum(),
        'mean_area': grouped[area_col].mean(),
    }).reset_index(drop=True)

    # Calculate moment magnitude: Mw = (2/3) * log10(M0) - 6.07
    # M0 should be in N⋅m (dyne⋅cm = 1e-7 N⋅m)
    df_moment['Mw'] = (2.0/3.0) * np.log10(df_moment['total_M0']) - 6.07

    # Additional derived quantities
    df_moment['mean_M0_per_element'] = df_moment['total_M0'] / df_moment['count']

    # Sort by total moment (descending)
    df_moment = df_moment.sort_values('total_M0', ascending=False).reset_index(drop=True)

    # Print summary
    print("\n" + "="*70)
    print("SEISMIC MOMENT SUMMARY BY FAULT SEGMENT")
    print("="*70)
    print(f"\nNumber of fault segments: {len(df_moment)}")
    print(f"Total seismic moment: {df_moment['total_M0'].sum():.3e} N⋅m")
    print(f"Equivalent moment magnitude (total): Mw = {(2.0/3.0) * np.log10(df_moment['total_M0'].sum()) - 6.07:.2f}")
    print(f"\nMoment magnitude range: [{df_moment['Mw'].min():.2f}, {df_moment['Mw'].max():.2f}]")
    print(f"Mean shear modulus: {df_moment['mean_mu'].mean():.3e} Pa")

    # Display top segments
    print("\nTop 5 fault segments by seismic moment:")
    print(df_moment[['fault-tag', 'total_M0', 'Mw', 'mean_slip', 'total_area', 'count']].head())

    print("\n" + "="*70)

    return df_moment


def plot_seismic_moment_summary(
    df_moment: pd.DataFrame,
    plot_type: str = 'bar',
    figsize: tuple = (12, 8),
    save_path: str = None
) -> tuple:
    """
    Plot seismic moment summary by fault segment.

    Parameters
    ----------
    df_moment : pd.DataFrame
        Output from calculate_seismic_moment_by_fault()
    plot_type : str
        Type of plot: 'bar' or 'scatter' (default: 'bar')
    figsize : tuple
        Figure size (default: (12, 8))
    save_path : str, optional
        Path to save figure (default: None)

    Returns
    -------
    fig, axes : tuple
        Matplotlib figure and axes objects
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Plot 1: Total moment by fault-tag
    ax = axes[0, 0]
    if plot_type == 'bar':
        ax.bar(range(len(df_moment)), df_moment['total_M0'],
               color='steelblue', alpha=0.7)
        ax.set_xticks(range(len(df_moment)))
        ax.set_xticklabels(df_moment['fault-tag'], rotation=45, ha='right')
    else:
        ax.scatter(range(len(df_moment)), df_moment['total_M0'],
                  s=100, c='steelblue', alpha=0.7)

    ax.set_ylabel('Total Seismic Moment (N⋅m)')
    ax.set_title('Total M₀ by Fault Segment')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    # Plot 2: Moment magnitude
    ax = axes[0, 1]
    if plot_type == 'bar':
        ax.bar(range(len(df_moment)), df_moment['Mw'],
               color='tomato', alpha=0.7)
        ax.set_xticks(range(len(df_moment)))
        ax.set_xticklabels(df_moment['fault-tag'], rotation=45, ha='right')
    else:
        ax.scatter(range(len(df_moment)), df_moment['Mw'],
                  s=100, c='tomato', alpha=0.7)

    ax.set_ylabel('Moment Magnitude (Mw)')
    ax.set_title('Moment Magnitude by Fault Segment')
    ax.grid(True, alpha=0.3)

    # Plot 3: Moment vs slip
    ax = axes[1, 0]
    scatter = ax.scatter(df_moment['mean_slip'], df_moment['total_M0'],
                        s=df_moment['total_area']/1e6,  # Scale by area
                        c=df_moment['mean_mu'],
                        cmap='viridis', alpha=0.7)

    # Add fault-tag labels
    for idx, row in df_moment.iterrows():
        ax.annotate(str(row['fault-tag']),
                   (row['mean_slip'], row['total_M0']),
                   fontsize=8, alpha=0.7)

    ax.set_xlabel('Mean Slip (m)')
    ax.set_ylabel('Total Seismic Moment (N⋅m)')
    ax.set_title('M₀ vs Slip (size=area, color=μ)')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Mean μ (Pa)')

    # Plot 4: Summary statistics table
    ax = axes[1, 1]
    ax.axis('off')

    # Create summary text
    summary_text = f"""
SUMMARY STATISTICS
{'='*40}

Number of segments: {len(df_moment)}

Total M₀: {df_moment['total_M0'].sum():.3e} N⋅m
Total Mw: {(2.0/3.0) * np.log10(df_moment['total_M0'].sum()) - 6.07:.2f}

Mw range: [{df_moment['Mw'].min():.2f}, {df_moment['Mw'].max():.2f}]

Mean μ: {df_moment['mean_mu'].mean():.3e} Pa
Mean slip: {df_moment['mean_slip'].mean():.2f} m
Total area: {df_moment['total_area'].sum():.3e} m²

Top 3 segments (by M₀):
"""

    for i in range(min(3, len(df_moment))):
        row = df_moment.iloc[i]
        summary_text += f"\n{i+1}. Fault-{row['fault-tag']}: "
        summary_text += f"M₀={row['total_M0']:.2e}, Mw={row['Mw']:.2f}"

    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
           fontsize=10, verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()

    # if save_path:
    #     plt.savefig(save_path, dpi=150, bbox_inches='tight')
    #     print(f"\nPlot saved to: {save_path}")

    return fig, axes