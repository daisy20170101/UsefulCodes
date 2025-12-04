"""
DataFrame downsampling functions for large ground motion datasets.

This module provides various methods to reduce the size of large DataFrames
while preserving important spatial and statistical properties.
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Tuple


def random_sample(df: pd.DataFrame, frac: float = 0.1, random_state: int = 42) -> pd.DataFrame:
    """
    Random downsampling - keep random fraction of data.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    frac : float, optional
        Fraction of data to keep (default: 0.1 for 10%)
    random_state : int, optional
        Random seed for reproducibility (default: 42)

    Returns
    -------
    pd.DataFrame
        Randomly sampled DataFrame

    Examples
    --------
    >>> # Keep 10% of data
    >>> df_sampled = random_sample(df, frac=0.1)
    >>>
    >>> # Keep 5% of data
    >>> df_sampled = random_sample(df, frac=0.05)
    """
    return df.sample(frac=frac, random_state=random_state).reset_index(drop=True)


def systematic_sample(df: pd.DataFrame, every_nth: int = 10) -> pd.DataFrame:
    """
    Systematic downsampling - keep every Nth row.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    every_nth : int, optional
        Keep every Nth row (default: 10)

    Returns
    -------
    pd.DataFrame
        Systematically sampled DataFrame

    Examples
    --------
    >>> # Keep every 10th row
    >>> df_sampled = systematic_sample(df, every_nth=10)
    >>>
    >>> # Keep every 100th row
    >>> df_sampled = systematic_sample(df, every_nth=100)
    """
    return df.iloc[::every_nth].reset_index(drop=True)


def spatial_grid_downsample(
    df: pd.DataFrame,
    grid_spacing_km: float = 1.0,
    x_col: str = 'x (km)',
    y_col: str = 'y (km)',
    agg_method: str = 'max',
    agg_column: str = 'PGA (%g)'
) -> pd.DataFrame:
    """
    Downsample by spatial grid - keep one point per grid cell.

    Divides the spatial domain into grid cells and keeps one representative
    point per cell. The representative can be chosen by max/min/mean of
    a specified column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with spatial coordinates
    grid_spacing_km : float, optional
        Grid spacing in kilometers (default: 1.0)
    x_col : str, optional
        Name of x-coordinate column (default: 'x (km)')
    y_col : str, optional
        Name of y-coordinate column (default: 'y (km)')
    agg_method : str, optional
        Aggregation method: 'max', 'min', 'mean', 'first' (default: 'max')
        Determines which point to keep in each grid cell
    agg_column : str, optional
        Column to use for aggregation (default: 'PGA (%g)')

    Returns
    -------
    pd.DataFrame
        Spatially downsampled DataFrame

    Examples
    --------
    >>> # 2 km grid, keep max PGA point in each cell
    >>> df_sampled = spatial_grid_downsample(df, grid_spacing_km=2.0)
    >>>
    >>> # 5 km grid, keep max PGV point in each cell
    >>> df_sampled = spatial_grid_downsample(
    ...     df,
    ...     grid_spacing_km=5.0,
    ...     agg_method='max',
    ...     agg_column='PGV (cm/s)'
    ... )
    """
    df_copy = df.copy()

    # Create grid indices
    df_copy['grid_x'] = (df_copy[x_col] / grid_spacing_km).astype(int)
    df_copy['grid_y'] = (df_copy[y_col] / grid_spacing_km).astype(int)

    # Group by grid cell and select representative point
    if agg_method == 'max':
        idx = df_copy.groupby(['grid_x', 'grid_y'])[agg_column].idxmax()
    elif agg_method == 'min':
        idx = df_copy.groupby(['grid_x', 'grid_y'])[agg_column].idxmin()
    elif agg_method == 'mean':
        # For mean, take the point closest to cell center
        df_copy['dist_to_center'] = np.sqrt(
            (df_copy[x_col] - (df_copy['grid_x'] + 0.5) * grid_spacing_km)**2 +
            (df_copy[y_col] - (df_copy['grid_y'] + 0.5) * grid_spacing_km)**2
        )
        idx = df_copy.groupby(['grid_x', 'grid_y'])['dist_to_center'].idxmin()
    elif agg_method == 'first':
        idx = df_copy.groupby(['grid_x', 'grid_y']).apply(lambda x: x.index[0])
    else:
        raise ValueError(f"Unknown agg_method: {agg_method}. Use 'max', 'min', 'mean', or 'first'")

    df_sampled = df_copy.loc[idx]

    # Remove temporary columns
    cols_to_drop = ['grid_x', 'grid_y']
    if 'dist_to_center' in df_sampled.columns:
        cols_to_drop.append('dist_to_center')
    df_sampled = df_sampled.drop(columns=cols_to_drop)

    return df_sampled.reset_index(drop=True)


def distance_filter(
    df: pd.DataFrame,
    max_distance: float,
    distance_col: str = 'r_rup (km)'
) -> pd.DataFrame:
    """
    Filter data by maximum distance.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    max_distance : float
        Maximum distance to keep (km)
    distance_col : str, optional
        Name of distance column (default: 'r_rup (km)')

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame

    Examples
    --------
    >>> # Keep only sites within 100 km
    >>> df_sampled = distance_filter(df, max_distance=100)
    >>>
    >>> # Keep sites within 50 km using r_jb
    >>> df_sampled = distance_filter(df, max_distance=50, distance_col='r_jb (km)')
    """
    return df[df[distance_col] <= max_distance].reset_index(drop=True)


def distance_stratified_sample(
    df: pd.DataFrame,
    near_field_km: float = 20,
    near_frac: float = 0.5,
    far_frac: float = 0.1,
    distance_col: str = 'r_rup (km)',
    random_state: int = 42
) -> pd.DataFrame:
    """
    Sample with different fractions for near-field vs far-field.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    near_field_km : float, optional
        Distance threshold for near-field (default: 20 km)
    near_frac : float, optional
        Fraction to keep in near-field (default: 0.5 for 50%)
    far_frac : float, optional
        Fraction to keep in far-field (default: 0.1 for 10%)
    distance_col : str, optional
        Name of distance column (default: 'r_rup (km)')
    random_state : int, optional
        Random seed (default: 42)

    Returns
    -------
    pd.DataFrame
        Stratified sample

    Examples
    --------
    >>> # Keep 50% near-field (<20km), 10% far-field
    >>> df_sampled = distance_stratified_sample(df)
    >>>
    >>> # Custom thresholds
    >>> df_sampled = distance_stratified_sample(
    ...     df,
    ...     near_field_km=30,
    ...     near_frac=0.7,
    ...     far_frac=0.05
    ... )
    """
    near_field = df[df[distance_col] < near_field_km].sample(
        frac=near_frac, random_state=random_state
    )
    far_field = df[df[distance_col] >= near_field_km].sample(
        frac=far_frac, random_state=random_state
    )

    return pd.concat([near_field, far_field]).reset_index(drop=True)


def threshold_filter(
    df: pd.DataFrame,
    threshold: float,
    column: str = 'PGA (%g)',
    operator: str = '>='
) -> pd.DataFrame:
    """
    Filter data by threshold on a specified column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    threshold : float
        Threshold value
    column : str, optional
        Column to filter on (default: 'PGA (%g)')
    operator : str, optional
        Comparison operator: '>=', '>', '<=', '<' (default: '>=')

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame

    Examples
    --------
    >>> # Keep only high PGA sites (>= 1%g)
    >>> df_sampled = threshold_filter(df, threshold=1.0, column='PGA (%g)')
    >>>
    >>> # Keep low PGV sites (< 10 cm/s)
    >>> df_sampled = threshold_filter(
    ...     df,
    ...     threshold=10.0,
    ...     column='PGV (cm/s)',
    ...     operator='<'
    ... )
    """
    if operator == '>=':
        return df[df[column] >= threshold].reset_index(drop=True)
    elif operator == '>':
        return df[df[column] > threshold].reset_index(drop=True)
    elif operator == '<=':
        return df[df[column] <= threshold].reset_index(drop=True)
    elif operator == '<':
        return df[df[column] < threshold].reset_index(drop=True)
    else:
        raise ValueError(f"Unknown operator: {operator}. Use '>=', '>', '<=', or '<'")


def percentile_filter(
    df: pd.DataFrame,
    top_percent: float = 0.2,
    column: str = 'PGA (%g)'
) -> pd.DataFrame:
    """
    Keep only top N% of data by specified column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    top_percent : float, optional
        Fraction to keep (default: 0.2 for top 20%)
    column : str, optional
        Column to rank by (default: 'PGA (%g)')

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with top percentile

    Examples
    --------
    >>> # Keep top 20% by PGA
    >>> df_sampled = percentile_filter(df, top_percent=0.2)
    >>>
    >>> # Keep top 10% by PGV
    >>> df_sampled = percentile_filter(
    ...     df,
    ...     top_percent=0.1,
    ...     column='PGV (cm/s)'
    ... )
    """
    threshold = df[column].quantile(1 - top_percent)
    return df[df[column] >= threshold].reset_index(drop=True)


def distance_binned_sample(
    df: pd.DataFrame,
    n_samples_per_bin: int = 100,
    distance_bins: Optional[List[float]] = None,
    distance_col: str = 'r_rup (km)',
    random_state: int = 42
) -> pd.DataFrame:
    """
    Sample uniformly across distance bins.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    n_samples_per_bin : int, optional
        Number of samples per distance bin (default: 100)
    distance_bins : list of float, optional
        Distance bin edges in km (default: [0, 10, 20, 50, 100, 200, inf])
    distance_col : str, optional
        Name of distance column (default: 'r_rup (km)')
    random_state : int, optional
        Random seed (default: 42)

    Returns
    -------
    pd.DataFrame
        Stratified sample with equal representation across distance bins

    Examples
    --------
    >>> # Default bins with 100 samples each
    >>> df_sampled = distance_binned_sample(df, n_samples_per_bin=100)
    >>>
    >>> # Custom distance bins
    >>> df_sampled = distance_binned_sample(
    ...     df,
    ...     n_samples_per_bin=200,
    ...     distance_bins=[0, 5, 15, 30, 60, 100, np.inf]
    ... )
    """
    if distance_bins is None:
        distance_bins = [0, 10, 20, 50, 100, 200, np.inf]

    df_copy = df.copy()
    df_copy['distance_bin'] = pd.cut(df_copy[distance_col], bins=distance_bins)

    # Sample from each bin
    sampled_dfs = []
    for bin_name in df_copy['distance_bin'].cat.categories:
        bin_data = df_copy[df_copy['distance_bin'] == bin_name]
        if len(bin_data) > 0:
            n = min(n_samples_per_bin, len(bin_data))
            sampled_dfs.append(bin_data.sample(n=n, random_state=random_state))

    df_sampled = pd.concat(sampled_dfs)
    df_sampled = df_sampled.drop(columns=['distance_bin'])

    return df_sampled.reset_index(drop=True)


def smart_downsample(
    df: pd.DataFrame,
    near_field_km: float = 30,
    grid_near: float = 0.5,
    grid_far: float = 5.0,
    x_col: str = 'x (km)',
    y_col: str = 'y (km)',
    distance_col: str = 'r_rup (km)',
    agg_column: str = 'PGA (%g)'
) -> pd.DataFrame:
    """
    Smart downsampling with distance-dependent grid spacing.

    Uses fine spatial grid for near-field (important for engineering)
    and coarse grid for far-field (smoother variations, less critical).

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    near_field_km : float, optional
        Distance threshold for near-field (default: 30 km)
    grid_near : float, optional
        Grid spacing for near-field in km (default: 0.5 km)
    grid_far : float, optional
        Grid spacing for far-field in km (default: 5.0 km)
    x_col : str, optional
        Name of x-coordinate column (default: 'x (km)')
    y_col : str, optional
        Name of y-coordinate column (default: 'y (km)')
    distance_col : str, optional
        Name of distance column (default: 'r_rup (km)')
    agg_column : str, optional
        Column to use for selecting representative point (default: 'PGA (%g)')

    Returns
    -------
    pd.DataFrame
        Intelligently downsampled DataFrame

    Examples
    --------
    >>> # Default: 0.5km near-field, 5km far-field
    >>> df_sampled = smart_downsample(df)
    >>>
    >>> # Custom settings for large region
    >>> df_sampled = smart_downsample(
    ...     df,
    ...     near_field_km=50,
    ...     grid_near=1.0,
    ...     grid_far=10.0
    ... )

    Notes
    -----
    This method is recommended for most seismic analyses as it:
    - Preserves detail in near-field (most important for engineering)
    - Reduces far-field points (less critical, smoother variations)
    - Maintains spatial coverage across the domain
    - Typically achieves 90-95% file size reduction
    """
    # Split near and far field
    df_near = df[df[distance_col] < near_field_km].copy()
    df_far = df[df[distance_col] >= near_field_km].copy()

    # Downsample each with different grid spacing
    if len(df_near) > 0:
        df_near_sampled = spatial_grid_downsample(
            df_near,
            grid_spacing_km=grid_near,
            x_col=x_col,
            y_col=y_col,
            agg_column=agg_column
        )
    else:
        df_near_sampled = df_near

    if len(df_far) > 0:
        df_far_sampled = spatial_grid_downsample(
            df_far,
            grid_spacing_km=grid_far,
            x_col=x_col,
            y_col=y_col,
            agg_column=agg_column
        )
    else:
        df_far_sampled = df_far

    # Combine
    df_sampled = pd.concat([df_near_sampled, df_far_sampled])

    return df_sampled.reset_index(drop=True)


def combined_downsample(
    df: pd.DataFrame,
    max_distance: Optional[float] = None,
    min_threshold: Optional[float] = None,
    threshold_column: str = 'PGA (%g)',
    grid_spacing: Optional[float] = None,
    random_frac: Optional[float] = None,
    random_state: int = 42
) -> pd.DataFrame:
    """
    Combined downsampling with multiple filters.

    Apply multiple downsampling methods in sequence:
    1. Distance filter (if specified)
    2. Threshold filter (if specified)
    3. Spatial grid (if specified)
    4. Random sampling (if specified)

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    max_distance : float, optional
        Maximum distance to keep (km)
    min_threshold : float, optional
        Minimum value threshold for threshold_column
    threshold_column : str, optional
        Column for threshold filtering (default: 'PGA (%g)')
    grid_spacing : float, optional
        Spatial grid spacing (km)
    random_frac : float, optional
        Final random sampling fraction
    random_state : int, optional
        Random seed (default: 42)

    Returns
    -------
    pd.DataFrame
        Downsampled DataFrame

    Examples
    --------
    >>> # Distance + threshold + grid
    >>> df_sampled = combined_downsample(
    ...     df,
    ...     max_distance=100,
    ...     min_threshold=0.5,
    ...     grid_spacing=2.0
    ... )
    >>>
    >>> # All filters
    >>> df_sampled = combined_downsample(
    ...     df,
    ...     max_distance=150,
    ...     min_threshold=1.0,
    ...     threshold_column='PGA (%g)',
    ...     grid_spacing=3.0,
    ...     random_frac=0.5
    ... )
    """
    df_result = df.copy()

    # Apply filters in sequence
    if max_distance is not None:
        df_result = distance_filter(df_result, max_distance=max_distance)

    if min_threshold is not None:
        df_result = threshold_filter(
            df_result,
            threshold=min_threshold,
            column=threshold_column,
            operator='>='
        )

    if grid_spacing is not None:
        df_result = spatial_grid_downsample(df_result, grid_spacing_km=grid_spacing)

    if random_frac is not None:
        df_result = random_sample(df_result, frac=random_frac, random_state=random_state)

    return df_result


def print_downsample_summary(df_original: pd.DataFrame, df_sampled: pd.DataFrame):
    """
    Print summary statistics comparing original and downsampled DataFrames.

    Parameters
    ----------
    df_original : pd.DataFrame
        Original DataFrame
    df_sampled : pd.DataFrame
        Downsampled DataFrame

    Examples
    --------
    >>> df_sampled = smart_downsample(df)
    >>> print_downsample_summary(df, df_sampled)
    """
    print("="*60)
    print("Downsampling Summary")
    print("="*60)
    print(f"Original size:    {len(df_original):,} points")
    print(f"Downsampled size: {len(df_sampled):,} points")
    print(f"Reduction:        {(1 - len(df_sampled)/len(df_original))*100:.1f}%")
    print(f"Remaining:        {len(df_sampled)/len(df_original)*100:.1f}%")

    # Compare statistics for key columns
    common_cols = ['PGA (%g)', 'PGV (cm/s)', 'r_rup (km)']
    available_cols = [col for col in common_cols if col in df_original.columns]

    if available_cols:
        print("\nStatistical Comparison:")
        print("-"*60)
        for col in available_cols:
            orig_mean = df_original[col].mean()
            samp_mean = df_sampled[col].mean()
            orig_max = df_original[col].max()
            samp_max = df_sampled[col].max()

            print(f"\n{col}:")
            print(f"  Mean: {orig_mean:.3f} → {samp_mean:.3f} "
                  f"({(samp_mean/orig_mean-1)*100:+.1f}%)")
            print(f"  Max:  {orig_max:.3f} → {samp_max:.3f} "
                  f"({(samp_max/orig_max-1)*100:+.1f}%)")

    print("="*60)


# Example usage
if __name__ == "__main__":
    print("Ground Motion DataFrame Downsampling Functions")
    print("="*60)
    print("\nAvailable functions:")
    print("  1. random_sample() - Random sampling")
    print("  2. systematic_sample() - Every Nth row")
    print("  3. spatial_grid_downsample() - Spatial grid")
    print("  4. distance_filter() - Distance threshold")
    print("  5. distance_stratified_sample() - Near/far field stratification")
    print("  6. threshold_filter() - Value threshold")
    print("  7. percentile_filter() - Top N% by value")
    print("  8. distance_binned_sample() - Uniform across distance bins")
    print("  9. smart_downsample() - Distance-dependent grid (RECOMMENDED)")
    print(" 10. combined_downsample() - Multiple filters")
    print(" 11. print_downsample_summary() - Compare original vs downsampled")
    print("\nImport with:")
    print("  from GMfunc.func_desample import smart_downsample")
