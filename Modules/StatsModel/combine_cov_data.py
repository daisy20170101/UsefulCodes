"""
Utility functions for combining covariance analysis data from multiple CSV files.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional, Tuple


def extract_appendix_from_filename(filename: str, prefix: str = 'stress_') -> Optional[str]:
    """
    Extract the model name (appendix) after a specified prefix from a filename.

    Handles new filename format: stress_modelname_final.csv or stress_modelname_t0.csv
    and legacy format: stress_modelname.csv

    Parameters
    ----------
    filename : str
        The filename to extract from (e.g., 'stress_model1_final.csv')
    prefix : str
        The prefix string to search for (default: 'stress_')

    Returns
    -------
    str or None
        The model name after the prefix (e.g., 'model1'), or None if prefix not found

    Example
    -------
    >>> extract_appendix_from_filename('stress_A1B2_final.csv', 'stress_')
    'A1B2'
    >>> extract_appendix_from_filename('stress_A1B2_t0.csv', 'stress_')
    'A1B2'
    >>> extract_appendix_from_filename('stress_A1B2.csv', 'stress_')
    'A1B2'
    """
    # Remove .csv extension if present
    name = filename.replace('.csv', '')

    # Find the prefix
    if prefix in name:
        idx = name.find(prefix) + len(prefix)
        appendix = name[idx:]

        # Remove _final or _t0 suffix if present
        if appendix.endswith('_final'):
            appendix = appendix[:-6]  # Remove '_final'
        elif appendix.endswith('_t0'):
            appendix = appendix[:-3]  # Remove '_t0'

        return appendix
    else:
        return None


def list_appendices_from_folder(
    folder_path: str,
    prefix: str = 'stress_',
    suffix: str = '.csv'
) -> List[str]:
    """
    List all unique appendices after a prefix from CSV filenames in a folder.

    Parameters
    ----------
    folder_path : str
        Path to the folder containing CSV files
    prefix : str
        The prefix string to search for (default: 'stress_')
    suffix : str
        File extension to filter (default: '.csv')

    Returns
    -------
    list of str
        List of unique appendices found

    Example
    -------
    >>> appendices = list_appendices_from_folder('/cov', prefix='stress_')
    >>> print(appendices)
    ['model1', 'model2', 'model3']
    """
    folder = Path(folder_path)

    if not folder.exists():
        raise FileNotFoundError(f"Folder not found: {folder_path}")

    appendices = []

    # Iterate through all files in the folder
    for file_path in folder.glob(f'*{suffix}'):
        filename = file_path.name

        # Extract appendix
        appendix = extract_appendix_from_filename(filename, prefix)

        if appendix is not None:
            appendices.append(appendix)

    # Remove duplicates and sort
    appendices = sorted(list(set(appendices)))

    print(f"Found {len(appendices)} unique appendices: {appendices}")

    return appendices


def combine_csv_files(
    folder_path: str,
    columns_list: List[str],
    prefix: str = 'stress_',
    suffix: str = '.csv',
    add_source_column: bool = True,
    source_column_name: str = 'model',
    stress_drop: bool = True,
    group_by_fault: bool = False,
    groupby_column: str = 'fault-tag',
    agg_func: str = 'mean',
    include_count: bool = True,
    add_seg_tag: bool = True,
    seg_tag_column: str = 'seg-tag',
) -> pd.DataFrame:
    """
    Combine multiple CSV files from a folder into a single DataFrame.

    Parameters
    ----------
    folder_path : str
        Path to the folder containing CSV files
    columns_list : list of str
        List of column names to extract from each CSV file
    prefix : str
        The prefix string to search for in filenames (default: 'stress_')
    suffix : str
        File extension to filter (default: '.csv')
    add_source_column : bool
        Whether to add a column indicating the source file (default: True)
    source_column_name : str
        Name for the source column (default: 'model')
    stress_drop : bool
        Whether to calculate stress drop from t0 and final files (default: True)
        Also calculates T0 and Pf from stress components.
    group_by_fault : bool
        If True, group by fault-tag and return aggregated values instead of
        all individual elements. Also creates composite fault-tags as
        'modelname_faulttag' to keep each model's faults distinct (default: False)
    groupby_column : str
        Column name to group by when group_by_fault=True (default: 'fault-tag')
        This column will be modified to 'modelname_originalvalue' format
    agg_func : str
        Aggregation function when group_by_fault=True: 'mean', 'median', etc.
        (default: 'mean')
    include_count : bool
        If True and group_by_fault=True, include count column (default: True)
    add_seg_tag : bool
        If True, add a 'seg-tag' column containing the ORIGINAL fault-tag values
        before any modifications (e.g., before adding model prefix) (default: True)
    seg_tag_column : str
        Name for the segment tag column (default: 'seg-tag')

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with data from all CSV files.
        If group_by_fault=True, returns grouped/aggregated data.

    Example
    -------
    >>> # Get all individual elements
    >>> columns = ['fault-tag', 'rake', 'ASl', 'Vr']
    >>> df_all = combine_csv_files('/cov', columns, group_by_fault=False)
    >>> print(df_all.shape)
    (100000, 5)  # Many rows

    >>> # Get mean per fault-tag (each model's faults kept separate)
    >>> df_mean = combine_csv_files('/cov', columns, group_by_fault=True)
    >>> print(df_mean.shape)
    (150, 5)  # One row per model_faulttag combination
    >>> print(df_mean['fault-tag'].head())
    # 'modelA_1', 'modelA_2', 'modelA_3', 'modelB_1', 'modelB_2', ...

    Notes
    -----
    When group_by_fault=True:
    - Original fault-tag '1' from model 'ABC' becomes 'ABC_1'
    - Original fault-tag '2' from model 'XYZ' becomes 'XYZ_2'
    - This ensures each model's faults are kept distinct
    - If model 'ABC' and model 'XYZ' both have fault-tag '1', they won't be merged

    When add_seg_tag=True:
    - A 'seg-tag' column is added containing the ORIGINAL fault-tag values
    - This preserves the fault-tag BEFORE it's modified by group_by_fault
    - Useful when you want both: original fault IDs AND model-prefixed fault IDs
    - Example: If original fault-tag is '1', seg-tag will be '1' even if
      fault-tag gets modified to 'ABC_1' by group_by_fault

    Column differences:
    - 'model': Just the model name (e.g., 'ABC')
    - 'fault-tag': Fault identifier, modified to 'modelname_faultid' if group_by_fault=True
      (e.g., original '1' becomes 'ABC_1')
    - 'seg-tag': Original fault-tag value before any modifications (e.g., '1', '2', '3')

    Example with group_by_fault=True:
    - Original data has fault-tag: 1, 1, 2, 2, 3
    - After combine_csv_files():
      - fault-tag: ABC_1, ABC_1, ABC_2, ABC_2, ABC_3  (modified)
      - seg-tag:   1,     1,     2,     2,     3       (original)

    IMPORTANT - Temporal Differences in Variables:
    - PSR (Peak Slip Rate): MAXIMUM slip rate magnitude over ENTIRE time series
    - SRs, SRd (Slip Rate components): Values at FINAL time step only
    - Therefore: PSR >> sqrt(SRs² + SRd²) at final time step
    - Example: PSR=26 m/s (peak during rupture) vs SRs,SRd~0.02 m/s (near zero at end)
    - PSR must be extracted from ParaView during data extraction, not calculated from final SRs/SRd
    """
    folder = Path(folder_path)

    if not folder.exists():
        raise FileNotFoundError(f"Folder not found: {folder_path}")

    # Get list of appendices
    appendices = list_appendices_from_folder(folder_path, prefix, suffix)

    if len(appendices) == 0:
        raise ValueError(f"No files found with prefix '{prefix}' in {folder_path}")

    # List to store individual DataFrames
    df_list = []

    # Process each file
    for appendix in appendices:
        filename = f"{prefix}{appendix}_final{suffix}"
        file_path = folder / filename

        if not file_path.exists():
            print(f"Warning: File not found: {file_path}")
            continue

        try:
            # Read CSV file
            df = pd.read_csv(file_path)

            if stress_drop:
                filename0= f"{prefix}{appendix}_t0{suffix}"
                df0= pd.read_csv(folder/filename0)

            # Check if all required columns exist
            missing_cols = [col for col in columns_list if col not in df.columns]
            if missing_cols:
                print(f"Warning: Missing columns {missing_cols} in {filename}, skipping...")
                continue

            # Select only the required columns
            df_subset = df[columns_list].copy()

            # Add source column if requested
            if add_source_column:
                df_subset[source_column_name] = appendix

            # Add seg-tag column (copy of original fault-tag before modification)
            if add_seg_tag:
                if groupby_column in df_subset.columns:
                    # Store original fault-tag value in seg-tag
                    df_subset[seg_tag_column] = df_subset[groupby_column].copy()
                else:
                    # If fault-tag doesn't exist, use row index
                    df_subset[seg_tag_column] = np.arange(len(df_subset))
                    print(f"Warning: '{groupby_column}' not found in {filename}, using row index for seg-tag")

            # Create composite fault-tag if groupby_column exists and group_by_fault is True
            if group_by_fault and groupby_column in df_subset.columns:
                # Create new column: modelname_faulttag
                df_subset[groupby_column] = appendix + '_' + df_subset[groupby_column].astype(str)

            if stress_drop:
                df_subset['sdrop'] = np.sqrt(df0['Td0']**2+ df0['Ts0']**2) - np.sqrt(df_subset['Td0']**2+ df_subset['Ts0']**2)

                df_subset['T0']= np.sqrt(df0['Td0']**2+ df0['Ts0']**2)
                df_subset['Pf'] = np.abs(df0['Pn0'])

            df_list.append(df_subset)

            print(f"Loaded {len(df_subset)} rows from {filename}")

        except Exception as e:
            print(f"Error reading {filename}: {e}")
            continue

    if len(df_list) == 0:
        raise ValueError("No valid data could be loaded from any files")

    # Combine all DataFrames
    df_all = pd.concat(df_list, axis=0, ignore_index=True)

    print(f"\nCombined {len(df_list)} files into DataFrame with shape: {df_all.shape}")
    print(f"Columns: {df_all.columns.tolist()}")

    # Group by fault-tag if requested
    if group_by_fault:
        if groupby_column not in df_all.columns:
            raise ValueError(f"Groupby column '{groupby_column}' not found in combined data. "
                           f"Make sure to include it in columns_list!")

        print(f"\nGrouping by '{groupby_column}' with agg_func='{agg_func}'...")

        # Get numeric columns to aggregate (exclude groupby_column and source_column)
        exclude_cols = [groupby_column]
        if add_source_column:
            exclude_cols.append(source_column_name)

        numeric_cols = [col for col in df_all.columns
                       if col not in exclude_cols
                       and pd.api.types.is_numeric_dtype(df_all[col])]

        if len(numeric_cols) == 0:
            raise ValueError("No numeric columns to aggregate!")

        # Prepare aggregation
        if isinstance(agg_func, str):
            agg_dict = {col: agg_func for col in numeric_cols}
        else:
            agg_dict = agg_func

        # Group and aggregate
        df_grouped = df_all.groupby(groupby_column)[numeric_cols].agg(agg_dict)

        # Add count if requested
        if include_count:
            count_series = df_all.groupby(groupby_column).size()
            df_grouped.insert(0, 'count', count_series)

        # Reset index
        df_grouped = df_grouped.reset_index()

        print(f"Grouped data shape: {df_grouped.shape}")
        print(f"Number of unique {groupby_column}: {len(df_grouped)}")

        return df_grouped
    else:
        return df_all


def combine_csv_files_advanced(
    folder_path: str,
    columns_list: List[str],
    prefix: str = 'stress_',
    suffix: str = '.csv',
    filter_appendices: Optional[List[str]] = None,
    add_source_column: bool = True,
    source_column_name: str = 'model',
    dropna: bool = False
) -> Tuple[pd.DataFrame, dict]:
    """
    Advanced version: Combine CSV files and return summary statistics.

    Parameters
    ----------
    folder_path : str
        Path to the folder containing CSV files
    columns_list : list of str
        List of column names to extract from each CSV file
    prefix : str
        The prefix string to search for in filenames (default: 'stress_')
    suffix : str
        File extension to filter (default: '.csv')
    filter_appendices : list of str, optional
        Only process files with these appendices. If None, processes all.
    add_source_column : bool
        Whether to add a column indicating the source file (default: True)
    source_column_name : str
        Name for the source column (default: 'model')
    dropna : bool
        Whether to drop rows with NaN values (default: False)

    Returns
    -------
    df_all : pd.DataFrame
        Combined DataFrame with data from all CSV files
    summary : dict
        Dictionary with summary information about each source file

    Example
    -------
    >>> columns = ['stress', 'slip', 'vr']
    >>> df_all, summary = combine_csv_files_advanced('/cov', columns, prefix='stress_')
    >>> print(summary['model1']['n_rows'])
    2500
    """
    folder = Path(folder_path)

    if not folder.exists():
        raise FileNotFoundError(f"Folder not found: {folder_path}")

    # Get list of appendices
    all_appendices = list_appendices_from_folder(folder_path, prefix, suffix)

    # Filter if specified
    if filter_appendices is not None:
        appendices = [a for a in all_appendices if a in filter_appendices]
        print(f"Filtered to {len(appendices)} appendices: {appendices}")
    else:
        appendices = all_appendices

    if len(appendices) == 0:
        raise ValueError("No valid appendices to process")

    # List to store individual DataFrames
    df_list = []
    summary = {}

    # Process each file
    for appendix in appendices:
        filename = f"{prefix}{appendix}{suffix}"
        file_path = folder / filename

        if not file_path.exists():
            print(f"Warning: File not found: {file_path}")
            continue

        try:
            # Read CSV file
            df = pd.read_csv(file_path)

            # Check if all required columns exist
            missing_cols = [col for col in columns_list if col not in df.columns]
            if missing_cols:
                print(f"Warning: Missing columns {missing_cols} in {filename}, skipping...")
                summary[appendix] = {'status': 'missing_columns', 'missing': missing_cols}
                continue

            # Select only the required columns
            df_subset = df[columns_list].copy()

            # Store summary before any filtering
            summary[appendix] = {
                'status': 'success',
                'n_rows': len(df_subset),
                'n_cols': len(columns_list),
                'file_path': str(file_path)
            }

            # Drop NaN if requested
            if dropna:
                n_before = len(df_subset)
                df_subset = df_subset.dropna()
                n_after = len(df_subset)
                summary[appendix]['n_dropped_nan'] = n_before - n_after

            # Add source column if requested
            if add_source_column:
                df_subset[source_column_name] = appendix

            df_list.append(df_subset)

            print(f"Loaded {len(df_subset)} rows from {filename}")

        except Exception as e:
            print(f"Error reading {filename}: {e}")
            summary[appendix] = {'status': 'error', 'error': str(e)}
            continue

    if len(df_list) == 0:
        raise ValueError("No valid data could be loaded from any files")

    # Combine all DataFrames
    df_all = pd.concat(df_list, axis=0, ignore_index=True)

    print(f"\n{'='*60}")
    print(f"Combined {len(df_list)} files into DataFrame with shape: {df_all.shape}")
    print(f"Columns: {df_all.columns.tolist()}")
    print(f"{'='*60}")

    return df_all, summary


def group_by_fault_tag(df: pd.DataFrame,
                       variable_list: List[str],
                       groupby_column: str = 'fault-tag',
                       agg_func: str = 'mean',
                       include_count: bool = True) -> pd.DataFrame:
    """
    Group DataFrame by fault-tag and calculate statistics for specified variables.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing fault-tag and variables to aggregate
    variable_list : list of str
        List of variable names to calculate statistics for
    groupby_column : str
        Column name to group by (default: 'fault-tag')
    agg_func : str or dict
        Aggregation function: 'mean', 'median', 'std', 'min', 'max', 'sum'
        or a dictionary mapping variable names to functions
        (default: 'mean')
    include_count : bool
        If True, include count of samples per group (default: True)

    Returns
    -------
    pd.DataFrame
        Grouped DataFrame with statistics for each fault-tag

    Example
    -------
    >>> # Group by fault-tag and calculate mean
    >>> df_grouped = group_by_fault_tag(
    ...     ssTable1,
    ...     variable_list=['rake', 'ASl', 'Vr', 'PSR', 'T0', 'Pf']
    ... )

    >>> # Group with multiple statistics
    >>> df_grouped = group_by_fault_tag(
    ...     ssTable1,
    ...     variable_list=['rake', 'ASl', 'Vr'],
    ...     agg_func={'rake': 'mean', 'ASl': 'median', 'Vr': 'std'}
    ... )

    >>> # Group by a different column
    >>> df_grouped = group_by_fault_tag(
    ...     ssTable1,
    ...     variable_list=['T0', 'Pf'],
    ...     groupby_column='model',
    ...     agg_func='mean'
    ... )
    """
    # Check if groupby column exists
    if groupby_column not in df.columns:
        raise ValueError(f"Groupby column '{groupby_column}' not found in DataFrame. "
                        f"Available columns: {df.columns.tolist()}")

    # Check which variables exist in the DataFrame
    missing_vars = [var for var in variable_list if var not in df.columns]
    if missing_vars:
        print(f"Warning: Variables {missing_vars} not found in DataFrame, skipping...")

    available_vars = [var for var in variable_list if var in df.columns]

    if len(available_vars) == 0:
        raise ValueError(f"None of the requested variables found in DataFrame. "
                        f"Requested: {variable_list}, Available: {df.columns.tolist()}")

    # Prepare aggregation
    if isinstance(agg_func, str):
        # Single function for all variables
        agg_dict = {var: agg_func for var in available_vars}
    elif isinstance(agg_func, dict):
        # Different function for each variable
        agg_dict = {var: agg_func.get(var, 'mean')
                   for var in available_vars if var in agg_func or var in available_vars}
    else:
        raise ValueError(f"agg_func must be str or dict, got {type(agg_func)}")

    # Group and aggregate
    grouped = df.groupby(groupby_column)[available_vars].agg(agg_dict)

    # Add count if requested
    if include_count:
        count_series = df.groupby(groupby_column).size()
        grouped.insert(0, 'count', count_series)

    # Reset index to make fault-tag a column
    grouped = grouped.reset_index()

    # Print summary
    print("="*70)
    print(f"GROUPING SUMMARY: {agg_func.upper() if isinstance(agg_func, str) else 'CUSTOM'}")
    print("="*70)
    print(f"Grouped by: {groupby_column}")
    print(f"Number of groups: {len(grouped)}")
    print(f"Variables aggregated: {available_vars}")
    print(f"Total rows in input: {len(df)}")
    if include_count:
        print(f"Average rows per group: {df.groupby(groupby_column).size().mean():.1f}")
    print("="*70 + "\n")

    return grouped


def group_by_fault_tag_multi_stats(df: pd.DataFrame,
                                    variable_list: List[str],
                                    groupby_column: str = 'fault-tag',
                                    stats: List[str] = ['mean', 'std', 'count']) -> pd.DataFrame:
    """
    Group by fault-tag and calculate multiple statistics for each variable.

    This creates a multi-level column structure with variables and statistics.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    variable_list : list of str
        Variables to aggregate
    groupby_column : str
        Column to group by (default: 'fault-tag')
    stats : list of str
        List of statistics to calculate: 'mean', 'std', 'min', 'max', 'median', 'count'
        (default: ['mean', 'std', 'count'])

    Returns
    -------
    pd.DataFrame
        DataFrame with multi-level columns: (variable, statistic)

    Example
    -------
    >>> # Calculate mean, std, and count for each variable
    >>> df_stats = group_by_fault_tag_multi_stats(
    ...     ssTable1,
    ...     variable_list=['rake', 'ASl', 'Vr', 'PSR', 'T0', 'Pf'],
    ...     stats=['mean', 'std', 'min', 'max', 'count']
    ... )
    >>>
    >>> # Access specific statistic
    >>> print(df_stats['rake']['mean'])
    >>> print(df_stats['Vr']['std'])
    """
    # Check if groupby column exists
    if groupby_column not in df.columns:
        raise ValueError(f"Groupby column '{groupby_column}' not found in DataFrame")

    # Check which variables exist
    available_vars = [var for var in variable_list if var in df.columns]
    missing_vars = [var for var in variable_list if var not in df.columns]

    if missing_vars:
        print(f"Warning: Variables {missing_vars} not found, skipping...")

    if len(available_vars) == 0:
        raise ValueError("None of the requested variables found in DataFrame")

    # Create aggregation dictionary
    agg_dict = {var: stats for var in available_vars}

    # Group and aggregate
    grouped = df.groupby(groupby_column)[available_vars].agg(agg_dict)

    # Reset index
    grouped = grouped.reset_index()

    # Print summary
    print("="*70)
    print("MULTI-STATISTIC GROUPING SUMMARY")
    print("="*70)
    print(f"Grouped by: {groupby_column}")
    print(f"Number of groups: {len(grouped)}")
    print(f"Variables: {available_vars}")
    print(f"Statistics: {stats}")
    print("="*70 + "\n")

    return grouped


def plot_grouped_stats(df_stats: pd.DataFrame,
                       variables: Optional[List[str]] = None,
                       groupby_column: str = 'fault-tag',
                       stats_to_plot: List[str] = ['mean'],
                       plot_type: str = 'bar',
                       figsize: Optional[tuple] = None,
                       colors: Optional[List[str]] = None,
                       show_errorbar: bool = True,
                       errorbar_stat: str = 'std',
                       pregroup_by_segtag: bool = False,
                       segtag_column: str = 'seg-tag',
                       segtag_agg: str = 'mean') -> tuple:
    """
    Plot grouped statistics from group_by_fault_tag or group_by_fault_tag_multi_stats.

    Parameters
    ----------
    df_stats : pd.DataFrame
        Output from group_by_fault_tag() or group_by_fault_tag_multi_stats()
    variables : list of str, optional
        Variables to plot. If None, plots all numeric columns except groupby_column
    groupby_column : str
        Name of the grouping column (default: 'fault-tag')
    stats_to_plot : list of str
        Statistics to plot if multi-level columns (default: ['mean'])
        Can be ['mean'], ['mean', 'std'], etc.
    plot_type : str
        'bar', 'barh' (horizontal), 'line', 'heatmap', or 'box'
    figsize : tuple, optional
        Figure size (width, height). Auto-calculated if None.
    colors : list of str, optional
        Colors for each variable. Auto-generated if None.
    show_errorbar : bool
        If True and 'std' available, show error bars (default: True)
    errorbar_stat : str
        Statistic to use for error bars: 'std' or 'sem' (default: 'std')
    pregroup_by_segtag : bool
        If True, first groups data by 'seg-tag' column before plotting.
        This aggregates records with the same original fault-tag (default: False)
    segtag_column : str
        Name of the seg-tag column to group by (default: 'seg-tag')
    segtag_agg : str
        Aggregation function when grouping by seg-tag: 'mean', 'median', etc.
        (default: 'mean')

    Returns
    -------
    fig, axes : matplotlib figure and axes

    Example
    -------
    >>> # Simple mean plot
    >>> df_stats = group_by_fault_tag(ssTable1, ['rake', 'ASl', 'Vr'])
    >>> fig, axes = plot_grouped_stats(df_stats, plot_type='bar')

    >>> # Multi-stat with error bars
    >>> df_stats = group_by_fault_tag_multi_stats(ssTable1, ['rake', 'ASl', 'Vr'])
    >>> fig, axes = plot_grouped_stats(df_stats, stats_to_plot=['mean'],
    ...                                 show_errorbar=True)

    >>> # Group by seg-tag first (aggregate by original fault-tag)
    >>> fig, axes = plot_grouped_stats(df_stats, pregroup_by_segtag=True)

    >>> # Heatmap of all statistics
    >>> fig, axes = plot_grouped_stats(df_stats, plot_type='heatmap')
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Pre-group by seg-tag if requested
    if pregroup_by_segtag:
        if segtag_column not in df_stats.columns:
            raise ValueError(f"Column '{segtag_column}' not found in DataFrame. "
                           f"Available columns: {df_stats.columns.tolist()}")

        print(f"Pre-grouping by '{segtag_column}' using '{segtag_agg}' aggregation...")

        # Get numeric columns to aggregate (exclude seg-tag and groupby columns)
        numeric_cols = [col for col in df_stats.columns
                       if col not in [segtag_column, groupby_column, 'count']
                       and pd.api.types.is_numeric_dtype(df_stats[col])]

        if len(numeric_cols) == 0:
            raise ValueError("No numeric columns to aggregate after pre-grouping!")

        # Group by seg-tag and aggregate
        df_stats = df_stats.groupby(segtag_column)[numeric_cols].agg(segtag_agg).reset_index()

        # Use seg-tag as the new groupby column for plotting
        groupby_column = segtag_column

        print(f"After pre-grouping: {len(df_stats)} unique {segtag_column} values")

    # Detect if multi-level columns
    is_multilevel = isinstance(df_stats.columns, pd.MultiIndex)

    # Get groupby values
    if groupby_column not in df_stats.columns:
        # For multi-level, groupby_column might be in level 0
        if is_multilevel and groupby_column in df_stats.columns.get_level_values(0):
            group_values = df_stats[groupby_column].iloc[:, 0].values
        else:
            raise ValueError(f"Groupby column '{groupby_column}' not found")
    else:
        group_values = df_stats[groupby_column].values

    # Determine variables to plot
    if variables is None:
        if is_multilevel:
            # Get all level-0 columns except groupby_column
            variables = [col for col in df_stats.columns.get_level_values(0).unique()
                        if col != groupby_column and col != 'count']
        else:
            # Get all numeric columns except groupby_column and count
            variables = [col for col in df_stats.columns
                        if col not in [groupby_column, 'count']
                        and pd.api.types.is_numeric_dtype(df_stats[col])]

    if len(variables) == 0:
        raise ValueError("No variables to plot")

    # Set default colors
    if colors is None:
        colors = plt.cm.tab10(np.linspace(0, 1, len(variables)))

    # Calculate figure size
    # Note: For bar/line plots, figsize is calculated in _plot_bars_or_lines
    # to account for 3x2 grid layout
    if figsize is None and plot_type == 'heatmap':
        figsize = (min(12, len(variables) * 1.5), min(10, len(group_values) * 0.5))

    # Create figure based on plot type
    if plot_type == 'heatmap':
        return _plot_heatmap(df_stats, variables, group_values, groupby_column,
                            stats_to_plot, is_multilevel, figsize)
    else:
        return _plot_bars_or_lines(df_stats, variables, group_values, groupby_column,
                                   stats_to_plot, is_multilevel, plot_type, figsize,
                                   colors, show_errorbar, errorbar_stat)


def _plot_bars_or_lines(df_stats, variables, group_values, groupby_column,
                        stats_to_plot, is_multilevel, plot_type, figsize,
                        colors, show_errorbar, errorbar_stat):
    """Helper function for bar and line plots."""
    import matplotlib.pyplot as plt

    n_vars = len(variables)

    # Create 3x2 grid layout
    ncols = 3
    nrows = int(np.ceil(n_vars / ncols))

    # Adjust figsize for grid layout
    if figsize is None:
        figsize = (14, 4.5 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    # Hide extra subplots if n_vars < 6
    for idx in range(n_vars, nrows * ncols):
        axes[idx].axis('off')

    for i, var in enumerate(variables):
        ax = axes[i]

        # Extract data
        if is_multilevel:
            # Multi-level columns: (variable, statistic)
            if var not in df_stats.columns.get_level_values(0):
                continue

            # Get mean values
            stat_main = stats_to_plot[0]
            if (var, stat_main) in df_stats.columns:
                y_values = df_stats[(var, stat_main)].values
            else:
                print(f"Warning: ({var}, {stat_main}) not found, skipping...")
                continue

            # Get error values if requested
            yerr = None
            if show_errorbar and errorbar_stat in df_stats[var].columns:
                yerr = df_stats[(var, errorbar_stat)].values

        else:
            # Simple columns
            if var not in df_stats.columns:
                continue
            y_values = df_stats[var].values
            yerr = None

        # Plot
        x_pos = np.arange(len(group_values))

        if plot_type == 'bar':
            ax.bar(x_pos, y_values, color=colors[i], alpha=0.7,
                  yerr=yerr, capsize=5, label=var)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(group_values, rotation=45, ha='right')

        elif plot_type == 'barh':
            ax.barh(x_pos, y_values, color=colors[i], alpha=0.7,
                   xerr=yerr, capsize=5, label=var)
            ax.set_yticks(x_pos)
            ax.set_yticklabels(group_values)

        elif plot_type == 'line':
            if yerr is not None:
                ax.errorbar(x_pos, y_values, yerr=yerr, marker='o',
                           linestyle='-', linewidth=2, markersize=8,
                           color=colors[i], capsize=5, label=var)
            else:
                ax.plot(x_pos, y_values, marker='o', linestyle='-',
                       linewidth=2, markersize=8, color=colors[i], label=var)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(group_values, rotation=45, ha='right')

        # Labels and formatting
        ax.set_title(var, fontsize=13, fontweight='bold')
        ax.set_xlabel(groupby_column, fontsize=11)
        ax.set_ylabel(stats_to_plot[0].capitalize(), fontsize=11)
        ax.grid(True, alpha=0.3, axis='y' if plot_type != 'barh' else 'x')

    plt.tight_layout()
    return fig, axes


def _plot_heatmap(df_stats, variables, group_values, groupby_column,
                 stats_to_plot, is_multilevel, figsize):
    """Helper function for heatmap plot."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Prepare data for heatmap
    if is_multilevel:
        # For multi-level, create separate heatmap for each statistic
        n_stats = len(stats_to_plot)
        fig, axes = plt.subplots(1, n_stats, figsize=figsize, squeeze=False)
        axes = axes.flatten()

        for i, stat in enumerate(stats_to_plot):
            # Extract data for this statistic
            data_list = []
            for var in variables:
                if (var, stat) in df_stats.columns:
                    data_list.append(df_stats[(var, stat)].values)
                else:
                    data_list.append(np.full(len(group_values), np.nan))

            data_matrix = np.column_stack(data_list).T

            # Create heatmap
            sns.heatmap(data_matrix, annot=True, fmt='.2f', cmap='RdYlBu_r',
                       xticklabels=group_values, yticklabels=variables,
                       ax=axes[i], cbar_kws={'label': stat})

            axes[i].set_title(f'{stat.capitalize()} by {groupby_column}',
                            fontsize=13, fontweight='bold')
            axes[i].set_xlabel(groupby_column, fontsize=11)
            axes[i].set_ylabel('Variables', fontsize=11)

    else:
        # Simple heatmap for single statistic
        fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Extract data
        data_list = [df_stats[var].values for var in variables if var in df_stats.columns]
        data_matrix = np.column_stack(data_list).T

        # Create heatmap
        sns.heatmap(data_matrix, annot=True, fmt='.2f', cmap='RdYlBu_r',
                   xticklabels=group_values, yticklabels=variables,
                   ax=ax, cbar_kws={'label': 'Value'})

        ax.set_title(f'Statistics by {groupby_column}',
                    fontsize=13, fontweight='bold')
        ax.set_xlabel(groupby_column, fontsize=11)
        ax.set_ylabel('Variables', fontsize=11)

        axes = [ax]

    plt.tight_layout()
    return fig, axes


def plot_grouped_comparison(df_stats: pd.DataFrame,
                            variables: List[str],
                            groupby_column: str = 'fault-tag',
                            normalize: bool = False,
                            pregroup_by_segtag: bool = False,
                            segtag_column: str = 'seg-tag',
                            segtag_agg: str = 'mean') -> tuple:
    """
    Create a comprehensive comparison plot for grouped statistics.

    Shows all variables in one plot for easy comparison across groups.

    Parameters
    ----------
    df_stats : pd.DataFrame
        Output from group_by_fault_tag() with 'mean' statistic
    variables : list of str
        Variables to compare
    groupby_column : str
        Grouping column name (default: 'fault-tag')
    normalize : bool
        If True, normalize each variable to [0, 1] for comparison (default: False)
    pregroup_by_segtag : bool
        If True, first groups data by 'seg-tag' column before plotting (default: False)
    segtag_column : str
        Name of the seg-tag column to group by (default: 'seg-tag')
    segtag_agg : str
        Aggregation function when grouping by seg-tag (default: 'mean')

    Returns
    -------
    fig, ax : matplotlib figure and axis

    Example
    -------
    >>> df_stats = group_by_fault_tag(ssTable1, ['rake', 'ASl', 'Vr', 'PSR'])
    >>> fig, ax = plot_grouped_comparison(df_stats, ['rake', 'ASl', 'Vr'],
    ...                                    normalize=True)

    >>> # Group by seg-tag first
    >>> fig, ax = plot_grouped_comparison(df_stats, ['rake', 'ASl', 'Vr'],
    ...                                    pregroup_by_segtag=True)
    """
    import matplotlib.pyplot as plt

    # Pre-group by seg-tag if requested
    if pregroup_by_segtag:
        if segtag_column not in df_stats.columns:
            raise ValueError(f"Column '{segtag_column}' not found in DataFrame. "
                           f"Available columns: {df_stats.columns.tolist()}")

        print(f"Pre-grouping by '{segtag_column}' using '{segtag_agg}' aggregation...")

        # Get numeric columns to aggregate (exclude seg-tag and groupby columns)
        numeric_cols = [col for col in df_stats.columns
                       if col not in [segtag_column, groupby_column, 'count']
                       and col in variables
                       and pd.api.types.is_numeric_dtype(df_stats[col])]

        if len(numeric_cols) == 0:
            raise ValueError("No numeric columns to aggregate after pre-grouping!")

        # Group by seg-tag and aggregate
        df_stats = df_stats.groupby(segtag_column)[numeric_cols].agg(segtag_agg).reset_index()

        # Use seg-tag as the new groupby column for plotting
        groupby_column = segtag_column

        print(f"After pre-grouping: {len(df_stats)} unique {segtag_column} values")

    # Detect if multi-level
    is_multilevel = isinstance(df_stats.columns, pd.MultiIndex)

    # Get group values
    if groupby_column not in df_stats.columns:
        if is_multilevel and groupby_column in df_stats.columns.get_level_values(0):
            group_values = df_stats[groupby_column].iloc[:, 0].values
        else:
            raise ValueError(f"Groupby column '{groupby_column}' not found")
    else:
        group_values = df_stats[groupby_column].values

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    x_pos = np.arange(len(group_values))
    width = 0.8 / len(variables)

    # Plot each variable
    for i, var in enumerate(variables):
        # Extract values
        if is_multilevel:
            if (var, 'mean') in df_stats.columns:
                values = df_stats[(var, 'mean')].values
            else:
                print(f"Warning: ({var}, 'mean') not found, skipping...")
                continue
        else:
            if var in df_stats.columns:
                values = df_stats[var].values
            else:
                print(f"Warning: {var} not found, skipping...")
                continue

        # Normalize if requested
        if normalize:
            values = (values - values.min()) / (values.max() - values.min() + 1e-10)

        # Plot
        offset = (i - len(variables)/2) * width + width/2
        ax.bar(x_pos + offset, values, width, label=var, alpha=0.7)

    # Formatting
    ax.set_xlabel(groupby_column, fontsize=12, fontweight='bold')
    ax.set_ylabel('Normalized Value' if normalize else 'Mean Value',
                 fontsize=12, fontweight='bold')
    ax.set_title(f'Variable Comparison by {groupby_column}',
                fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(group_values, rotation=45, ha='right')
    ax.legend(loc='best', frameon=True, shadow=True)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    return fig, ax


def explain_psr_variability():
    """
    Explain why PSR (Peak Slip Rate) often shows high standard deviation.

    Returns
    -------
    str
        Detailed explanation of PSR variability

    Notes
    -----
    PSR is one of the most variable parameters in earthquake rupture dynamics.
    This is both physically meaningful and expected.

    Example
    -------
    >>> from combine_cov_data import explain_psr_variability
    >>> explain_psr_variability()
    """
    explanation = """
    ╔════════════════════════════════════════════════════════════════════╗
    ║       WHY PSR (Peak Slip Rate) HAS HIGH STANDARD DEVIATION        ║
    ╚════════════════════════════════════════════════════════════════════╝

    PSR typically shows MUCH higher variability than other rupture parameters.
    This is NORMAL and EXPECTED for several reasons:

    1. PHYSICAL REASONS (Earthquake Mechanics)
    ═══════════════════════════════════════════

    a) Instantaneous Peak vs. Integrated Quantities:
       • PSR = MAXIMUM slip rate at a SINGLE moment in time
       • ASl (accumulated slip) = INTEGRATED over entire rupture
       • PSR is like "peak speed" vs ASl is like "total distance"
       • Peak values are inherently more variable than averages

    b) Highly Sensitive to Local Conditions:
       • Stress heterogeneity → PSR varies by orders of magnitude
       • Small changes in stress drop → Large changes in PSR
       • Fault geometry variations → PSR spikes at kinks/steps
       • Rupture velocity changes → PSR varies accordingly

    c) Near-field vs Far-field Effects:
       • Near rupture nucleation: Very high PSR
       • Far from nucleation: Lower PSR
       • Stopping phases: PSR can spike again
       • Healing fronts: PSR drops rapidly

    d) Super-shear vs Sub-shear Rupture:
       • Super-shear rupture: PSR can be 10x higher
       • Mode switching during rupture → PSR jumps
       • Mach cones create PSR peaks

    2. NUMERICAL/DATA REASONS
    ═══════════════════════════

    a) Time Derivative Amplifies Noise:
       • PSR = d(slip)/dt at peak
       • Derivatives amplify high-frequency variations
       • Numerical discretization affects peak detection

    b) Temporal Resolution:
       • True peak might occur between time steps
       • Under-sampled → underestimate PSR
       • Interpolation errors → PSR uncertainty

    c) Spatial Averaging:
       • Cell size affects measured PSR
       • Fine mesh → higher PSR peaks captured
       • Coarse mesh → peaks smoothed out

    3. TYPICAL VARIABILITY RANGES
    ══════════════════════════════

    For earthquake rupture parameters:

    ┌─────────────┬───────────────┬──────────────────────────┐
    │ Parameter   │ Typical CV*   │ Variability Level        │
    ├─────────────┼───────────────┼──────────────────────────┤
    │ rake        │ 10-30%        │ LOW  (geometric)         │
    │ ASl         │ 30-60%        │ MODERATE (integrated)    │
    │ Vr          │ 20-40%        │ MODERATE (stable)        │
    │ T0, Pf      │ 30-50%        │ MODERATE (initial state) │
    │ PSR         │ 80-200%+      │ HIGH (instantaneous)     │
    │ sdrop       │ 50-100%       │ MODERATE-HIGH            │
    └─────────────┴───────────────┴──────────────────────────┘

    *CV = Coefficient of Variation (std/mean × 100%)

    4. WHAT THIS MEANS FOR YOUR ANALYSIS
    ═════════════════════════════════════

    ✓ High PSR variability is PHYSICALLY REALISTIC
    ✓ It reflects the complex, transient nature of rupture
    ✓ PSR contains important information about rupture dynamics

    ✗ DON'T remove or ignore high PSR variability
    ✗ DON'T treat PSR outliers as "errors" automatically

    → DO investigate PSR patterns across faults
    → DO check for correlations with stress, Vr, geometry
    → DO consider log-scale for PSR in some analyses

    5. DIAGNOSTIC CHECKS
    ════════════════════

    If you're concerned about PSR variability:

    a) Check distribution shape:
       • Log-normal? → Expected for PSR
       • Bimodal? → May indicate rupture modes
       • Outliers? → Check for numerical artifacts

    b) Compare with other variables:
       • High PSR with high Vr? → Super-shear rupture
       • High PSR with low ASl? → Short duration pulse
       • PSR vs stress_drop correlation? → Expected

    c) Spatial patterns:
       • PSR peaks near nucleation? → Normal
       • PSR peaks at barriers? → Expected
       • Random PSR distribution? → May indicate issues

    6. WHEN TO BE CONCERNED
    ═══════════════════════

    PSR variability is concerning if:

    ⚠ PSR shows NO correlation with any other parameter
    ⚠ PSR values are physically unreasonable (>100 m/s)
    ⚠ PSR has NaN or negative values
    ⚠ PSR std is >5x the mean (extreme even for PSR)

    Otherwise, accept high PSR variability as a feature of
    earthquake rupture dynamics, not a bug!

    ════════════════════════════════════════════════════════════════════
    """
    print(explanation)


# Example usage
if __name__ == "__main__":
    # Example: List appendices
    folder = "/Users/DuoL/Documents/NSHM/cov"
    appendices = list_appendices_from_folder(folder, prefix='stress_')
    print(f"Found appendices: {appendices}")

    # Example: Combine files
    columns = ['stress', 'slip', 'vr', 'rake_deg', 's2ratio', 'xnuc']
    df_all = combine_csv_files(folder, columns, prefix='stress_')
    print(df_all.head())
    print(df_all.info())
