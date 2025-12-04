"""
Download GeoNet Strong Motion Summary Files

This module provides functions to download geocsv summary files from GeoNet's
strong motion data repository.

GeoNet URL structure:
https://data.geonet.org.nz/seismic-products/strong-motion/geocsv/summary-latest/id_prefix%3D0000/
"""

import urllib.request
import urllib.error
import pandas as pd
import io
import os
from typing import Optional, Union


def download_geonet_summary(
    filename: str,
    save_path: Optional[str] = None,
    id_prefix: str = "0000",
    return_dataframe: bool = True
) -> Union[pd.DataFrame, str, None]:
    """
    Download a GeoNet strong motion summary file from the data server.

    Parameters
    ----------
    filename : str
        Name of the summary file to download (e.g., '2024p587797.summary.geocsv')
    save_path : str, optional
        Local path to save the downloaded file. If None, file is not saved to disk.
    id_prefix : str, optional
        ID prefix for the URL path (default: '0000')
    return_dataframe : bool, optional
        If True, return data as pandas DataFrame. If False, return raw text (default: True)

    Returns
    -------
    pd.DataFrame or str or None
        - If return_dataframe=True: pandas DataFrame with the geocsv data
        - If return_dataframe=False: raw text content of the file
        - None if download fails

    Examples
    --------
    >>> # Download and return as DataFrame
    >>> df = download_geonet_summary('2024p587797.summary.geocsv')

    >>> # Download and save to file
    >>> df = download_geonet_summary('2024p587797.summary.geocsv',
    ...                               save_path='./data/2024p587797.summary.geocsv')

    >>> # Download as raw text
    >>> text = download_geonet_summary('2024p587797.summary.geocsv',
    ...                                 return_dataframe=False)

    Notes
    -----
    - GeoCSV files contain header lines starting with '#'
    - Data section starts after the header
    - Typical columns: Station, Network, Location, Epicentral Distance, PGA, PGV, SA, etc.
    - Connection timeout is set to 30 seconds
    """

    # Construct the full URL
    base_url = "https://data.geonet.org.nz/seismic-products/strong-motion/geocsv/summary-latest"
    url = f"{base_url}/id_prefix%3D{id_prefix}/{filename}"

    print(f"Downloading from: {url}")

    try:
        # Download the file with timeout
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read().decode('utf-8')
            print(f"Successfully downloaded {filename} ({len(content)} bytes)")

            # Save to file if path is provided
            if save_path:
                os.makedirs(os.path.dirname(save_path), exist_ok=True)
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write(content)
                print(f"Saved to: {save_path}")

            # Return as DataFrame or raw text
            if return_dataframe:
                # Parse GeoCSV format
                # Skip header lines starting with '#'
                lines = content.split('\n')

                # Find where data starts (after header lines)
                data_start_idx = 0
                for i, line in enumerate(lines):
                    if not line.startswith('#') and line.strip():
                        data_start_idx = i
                        break

                # Read into DataFrame
                data_content = '\n'.join(lines[data_start_idx:])
                df = pd.read_csv(io.StringIO(data_content), sep=',')

                print(f"Loaded DataFrame with {len(df)} rows and {len(df.columns)} columns")
                print(f"Columns: {list(df.columns)}")

                return df
            else:
                return content

    except urllib.error.HTTPError as e:
        print(f"HTTP Error {e.code}: {e.reason}")
        print(f"Failed to download {filename}")
        print(f"URL: {url}")
        return None

    except urllib.error.URLError as e:
        print(f"URL Error: {e.reason}")
        print(f"Check your internet connection or URL validity")
        return None

    except Exception as e:
        print(f"Unexpected error: {type(e).__name__}: {e}")
        return None


def download_multiple_summaries(
    filenames: list[str],
    save_dir: Optional[str] = None,
    id_prefix: str = "0000",
    return_dataframes: bool = True
) -> dict:
    """
    Download multiple GeoNet summary files.

    Parameters
    ----------
    filenames : list of str
        List of summary filenames to download
    save_dir : str, optional
        Directory to save downloaded files. If None, files are not saved.
    id_prefix : str, optional
        ID prefix for the URL path (default: '0000')
    return_dataframes : bool, optional
        If True, return DataFrames. If False, return raw text (default: True)

    Returns
    -------
    dict
        Dictionary mapping filenames to their downloaded content (DataFrame or str)
        Failed downloads will have None values

    Examples
    --------
    >>> files = ['2024p587797.summary.geocsv', '2024p587798.summary.geocsv']
    >>> results = download_multiple_summaries(files, save_dir='./data/')
    """

    results = {}

    for filename in filenames:
        print(f"\n{'='*60}")
        print(f"Processing: {filename}")
        print(f"{'='*60}")

        # Construct save path if save_dir is provided
        save_path = None
        if save_dir:
            save_path = os.path.join(save_dir, filename)

        # Download the file
        result = download_geonet_summary(
            filename=filename,
            save_path=save_path,
            id_prefix=id_prefix,
            return_dataframe=return_dataframes
        )

        results[filename] = result

    # Summary
    print(f"\n{'='*60}")
    print(f"Download Summary:")
    print(f"{'='*60}")
    successful = sum(1 for v in results.values() if v is not None)
    print(f"Total files: {len(filenames)}")
    print(f"Successful: {successful}")
    print(f"Failed: {len(filenames) - successful}")

    return results


def get_event_id_from_filename(filename: str) -> str:
    """
    Extract event ID from GeoNet summary filename.

    Parameters
    ----------
    filename : str
        GeoNet summary filename (e.g., '2024p587797.summary.geocsv')

    Returns
    -------
    str
        Event ID (e.g., '2024p587797')

    Examples
    --------
    >>> get_event_id_from_filename('2024p587797.summary.geocsv')
    '2024p587797'
    """
    return filename.split('.')[0]


# Example usage
if __name__ == "__main__":
    # Example 1: Download single file as DataFrame
    print("Example 1: Download single file")
    df = download_geonet_summary('2024p587797.summary.geocsv')
    if df is not None:
        print("\nFirst few rows:")
        print(df.head())

    # Example 2: Download and save to file
    print("\n" + "="*60)
    print("Example 2: Download and save")
    df = download_geonet_summary(
        '2024p587797.summary.geocsv',
        save_path='./geonet_data/2024p587797.summary.geocsv'
    )

    # Example 3: Download multiple files
    print("\n" + "="*60)
    print("Example 3: Download multiple files")
    files = [
        '2024p587797.summary.geocsv',
        '2024p587798.summary.geocsv'
    ]
    results = download_multiple_summaries(files, save_dir='./geonet_data/')
