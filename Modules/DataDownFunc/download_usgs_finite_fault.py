"""
Download USGS Finite Fault Model data for earthquakes
Includes both finite fault model maps and moment rate function files
"""

import requests
import os
from pathlib import Path
import zipfile

def get_product_metadata(event_id, product_type='finite-fault'):
    """
    Get product metadata to find the correct timestamp and code

    Parameters:
    -----------
    event_id : str
        USGS event ID
    product_type : str
        Product type (e.g., 'finite-fault', 'shakemap')

    Returns:
    --------
    dict : Product metadata including timestamp and source code
    """

    # Try to get product info from event geojson
    url = f"https://earthquake.usgs.gov/earthquakes/feed/v1.0/detail/{event_id}.geojson"

    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            products = data.get('properties', {}).get('products', {})

            if product_type in products:
                # Get the most recent product (first in list)
                product_list = products[product_type]
                if product_list:
                    product = product_list[0]
                    return {
                        'code': product.get('code'),
                        'updateTime': product.get('updateTime'),
                        'source': product.get('source', 'us'),
                        'properties': product.get('properties', {}),
                        'contents': product.get('contents', {})
                    }
    except Exception as e:
        print(f"[INFO] Could not retrieve product metadata: {e}")

    return None


def download_usgs_finite_fault(event_id, output_dir=".", source_code='us', version='1'):
    """
    Download USGS Finite Fault Model data for a given earthquake

    Parameters:
    -----------
    event_id : str
        USGS event ID (e.g., 'us6000jllz' for 2023 Turkey earthquake)
        Can include version suffix like 'us6000jllz_1' or just 'us6000jllz'
    output_dir : str
        Directory to save downloaded files
    source_code : str
        Source code (usually 'us'), will be auto-detected if possible
    version : str or int
        Version number (default: 1)

    Returns:
    --------
    dict : Dictionary with paths to downloaded files
    """

    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Parse event_id to handle cases like 'us6000jllz_1'
    if '_' in event_id:
        base_event_id, version = event_id.rsplit('_', 1)
    else:
        base_event_id = event_id

    # Try to get product metadata to find the correct URL structure
    print(f"[INFO] Checking product metadata for event: {base_event_id}")
    metadata = get_product_metadata(base_event_id, 'finite-fault')

    if metadata:
        print(f"[OK] Found finite-fault product from source: {metadata['source']}")
        source_code = metadata['source']
        update_time = metadata['updateTime']

        # Extract product code which might have version suffix
        product_code = metadata['code']
        print(f"[INFO] Product code: {product_code}")
        print(f"[INFO] Update time: {update_time}")

        # Build base URL using the actual product code and timestamp
        base_url = f"https://earthquake.usgs.gov/product/finite-fault/{product_code}/{source_code}/{update_time}"

        # List available files from metadata
        if 'contents' in metadata and metadata['contents']:
            print(f"[INFO] Available files in product:")
            available_files = {}
            for file_path, file_info in metadata['contents'].items():
                filename = os.path.basename(file_path)
                available_files[filename] = file_info.get('url')
                print(f"  - {filename}")
    else:
        # Fallback to standard URL structure
        print(f"[WARNING] Could not retrieve metadata, using standard URL structure")
        product_code = f"{base_event_id}_{version}" if version != '1' else base_event_id
        base_url = f"https://earthquake.usgs.gov/product/finite-fault/{product_code}/{source_code}"
        available_files = {}

    # Common file extensions and their descriptions
    files_to_download = {
        # Archive files (PRIORITY - contain complete model data)
        'fits.zip': 'Complete Finite Fault Model Archive (fits)',
        'FFM.zip': 'Finite Fault Model Archive',
        'inversion.zip': 'Inversion Data Archive',

        # Moment rate function files
        'moment_rate.mr': 'Moment Rate Function Data',
        'moment_rate.png': 'Moment Rate Function Plot',
        'moment_rate.pdf': 'Moment Rate Function PDF',
        'mr.dat': 'Moment Rate Function Data (alternative)',
        'mrtime.png': 'Moment Rate Function Plot (alternative)',
        'mrtime.pdf': 'Moment Rate Function PDF (alternative)',
        'MRTIME.OUT': 'Moment Rate Time Output',

        # Finite fault model maps and figures
        'basemap.png': 'Finite Fault Model Map (basemap)',
        'basemap.pdf': 'Finite Fault Model Map PDF (basemap)',
        'basemap.ps': 'Finite Fault Model Map PostScript (basemap)',
        'map.png': 'Finite Fault Model Map',
        'map.pdf': 'Finite Fault Model Map PDF',
        'slip.png': 'Slip Distribution Map',
        'slip.pdf': 'Slip Distribution Map PDF',
        'slip1.png': 'Slip Distribution Map (version 1)',
        'slip2.png': 'Slip Distribution Map (version 2)',
        'slip3.png': 'Slip Distribution Map (version 3)',
        'coulomb.png': 'Coulomb Stress Map',
        'coulomb.pdf': 'Coulomb Stress Map PDF',
        'cross_section.png': 'Cross Section',
        'cross_section.pdf': 'Cross Section PDF',

        # Model parameters and data files
        'basic_inversion.param': 'Inversion Parameters',
        'FFM.param': 'Finite Fault Model Parameters',
        'INVERSION.INP': 'Inversion Input File',
        'VELMOD.INP': 'Velocity Model Input',
        'basic_inversion.fsp': 'FSP Format Model',
        'complete_inversion.fsp': 'Complete Inversion FSP',

        # Fault geometry and slip files
        'fault.param': 'Fault Parameters',
        'fault.txt': 'Fault Geometry',
        'slip_dist.txt': 'Slip Distribution Data',
        'FFM.dat': 'Finite Fault Model Data',
        'FFM.fsp': 'Finite Fault Model FSP Format',
        'Solucion.txt': 'Solution File',
        'complete_model.fsp': 'Complete Model FSP',

        # CMT solution
        'CMTSOLUTION': 'CMT Solution File',
        'CMTSOLUTION.txt': 'CMT Solution Text',

        # Waveform comparison files
        'wave.pdf': 'Waveform Comparison PDF',
        'wave.png': 'Waveform Comparison PNG',
        'tele.pdf': 'Teleseismic Waveforms PDF',
        'tele.png': 'Teleseismic Waveforms PNG',
        'surf.pdf': 'Surface Waves PDF',
        'surf.png': 'Surface Waves PNG',
        'body.pdf': 'Body Waves PDF',
        'body.png': 'Body Waves PNG',
        'strong.pdf': 'Strong Motion PDF',
        'strong.png': 'Strong Motion PNG',

        # Additional data files
        'synm.str': 'Synthetic Strong Motion',
        'data.txt': 'Data File',
        'stations.txt': 'Stations List',
        'hypocenter.txt': 'Hypocenter Information',
    }

    downloaded_files = {}
    failed_files = []

    print(f"\nDownloading USGS Finite Fault Model data for event: {base_event_id}")
    print(f"Output directory: {output_path.absolute()}\n")

    for filename, description in files_to_download.items():
        # Use direct URL from metadata if available
        if filename in available_files and available_files[filename]:
            url = available_files[filename]
        else:
            # Construct URL
            url = f"{base_url}/{filename}"

        output_file = output_path / filename

        try:
            print(f"Trying to download: {description} ({filename})...", end=" ")
            response = requests.get(url, timeout=30)

            if response.status_code == 200:
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                downloaded_files[filename] = str(output_file)
                print(f"✓ Downloaded ({len(response.content)/1024:.1f} KB)")
            else:
                print(f"✗ Not available (HTTP {response.status_code})")
                failed_files.append(filename)

        except requests.exceptions.RequestException as e:
            print(f"✗ Error: {e}")
            failed_files.append(filename)

    print(f"\n{'='*70}")
    print(f"Download Summary:")
    print(f"  Successfully downloaded: {len(downloaded_files)} files")
    print(f"  Not available/failed: {len(failed_files)} files")
    print(f"{'='*70}\n")

    if downloaded_files:
        print("Downloaded files:")
        for filename in downloaded_files.keys():
            print(f"  - {filename}")

    # Extract zip files if any were downloaded
    extract_zip_files(downloaded_files, output_path)

    return downloaded_files


def extract_zip_files(downloaded_files, output_dir):
    """
    Extract downloaded zip files and list their contents

    Parameters:
    -----------
    downloaded_files : dict
        Dictionary of downloaded files
    output_dir : Path
        Output directory
    """

    zip_files = [f for f in downloaded_files.keys() if f.endswith('.zip')]

    if not zip_files:
        return

    print(f"\n{'='*70}")
    print(f"Extracting ZIP archives:")
    print(f"{'='*70}\n")

    for zip_filename in zip_files:
        zip_path = Path(downloaded_files[zip_filename])

        try:
            # Create extraction directory (remove .zip extension)
            extract_dir = output_dir / zip_filename.replace('.zip', '')
            extract_dir.mkdir(exist_ok=True)

            print(f"Extracting {zip_filename} to {extract_dir.name}/")

            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)

                # List contents
                file_list = zip_ref.namelist()
                print(f"  Extracted {len(file_list)} files:")
                for filename in file_list[:10]:  # Show first 10 files
                    print(f"    - {filename}")
                if len(file_list) > 10:
                    print(f"    ... and {len(file_list) - 10} more files")

            print(f"  ✓ Successfully extracted to: {extract_dir}\n")

        except Exception as e:
            print(f"  ✗ Error extracting {zip_filename}: {e}\n")


def download_usgs_shakemap(event_id, output_dir="."):
    """
    Download USGS ShakeMap data for a given earthquake

    Parameters:
    -----------
    event_id : str
        USGS event ID
    output_dir : str
        Directory to save downloaded files
    """

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    base_url = f"https://earthquake.usgs.gov/product/shakemap/{event_id}/us"

    shakemap_files = {
        'intensity.jpg': 'ShakeMap Intensity Map',
        'pga.jpg': 'Peak Ground Acceleration Map',
        'pgv.jpg': 'Peak Ground Velocity Map',
        'psa0p3.jpg': 'Spectral Acceleration 0.3s Map',
        'psa1p0.jpg': 'Spectral Acceleration 1.0s Map',
        'psa3p0.jpg': 'Spectral Acceleration 3.0s Map',
        'grid.xml': 'ShakeMap Grid Data (XML)',
        'stationlist.txt': 'Station List',
        'cont_mi.json': 'Intensity Contours (JSON)',
        'cont_pga.json': 'PGA Contours (JSON)',
        'cont_pgv.json': 'PGV Contours (JSON)',
    }

    downloaded_files = {}

    print(f"\nDownloading USGS ShakeMap data for event: {event_id}\n")

    for filename, description in shakemap_files.items():
        url = f"{base_url}/1/{filename}"
        output_file = output_path / filename

        try:
            print(f"Trying to download: {description} ({filename})...", end=" ")
            response = requests.get(url, timeout=30)

            if response.status_code == 200:
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                downloaded_files[filename] = str(output_file)
                print(f"✓ Downloaded ({len(response.content)/1024:.1f} KB)")
            else:
                print(f"✗ Not available (HTTP {response.status_code})")

        except requests.exceptions.RequestException as e:
            print(f"✗ Error: {e}")

    return downloaded_files


def get_event_info(event_id):
    """
    Get basic event information from USGS

    Parameters:
    -----------
    event_id : str
        USGS event ID

    Returns:
    --------
    dict : Event information
    """

    url = f"https://earthquake.usgs.gov/fdsnws/event/1/query?eventid={event_id}&format=geojson"

    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            properties = data['properties']
            geometry = data['geometry']

            info = {
                'title': properties.get('title'),
                'magnitude': properties.get('mag'),
                'time': properties.get('time'),
                'place': properties.get('place'),
                'longitude': geometry['coordinates'][0],
                'latitude': geometry['coordinates'][1],
                'depth': geometry['coordinates'][2],
                'url': properties.get('url')
            }

            print("\nEvent Information:")
            print(f"  Title: {info['title']}")
            print(f"  Magnitude: {info['magnitude']}")
            print(f"  Location: {info['latitude']:.4f}°N, {info['longitude']:.4f}°E")
            print(f"  Depth: {info['depth']:.1f} km")
            print(f"  URL: {info['url']}\n")

            return info
        else:
            print(f"Could not retrieve event information (HTTP {response.status_code})")
            return None

    except Exception as e:
        print(f"Error retrieving event info: {e}")
        return None


# ============================================================================
# Example usage
# ============================================================================

if __name__ == "__main__":
    import sys

    # Example: 2023 Turkey-Syria earthquake (M7.8)
    # URL: https://earthquake.usgs.gov/product/finite-fault/us6000jllz_1/us/1676951251912/moment_rate.mr
    event_id = "us6000jllz"  # or "us6000jllz_1" to specify version
    output_dir = "./usgs_data_turkey"

    # Get event information
    print("="*70)
    event_info = get_event_info(event_id)

    # Download finite fault model data
    finite_fault_files = download_usgs_finite_fault(event_id, output_dir)

    # Download shakemap data (optional)
    # shakemap_files = download_usgs_shakemap(event_id, output_dir)

    print("\n" + "="*70)
    print("Download complete!")
    print("="*70)


# ============================================================================
# Other examples of USGS event IDs:
# ============================================================================
"""
# 2023 Turkey-Syria (M7.8)
event_id = "us6000jllz"
files = download_usgs_finite_fault(event_id, "./usgs_data_turkey")

# 2016 Kaikoura, New Zealand (M7.8)
event_id = "us10006g7d"
files = download_usgs_finite_fault(event_id, "./usgs_data_kaikoura")

# 2015 Nepal (M7.8)
event_id = "us20002926"
files = download_usgs_finite_fault(event_id, "./usgs_data_nepal")

# 2019 Ridgecrest, California (M7.1)
event_id = "ci38457511"
files = download_usgs_finite_fault(event_id, "./usgs_data_ridgecrest")

# 2011 Tohoku, Japan (M9.1)
event_id = "official20110311054624120_30"
files = download_usgs_finite_fault(event_id, "./usgs_data_tohoku")

# Command line usage:
# python download_usgs_finite_fault.py us6000jllz ./output_directory
"""
