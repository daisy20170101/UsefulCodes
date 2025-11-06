#!/usr/bin/env python3
"""
Download broadband seismic data from IRIS for New Zealand earthquake
Event: 2025-11-06T08:09:54Z

This script uses ObsPy to download waveform data from IRIS FDSN web services
for the New Zealand earthquake event.

Requirements:
    pip install obspy

Usage:
    python download_nz_earthquake_data.py
"""

import os
from datetime import datetime, timedelta
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.inventory import Inventory


def download_earthquake_data(
    event_time_str="2025-11-06T08:09:54Z",
    output_dir="./nz_earthquake_data",
    network="*",
    station="*",
    location="*",
    channel="BH*",  # Broadband, high-gain seismometers
    time_before=300,  # seconds before event
    time_after=3600,  # seconds after event (1 hour)
    max_radius=20.0,  # degrees from epicenter
    min_radius=0.0,
    latitude=-37.5,  # Approximate New Zealand location
    longitude=179.0,
):
    """
    Download broadband seismic data from IRIS for specified earthquake

    Parameters:
    -----------
    event_time_str : str
        Event origin time in ISO format
    output_dir : str
        Directory to save downloaded data
    network : str
        Network code(s), wildcards allowed (e.g., 'IU', '*')
    station : str
        Station code(s), wildcards allowed
    location : str
        Location code(s), wildcards allowed
    channel : str
        Channel code(s) - BH* for broadband high-gain
        - BHZ: Vertical component
        - BHN: North component
        - BHE: East component
    time_before : float
        Seconds before event to download
    time_after : float
        Seconds after event to download
    max_radius : float
        Maximum distance from epicenter in degrees
    min_radius : float
        Minimum distance from epicenter in degrees
    latitude : float
        Event latitude (will be updated from catalog if available)
    longitude : float
        Event longitude (will be updated from catalog if available)

    Returns:
    --------
    stream : obspy.Stream
        Downloaded waveform data
    inventory : obspy.Inventory
        Station metadata
    event : obspy.core.event.Event
        Event information from catalog
    """

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Initialize IRIS client
    print("Connecting to IRIS FDSN web services...")
    client = Client("IRIS")

    # Convert event time to UTCDateTime
    event_time = UTCDateTime(event_time_str)
    print(f"Event time: {event_time}")

    # Define time window
    start_time = event_time - time_before
    end_time = event_time + time_after

    print(f"Time window: {start_time} to {end_time}")
    print(f"Duration: {time_before + time_after} seconds ({(time_before + time_after)/60:.1f} minutes)")

    # Try to get event information from catalog
    event = None
    try:
        print("\nSearching for event in USGS catalog...")
        catalog = client.get_events(
            starttime=event_time - 60,  # Search 1 minute before
            endtime=event_time + 60,    # Search 1 minute after
            minlatitude=-50,
            maxlatitude=-30,
            minlongitude=165,
            maxlongitude=-175,
        )

        if len(catalog) > 0:
            event = catalog[0]
            origin = event.preferred_origin() or event.origins[0]
            magnitude = event.preferred_magnitude() or event.magnitudes[0]

            latitude = origin.latitude
            longitude = origin.longitude
            depth = origin.depth / 1000.0  # Convert to km

            print(f"✓ Found event in catalog:")
            print(f"  Origin time: {origin.time}")
            print(f"  Latitude: {latitude:.4f}°")
            print(f"  Longitude: {longitude:.4f}°")
            print(f"  Depth: {depth:.2f} km")
            print(f"  Magnitude: {magnitude.mag} {magnitude.magnitude_type}")

            # Save event information
            event_file = os.path.join(output_dir, "event_info.xml")
            catalog.write(event_file, format="QUAKEML")
            print(f"✓ Saved event info to: {event_file}")
        else:
            print("Warning: Event not found in catalog, using provided coordinates")

    except Exception as e:
        print(f"Warning: Could not retrieve event from catalog: {e}")
        print("Using provided coordinates")

    # Download waveform data
    print(f"\nDownloading waveform data...")
    print(f"  Network: {network}")
    print(f"  Station: {station}")
    print(f"  Location: {location}")
    print(f"  Channel: {channel}")
    print(f"  Radius: {min_radius}° to {max_radius}° from epicenter")

    try:
        stream = client.get_waveforms(
            network=network,
            station=station,
            location=location,
            channel=channel,
            starttime=start_time,
            endtime=end_time,
            latitude=latitude,
            longitude=longitude,
            minradius=min_radius,
            maxradius=max_radius,
        )

        print(f"✓ Downloaded {len(stream)} traces")
        print(f"\nStream summary:")
        print(stream)

        # Save waveforms
        waveform_file = os.path.join(output_dir, "waveforms.mseed")
        stream.write(waveform_file, format="MSEED")
        print(f"✓ Saved waveforms to: {waveform_file}")

    except Exception as e:
        print(f"Error downloading waveforms: {e}")
        stream = None

    # Download station metadata
    print(f"\nDownloading station metadata...")

    try:
        inventory = client.get_stations(
            network=network,
            station=station,
            location=location,
            channel=channel,
            starttime=start_time,
            endtime=end_time,
            latitude=latitude,
            longitude=longitude,
            minradius=min_radius,
            maxradius=max_radius,
            level="response",  # Include instrument response
        )

        print(f"✓ Downloaded metadata for {len(inventory.get_contents()['stations'])} stations")
        print(f"\nInventory summary:")
        print(inventory)

        # Save inventory
        inventory_file = os.path.join(output_dir, "station_inventory.xml")
        inventory.write(inventory_file, format="STATIONXML")
        print(f"✓ Saved inventory to: {inventory_file}")

    except Exception as e:
        print(f"Error downloading inventory: {e}")
        inventory = None

    print("\n" + "="*70)
    print("DOWNLOAD COMPLETE")
    print("="*70)
    print(f"Output directory: {output_dir}")
    print(f"Files created:")
    if event:
        print(f"  - event_info.xml (event metadata)")
    if stream:
        print(f"  - waveforms.mseed ({len(stream)} traces)")
    if inventory:
        print(f"  - station_inventory.xml (station metadata)")

    return stream, inventory, event


def plot_waveforms(stream, event=None, output_dir="./nz_earthquake_data"):
    """
    Create basic plots of downloaded waveforms

    Parameters:
    -----------
    stream : obspy.Stream
        Waveform data to plot
    event : obspy.core.event.Event, optional
        Event information for annotation
    output_dir : str
        Directory to save plots
    """

    if stream is None or len(stream) == 0:
        print("No waveforms to plot")
        return

    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt

        # Plot all traces
        print("\nGenerating waveform plot...")
        fig = stream.plot(handle=True, show=False)
        plot_file = os.path.join(output_dir, "waveforms_all.png")
        fig.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"✓ Saved plot to: {plot_file}")

        # Plot first 3 components if available
        if len(stream) >= 3:
            print("Generating 3-component plot...")
            st_subset = stream[:3]
            fig = st_subset.plot(handle=True, show=False)
            plot_file = os.path.join(output_dir, "waveforms_3comp.png")
            fig.savefig(plot_file, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"✓ Saved plot to: {plot_file}")

    except Exception as e:
        print(f"Warning: Could not create plots: {e}")


def main():
    """Main function"""

    print("="*70)
    print("NEW ZEALAND EARTHQUAKE BROADBAND DATA DOWNLOAD")
    print("="*70)
    print("Event: 2025-11-06T08:09:54Z")
    print("Data source: IRIS FDSN")
    print("="*70)
    print()

    # Download data for the New Zealand earthquake
    stream, inventory, event = download_earthquake_data(
        event_time_str="2025-11-06T08:09:54Z",
        output_dir="./nz_earthquake_2025-11-06",
        network="IU,II,G",  # Global networks
        station="*",
        location="*",
        channel="BH*",  # Broadband high-gain
        time_before=300,  # 5 minutes before
        time_after=3600,  # 1 hour after
        max_radius=20.0,  # Within 20 degrees
        min_radius=0.0,
    )

    # Create plots if matplotlib is available
    if stream:
        try:
            plot_waveforms(stream, event, output_dir="./nz_earthquake_2025-11-06")
        except ImportError:
            print("\nNote: Install matplotlib to generate plots:")
            print("  pip install matplotlib")

    print("\nTo process this data further, use ObsPy:")
    print("  from obspy import read, read_inventory")
    print("  stream = read('nz_earthquake_2025-11-06/waveforms.mseed')")
    print("  inventory = read_inventory('nz_earthquake_2025-11-06/station_inventory.xml')")
    print("\nExample processing:")
    print("  stream.filter('bandpass', freqmin=0.01, freqmax=1.0)")
    print("  stream.detrend('linear')")
    print("  stream.taper(max_percentage=0.05)")


if __name__ == "__main__":
    main()
