#!/usr/bin/env pvpython
"""
Generate snapshots of SRs (slip rate strike) and SRd (slip rate dip) from fault output
Uses ParaView state file for initial setup

Usage:
    pvpython generate_fault_snapshots.py [options]

Requirements:
    - ParaView with pvpython
    - State file: /Users/DuoL/Documents/NSHM/Central/Paraviews/rakeNew.pvsm
"""

import sys
import os
import argparse
from paraview.simple import *

# Configuration
DEFAULT_STATE_FILE = "/Users/DuoL/Documents/NSHM/Central/Paraviews/rakeNew.pvsm"
DEFAULT_OUTPUT_DIR = "./snapshots"
DEFAULT_RESOLUTION = (1920, 1080)

class FaultSnapshotGenerator:
    """Generate snapshots of fault slip rate components"""

    def __init__(self, state_file=None, output_dir=None, resolution=None):
        """
        Initialize snapshot generator

        Parameters:
        -----------
        state_file : str
            Path to ParaView state file (.pvsm)
        output_dir : str
            Output directory for snapshots
        resolution : tuple
            Image resolution (width, height)
        """
        self.state_file = state_file or DEFAULT_STATE_FILE
        self.output_dir = output_dir or DEFAULT_OUTPUT_DIR
        self.resolution = resolution or DEFAULT_RESOLUTION

        # Create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print(f"Created output directory: {self.output_dir}")

    def load_state(self):
        """Load ParaView state file"""
        if not os.path.exists(self.state_file):
            raise FileNotFoundError(f"State file not found: {self.state_file}")

        print(f"Loading state file: {self.state_file}")
        try:
            LoadState(self.state_file)
            print("State file loaded successfully")
            return True
        except Exception as e:
            print(f"Error loading state file: {e}")
            return False

    def get_active_view(self):
        """Get the active render view"""
        renderView = GetActiveViewOrCreate('RenderView')
        return renderView

    def set_color_by_field(self, source, field_name):
        """
        Set the color mapping to a specific field

        Parameters:
        -----------
        source : ParaView source
            Data source
        field_name : str
            Field name to color by (e.g., 'SRs', 'SRd')
        """
        # Get display properties
        display = GetDisplayProperties(source)

        # Color by the specified field
        ColorBy(display, ('POINTS', field_name))

        # Get color transfer function
        lut = GetColorTransferFunction(field_name)

        # Rescale to data range
        display.RescaleTransferFunctionToDataRange(True, False)

        # Update the view
        renderView = self.get_active_view()
        renderView.Update()

        print(f"Colored by field: {field_name}")

        return lut

    def setup_color_bar(self, field_name, lut):
        """
        Setup and configure color bar

        Parameters:
        -----------
        field_name : str
            Field name for the color bar label
        lut : ColorTransferFunction
            Lookup table
        """
        renderView = self.get_active_view()

        # Get or create scalar bar
        scalarBar = GetScalarBar(lut, renderView)

        # Configure scalar bar
        scalarBar.Title = field_name
        scalarBar.ComponentTitle = ''
        scalarBar.TitleFontSize = 16
        scalarBar.LabelFontSize = 14
        scalarBar.Visibility = 1

        # Position scalar bar
        scalarBar.Position = [0.85, 0.1]
        scalarBar.ScalarBarLength = 0.8

        print(f"Color bar configured for: {field_name}")

    def save_snapshot(self, output_filename, field_name=""):
        """
        Save current view as snapshot

        Parameters:
        -----------
        output_filename : str
            Output filename (without extension)
        field_name : str
            Field name (for filename)
        """
        renderView = self.get_active_view()

        # Set view size
        renderView.ViewSize = self.resolution

        # Generate full output path
        output_path = os.path.join(self.output_dir, f"{output_filename}.png")

        # Save screenshot
        SaveScreenshot(output_path, renderView,
                      ImageResolution=self.resolution,
                      TransparentBackground=0)

        print(f"Saved snapshot: {output_path}")
        return output_path

    def generate_SRs_snapshot(self, source=None, output_name="fault_SRs"):
        """
        Generate snapshot of SRs (slip rate strike component)

        Parameters:
        -----------
        source : ParaView source
            Data source (if None, uses active source)
        output_name : str
            Output filename prefix
        """
        print("\n" + "="*60)
        print("Generating SRs (Slip Rate Strike) Snapshot")
        print("="*60)

        if source is None:
            source = GetActiveSource()

        if source is None:
            print("Error: No active source found")
            return None

        # Set color by SRs field
        lut = self.set_color_by_field(source, 'SRs')

        # Setup color bar
        self.setup_color_bar('SRs (m/s)', lut)

        # Save snapshot
        output_path = self.save_snapshot(output_name, 'SRs')

        return output_path

    def generate_SRd_snapshot(self, source=None, output_name="fault_SRd"):
        """
        Generate snapshot of SRd (slip rate dip component)

        Parameters:
        -----------
        source : ParaView source
            Data source (if None, uses active source)
        output_name : str
            Output filename prefix
        """
        print("\n" + "="*60)
        print("Generating SRd (Slip Rate Dip) Snapshot")
        print("="*60)

        if source is None:
            source = GetActiveSource()

        if source is None:
            print("Error: No active source found")
            return None

        # Set color by SRd field
        lut = self.set_color_by_field(source, 'SRd')

        # Setup color bar
        self.setup_color_bar('SRd (m/s)', lut)

        # Save snapshot
        output_path = self.save_snapshot(output_name, 'SRd')

        return output_path

    def generate_all_snapshots(self):
        """Generate all snapshots (SRs and SRd)"""
        print("\n" + "="*60)
        print("FAULT SNAPSHOT GENERATOR")
        print("="*60)
        print(f"State file: {self.state_file}")
        print(f"Output directory: {self.output_dir}")
        print(f"Resolution: {self.resolution[0]} x {self.resolution[1]}")
        print("="*60)

        # Load state file
        if not self.load_state():
            return False

        # Get active source
        source = GetActiveSource()

        if source is None:
            print("Error: No data source found in state file")
            return False

        # Generate snapshots
        try:
            srs_path = self.generate_SRs_snapshot(source)
            srd_path = self.generate_SRd_snapshot(source)

            print("\n" + "="*60)
            print("SNAPSHOTS GENERATED SUCCESSFULLY")
            print("="*60)
            print(f"SRs snapshot: {srs_path}")
            print(f"SRd snapshot: {srd_path}")
            print("="*60)

            return True

        except Exception as e:
            print(f"\nError generating snapshots: {e}")
            import traceback
            traceback.print_exc()
            return False


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate fault slip rate snapshots using ParaView',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default state file and settings
  pvpython generate_fault_snapshots.py

  # Specify custom state file
  pvpython generate_fault_snapshots.py -s /path/to/state.pvsm

  # Specify output directory and resolution
  pvpython generate_fault_snapshots.py -o ./output -r 3840 2160

  # Generate only SRs snapshot
  pvpython generate_fault_snapshots.py --field SRs
        """
    )

    parser.add_argument('-s', '--state-file',
                       default=DEFAULT_STATE_FILE,
                       help=f'Path to ParaView state file (default: {DEFAULT_STATE_FILE})')

    parser.add_argument('-o', '--output-dir',
                       default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory for snapshots (default: {DEFAULT_OUTPUT_DIR})')

    parser.add_argument('-r', '--resolution',
                       nargs=2, type=int,
                       default=DEFAULT_RESOLUTION,
                       metavar=('WIDTH', 'HEIGHT'),
                       help=f'Image resolution (default: {DEFAULT_RESOLUTION[0]} {DEFAULT_RESOLUTION[1]})')

    parser.add_argument('-f', '--field',
                       choices=['SRs', 'SRd', 'all'],
                       default='all',
                       help='Field to generate snapshot for (default: all)')

    parser.add_argument('--output-name',
                       default='fault',
                       help='Output filename prefix (default: fault)')

    return parser.parse_args()


def main():
    """Main function"""
    args = parse_arguments()

    # Create generator
    generator = FaultSnapshotGenerator(
        state_file=args.state_file,
        output_dir=args.output_dir,
        resolution=tuple(args.resolution)
    )

    # Load state
    if not generator.load_state():
        sys.exit(1)

    # Get active source
    source = GetActiveSource()

    if source is None:
        print("Error: No data source found in state file")
        sys.exit(1)

    # Generate requested snapshots
    try:
        if args.field == 'all':
            generator.generate_SRs_snapshot(source, f"{args.output_name}_SRs")
            generator.generate_SRd_snapshot(source, f"{args.output_name}_SRd")
        elif args.field == 'SRs':
            generator.generate_SRs_snapshot(source, f"{args.output_name}_SRs")
        elif args.field == 'SRd':
            generator.generate_SRd_snapshot(source, f"{args.output_name}_SRd")

        print("\n" + "="*60)
        print("SUCCESS!")
        print("="*60)

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
