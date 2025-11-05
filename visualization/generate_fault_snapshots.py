#!/usr/bin/env pvpython
"""
Generate snapshots of SRs (slip rate strike) and SRd (slip rate dip) from fault output
Can load data from XDMF file or ParaView state file

Usage:
    pvpython generate_fault_snapshots.py -i input.xdmf [options]
    pvpython generate_fault_snapshots.py -s state.pvsm [options]

Requirements:
    - ParaView with pvpython
    - XDMF data file with SRs and SRd fields, or
    - State file: rakeNew.pvsm (in visualization folder or custom path)
"""

import sys
import os
import argparse
from paraview.simple import *

# Configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_STATE_FILE = os.path.join(SCRIPT_DIR, "rakeNew.pvsm")
FALLBACK_STATE_FILE = "/Users/DuoL/Documents/NSHM/Central/Paraviews/rakeNew.pvsm"
DEFAULT_OUTPUT_DIR = "./snapshots"
DEFAULT_RESOLUTION = (1920, 1080)

class FaultSnapshotGenerator:
    """Generate snapshots of fault slip rate components"""

    def __init__(self, data_file=None, state_file=None, output_dir=None, resolution=None):
        """
        Initialize snapshot generator

        Parameters:
        -----------
        data_file : str
            Path to XDMF data file (optional, can use state file instead)
        state_file : str
            Path to ParaView state file (.pvsm)
        output_dir : str
            Output directory for snapshots
        resolution : tuple
            Image resolution (width, height)
        """
        self.data_file = data_file
        self.state_file = state_file
        self.output_dir = output_dir or DEFAULT_OUTPUT_DIR
        self.resolution = resolution or DEFAULT_RESOLUTION

        # Determine state file to use
        if not self.state_file:
            if os.path.exists(DEFAULT_STATE_FILE):
                self.state_file = DEFAULT_STATE_FILE
            elif os.path.exists(FALLBACK_STATE_FILE):
                self.state_file = FALLBACK_STATE_FILE

        # Create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print(f"Created output directory: {self.output_dir}")

    def load_xdmf_file(self):
        """
        Load XDMF data file

        Returns:
        --------
        source : ParaView source object
        """
        if not self.data_file:
            raise ValueError("No data file specified")

        if not os.path.exists(self.data_file):
            raise FileNotFoundError(f"Data file not found: {self.data_file}")

        print(f"Loading XDMF file: {self.data_file}")

        try:
            # Load XDMF file
            source = XDMFReader(FileNames=[self.data_file])

            # Get render view
            renderView = GetActiveViewOrCreate('RenderView')

            # Show data
            display = Show(source, renderView)

            # Reset camera
            renderView.ResetCamera()

            # Apply state file settings if available
            if self.state_file and os.path.exists(self.state_file):
                print(f"Applying settings from state file: {self.state_file}")
                # Note: We can't fully load state with a different data source,
                # but we could load camera settings, etc.
                # For now, we'll just use default view settings

            print("XDMF file loaded successfully")
            return source

        except Exception as e:
            print(f"Error loading XDMF file: {e}")
            import traceback
            traceback.print_exc()
            return None

    def load_state(self):
        """Load ParaView state file"""
        if not self.state_file:
            print("No state file specified")
            return False

        if not os.path.exists(self.state_file):
            print(f"Warning: State file not found: {self.state_file}")
            return False

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
  # Load from XDMF file
  pvpython generate_fault_snapshots.py -i fault_output.xdmf

  # Load from XDMF file with custom settings
  pvpython generate_fault_snapshots.py -i fault_output.xdmf -s rakeNew.pvsm

  # Use state file only (legacy mode)
  pvpython generate_fault_snapshots.py -s /path/to/state.pvsm

  # High resolution output from XDMF
  pvpython generate_fault_snapshots.py -i data.xdmf -r 3840 2160

  # Generate only SRs snapshot
  pvpython generate_fault_snapshots.py -i data.xdmf --field SRs
        """
    )

    parser.add_argument('-i', '--input', '--data-file',
                       dest='data_file',
                       help='Path to XDMF data file')

    parser.add_argument('-s', '--state-file',
                       help='Path to ParaView state file (default: ./visualization/rakeNew.pvsm)')

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

    print("="*70)
    print("FAULT SNAPSHOT GENERATOR")
    print("="*70)

    # Create generator
    generator = FaultSnapshotGenerator(
        data_file=args.data_file,
        state_file=args.state_file,
        output_dir=args.output_dir,
        resolution=tuple(args.resolution)
    )

    # Determine which mode to use
    source = None

    if args.data_file:
        # Mode 1: Load XDMF file directly
        print(f"Input XDMF file: {args.data_file}")
        source = generator.load_xdmf_file()

        if source is None:
            print("Error: Failed to load XDMF file")
            sys.exit(1)

    elif args.state_file or generator.state_file:
        # Mode 2: Load from state file (legacy mode)
        print(f"State file: {generator.state_file or args.state_file}")

        if not generator.load_state():
            print("Error: Failed to load state file")
            sys.exit(1)

        # Get active source from state
        source = GetActiveSource()

        if source is None:
            print("Error: No data source found in state file")
            sys.exit(1)

    else:
        print("Error: Either --input (XDMF file) or --state-file must be specified")
        sys.exit(1)

    print(f"Output directory: {args.output_dir}")
    print(f"Resolution: {args.resolution[0]} x {args.resolution[1]}")
    print(f"Fields to generate: {args.field}")
    print("="*70)

    # Generate requested snapshots
    try:
        if args.field == 'all':
            generator.generate_SRs_snapshot(source, f"{args.output_name}_SRs")
            generator.generate_SRd_snapshot(source, f"{args.output_name}_SRd")
        elif args.field == 'SRs':
            generator.generate_SRs_snapshot(source, f"{args.output_name}_SRs")
        elif args.field == 'SRd':
            generator.generate_SRd_snapshot(source, f"{args.output_name}_SRd")

        print("\n" + "="*70)
        print("SUCCESS! Snapshots generated successfully")
        print("="*70)

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
