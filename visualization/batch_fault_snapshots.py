#!/usr/bin/env pvpython
"""
Batch process fault outputs to generate SRs and SRd snapshots
for multiple timesteps or data files

Usage:
    pvpython batch_fault_snapshots.py --data-dir /path/to/data --output-dir ./snapshots
"""

import sys
import os
import argparse
import glob
from paraview.simple import *

# Configuration
DEFAULT_STATE_FILE = "/Users/DuoL/Documents/NSHM/Central/Paraviews/rakeNew.pvsm"
DEFAULT_OUTPUT_DIR = "./snapshots_batch"
DEFAULT_RESOLUTION = (1920, 1080)

class BatchFaultProcessor:
    """Batch process fault data files to generate snapshots"""

    def __init__(self, state_file=None, output_dir=None, resolution=None):
        self.state_file = state_file or DEFAULT_STATE_FILE
        self.output_dir = output_dir or DEFAULT_OUTPUT_DIR
        self.resolution = resolution or DEFAULT_RESOLUTION

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print(f"Created output directory: {self.output_dir}")

    def find_data_files(self, data_dir, pattern="*.xdmf"):
        """
        Find all data files matching pattern

        Parameters:
        -----------
        data_dir : str
            Directory containing data files
        pattern : str
            File pattern to match (default: *.xdmf)

        Returns:
        --------
        list : List of file paths
        """
        search_path = os.path.join(data_dir, pattern)
        files = sorted(glob.glob(search_path))
        print(f"Found {len(files)} data files matching '{pattern}'")
        return files

    def load_data_file(self, data_file):
        """
        Load a data file

        Parameters:
        -----------
        data_file : str
            Path to data file

        Returns:
        --------
        source : ParaView source object
        """
        print(f"\nLoading data file: {data_file}")

        # Determine file type and load accordingly
        ext = os.path.splitext(data_file)[1].lower()

        if ext == '.xdmf':
            source = XDMFReader(FileNames=[data_file])
        elif ext == '.vtk':
            source = LegacyVTKReader(FileNames=[data_file])
        elif ext == '.vtu':
            source = XMLUnstructuredGridReader(FileNames=[data_file])
        elif ext == '.pvd':
            source = PVDReader(FileName=data_file)
        else:
            print(f"Warning: Unknown file type '{ext}', attempting generic reader")
            source = OpenDataFile(data_file)

        return source

    def generate_snapshots_for_file(self, data_file, timestep=None):
        """
        Generate SRs and SRd snapshots for a single data file

        Parameters:
        -----------
        data_file : str
            Path to data file
        timestep : int
            Timestep number (for filename)
        """
        # Load data
        source = self.load_data_file(data_file)

        if source is None:
            print(f"Error: Could not load data file: {data_file}")
            return False

        # Get render view
        renderView = GetActiveViewOrCreate('RenderView')
        renderView.ViewSize = self.resolution

        # Get display properties
        display = Show(source, renderView)

        # Base filename
        base_name = os.path.splitext(os.path.basename(data_file))[0]
        if timestep is not None:
            base_name = f"{base_name}_t{timestep:05d}"

        # Generate SRs snapshot
        print(f"Generating SRs snapshot...")
        self.generate_field_snapshot(source, display, renderView,
                                     'SRs', f"{base_name}_SRs")

        # Generate SRd snapshot
        print(f"Generating SRd snapshot...")
        self.generate_field_snapshot(source, display, renderView,
                                     'SRd', f"{base_name}_SRd")

        # Clean up
        Delete(source)

        return True

    def generate_field_snapshot(self, source, display, renderView,
                                field_name, output_name):
        """
        Generate snapshot for a specific field

        Parameters:
        -----------
        source : ParaView source
        display : Display properties
        renderView : Render view
        field_name : str
            Field name ('SRs' or 'SRd')
        output_name : str
            Output filename (without extension)
        """
        # Color by field
        ColorBy(display, ('POINTS', field_name))

        # Get color transfer function
        lut = GetColorTransferFunction(field_name)

        # Rescale to data range
        display.RescaleTransferFunctionToDataRange(True, False)

        # Setup color bar
        scalarBar = GetScalarBar(lut, renderView)
        scalarBar.Title = field_name
        scalarBar.ComponentTitle = 'm/s'
        scalarBar.TitleFontSize = 16
        scalarBar.LabelFontSize = 14
        scalarBar.Visibility = 1
        scalarBar.Position = [0.85, 0.1]
        scalarBar.ScalarBarLength = 0.8

        # Update view
        renderView.Update()

        # Save snapshot
        output_path = os.path.join(self.output_dir, f"{output_name}.png")
        SaveScreenshot(output_path, renderView,
                      ImageResolution=self.resolution,
                      TransparentBackground=0)

        print(f"  Saved: {output_path}")

    def process_directory(self, data_dir, pattern="*.xdmf"):
        """
        Process all data files in a directory

        Parameters:
        -----------
        data_dir : str
            Directory containing data files
        pattern : str
            File pattern to match
        """
        print("="*70)
        print("BATCH FAULT SNAPSHOT PROCESSOR")
        print("="*70)
        print(f"Data directory: {data_dir}")
        print(f"File pattern: {pattern}")
        print(f"Output directory: {self.output_dir}")
        print(f"Resolution: {self.resolution[0]} x {self.resolution[1]}")
        print("="*70)

        # Find data files
        data_files = self.find_data_files(data_dir, pattern)

        if not data_files:
            print(f"No data files found in {data_dir}")
            return False

        # Process each file
        success_count = 0
        for i, data_file in enumerate(data_files, 1):
            print(f"\n[{i}/{len(data_files)}] Processing: {os.path.basename(data_file)}")

            try:
                if self.generate_snapshots_for_file(data_file, timestep=i):
                    success_count += 1
            except Exception as e:
                print(f"Error processing {data_file}: {e}")
                import traceback
                traceback.print_exc()

        print("\n" + "="*70)
        print(f"COMPLETED: {success_count}/{len(data_files)} files processed successfully")
        print("="*70)

        return True

    def process_timesteps_from_state(self):
        """
        Process timesteps from loaded ParaView state file
        """
        print("="*70)
        print("PROCESSING TIMESTEPS FROM STATE FILE")
        print("="*70)

        # Load state file
        if not os.path.exists(self.state_file):
            print(f"Error: State file not found: {self.state_file}")
            return False

        print(f"Loading state file: {self.state_file}")
        LoadState(self.state_file)

        # Get active source
        source = GetActiveSource()
        if source is None:
            print("Error: No active source in state file")
            return False

        # Get render view
        renderView = GetActiveViewOrCreate('RenderView')
        renderView.ViewSize = self.resolution

        # Get display
        display = GetDisplayProperties(source)

        # Get time information
        animationScene = GetAnimationScene()
        timeKeeper = animationScene.TimeKeeper
        timesteps = timeKeeper.TimestepValues

        print(f"Found {len(timesteps)} timesteps")

        # Process each timestep
        for i, time in enumerate(timesteps):
            print(f"\n[{i+1}/{len(timesteps)}] Processing timestep {i}: t = {time}")

            # Set time
            animationScene.AnimationTime = time
            renderView.Update()

            # Generate snapshots
            self.generate_field_snapshot(source, display, renderView,
                                        'SRs', f"fault_t{i:05d}_SRs")
            self.generate_field_snapshot(source, display, renderView,
                                        'SRd', f"fault_t{i:05d}_SRd")

        print("\n" + "="*70)
        print(f"COMPLETED: Processed {len(timesteps)} timesteps")
        print("="*70)

        return True


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Batch process fault data to generate snapshots',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all XDMF files in a directory
  pvpython batch_fault_snapshots.py -d /path/to/data -o ./snapshots

  # Process VTK files with custom pattern
  pvpython batch_fault_snapshots.py -d /path/to/data -p "fault_*.vtk"

  # Process timesteps from state file
  pvpython batch_fault_snapshots.py --state-mode -s /path/to/state.pvsm

  # High resolution output
  pvpython batch_fault_snapshots.py -d /path/to/data -r 3840 2160
        """
    )

    parser.add_argument('-d', '--data-dir',
                       help='Directory containing data files')

    parser.add_argument('-p', '--pattern',
                       default='*.xdmf',
                       help='File pattern to match (default: *.xdmf)')

    parser.add_argument('-o', '--output-dir',
                       default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')

    parser.add_argument('-r', '--resolution',
                       nargs=2, type=int,
                       default=DEFAULT_RESOLUTION,
                       metavar=('WIDTH', 'HEIGHT'),
                       help=f'Image resolution (default: {DEFAULT_RESOLUTION[0]} {DEFAULT_RESOLUTION[1]})')

    parser.add_argument('-s', '--state-file',
                       default=DEFAULT_STATE_FILE,
                       help=f'ParaView state file (default: {DEFAULT_STATE_FILE})')

    parser.add_argument('--state-mode',
                       action='store_true',
                       help='Process timesteps from state file instead of data directory')

    return parser.parse_args()


def main():
    """Main function"""
    args = parse_arguments()

    # Create processor
    processor = BatchFaultProcessor(
        state_file=args.state_file,
        output_dir=args.output_dir,
        resolution=tuple(args.resolution)
    )

    # Process based on mode
    if args.state_mode:
        # Process timesteps from state file
        success = processor.process_timesteps_from_state()
    else:
        # Process data directory
        if not args.data_dir:
            print("Error: --data-dir is required when not in --state-mode")
            sys.exit(1)

        if not os.path.exists(args.data_dir):
            print(f"Error: Data directory not found: {args.data_dir}")
            sys.exit(1)

        success = processor.process_directory(args.data_dir, args.pattern)

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
