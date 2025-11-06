#!/usr/bin/env pvpython
"""
Generate snapshot from ParaView state file with custom data and parameters

Loads vrNew2.pvsm state file and allows customization of:
- Fault data file
- Timestep to visualize
- Variable name to color by
- Color scale range (vmin, vmax)
- Automatic camera adjustment to center and zoom on fault (150% view)

The camera automatically centers on the fault geometry and zooms to show
150% of the fault bounds, giving nice composition with margins.

Usage:
    pvpython snapshot_from_state.py -i fault.xdmf -t 10 -v SRs
    pvpython snapshot_from_state.py -i fault.xdmf -t 0 -v SRd -o output.png
    pvpython snapshot_from_state.py -i fault.xdmf -t 10 -v SRs --vmin 0 --vmax 10
"""

import sys
import os
import argparse
from paraview.simple import *

# Configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_STATE_FILE = os.path.join(SCRIPT_DIR, "vrNew2.pvsm")
FALLBACK_STATE_FILE = "/Users/DuoL/Documents/GitHub/UsefulCodes/visualization/vrNew2.pvsm"
DEFAULT_OUTPUT = "./snapshot.png"
DEFAULT_RESOLUTION = (1920, 1080)


class StateFileSnapshotGenerator:
    """Generate snapshot from ParaView state file with custom parameters"""

    def __init__(self, state_file, data_file, timestep, variable, output_file, resolution, vmin=None, vmax=None):
        """
        Initialize generator

        Parameters:
        -----------
        state_file : str
            Path to ParaView state file (.pvsm)
        data_file : str
            Path to fault data file (XDMF, VTK, etc.)
        timestep : int
            Timestep index to visualize
        variable : str
            Variable name to color by (e.g., 'SRs', 'SRd', 'Vs', etc.)
        output_file : str
            Output file path for snapshot
        resolution : tuple
            Image resolution (width, height)
        vmin : float, optional
            Minimum value for color scale (default: auto from data)
        vmax : float, optional
            Maximum value for color scale (default: auto from data)
        """
        self.state_file = state_file
        self.data_file = data_file
        self.timestep = timestep
        self.variable = variable
        self.output_file = output_file
        self.resolution = resolution
        self.vmin = vmin
        self.vmax = vmax

        # Validate inputs
        if not os.path.exists(self.state_file):
            raise FileNotFoundError(f"State file not found: {self.state_file}")

        if not os.path.exists(self.data_file):
            raise FileNotFoundError(f"Data file not found: {self.data_file}")

        # Validate vmin/vmax
        if self.vmin is not None and self.vmax is not None:
            if self.vmin >= self.vmax:
                raise ValueError(f"vmin ({self.vmin}) must be less than vmax ({self.vmax})")

    def load_state_and_replace_data(self):
        """
        Load state file and replace data source with user's file

        Returns:
        --------
        source : ParaView source object
        """
        print(f"Loading state file: {self.state_file}")

        try:
            # Load the state file
            LoadState(self.state_file)
            print("✓ State file loaded successfully")

            # Get the active source from the state file
            original_source = GetActiveSource()

            if original_source is None:
                print("Warning: No active source in state file")
                print("Attempting to get first available source...")
                sources = GetSources()
                if sources:
                    # Get first source
                    original_source = list(sources.values())[0]
                    SetActiveSource(original_source)

            if original_source:
                print(f"✓ Found original source: {original_source}")

            # Load new data file
            print(f"\nLoading data file: {self.data_file}")
            new_source = self.load_data_file()

            if new_source is None:
                raise RuntimeError("Failed to load data file")

            # Get render view
            renderView = GetActiveViewOrCreate('RenderView')

            # Replace the original source's display with new source
            if original_source:
                # Get original display properties to preserve settings
                original_display = GetDisplayProperties(original_source, renderView)

                # Hide original source
                Hide(original_source, renderView)

                # Show new source
                new_display = Show(new_source, renderView)

                # Copy some display properties from original if available
                try:
                    if hasattr(original_display, 'Representation'):
                        new_display.Representation = original_display.Representation
                except:
                    pass

            else:
                # No original source, just show the new one
                Show(new_source, renderView)

            # Set as active source
            SetActiveSource(new_source)

            # Reset camera to fit new data
            renderView.ResetCamera()

            print("✓ Data source replaced successfully")
            return new_source

        except Exception as e:
            print(f"Error loading state and replacing data: {e}")
            import traceback
            traceback.print_exc()
            return None

    def load_data_file(self):
        """
        Load data file based on extension

        Returns:
        --------
        source : ParaView source object
        """
        ext = os.path.splitext(self.data_file)[1].lower()

        try:
            if ext == '.xdmf':
                source = XDMFReader(FileNames=[self.data_file])
                print("✓ Loaded XDMF file")
            elif ext == '.vtk':
                source = LegacyVTKReader(FileNames=[self.data_file])
                print("✓ Loaded VTK file")
            elif ext == '.vtu':
                source = XMLUnstructuredGridReader(FileNames=[self.data_file])
                print("✓ Loaded VTU file")
            elif ext == '.pvd':
                source = PVDReader(FileName=self.data_file)
                print("✓ Loaded PVD file")
            else:
                print(f"Unknown file type '{ext}', trying generic reader")
                source = OpenDataFile(self.data_file)

            return source

        except Exception as e:
            print(f"Error loading data file: {e}")
            import traceback
            traceback.print_exc()
            return None

    def set_timestep(self, source):
        """
        Set the visualization to a specific timestep

        Parameters:
        -----------
        source : ParaView source
            Data source
        """
        try:
            # Get available timesteps
            animationScene = GetAnimationScene()
            timeKeeper = animationScene.TimeKeeper
            timesteps = timeKeeper.TimestepValues

            if not timesteps or len(timesteps) == 0:
                print("Warning: No timesteps found in data")
                return

            print(f"\nAvailable timesteps: {len(timesteps)} timesteps")

            # Validate timestep index
            if self.timestep < 0 or self.timestep >= len(timesteps):
                print(f"Warning: Timestep index {self.timestep} out of range [0, {len(timesteps)-1}]")
                print(f"Using timestep 0 instead")
                self.timestep = 0

            # Get the time value for this timestep
            time_value = timesteps[self.timestep]
            print(f"Setting timestep {self.timestep}: t = {time_value}")

            # Set animation time
            animationScene.AnimationTime = time_value

            # CRITICAL: Force pipeline to update with new timestep data
            source.UpdatePipeline(time_value)

            # Update view
            renderView = GetActiveViewOrCreate('RenderView')
            renderView.Update()

            # Force render
            Render(renderView)

            print(f"✓ Timestep set to {self.timestep} (t = {time_value})")

        except Exception as e:
            print(f"Error setting timestep: {e}")
            import traceback
            traceback.print_exc()

    def set_coloring(self, source):
        """
        Set the coloring to the specified variable

        Parameters:
        -----------
        source : ParaView source
            Data source
        """
        try:
            print(f"\nSetting color mapping to variable: {self.variable}")

            # Get display properties
            renderView = GetActiveViewOrCreate('RenderView')
            display = GetDisplayProperties(source, renderView)

            # Check if variable exists in the data
            point_arrays = source.PointData
            cell_arrays = source.CellData

            variable_found = False
            variable_association = None

            # Check point data
            if point_arrays:
                for i in range(len(point_arrays)):
                    if point_arrays[i].Name == self.variable:
                        variable_found = True
                        variable_association = 'POINTS'
                        print(f"✓ Found '{self.variable}' in Point Data")
                        break

            # Check cell data if not found in points
            if not variable_found and cell_arrays:
                for i in range(len(cell_arrays)):
                    if cell_arrays[i].Name == self.variable:
                        variable_found = True
                        variable_association = 'CELLS'
                        print(f"✓ Found '{self.variable}' in Cell Data")
                        break

            if not variable_found:
                print(f"Warning: Variable '{self.variable}' not found in data")
                print("\nAvailable Point Data arrays:")
                if point_arrays:
                    for i in range(len(point_arrays)):
                        print(f"  - {point_arrays[i].Name}")
                print("\nAvailable Cell Data arrays:")
                if cell_arrays:
                    for i in range(len(cell_arrays)):
                        print(f"  - {cell_arrays[i].Name}")
                raise ValueError(f"Variable '{self.variable}' not found in dataset")

            # Color by the specified variable
            ColorBy(display, (variable_association, self.variable))

            # Get color transfer function
            lut = GetColorTransferFunction(self.variable)

            # Set color range
            if self.vmin is not None and self.vmax is not None:
                # Use custom range
                print(f"  Setting custom color range: [{self.vmin}, {self.vmax}]")
                lut.RescaleTransferFunction(self.vmin, self.vmax)
                display.SetScalarBarVisibility(renderView, True)
            elif self.vmin is not None or self.vmax is not None:
                # Only one bound specified - get data range and use it for the other
                display.RescaleTransferFunctionToDataRange(True, False)
                data_range = lut.GetDataRange()
                actual_vmin = self.vmin if self.vmin is not None else data_range[0]
                actual_vmax = self.vmax if self.vmax is not None else data_range[1]
                print(f"  Setting partial custom color range: [{actual_vmin}, {actual_vmax}]")
                lut.RescaleTransferFunction(actual_vmin, actual_vmax)
                display.SetScalarBarVisibility(renderView, True)
            else:
                # Automatic range from data
                print(f"  Using automatic color range from data")
                display.RescaleTransferFunctionToDataRange(True, False)

            # Setup color bar
            scalarBar = GetScalarBar(lut, renderView)
            scalarBar.Title = self.variable
            scalarBar.ComponentTitle = ''
            scalarBar.Visibility = 1

            # Update view
            renderView.Update()

            print(f"✓ Color mapping applied to '{self.variable}'")

        except Exception as e:
            print(f"Error setting coloring: {e}")
            import traceback
            traceback.print_exc()
            raise

    def adjust_camera_to_fault(self, source):
        """
        Adjust camera to center and zoom on fault geometry with 150% outline

        Parameters:
        -----------
        source : ParaView source
            Data source to zoom to
        """
        try:
            print(f"\nAdjusting camera view to fault geometry...")

            # Get render view
            renderView = GetActiveViewOrCreate('RenderView')

            # Get data bounds [xmin, xmax, ymin, ymax, zmin, zmax]
            dataInfo = source.GetDataInformation()
            bounds = dataInfo.GetBounds()

            if not bounds or len(bounds) != 6:
                print("Warning: Could not get data bounds, using default camera")
                return

            xmin, xmax, ymin, ymax, zmin, zmax = bounds

            # Calculate center of geometry
            center_x = (xmin + xmax) / 2.0
            center_y = (ymin + ymax) / 2.0
            center_z = (zmin + zmax) / 2.0

            # Calculate size of geometry
            size_x = xmax - xmin
            size_y = ymax - ymin
            size_z = zmax - zmin

            print(f"  Fault bounds: X=[{xmin:.1f}, {xmax:.1f}], Y=[{ymin:.1f}, {ymax:.1f}], Z=[{zmin:.1f}, {zmax:.1f}]")
            print(f"  Fault center: ({center_x:.1f}, {center_y:.1f}, {center_z:.1f})")
            print(f"  Fault size: {size_x:.1f} x {size_y:.1f} x {size_z:.1f}")

            # Expand bounds by 150% (multiply each dimension by 1.5)
            # This means the geometry takes up 2/3 of view with nice margins
            expansion_factor = 1.5
            expanded_size_x = size_x * expansion_factor
            expanded_size_y = size_y * expansion_factor
            expanded_size_z = size_z * expansion_factor

            # Calculate expanded bounds centered on geometry
            expanded_xmin = center_x - expanded_size_x / 2.0
            expanded_xmax = center_x + expanded_size_x / 2.0
            expanded_ymin = center_y - expanded_size_y / 2.0
            expanded_ymax = center_y + expanded_size_y / 2.0
            expanded_zmin = center_z - expanded_size_z / 2.0
            expanded_zmax = center_z + expanded_size_z / 2.0

            expanded_bounds = [expanded_xmin, expanded_xmax,
                             expanded_ymin, expanded_ymax,
                             expanded_zmin, expanded_zmax]

            print(f"  Expanded bounds (150%): X=[{expanded_xmin:.1f}, {expanded_xmax:.1f}], Y=[{expanded_ymin:.1f}, {expanded_ymax:.1f}], Z=[{expanded_zmin:.1f}, {expanded_zmax:.1f}]")

            # Reset camera to expanded bounds
            renderView.ResetCamera(expanded_bounds)

            # Set focal point to center of geometry
            renderView.CameraFocalPoint = [center_x, center_y, center_z]

            # Update view
            renderView.Update()
            Render(renderView)

            print(f"✓ Camera adjusted: centered on fault with 150% view")

        except Exception as e:
            print(f"Warning: Error adjusting camera: {e}")
            print("Continuing with default camera view...")
            import traceback
            traceback.print_exc()

    def save_snapshot(self):
        """
        Save the current view as a snapshot
        """
        try:
            print(f"\nSaving snapshot to: {self.output_file}")

            # Get render view
            renderView = GetActiveViewOrCreate('RenderView')

            # Set view size
            renderView.ViewSize = self.resolution

            # Create output directory if needed
            output_dir = os.path.dirname(self.output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)
                print(f"Created output directory: {output_dir}")

            # Save screenshot
            SaveScreenshot(self.output_file, renderView,
                          ImageResolution=self.resolution,
                          TransparentBackground=0)

            print(f"✓ Snapshot saved successfully")
            print(f"  File: {self.output_file}")
            print(f"  Resolution: {self.resolution[0]} x {self.resolution[1]}")

        except Exception as e:
            print(f"Error saving snapshot: {e}")
            import traceback
            traceback.print_exc()
            raise

    def generate(self):
        """
        Main workflow: load state, set parameters, save snapshot
        """
        print("="*70)
        print("PARAVIEW STATE-BASED SNAPSHOT GENERATOR")
        print("="*70)
        print(f"State file: {self.state_file}")
        print(f"Data file: {self.data_file}")
        print(f"Timestep: {self.timestep}")
        print(f"Variable: {self.variable}")
        if self.vmin is not None or self.vmax is not None:
            vmin_str = str(self.vmin) if self.vmin is not None else "auto"
            vmax_str = str(self.vmax) if self.vmax is not None else "auto"
            print(f"Color range: [{vmin_str}, {vmax_str}]")
        print(f"Output: {self.output_file}")
        print("="*70)

        try:
            # Load state and replace data source
            source = self.load_state_and_replace_data()

            if source is None:
                raise RuntimeError("Failed to load and setup data source")

            # Set timestep
            self.set_timestep(source)

            # Set variable coloring
            self.set_coloring(source)

            # Adjust camera to center and zoom on fault (150% view)
            self.adjust_camera_to_fault(source)

            # Save snapshot
            self.save_snapshot()

            print("\n" + "="*70)
            print("SUCCESS!")
            print("="*70)
            return True

        except Exception as e:
            print("\n" + "="*70)
            print("FAILED!")
            print("="*70)
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
            return False


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate snapshot from ParaView state file with custom parameters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate snapshot of SRs at timestep 10
  pvpython snapshot_from_state.py -i fault.xdmf -t 10 -v SRs

  # Generate snapshot of SRd at timestep 0
  pvpython snapshot_from_state.py -i fault.xdmf -t 0 -v SRd -o output.png

  # Generate snapshot with custom state file
  pvpython snapshot_from_state.py -i fault.xdmf -t 5 -v Vs -s custom.pvsm

  # High resolution output
  pvpython snapshot_from_state.py -i fault.xdmf -t 10 -v SRs -r 3840 2160

  # Custom color scale range
  pvpython snapshot_from_state.py -i fault.xdmf -t 10 -v SRs --vmin 0 --vmax 10

  # Set only maximum value (minimum from data)
  pvpython snapshot_from_state.py -i fault.xdmf -t 10 -v SRd --vmax 5.0
        """
    )

    parser.add_argument('-i', '--input', '--data-file',
                       dest='data_file',
                       required=True,
                       help='Path to fault data file (XDMF, VTK, VTU, PVD)')

    parser.add_argument('-t', '--timestep',
                       type=int,
                       required=True,
                       help='Timestep index to visualize (0-based)')

    parser.add_argument('-v', '--variable',
                       required=True,
                       help='Variable name to color by (e.g., SRs, SRd, Vs)')

    parser.add_argument('-s', '--state-file',
                       default=None,
                       help='Path to ParaView state file (default: vrNew2.pvsm in script dir)')

    parser.add_argument('-o', '--output',
                       default=DEFAULT_OUTPUT,
                       help=f'Output snapshot file path (default: {DEFAULT_OUTPUT})')

    parser.add_argument('-r', '--resolution',
                       nargs=2, type=int,
                       default=DEFAULT_RESOLUTION,
                       metavar=('WIDTH', 'HEIGHT'),
                       help=f'Image resolution (default: {DEFAULT_RESOLUTION[0]} {DEFAULT_RESOLUTION[1]})')

    parser.add_argument('--vmin',
                       type=float,
                       default=None,
                       help='Minimum value for color scale (default: auto from data)')

    parser.add_argument('--vmax',
                       type=float,
                       default=None,
                       help='Maximum value for color scale (default: auto from data)')

    return parser.parse_args()


def main():
    """Main function"""
    args = parse_arguments()

    # Determine state file to use
    state_file = args.state_file
    if not state_file:
        if os.path.exists(DEFAULT_STATE_FILE):
            state_file = DEFAULT_STATE_FILE
        elif os.path.exists(FALLBACK_STATE_FILE):
            state_file = FALLBACK_STATE_FILE
        else:
            print(f"Error: State file not found")
            print(f"  Looked for: {DEFAULT_STATE_FILE}")
            print(f"  Looked for: {FALLBACK_STATE_FILE}")
            print(f"  Please specify with -s/--state-file")
            sys.exit(1)

    # Create generator
    generator = StateFileSnapshotGenerator(
        state_file=state_file,
        data_file=args.data_file,
        timestep=args.timestep,
        variable=args.variable,
        output_file=args.output,
        resolution=tuple(args.resolution),
        vmin=args.vmin,
        vmax=args.vmax
    )

    # Generate snapshot
    success = generator.generate()

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
