#!/usr/bin/env pvpython
"""
Diagnose ParaView state file and check for data sources

Usage:
    pvpython check_state_file.py <state_file.pvsm>
"""

import sys
import os
from paraview.simple import *

def check_state_file(state_file):
    """Check state file and list available sources"""

    print("="*70)
    print("ParaView State File Diagnostic")
    print("="*70)

    if not os.path.exists(state_file):
        print(f"ERROR: State file not found: {state_file}")
        return False

    print(f"\nState file: {state_file}")
    print(f"File size: {os.path.getsize(state_file)} bytes")

    # Try to load state
    print("\nAttempting to load state file...")
    try:
        LoadState(state_file)
        print("✓ State file loaded successfully")
    except Exception as e:
        print(f"✗ Failed to load state file: {e}")
        return False

    # Check for data sources
    print("\n" + "-"*70)
    print("Data Sources:")
    print("-"*70)

    sources = GetSources()

    if not sources:
        print("✗ No data sources found in state file")
        print("\nThis state file cannot be used alone.")
        print("You need to either:")
        print("  1. Use an XDMF file as input: -i your_data.xdmf")
        print("  2. Load a state file that already contains data")
        return False

    print(f"✓ Found {len(sources)} data source(s):\n")

    for (name, id), source in sources.items():
        print(f"  Source: {name}")
        print(f"    ID: {id}")
        print(f"    Type: {type(source).__name__}")

        # Try to get data info
        try:
            if hasattr(source, 'PointData'):
                point_arrays = source.PointData
                if point_arrays:
                    print(f"    Point Data Arrays:")
                    for i in range(len(point_arrays)):
                        array = point_arrays[i]
                        print(f"      - {array.Name}")
        except:
            pass

        print()

    # Check for active source
    active_source = GetActiveSource()
    if active_source:
        print(f"✓ Active source: {active_source}")
    else:
        print("✗ No active source")

    # Check for render view
    print("\n" + "-"*70)
    print("Render View:")
    print("-"*70)

    try:
        render_view = GetActiveViewOrCreate('RenderView')
        print(f"✓ Render view available")
        print(f"  View size: {render_view.ViewSize}")
    except Exception as e:
        print(f"✗ No render view: {e}")

    print("\n" + "="*70)
    print("DIAGNOSIS COMPLETE")
    print("="*70)

    return True

def main():
    if len(sys.argv) < 2:
        print("Usage: pvpython check_state_file.py <state_file.pvsm>")
        print("\nExample:")
        print("  pvpython check_state_file.py rakeNew.pvsm")
        sys.exit(1)

    state_file = sys.argv[1]
    success = check_state_file(state_file)

    if not success:
        print("\n" + "="*70)
        print("RECOMMENDED WORKFLOW")
        print("="*70)
        print("\nInstead of using state file only, use XDMF input:\n")
        print("  pvpython generate_fault_snapshots.py -i your_data.xdmf\n")
        print("The state file is optional and only used for camera/view settings.")
        print("="*70)
        sys.exit(1)

if __name__ == '__main__':
    main()
