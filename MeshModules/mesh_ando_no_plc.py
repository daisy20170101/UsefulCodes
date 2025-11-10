#!/usr/bin/env python3
"""
Enhanced version with larger buffers to avoid PLC errors

Usage: python mesh_ando_no_plc.py fault.geo
"""

from mesh_ando_robust import generate_Ando_fault_robust

# Wrapper with larger buffers
def generate_without_plc_errors(geo_file, **kwargs):
    """
    Generate mesh with enhanced buffers to avoid PLC errors
    """
    # Default parameters optimized to avoid intersections
    defaults = {
        'lc': kwargs.get('lc', 2000),
        'lc_coarse': kwargs.get('lc_coarse', 10000),
        'box_length': kwargs.get('box_length', 150000),  # Increased from 100k
        'box_width': kwargs.get('box_width', 150000),    # Increased from 100k
        'box_depth': kwargs.get('box_depth', 120000),    # Increased from 100k
        'output_mesh': kwargs.get('output_mesh', 'ando_fault_mesh.msh'),
        'show_gui': kwargs.get('show_gui', False),
        'verbose': kwargs.get('verbose', True)
    }

    print("Using enhanced buffers to avoid PLC errors:")
    print(f"  Box padding: {defaults['box_length']/1e3:.0f} km × {defaults['box_width']/1e3:.0f} km")
    print(f"  Box depth: {defaults['box_depth']/1e3:.0f} km")

    return generate_Ando_fault_robust(geo_file, **defaults)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python mesh_ando_no_plc.py <geo_file>")
        sys.exit(1)

    result = generate_without_plc_errors(
        sys.argv[1],
        show_gui='--gui' in sys.argv
    )

    print(f"\n✓ Success! Generated {result['num_nodes']:,} nodes")
