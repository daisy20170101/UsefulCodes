"""
Convenience script to run the full segmented meshing workflow:
1. Extract and mesh individual fault segments from .geo file
2. Merge segment meshes into final 3D mesh
"""

import os
import sys
from mesh_fault_segments import mesh_individual_fault_segments
from mesh_ando_from_segments import generate_mesh_from_fault_segments


def run_segmented_meshing(geo_file,
                          segment_dir='fault_segments',
                          lc_fault=2000,
                          lc_coarse=10000,
                          box_length=100e3,
                          box_width=100e3,
                          box_depth=100e3,
                          output_mesh='ando_fault_mesh_merged.msh',
                          remove_criteria=None):
    """
    Complete workflow: segment -> mesh -> merge

    Parameters:
    -----------
    geo_file : str
        Path to .geo file with fault geometry
    segment_dir : str
        Directory to store individual segment meshes
    lc_fault : float
        Mesh size on fault surfaces (meters)
    lc_coarse : float
        Mesh size away from faults (meters)
    box_length : float
        Box padding in X direction (meters)
    box_width : float
        Box padding in Y direction (meters)
    box_depth : float
        Box depth (meters)
    output_mesh : str
        Final output mesh filename
    remove_criteria : dict, optional
        Criteria for removing surfaces from .geo file

    Returns:
    --------
    output_mesh : str
        Path to final mesh file
    segment_info : list
        Information about each fault segment
    """

    print("=" * 80)
    print(" SEGMENTED FAULT MESHING WORKFLOW")
    print("=" * 80)
    print(f"Input .geo file: {geo_file}")
    print(f"Segment directory: {segment_dir}")
    print(f"Fault mesh size: {lc_fault} m")
    print(f"Coarse mesh size: {lc_coarse} m")
    print(f"Output mesh: {output_mesh}")
    print()

    # Step 1: Extract and mesh individual fault segments
    print("\n" + "=" * 80)
    print(" STEP 1: EXTRACTING AND MESHING INDIVIDUAL FAULT SEGMENTS")
    print("=" * 80 + "\n")

    segment_files, segment_info = mesh_individual_fault_segments(
        geo_file=geo_file,
        output_dir=segment_dir,
        lc=lc_fault,
        remove_criteria=remove_criteria
    )

    if not segment_files:
        raise ValueError("No fault segments were successfully meshed!")

    print(f"\n✓ Step 1 complete: {len(segment_files)} fault segments meshed")

    # Step 2: Merge segments into 3D mesh
    print("\n" + "=" * 80)
    print(" STEP 2: MERGING SEGMENTS INTO 3D MESH")
    print("=" * 80 + "\n")

    final_mesh = generate_mesh_from_fault_segments(
        segment_files=segment_files,
        lc_coarse=lc_coarse,
        box_length=box_length,
        box_width=box_width,
        box_depth=box_depth,
        output_mesh=output_mesh,
        use_distance_field=False
    )

    print("\n" + "=" * 80)
    print(" WORKFLOW COMPLETE!")
    print("=" * 80)
    print(f"✓ Meshed {len(segment_files)} fault segments")
    print(f"✓ Generated 3D mesh: {final_mesh}")
    print()

    return final_mesh, segment_info


if __name__ == "__main__":
    # Example usage
    import argparse

    parser = argparse.ArgumentParser(description="Segment and mesh fault geometry")
    parser.add_argument("geo_file", help="Path to .geo file")
    parser.add_argument("--segment-dir", default="fault_segments", help="Directory for segment meshes")
    parser.add_argument("--lc-fault", type=float, default=2000, help="Fault mesh size (m)")
    parser.add_argument("--lc-coarse", type=float, default=10000, help="Coarse mesh size (m)")
    parser.add_argument("--box-length", type=float, default=100e3, help="Box padding X (m)")
    parser.add_argument("--box-width", type=float, default=100e3, help="Box padding Y (m)")
    parser.add_argument("--box-depth", type=float, default=100e3, help="Box depth (m)")
    parser.add_argument("--output", default="ando_fault_mesh_merged.msh", help="Output mesh file")

    args = parser.parse_args()

    # Run the workflow
    final_mesh, segment_info = run_segmented_meshing(
        geo_file=args.geo_file,
        segment_dir=args.segment_dir,
        lc_fault=args.lc_fault,
        lc_coarse=args.lc_coarse,
        box_length=args.box_length,
        box_width=args.box_width,
        box_depth=args.box_depth,
        output_mesh=args.output,
        remove_criteria={
            'z_near_zero_tol': 1000.0,
            'z_below': -14000.0
        }
    )

    print("\nSegment details:")
    for info in segment_info:
        print(f"  {os.path.basename(info['file'])}: {info['num_elements']} elements, tag={info['tag']}")
