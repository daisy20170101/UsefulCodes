#!/usr/bin/env python3
"""
Check if fault surfaces can be embedded in a box volume

This checks for common issues that prevent successful embedding
"""

import gmsh
import sys


def check_embeddable(geo_file):
    """Check if faults from geo file can be embedded"""

    print("="*70)
    print("EMBEDDABILITY CHECK")
    print("="*70)

    gmsh.initialize()
    gmsh.model.add("check")

    # Load and scale
    gmsh.merge(geo_file)
    gmsh.model.geo.synchronize()

    all_entities = []
    for dim in range(4):
        all_entities.extend(gmsh.model.getEntities(dim))
    gmsh.model.geo.dilate(all_entities, 0, 0, 0, 1000, 1000, 1000)
    gmsh.model.geo.synchronize()

    # Filter surfaces
    surfaces = gmsh.model.getEntities(dim=2)
    fault_surfaces = []

    for surf_dim, surf_tag in surfaces:
        bbox = gmsh.model.getBoundingBox(surf_dim, surf_tag)
        center_z = (bbox[2] + bbox[5]) / 2

        if not (abs(center_z) < 1000 or center_z < -14000.0):
            fault_surfaces.append(surf_tag)

    print(f"\nFound {len(fault_surfaces)} fault surfaces to embed")

    if not fault_surfaces:
        print("✗ No fault surfaces!")
        gmsh.finalize()
        return False

    # Get overall bounds
    xmin, ymin, zmin = float('inf'), float('inf'), float('inf')
    xmax, ymax, zmax = float('-inf'), float('-inf'), float('-inf')

    print("\nFault surface bounds:")
    for i, surf_tag in enumerate(fault_surfaces):
        bbox = gmsh.model.getBoundingBox(2, surf_tag)
        print(f"  Surface {surf_tag}: "
              f"X=[{bbox[0]:.0f}, {bbox[3]:.0f}], "
              f"Y=[{bbox[1]:.0f}, {bbox[4]:.0f}], "
              f"Z=[{bbox[2]:.0f}, {bbox[5]:.0f}]")

        xmin = min(xmin, bbox[0])
        ymin = min(ymin, bbox[1])
        zmin = min(zmin, bbox[2])
        xmax = max(xmax, bbox[3])
        ymax = max(ymax, bbox[4])
        zmax = max(zmax, bbox[5])

    print(f"\nCombined fault bounds:")
    print(f"  X: [{xmin:.0f}, {xmax:.0f}] ({(xmax-xmin)/1e3:.1f} km)")
    print(f"  Y: [{ymin:.0f}, {ymax:.0f}] ({(ymax-ymin)/1e3:.1f} km)")
    print(f"  Z: [{zmin:.0f}, {zmax:.0f}] ({(zmax-zmin)/1e3:.1f} km)")

    # Create a simple box
    print("\nCreating test box...")
    box_length = 150000
    box_width = 150000
    box_depth = 120000

    buffer_x = max(10000, (xmax - xmin) * 0.1)
    buffer_y = max(10000, (ymax - ymin) * 0.1)

    x0 = xmin - box_length/2 - buffer_x
    x1 = xmax + box_length/2 + buffer_x
    y0 = ymin - box_width/2 - buffer_y
    y1 = ymax + box_width/2 + buffer_y
    z0 = -box_depth
    z1 = 1000.0

    print(f"Box: X=[{x0:.0f}, {x1:.0f}], Y=[{y0:.0f}, {y1:.0f}], Z=[{z0:.0f}, {z1:.0f}]")

    # Check if any fault extends outside box
    issues = []

    for surf_tag in fault_surfaces:
        bbox = gmsh.model.getBoundingBox(2, surf_tag)

        outside = []
        if bbox[0] < x0: outside.append(f"X too low ({bbox[0]:.0f} < {x0:.0f})")
        if bbox[3] > x1: outside.append(f"X too high ({bbox[3]:.0f} > {x1:.0f})")
        if bbox[1] < y0: outside.append(f"Y too low ({bbox[1]:.0f} < {y0:.0f})")
        if bbox[4] > y1: outside.append(f"Y too high ({bbox[4]:.0f} > {y1:.0f})")
        if bbox[2] < z0: outside.append(f"Z too low ({bbox[2]:.0f} < {z0:.0f})")
        if bbox[5] > z1: outside.append(f"Z too high ({bbox[5]:.0f} > {z1:.0f})")

        if outside:
            issues.append((surf_tag, outside))
            print(f"\n✗ Surface {surf_tag} extends outside box:")
            for issue in outside:
                print(f"     {issue}")

    # Summary
    print("\n" + "="*70)
    if not issues:
        print("✓ All fault surfaces fit within box")
        print("  Embedding should work")
    else:
        print(f"✗ {len(issues)} surface(s) extend outside box")
        print("  This will cause embedding to fail!")
        print("\n  Solutions:")
        print("  1. Increase box size (--box-length, --box-width, --box-depth)")
        print("  2. Remove problematic surfaces manually")
        print("  3. Trim fault surfaces to box boundaries")

    print("="*70)

    gmsh.finalize()

    return len(issues) == 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_embeddable.py <geo_file>")
        sys.exit(1)

    success = check_embeddable(sys.argv[1])
    sys.exit(0 if success else 1)
