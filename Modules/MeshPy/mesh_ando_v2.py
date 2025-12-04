#!/usr/bin/env python3
import gmsh
import numpy as np


def mesh_fault_box(
    msh_file: str,
    fault_surface_tags,
    padding=(5e3, 5e3, 5e3),
    char_len=None,
    out_file="fault_box.msh",
):
    """
    Load a mesh/geometry file, take a list of surface tags representing faults,
    build a box around their center, embed the faults as internal boundaries,
    and generate a 3D tetrahedral mesh in that box.

    Parameters
    ----------
    msh_file : str
        Path to the input .msh (or any gmsh-readable) file that contains the fault surfaces.
        IMPORTANT: The file must have geometry information (e.g. saved with Mesh.SaveAll = 1).
    fault_surface_tags : list[int]
        List of geometric surface tags (dim=2) that represent the faults.
    padding : tuple[float, float, float]
        Extra half-width to add around the fault bounding box in x, y, z (meters).
    char_len : float or None
        Optional global characteristic length to control element size.
        If None, existing size settings are kept.
    out_file : str
        Name of the output .msh file with the 3D tetrahedral mesh.
    """

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    gmsh.merge(msh_file)

    # -------------------------------------------------------
    # 1. Compute bounding box and center of the fault surfaces
    # -------------------------------------------------------
    fault_surface_tags = list(fault_surface_tags)
    if not fault_surface_tags:
        raise ValueError("fault_surface_tags is empty.")

    xmin, ymin, zmin = np.inf, np.inf, np.inf
    xmax, ymax, zmax = -np.inf, -np.inf, -np.inf

    for s in fault_surface_tags:
        # Generic kernel-independent bounding box
        bxmin, bymin, bzmin, bxmax, bymax, bzmax = gmsh.model.getBoundingBox(2, s)
        xmin = min(xmin, bxmin)
        ymin = min(ymin, bymin)
        zmin = min(zmin, bzmin)
        xmax = max(xmax, bxmax)
        ymax = max(ymax, bymax)
        zmax = max(zmax, bzmax)

    if not np.isfinite([xmin, ymin, zmin, xmax, ymax, zmax]).all():
        gmsh.finalize()
        raise RuntimeError("Could not compute a finite bounding box for the fault surfaces.")

    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)

    px, py, pz = padding

    # Box extents
    bxmin = xmin - px
    bymin = ymin - py
    bzmin = zmin - pz
    bxmax = xmax + px
    bymax = ymax + py
    bzmax = zmax + pz

    Lx = bxmax - bxmin
    Ly = bymax - bymin
    Lz = bzmax - bzmin

    print(f"Fault center: ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    print(f"Box origin:  ({bxmin:.3f}, {bymin:.3f}, {bzmin:.3f})")
    print(f"Box size:    ({Lx:.3f}, {Ly:.3f}, {Lz:.3f})")

    # -------------------------------------------------------
    # 2. Create the box as a new OCC volume around the faults
    # -------------------------------------------------------
    # NOTE: We assume OCC kernel here for Booleans.
    box_tag = gmsh.model.occ.addBox(bxmin, bymin, bzmin, Lx, Ly, Lz)
    gmsh.model.occ.synchronize()

    # Physical group for faults (2D)
    pg_fault = gmsh.model.addPhysicalGroup(2, fault_surface_tags)
    gmsh.model.setPhysicalName(2, pg_fault, "fault_surfaces")

    # Physical group for box volume (3D)
    pg_box = gmsh.model.addPhysicalGroup(3, [box_tag])
    gmsh.model.setPhysicalName(3, pg_box, "fault_box_volume")

    # -------------------------------------------------------
    # 3. Embed the fault surfaces into the box as internal interface
    # -------------------------------------------------------
    gmsh.model.mesh.embed(2, fault_surface_tags, 3, box_tag)

    # Optional global size control
    if char_len is not None:
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", char_len)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", char_len)

    # -------------------------------------------------------
    # 4. Mesh: 2D (faults and box faces) then 3D tets
    #    SaveAll = 1 ensures internal 2D elements are kept
    # -------------------------------------------------------
    gmsh.model.mesh.clear()  # clear any existing mesh in the file
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(out_file)
    print(f"Written 3D mesh to: {out_file}")

    gmsh.finalize()


if __name__ == "__main__":
    # Example usage:
    #   - Input file: "fault_surfaces.msh"
    #   - Fault surfaces: 101, 103, 105
    #   - Padding: 5 km in each direction
    #   - Target size: 1000 m
    mesh_fault_box(
        msh_file="fault_surfaces.msh",
        fault_surface_tags=[101, 103, 105],
        padding=(5e3, 5e3, 5e3),
        char_len=1000.0,
        out_file="fault_box_3d.msh",
    )
