from MeshFunc.locate_topo_data import locate_modify_topo
import gmsh
import numpy as np


def generate_mesh_fault(MeshFile, stl_fault, stl_box, xc, yc, zmax, zmin, L, B, ncfile):
    """
    Generate a 3D mesh from a watertight box STL and a fault STL.
    xc, yc, zmin: center and the bottom of the domain
    L, B: length and width of the domain
    """

    gmsh.initialize()
    gmsh.model.add("KaikouraFault")
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

    # --- Step 1: merge the fault STL (surface only)
    gmsh.merge(stl_fault)
    gmsh.model.occ.synchronize()

    # Collect all current surfaces (faults)
    fault_surfaces = gmsh.model.getEntities(dim=2)
    fault_faces = [tag for (dim, tag) in fault_surfaces]

    # --- Step 2: merge the box STL (watertight volume)
    gmsh.merge(stl_box)
    gmsh.model.occ.synchronize()

    # Convert imported STL surfaces into a closed volume
    box_surfaces = gmsh.model.getEntities(dim=2)
    box_surface_tags = [s[1] for s in box_surfaces if s not in fault_surfaces]
    if box_surface_tags:
        loop = gmsh.model.occ.addSurfaceLoop(box_surface_tags)
        gmsh.model.occ.addVolume([loop])
        gmsh.model.occ.synchronize()
    # gmsh.model.geo.synchronize()

    # --- Step 3: distance-based refinement near faults
    print("Stage 2: Applying distance-based refinement...")
    resolution = 1000.0  # base element size
    r = 1000.0  # distance scaling

    distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance, "FacesList", fault_faces)
    gmsh.model.mesh.field.setNumber(distance, "Sampling", 50)

    threshold1 = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold1, "IField", distance)
    gmsh.model.mesh.field.setNumber(threshold1, "LcMin", resolution)
    gmsh.model.mesh.field.setNumber(threshold1, "LcMax", 40 * resolution)
    gmsh.model.mesh.field.setNumber(threshold1, "DistMin", 0.01 * r)
    gmsh.model.mesh.field.setNumber(threshold1, "DistMax", 15 * r)

    minimum = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [threshold1])
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum)

    # --- Step 4: generate mesh
    gmsh.model.mesh.generate(3)

    gmsh.model.occ.synchronize()

    # --- Step 5: classify surfaces for BCs
    faces = gmsh.model.getEntities(dim=2)
    print("All faces:", faces)
    gmsh.model.mesh.createGeometry()

    absorb_end = []
    fault_group = []

    for surface in faces:
        print(surface)
        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
        print(com)
        if np.allclose(com, [xc, yc, zmax], atol=1e3):
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], 101)
            print("Free surface:", surface[1])
            nodeTags, nodeCoord, nodeParams = gmsh.model.mesh.getNodes(2, surface[1])
            vxyz = nodeCoord.reshape((-1, 3))

        elif np.sqrt((com[0] - xc) ** 2 + (com[1] - yc) ** 2) > 0.75 * (L + B) / 4.0:
            absorb_end.append(surface[1])
        elif np.allclose(com, [xc, yc, zmin], atol=1e3):
            absorb_end.append(surface[1])
        else:
            fault_group.append()

    gmsh.model.occ.synchronize()

    print("Absorbing BC:", absorb_end)
    print("Dynamic BC:", fault_faces)
    print("fault: ", fault_group)

    if absorb_end:
        gmsh.model.addPhysicalGroup(2, absorb_end, 105)
    if fault_group:
        gmsh.model.addPhysicalGroup(2, fault_group, 103)

    gmsh.model.occ.synchronize()

    volumes = gmsh.model.occ.getEntities(dim=3)
    print(volumes)
    gmsh.model.addPhysicalGroup(3, [1], 1)

    gmsh.option.setNumber("Mesh.MshFileVersion", 4)
    gmsh.option.setNumber("Mesh.SaveAll", 1)

    # --- Step 6: save mesh
    gmsh.write(MeshFile + ".msh4")
    gmsh.finalize()
