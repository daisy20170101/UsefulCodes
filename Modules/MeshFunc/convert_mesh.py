# generate a xdmf file for Paraview
import meshio

# Load the .msh file
def convert_mesh(MeshFile):
    
    mesh = meshio.read(
    MeshFile+ '.msh2',  # string, os.PathLike, or a buffer/open file
    file_format="gmsh",  # optional if filename is a path; inferred from extension
    )

    meshio.write(
    MeshFile+".vtk",  # str, os.PathLike, or buffer/ open file
    mesh,
    file_format="vtk",  # optional if first argument is a path; inferred from extension
    )
    


# Load the .msh file
def convert_mesh_msh4(MeshFile):
    
    mesh = meshio.read(
    MeshFile+ '.msh4',  # string, os.PathLike, or a buffer/open file
    file_format="gmsh",  # optional if filename is a path; inferred from extension
    )

    meshio.write(
    MeshFile+".vtk",  # str, os.PathLike, or buffer/ open file
    mesh,
    file_format="vtk",  # optional if first argument is a path; inferred from extension
    )
    