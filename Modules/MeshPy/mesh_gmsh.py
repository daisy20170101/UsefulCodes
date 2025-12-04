
import gmsh
import numpy as np
from MeshFunc.locate_topo_data import *
from MeshPy.gmsh_roughness import add_roughness_to_surface

def generate_mesh_single_fault(MeshFile,stlFile1,domainBox,ncfile,rough=True,n_patches=10):
    
    lc = 20.0e3
    N = 200
    
    PI = np.pi/180


    gmsh.initialize()
    gmsh.model.add("HikWhaipapaFault")
    gmsh.option.setNumber("Mesh.Optimize", 1)  # Enables mesh optimization (default is off)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)  # Use Netgen mesh optimizer

    
    # Adjust the overall mesh size by modifying the characteristic length factor.
    gmsh.option.setNumber("Mesh.AngleToleranceFacetOverlap", 0.2)

    r = 2000.0
    resolution = 1000.0


    
    # Step 1: Load the STL surface fi le
    gmsh.merge(stlFile1)
 

    # Step 3: Create a compound entity containing both the box and the loaded STL surface
    # Step 4: Merge the compound entity
    gmsh.model.geo.synchronize()

    faces = gmsh.model.getEntities(2)
    print('fault faces:',faces)

    if rough:
        add_roughness_to_surface(
    		surface_tag=faces[0][1],
    		n_patches=n_patches,                      # Number of patches
    		radius_range=(200.0, 1000.0),     # 2-8 km patch radii
    		amplitude_factor=0.2,              # Roughness = 10% of radius
    		blend_exponent=2.0,                # Controls transition smoothness
    		random_seed=42                     # For reproducibility
	)

    gmsh.merge(domainBox)
    # gmsh.merge() # larger domain

    # asign physical group in gmsh
    gmsh.model.geo.synchronize()

    faces = gmsh.model.getEntities(dim=2)
    print('all faces:',faces)

    sl1 = gmsh.model.geo.addSurfaceLoop([1,2,3])
    v1 = gmsh.model.geo.addVolume([sl1])
    gmsh.model.geo.synchronize()
    
    gmsh.model.addPhysicalGroup(3, [v1], 1 )

    
    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(3)

    # asign physical group in gmsh
    faces = gmsh.model.getEntities(dim=2)
    print('all faces:',faces)

    absorb_end = [3]
    faults_end = [1]
    free_surf = [2]

    nodeTags, nodeCoord, nodeParams = gmsh.model.mesh.getNodes(2,free_surf[0])
    vxyz = nodeCoord.reshape((-1, 3)) 

    znew = locate_modify_topo(ncfile,vxyz[:,0],vxyz[:,1])
    
    for i in range(0, len(znew)): 
        if (znew[i]< vxyz[i,2]): # vxyz[i,2]
            gmsh.model.mesh.setNode(nodeTags[i], [nodeCoord[i*3],nodeCoord[i*3+1],vxyz[i,2]],[1,1])
        else:
            gmsh.model.mesh.setNode(nodeTags[i], [nodeCoord[i*3],nodeCoord[i*3+1],znew[i]],[1,1])

    print("dynamic BC:", faults_end)
    print("absorbing BC:", absorb_end)
            
    # add group 103 and 105 for test
    gmsh.model.addPhysicalGroup(2, absorb_end,105) 
    gmsh.model.addPhysicalGroup(2, faults_end,103) 
    gmsh.model.addPhysicalGroup(2, free_surf,101) 

    # Optionally, visualize the mesh
#     gmsh.fltk.run()
    gmsh.write(MeshFile+'.msh2') # type 2 Gmsh file  

    # Finalize Gmsh
    gmsh.finalize()

# generate three fault segemnts file



def generate_mesh_three_fault(MeshFile,stlFile1,stlFile2,stlFile3,domainBox,ncfile):
    
    lc = 20.0e3
    N = 200
    
    PI = np.pi/180


    gmsh.initialize()
    gmsh.model.add("HikWhaipapaFault")
    gmsh.option.setNumber("Mesh.Optimize", 1)  # Enables mesh optimization (default is off)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)  # Use Netgen mesh optimizer

    
    # Adjust the overall mesh size by modifying the characteristic length factor.
    gmsh.option.setNumber("Mesh.AngleToleranceFacetOverlap", 0.2)

    # Set the minimum volume threshold so that even very small tetrahedra are kept.
    # gmsh.option.setNumber("Mesh.NetgenMinVolume", 1e-12)

    r = 2000.0
    resolution = 1000.0

    # mesh_surf_size = 500.0
    
    # # Step 1: Create the box geometry
    # box_length = L
    # box_width = B
    # box_height = H + 3000.0
    
    # Step 1: Load the STL surface fi le
    gmsh.merge(stlFile1)
    gmsh.merge(stlFile2)
    gmsh.merge(stlFile3)

    # Step 3: Create a compound entity containing both the box and the loaded STL surface
    # Step 4: Merge the compound entity
    gmsh.model.geo.synchronize()

    s = gmsh.model.getEntities(2)
    print('s:',s)

    # nodeTags, nodeCoord, nodeParams = gmsh.model.mesh.getNodes(2,2)
    # vxyz = nodeCoord.reshape((-1, 3)) 
    
    # for i in range(0, len(vxyz[:,1])): 
    #     gmsh.model.mesh.setNode(nodeTags[i], [nodeCoord[i*3],nodeCoord[i*3+1],vxyz[i,2]-0.0],[1,1])


    gmsh.merge(domainBox)
    # gmsh.merge() # larger domain

    # asign physical group in gmsh
    gmsh.model.geo.synchronize()

    faces = gmsh.model.getEntities(dim=2)
    print('all faces:',faces)

    sl1 = gmsh.model.geo.addSurfaceLoop([1,2,3,4,5])
    v1 = gmsh.model.geo.addVolume([sl1])
    gmsh.model.geo.synchronize()
    
    gmsh.model.addPhysicalGroup(3, [v1], 1 )

    
    ## mesh size
    # distance = gmsh.model.mesh.field.add("Distance")
    # gmsh.model.mesh.field.setNumbers(distance, "FacesList", [1,2,3])
    
    # threshold = gmsh.model.mesh.field.add("Threshold")
    # gmsh.model.mesh.field.setNumber(threshold, "IField", distance)
    # gmsh.model.mesh.field.setNumber(threshold, "LcMin", 1.0*resolution)
    # gmsh.model.mesh.field.setNumber(threshold, "LcMax", 10*resolution)
    # gmsh.model.mesh.field.setNumber(threshold, "DistMin", r)
    # gmsh.model.mesh.field.setNumber(threshold, "DistMax", 2*r)
   
    # minimum = gmsh.model.mesh.field.add("Min")
    # gmsh.model.mesh.field.setNumbers(
    #     minimum, "FieldsList", [threshold])
    # gmsh.model.mesh.field.setAsBackgroundMesh(minimum)
    ## Generate the mesh
    
    # gmsh.model.mesh.setSize([(2,2)], mesh_surf_size)
    
    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(3)

    # asign physical group in gmsh
    faces = gmsh.model.getEntities(dim=2)
    print('all faces:',faces)

    absorb_end = [5]
    faults_end = [1]
    free_surf = [4]

    nodeTags, nodeCoord, nodeParams = gmsh.model.mesh.getNodes(2,4)
    vxyz = nodeCoord.reshape((-1, 3)) 

    znew = locate_modify_topo(ncfile,vxyz[:,0],vxyz[:,1])
    
    for i in range(0, len(znew)): 
        if (znew[i]< vxyz[i,2]): # vxyz[i,2]
            gmsh.model.mesh.setNode(nodeTags[i], [nodeCoord[i*3],nodeCoord[i*3+1],vxyz[i,2]],[1,1])
        else:
            gmsh.model.mesh.setNode(nodeTags[i], [nodeCoord[i*3],nodeCoord[i*3+1],znew[i]],[1,1])

    print("dynamic BC:", faults_end)
    print("absorbing BC:", absorb_end)
            
    # add group 103 and 105 for test
    gmsh.model.addPhysicalGroup(2, absorb_end,105) 
    gmsh.model.addPhysicalGroup(2, faults_end,103) 
    gmsh.model.addPhysicalGroup(2, [2],66) 
    gmsh.model.addPhysicalGroup(2, [3],65) 
    # gmsh.model.addPhysicalGroup(2, [4],68) 

    gmsh.model.addPhysicalGroup(2, free_surf,101) 

    # Optionally, visualize the mesh
#     gmsh.fltk.run()
    gmsh.write(MeshFile+'.msh2') # type 2 Gmsh file  

    # Finalize Gmsh
    gmsh.finalize()



# generate a xdmf file for Paraview
import meshio

# Load the .msh file
def convert_mesh(MeshFile, version='msh2'):

    if version=='msh2':
        
        mesh = meshio.read(
        MeshFile+ '.msh2',  # string, os.PathLike, or a buffer/open file
        file_format="gmsh",  # optional if filename is a path; inferred from extension
        )
    else:
        mesh = meshio.read(
        MeshFile+ '.msh4',  # string, os.PathLike, or a buffer/open file
        file_format="gmsh",  # optional if filename is a path; inferred from extension
        )
    meshio.write(
    MeshFile+".vtk",  # str, os.PathLike, or buffer/ open file
    mesh,
    file_format="vtk",  # optional if first argument is a path; inferred from extension
    )
    


from netCDF4 import Dataset
import pyproj


def locate_modify_topo(ncfile,xcoord,ycoord):
    
    topo,lon,lat = load_topo_data(ncfile)
    samples = len(lon)

    min_latitude= lat.min()
    min_longitude= lon.min()
    
    max_latitude= lat.max()
    max_longitude= lon.max()

    x = np.linspace((min_longitude),(max_longitude),samples)
    y = np.linspace((min_latitude),(max_latitude),samples)
    
    # print(x.shape,y.shape)
    # print(x,y)
    
    z_new = np.zeros(xcoord.shape)
    
    for k in range(z_new.size):
        indx = (np.abs(x-xcoord[k])).argmin()
        indy = (np.abs(y-ycoord[k])).argmin()
        z_new[k] = topo[indy,indx]

    return z_new


def load_topo_data(ncfile):
    nc = Dataset(ncfile,'r')
    lon= nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    topo = nc.variables['z'][:,:]
    lons,lats=np.meshgrid(lon,lat)

    samples = len(lon)
    
    min_latitude= lat.min()
    min_longitude= lon.min()
    
    max_latitude= lat.max()
    max_longitude= lon.max()
    
    mark_x = np.where( topo == -32768 )[0]
    mark_y = np.where( topo == -32768 )[1]
    # for x, y in zip(mark_x, mark_y) :
    #     slice = topo[max(0, x-1):x+1, max(0,y-1):y+1] # assuming a 5x5 square
    #     topo[x,y] = np.mean([i for i in slice.flatten() if i > 0])  # threshold is 0
    
    x_lon = np.linspace((min_longitude),(max_longitude),samples)
    y_lat = np.linspace((min_latitude),(max_latitude),samples)
    
    # UTM projection
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    # myproj = pyproj.Proj(proj='tmerc',lon_0=168, datum='WGS84')
    # myproj = pyproj.Proj(proj='utm',zone='59', datum='WGS84') # nzmg  
    myproj = pyproj.Proj(init='epsg:2193', datum='WGS84')
    
    xyz_map = pyproj.transform(lla, myproj, x_lon,y_lat,np.zeros(len(x_lon)), radians=False)
    x = xyz_map[0]
    y = xyz_map[1]
    # print(x,y)
    
    # Wellington city Epicenter coordinate
    lat_sou = -41.3
    lon_sou = 174.75 
    xyz_sou = pyproj.transform(lla, myproj, lon_sou,lat_sou, radians=False) # Epicenter in UTM domain
    print(xyz_sou)
    
    xmin = x[0] # Unit (m)
    xmax = x[-1]
    
    ymin = y[0]
    ymax = y[-1]

    return topo, x, y



# generate gmsh files with specific mesh size function
# import gmsh
# Helper function to return a node tag given two indices i and j:

# ext =5e3
# xmin,xmax,ymin,ymax =  xflt.min()-ext,xflt.max()+ext,yflt.min()-ext,yflt.max()+ext

def tag_fault(i,i0=100):
    '''tag fault surface'''
    return i+i0
        
def cal_mesh_lc(x0,y0,z0):
    ''' adaptive mesh based on distance to basin center'''
    xcc = (xmin+xmax)/2
    ycc = (ymin+ymax)/2
    zcc = -25e3
    lc = 5e3
    if z0 > -26000.0:
        if z0 > 0.0:
            lc =  500.0 # Densify the mesh for buildings
        else:
            if (np.sqrt((x0-xcc)**2 + (y0-ycc)**2 + (z0-zcc)**2)< 5e3):
                lc = 1000.0
            else:
                lc = 2000.0
        return lc
    else:
        return lc
        
def cal_mesh_size(entity_dim, entity_tag, x, y, z, lc):
    return cal_mesh_lc(x,y,z) 