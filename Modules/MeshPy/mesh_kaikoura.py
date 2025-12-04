# generate GMsh model
# need to instal and load GMSH
import gmsh
# Helper function to return a node tag given two indices i and j:

import numpy as np
import pandas as pd
import geopandas as gpd

def tag_fault(i,N):
    '''tag fault surface'''
    return N**2 +i 
        
def cal_mesh_lc(x0,y0,z0):
    ''' adaptive mesh based on distance to basin center'''
    xcc =1774453.12
    ycc =5451928.11
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

def generate_mesh_fault(
                        faultTable: gpd.GeoDataFrame,
                        additionTable: pd.DataFrame,
                        MeshFile:str = 'output',
                        add_fault=True, 
                        ncfile: str = 'topo.nc',
                        resolution=2000.0,
                        r_coarsen=1000.0,
                        ymin=5.2e6,
                        ymax=5.3e6,
                        xmin=1.7e6,
                        xmax=1.74e6,
                        zmin=-100.0e3,
                        ):
    
    '''create mesh file with fault segments, iteratively over fault segments, with any dip angle  '''

    # add faults s
    fwidth = 20000.0 
    f_ids = faultTable["Fault_ID"].unique()
    N = 100
    PI = np.pi/180

    gmsh.initialize()

    # mesh resolution
    r = r_coarsen
    lc = resolution

    # topo, xrange,yrange = load_topo_data(ncfile)
    
    print(xmin,xmax,ymin,ymax)
    
    # geometry of the domain around the fault
    L, B, H = (xmax-xmin)/1.5, (ymax-ymin)/1.7, -zmin 
    xc,yc,zc = (xmax+xmin)/2.0,(ymax+ymin)/2.0, -zmin

    x0, y0, z0 = xc-L/2, yc-B/2, zmin


    v1 = gmsh.model.occ.addBox(x0,y0,z0, L, B, H)
    
    gmsh.model.occ.synchronize()  # must add synchronize before addPhysicalGroup
    surface_tags_box = gmsh.model.getEntities(dim=2)
    print('box surfaces:',surface_tags_box)

    for k,ids in enumerate(f_ids):
        # print(k,ids)
        df = faultTable[faultTable['Fault_ID']==ids]

        fxyz  = df["geometry"].get_coordinates()
        dip  = df['Dip_pref'].to_numpy()[0]
        strik  = df['SR_pref'].to_numpy()[0]

        if (ids==547 or ids==546):
            dip= 50.0


        # if (df['Dip_dir']=='N' or df['Dip_dir']=='NW'):
        #     strik = df['SR_pref'].to_numpy()[0] + 180.0

        for id in range(len(fxyz['x'])-1): 
            # print(id)
            point1 = gmsh.model.occ.add_point(np.asarray(fxyz['x'])[id], np.asarray(fxyz['y'])[id],0.0,lc)
            point2 = gmsh.model.occ.add_point(np.asarray(fxyz['x'])[id+1], np.asarray(fxyz['y'])[id+1],0.0,lc)
    
            line1 = gmsh.model.occ.add_line(point1, point2)

            # dip direction SW or SE
            if (ids==586 ):

                face1 = gmsh.model.occ.extrude([(1,line1)], np.cos(dip*PI)*fwidth * (np.sin((strik+90)*PI)), np.cos(dip*PI)*fwidth * (np.cos((strik+90)*PI)),-fwidth*np.sin(dip*PI))  

            else:
                face1 = gmsh.model.occ.extrude([(1,line1)], np.cos(dip*PI)*fwidth * np.sin((strik+270)*PI), np.cos(dip*PI)*fwidth *  np.cos((strik+270)*PI),-fwidth*np.sin(dip*PI))  

            # face1 = gmsh.model.occ.extrude([(1,line1)],-np.cos(dip*PI)*fwidth * np.sin((strik+90)*PI),-np.cos(dip*PI)*fwidth *  np.cos((strik+90)*PI),-fwidth*np.sin(dip*PI))  
            # face1 = gmsh.model.occ.extrude([(1,line1)], -np.cos(dip*PI)*fwidth * np.cos((90-strik)*PI),np.cos(dip*PI)*fwidth *  np.sin((90-strik)*PI),-fwidth*np.sin(dip*PI))  
            gmsh.model.occ.synchronize()  # must add synchronize before addPhysicalGroup

           # Get the list of all surface tags
    
    if add_fault:

        fwidth=30000.0
        point1=gmsh.model.occ.add_point(additionTable['x'][0],additionTable['y'][0],0.0,lc)
        point2=gmsh.model.occ.add_point(additionTable['x'][1],additionTable['y'][1],0.0,lc)
        line1 = gmsh.model.occ.add_line(point1, point2)

        strike = additionTable['strike'].to_numpy()[0]
        dip = additionTable['dip'].to_numpy()[0]
        face1 = gmsh.model.occ.extrude([(1,line1)], np.cos(dip*PI)*fwidth * (np.sin((strik+270)*PI)), np.cos(dip*PI)*fwidth * (np.cos((strik+270)*PI)),
                                       -fwidth*np.sin(dip*PI))  
        gmsh.model.occ.synchronize()


    surface_tags_all = gmsh.model.getEntities(dim=2)
    surface_tags = surface_tags_all[6::]
    
    obstacles = []
    # Create a new physical group for surfaces
    for surface in surface_tags:
        obstacles.append(surface[1])
    
    # gmsh.model.addPhysicalGroup(2, obstacles,103) 
    print("fault surface:", obstacles)

    # gmsh.model.setPhysicalName(2, 103, "fault surface")
    gmsh.model.occ.synchronize()

    # for surface in surface_tags:
    #     v11 = gmsh.model.getEntities(3)
    #     print(v11,surface)
    #     gmsh.model.occ.fragment([(3,v11)],[surface],removeObject=True, removeTool=True )
    out, _ = gmsh.model.occ.fragment([(3, v1)], [(2, i) for i in obstacles])

    fault_surf = out[1::]
    print('fault surface:',fault_surf)
    gmsh.model.occ.synchronize()

    # Get all 3D entities (volumes) after fragmentation
    all_surfaces = gmsh.model.getEntities(2)
    print('all_surface:',all_surfaces)

    # Synchronize after removal
    gmsh.model.occ.synchronize()
    
    # Now filter fault surfaces (those that are at depth)
    surfaces_to_remove = []
    fault_surfaces =[]
    
    for surface in fault_surf:
        try:
            com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
            z_center = com[2]
            
            # Keep fault surfaces (they should be at depth)
            if z_center > 5.0:  # Fault surfaces are typically deep
                surfaces_to_remove.append(surface[1])

            else:
                if (z_center >-20000.0 and z_center<0.0):
                    fault_surfaces.append(surface[1])

        except Exception as e:
            print(f"Warning: Could not evaluate surface {surface[1]}: {e}")
    
    print('faces to remove:',surfaces_to_remove)

     # Remove selected surfaces
    for surface in surfaces_to_remove:
        gmsh.model.occ.remove([(2,surface)], recursive=True)
        
    gmsh.model.occ.synchronize()

    all_surfaces = gmsh.model.getEntities(2)
    print('left surfaces:',all_surfaces)
    print('fault surfaces:',fault_surfaces)
    
    remain_surfaces=[]
    for surface in all_surfaces:
        remain_surfaces.append(surface[1])

    # Fault surface trace:
    gmsh.model.occ.synchronize()


    # Distance-based refinement will be applied in Stage 2 after initial mesh generation
    
    # Two-stage mesh generation approach
    print("Stage 1: Generating initial rough mesh...")

    # Stage 2: Apply distance-based refinement
    print("Stage 2: Applying distance-based refinement...")
    
    # Create distance field for fault surfaces (now with existing geometry)
    distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance, "FacesList",fault_surfaces )
    gmsh.model.mesh.field.setNumbers(distance, "Sampling", [100])  # Higher sampling for better accuracy
    
    # Fine mesh very close to fault (0-300m)
    threshold1 = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold1, "IField", distance)
    gmsh.model.mesh.field.setNumber(threshold1, "LcMin", resolution)  # 1000m
    gmsh.model.mesh.field.setNumber(threshold1, "LcMax", 30*resolution)  # 2000m
    gmsh.model.mesh.field.setNumber(threshold1, "DistMin", 1*r)
    gmsh.model.mesh.field.setNumber(threshold1, "DistMax", 10*r)  # 300m
    
    # Combine all thresholds using minimum
    minimum = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [threshold1])
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum)
    
    # Set comprehensive mesh quality constraints
    print("Setting mesh quality constraints...")
    
    # Use the quality constraint function
    # set_mesh_quality_constraints(target_quality=0.7, min_quality=0.2)
    
    # Validate constraints
    # validate_mesh_quality_constraints()
    
    # Generate refined mesh:
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.model.occ.synchronize()

    # label all elements on face
    # faces = gmsh.model.getBoundary([(3, out[0][1])],combined=True, oriented=False, recursive=False)
    faces = gmsh.model.getEntities(dim=2)
    
    print('all faces:',faces)

    
    absorb_end = []
    faults_end = []

    for surface in faces:

        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])

        if np.allclose(com,[(xmax+xmin)/2.0, (ymax+ymin)/2.0, 0.0]):
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], 101)  # free-surfac boundary label

            print('free surface:',surface[1])
            nodeTags, nodeCoord, nodeParams = gmsh.model.mesh.getNodes(2,surface[1])
            vxyz = nodeCoord.reshape((-1, 3)) 
            
        elif np.sqrt((com[0]-xc)**2 + (com[1]-yc)**2) > 0.8 * (L+B)/4.0 : 
            absorb_end.append(surface[1])

        elif np.allclose(com,[(xmax+xmin)/2.0, (ymax+ymin)/2.0, zmin]):
            absorb_end.append(surface[1])
            
        else:
            faults_end.append(surface[1])
   
    # znew = locate_modify_topo(ncfile,vxyz[:,0],vxyz[:,1])

    # for i in range(0, len(znew)): 
    #     if (znew[i]>nodeCoord[i*3+2]):
    #         gmsh.model.mesh.setNode(nodeTags[i], [nodeCoord[i*3],nodeCoord[i*3+1],znew[i]], [nodeParams[i*2],nodeParams[i*2+1]])

    gmsh.model.occ.synchronize()

    print("dynamic BC:", faults_end)
    print("absorbing BC:", absorb_end)

    start = 165
    for isf, faultid in enumerate(faults_end[1::]):
        phys_tag = start + isf
        gmsh.model.addPhysicalGroup(2, [faultid], phys_tag)
        # gmsh.model.setPhysicalName(2, phys_tag, f"surf_{phys_tag}")


    gmsh.model.addPhysicalGroup(2, [faults_end[0]], 103)

    # add group 103 and 105 for test
    gmsh.model.addPhysicalGroup(2, absorb_end,105) 
    # gmsh.model.addPhysicalGroup(2, faults_end,103) 
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(3, [1],1) 
    
    gmsh.option.setNumber("Mesh.MshFileVersion", 4)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
            
    # Write mesh data:
    gmsh.write(MeshFile + ".msh4")
     
    # Creates  graphical user interface
    # if 'close' not in sys.argv:
    #     gmsh.fltk.run()
     
    # It finalize the Gmsh API
    gmsh.finalize()

def locate_modify_topo(ncfile,xcoord,ycoord):
    
    nc = Dataset(ncfile,'r')
    lon= nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    topo = nc.variables['z'][:,:]

    samples = len(lon)
    
    min_latitude= lat.min()
    min_longitude= lon.min()
    
    max_latitude= lat.max()
    max_longitude= lon.max()
    
    mark_x = np.where( topo == -32768 )[0]
    mark_y = np.where( topo == -32768 )[1]

    x_lon = np.linspace((min_longitude),(max_longitude),samples)
    y_lat = np.linspace((min_latitude),(max_latitude),samples)
    
    # UTM projection
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    myproj = pyproj.Proj(init='epsg:2193', datum='WGS84')
    
    x,y = pyproj.transform(lla, myproj, x_lon,y_lat, radians=False)

    # print(x.shape,y.shape)
    # print(x,y)
    
    z_new = np.zeros(xcoord.shape)
    
    for k in range(z_new.size):
        indx = (np.abs(x-xcoord[k])).argmin()
        indy = (np.abs(y-ycoord[k])).argmin()
        z_new[k] = topo[indy,indx]

    return z_new

from netCDF4 import Dataset
import pyproj


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
    
