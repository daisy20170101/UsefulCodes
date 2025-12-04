import seissolxdmf
import numpy as np
import pyproj
import matplotlib.tri as tri

def load_surf_data(xdmfFilename,variable):
                 
    sx = seissolxdmf.seissolxdmf(xdmfFilename)
    
    ndt = sx.ReadNdt()
    surfxyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    
    print(surfxyz.shape,connect.shape)
    
    # convert  # UTM projection
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    myproj = pyproj.Proj(init='epsg:2193', datum='WGS84')
    
    surf = pyproj.transform(myproj,lla,surfxyz[:,0],surfxyz[:,1], surfxyz[:,1], radians=False)
    
    # centers = (surfxyz[connect[:,0]] + surfxyz[connect[:,1]] + surfxyz[connect[:,2]])/3.
    # depi = np.sqrt( (centers[:,0]-origin[0])**2 + (centers[:,1]-origin[1])**2)
    
    # print(surf.shape)
    ############# load GMPEs data  ##############
    triang = tri.Triangulation(surf[0],surf[1],connect) # in longitude and latititude
    
    ##%%
    v1 = sx.ReadData(variable)
    # v2 = sx.ReadData('u2')
    # v3 = sx.ReadData('u3')
 
    
    return triang,v1