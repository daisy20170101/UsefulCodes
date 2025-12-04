import seissolxdmf
import pyproj
import matplotlib.tri as tri

def load_fault_data(xdmfFilename,variable):
                 
    sx = seissolxdmf.seissolxdmf(xdmfFilename)
    
    ndt = sx.ReadNdt()
    surfxyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    
    print(surfxyz.shape,connect.shape)
    
    # convert  # UTM projection
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    myproj = pyproj.Proj(init='epsg:2193', datum='WGS84')
    
    surf = pyproj.transform(myproj,lla,surfxyz[:,0],surfxyz[:,1], surfxyz[:,1], radians=False)
    
    # print(surf.shape)
    ############# load GMPEs data  ##############
    # triang = tri.Triangulation(surf[0],surf[1],connect) # in longitude and latititude
    triang = tri.Triangulation(surfxyz[:,0],surfxyz[:,1],connect)
    
    ##%%
    srs = sx.ReadData(variable)
    # srd = sx.ReadData('SRd')
    # slp = sx.ReadData('ASl')
    
    return  srs, triang