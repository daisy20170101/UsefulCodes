# This is an example to load Topographic date from Etopo
# load netcdf reader.

from netCDF4 import Dataset
import pyproj
import numpy as np


def load_topo_data(ncfile):
    '''Load topographic grd file (e.g. Gebco_01s.grd) to memory'''
    
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
    for x, y in zip(mark_x, mark_y) :
        slice = topo[max(0, x-1):x+1, max(0,y-1):y+1] # assuming a 5x5 square
        topo[x,y] = np.mean([i for i in slice.flatten() if i > 0])  # threshold is 0
    
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



def load_topo_data_cartesian(ncfile):
    nc = Dataset(ncfile,'r')
    lon= nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    topo = nc.variables['z'][:,:]
    lons,lats=np.meshgrid(lon,lat)

    samples = len(lon)
    
    min_latitude= -41.4
    min_longitude= 174.7
    
    max_latitude= -41.197
    max_longitude= 174.903
    
    mark_x = np.where( topo == -32768 )[0]
    mark_y = np.where( topo == -32768 )[1]
    for x, y in zip(mark_x, mark_y) :
        slice = topo[max(0, x-1):x+1, max(0,y-1):y+1] # assuming a 5x5 square
        topo[x,y] = np.mean([i for i in slice.flatten() if i > 0])  # threshold is 0
    
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