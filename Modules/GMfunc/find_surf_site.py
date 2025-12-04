# load seissol out put data and extract station
# load GM data and compare station-to-station
import seissolxdmf
import numpy as np
import scipy.spatial as spatial
import pyproj
import matplotlib.tri as tri


def find_surf_site(stafolder, stafile, xdmfFilename,modelname):
    
    staFile = open(stafolder  + stafile ,'r')
    
    sx = seissolxdmf.seissolxdmf(xdmfFilename)
    # Number of cells
    nElements = sx.ReadNElements()
    # Read time step
    dt = sx.ReadTimeStep()
    # Read number of time steps
    ndt = sx.ReadNdt()
    # load geometry array as a numpy array of shape ((nodes, 3))
    surfxyz= sx.ReadGeometry()
    # load connectivity array as a numpy array of shape ((nElements, 3 or 4))
    # The connectivity array gives for each cell a list of vertex ids.
    connect = sx.ReadConnect()
    # horizontal comp. and vertical
    # u = sx.ReadData('v1')
    # v = sx.ReadData('v2')
    # w = sx.ReadData('v3')
    # print('(ndt,nelemenet)= ', u.shape)
    # print('time interval=', dt)
    
    sitexyz = np.loadtxt(staFile)
    # coordinates convert if necessary
    # staxyz = pyproj.transform(myproj,lla,gmData.long,gmData.lat,radians=False)

    centers = (surfxyz[connect[:,0]] + surfxyz[connect[:,1]] + surfxyz[connect[:,2]])/3.
    Receiver = np.array(sitexyz)
    
    # Receiver = Receiver.transpose()
    
    # search for nearest points of stations
    tree = spatial.KDTree(centers)
    dist, ids = tree.query(Receiver)
    
    FidReceiversnew =  stafolder + modelname + '-'+ 'site_xyz.txt'
    fout = open(FidReceiversnew,'w')
    fout1 = open(stafolder + modelname + '-'+'site_number.txt','w')
    
    for k in range(sitexyz[:,0].size):
        #newrec = find_nearest_vector(centers, rec)
        newrec = centers[ids[k]]
        print(k,ids[k])
        fout.write("%f %f %f\n" %(newrec[0],newrec[1],newrec[2]))
        fout1.write("%d %f\n" %(ids[k],dist[k]))

        # data = np.array([u[:,ids[k]],v[:,ids[k]],w[:,ids[k]]])
        # np.savetxt(stafolder + modelname + '/sta'+ str(k)+'.txt',data.transpose())
        
    fout.close()
    fout1.close()



# load gpa and datas 
def get_indensity_table(stafolder, stafile, xdmfFilename,modelname,origin=(0,0)):

    '''create a table with  idensity data as a DataFrame'''
    
    staFile = open(stafolder  + stafile ,'r')
    
    sx = seissolxdmf.seissolxdmf(xdmfFilename)
    surfxyz= sx.ReadGeometry()
    connect = sx.ReadConnect()
    
    sitexyz = np.loadtxt(staFile)

    centers = (surfxyz[connect[:,0]] + surfxyz[connect[:,1]] + surfxyz[connect[:,2]])/3.

    Receiver = np.array(sitexyz)
        
    # search for nearest points of stations
    tree = spatial.KDTree(centers)
    dist, ids = tree.query(Receiver)

    print(f"Max index from KDTree query: {np.max(ids)}")
    print(f"Number of mesh centers: {len(centers)}")

    depi = np.sqrt( (centers[ids,0]-origin[0])**2 + (centers[ids,1]-origin[1])**2)

    print(f"Epicentral distance shape: {depi.shape}")

    ##%%
    pga = sx.ReadData('PGA')
    pgv = sx.ReadData('PGV')
    sa1 = sx.ReadData('SA01.000s')
    sa3 = sx.ReadData('SA03.000s')
    # sa0_3 = sx.ReadData('SA00.300s')

    print(f"Ground motion data shapes - PGA: {pga.shape}, PGV: {pgv.shape}, SA1.0: {sa1.shape}, SA3.0: {sa3.shape}")

    # Flatten arrays if they're 2D (e.g., shape (1, N) -> (N,))
    if pga.ndim > 1:
        pga = pga.flatten()
        pgv = pgv.flatten()
        sa1 = sa1.flatten()
        sa3 = sa3.flatten()
        # sa0_3 = sa0_3.flatten()
        print(f"Flattened to 1D - PGA: {pga.shape}, PGV: {pgv.shape}, SA1.0: {sa1.shape}, SA3.0: {sa3.shape}")

    # Check if indices are valid for the ground motion arrays
    if np.max(ids) >= len(pga):
        raise IndexError(f"Index {np.max(ids)} is out of bounds for ground motion arrays with size {len(pga)}. "
                        f"The mesh has {len(centers)} elements but ground motion data only has {len(pga)} elements. "
                        f"Ensure the XDMF file contains per-element data matching the mesh size.")

    # Extract data at receiver locations and ensure 1D arrays
    pga_sites = np.atleast_1d(pga[ids]).flatten()
    pgv_sites = np.atleast_1d(pgv[ids]).flatten()
    sa1_sites = np.atleast_1d(sa1[ids]).flatten()
    sa3_sites = np.atleast_1d(sa3[ids]).flatten()
    # sa0_3_sites = np.atleast_1d(sa0_3[ids]).flatten()
    depi_1d = np.atleast_1d(depi).flatten()

    # unit from m/s/s to %g 
    dataTable = {
        "distance (km)": depi_1d/1e3,
        "PGA (g)": pga_sites/9.8*100,
        "PGV (cm/s)": pgv_sites*100,
        "SA 3.0": sa3_sites/9.8*100,
        "SA 1.0": sa1_sites/9.8*100,
        # "SA 0.3": sa0_3_sites/9.8*100,
    }
    
    return dataTable


def get_intensity_table_vs30_baz(
    xdmfFilename: str,
    siteTable,
    source_coords,
    origin=(0, 0),
    x_col: str = 'x',
    y_col: str = 'y',
    z_col: str = 'z',
    vs30_col: str = 'Vs30',
    station_col: str = 'station'
):
    """
    Extract ground motion intensities at station locations with Vs30 and geometric metrics.

    Extracts PGA, PGV, SA values at specific station locations from XDMF output,
    computes distance metrics (r_rup, r_jb), azimuth angles, and includes Vs30
    from the site table.

    Parameters
    ----------
    xdmfFilename : str
        Path to the SeisSol XDMF ground motion output file
    siteTable : pd.DataFrame
        Station information table with columns:
        - x, y, z: station coordinates in meters (NZTM2000/EPSG:2193)
        - Vs30: time-averaged shear wave velocity (m/s)
        - station: station code/name (optional)
    source_coords : array-like
        Fault source coordinates as (x, y, z) in METERS
        - If 2D array (N x 3): multiple points defining fault surface
        - If 1D array (3,): single point (hypocenter)
    origin : tuple of float, optional
        (x, y) coordinates of epicenter in METERS (default: (0,0))
    x_col : str, optional
        Name of x-coordinate column in siteTable (default: 'x')
    y_col : str, optional
        Name of y-coordinate column in siteTable (default: 'y')
    z_col : str, optional
        Name of z-coordinate column in siteTable (default: 'z')
    vs30_col : str, optional
        Name of Vs30 column in siteTable (default: 'Vs30')
    station_col : str, optional
        Name of station code column (default: 'station')

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - station: station code (if available in siteTable)
        - x (m), y (m), z (m): station coordinates
        - lon, lat: geographic coordinates (degrees)
        - PGA (%g): Peak Ground Acceleration
        - PGV (cm/s): Peak Ground Velocity
        - SA_1.0 (%g): Spectral Acceleration at 1.0s
        - SA_3.0 (%g): Spectral Acceleration at 3.0s
        - r_epi (km): epicentral distance
        - r_rup (km): rupture distance
        - r_jb (km): Joyner-Boore distance
        - azimuth (deg): azimuth from origin to station
        - back_azimuth (deg): back-azimuth from station to origin
        - Vs30 (m/s): time-averaged shear wave velocity

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>>
    >>> # Create site table
    >>> siteTable = pd.DataFrame({
    ...     'station': ['WEL', 'SNZO', 'KNZ'],
    ...     'x': [1748616, 1750000, 1755000],
    ...     'y': [5427609, 5430000, 5435000],
    ...     'z': [0, 0, 0],
    ...     'Vs30': [300, 450, 600]
    ... })
    >>>
    >>> # Define source
    >>> hypocenter = np.array([1750000, 5428000, -10000])
    >>>
    >>> # Extract intensities
    >>> df = get_intensity_table_vs30_baz(
    ...     'surface-PGV.xdmf',
    ...     siteTable,
    ...     hypocenter,
    ...     origin=(1750000, 5428000)
    ... )
    >>>
    >>> print(df[['station', 'PGA (%g)', 'r_rup (km)', 'Vs30 (m/s)']].head())

    Notes
    -----
    - Nearest neighbor search finds closest mesh element to each station
    - All distance calculations done in kilometers
    - Azimuth convention: 0°=North, 90°=East, 180°=South, 270°=West
    - Requires siteTable with station coordinates and Vs30 values
    """
    import pandas as pd
    from obspy.geodetics import gps2dist_azimuth

    # Read XDMF data
    sx = seissolxdmf.seissolxdmf(xdmfFilename)
    surfxyz_m = sx.ReadGeometry()       # Vertex coordinates in meters
    connect = sx.ReadConnect()          # Element connectivity

    # Convert to km for distance calculations
    surfxyz = surfxyz_m / 1e3

    # Compute element centers in km
    centers = (surfxyz[connect[:,0]] + surfxyz[connect[:,1]] + surfxyz[connect[:,2]]) / 3.0

    # Extract station coordinates from siteTable
    station_coords_m = siteTable[[x_col, y_col, z_col]].values  # in meters
    station_coords = station_coords_m / 1e3  # convert to km

    # Find nearest mesh element for each station
    tree = spatial.KDTree(centers)
    dist, ids = tree.query(station_coords)

    # Load ground motion data
    pga = sx.ReadData('PGA')              # m/s²
    pgv = sx.ReadData('PGV')              # m/s
    sa1 = sx.ReadData('SA01.000s')        # m/s²
    sa3 = sx.ReadData('SA03.000s')        # m/s²

    # Flatten if needed
    if pga.ndim > 1:
        pga = pga.flatten()
        pgv = pgv.flatten()
        sa1 = sa1.flatten()
        sa3 = sa3.flatten()

    # Check indices validity
    if np.max(ids) >= len(pga):
        raise IndexError(
            f"Index {np.max(ids)} is out of bounds for ground motion arrays with size {len(pga)}. "
            f"The mesh has {len(centers)} elements but ground motion data only has {len(pga)} elements."
        )

    # Extract ground motion at station locations
    pga_sites = pga[ids]
    pgv_sites = pgv[ids]
    sa1_sites = sa1[ids]
    sa3_sites = sa3[ids]

    # Convert origin to km
    origin_km = (origin[0] / 1e3, origin[1] / 1e3)

    # Compute epicentral distance (r_epi)
    r_epi = np.sqrt(
        (station_coords[:,0] - origin_km[0])**2 +
        (station_coords[:,1] - origin_km[1])**2
    )

    # Ensure source_coords is 2D array and convert to km
    source_coords = np.atleast_2d(source_coords)
    if source_coords.shape[1] != 3:
        raise ValueError("source_coords must have shape (N, 3) with columns [x, y, z]")
    source_coords_km = source_coords / 1e3

    # Compute rupture distance (r_rup) - shortest 3D distance to fault
    r_rup = np.zeros(len(station_coords))
    for i, sta_coord in enumerate(station_coords):
        distances = np.sqrt(np.sum((source_coords_km - sta_coord)**2, axis=1))
        r_rup[i] = np.min(distances)

    # Compute Joyner-Boore distance (r_jb) - horizontal distance to fault projection
    source_xy = source_coords_km[:, :2]
    station_xy = station_coords[:, :2]

    r_jb = np.zeros(len(station_coords))
    for i, sta_xy in enumerate(station_xy):
        distances = np.sqrt(np.sum((source_xy - sta_xy)**2, axis=1))
        r_jb[i] = np.min(distances)

    # Convert station coordinates to lon/lat
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    myproj = pyproj.Proj(init='epsg:2193', datum='WGS84')

    station_lon, station_lat, _ = pyproj.transform(
        myproj, lla,
        station_coords_m[:,0],
        station_coords_m[:,1],
        station_coords_m[:,2],
        radians=False
    )

    # Convert origin to lon/lat for azimuth calculation
    origin_lon, origin_lat, _ = pyproj.transform(
        myproj, lla,
        origin[0], origin[1], 0.0,
        radians=False
    )

    # Compute azimuth and back-azimuth
    azimuth = np.zeros(len(station_coords))
    back_azimuth = np.zeros(len(station_coords))

    for i in range(len(station_coords)):
        _, az, baz = gps2dist_azimuth(
            origin_lat, origin_lon,
            station_lat[i], station_lon[i]
        )
        azimuth[i] = az
        back_azimuth[i] = baz

    # Build result DataFrame
    result = {
        # 'x (m)': station_coords_m[:,0],
        # 'y (m)': station_coords_m[:,1],
        # 'z (m)': station_coords_m[:,2],
        # 'lon': station_lon,
        # 'lat': station_lat,
        'PGA': pga_sites / 9.8,
        'PGV': pgv_sites,
        'pSA_1.0': sa1_sites / 9.8,
        'pSA_3.0': sa3_sites / 9.8,
        'r_epi': r_epi,
        'r_rup': r_rup,
        'r_jb': r_jb,
        'az': azimuth,
        'baz': back_azimuth,
        'Vs30': siteTable[vs30_col].values
    }

    # Add station code if available
    if station_col in siteTable.columns:
        result = {'sta': siteTable[station_col].values, **result}

    df_result = pd.DataFrame(result)

    return df_result


    ## load data and plot 

# import seissolxdmf
# import pyproj
# import matplotlib.tri as tri


# coastfile ='./Geometry/CoastNorth.txt.npy'
# coast =  np.load(coastfile)
# print(coast)

def load_surf_gm(xdmfFilename,origin=(0,0)):
    """
    Load surface ground motion data from SeisSol XDMF output file.

    Parameters
    ----------
    xdmfFilename : str
        Path to the SeisSol XDMF ground motion output file
    origin : tuple of float, optional
        (x, y) coordinates of the epicenter for distance calculation (default: (0,0))

    Returns
    -------
    pga : ndarray
        Peak Ground Acceleration at each mesh element (m/s²)
    pgv : ndarray
        Peak Ground Velocity at each mesh element (m/s)
    sa1 : ndarray
        Spectral Acceleration at 1.0s period (m/s²)
    sa3 : ndarray
        Spectral Acceleration at 3.0s period (m/s²)
    sa0_3 : ndarray
        Spectral Acceleration at 0.3s period (m/s²)
    triang : matplotlib.tri.Triangulation
        Triangulation object for plotting in lon/lat coordinates
    depi : ndarray
        Epicentral distance for each mesh element (m)
    """

    # Initialize XDMF reader
    sx = seissolxdmf.seissolxdmf(xdmfFilename)

    # Read mesh metadata
    ndt = sx.ReadNdt()                # Number of time steps
    surfxyz = sx.ReadGeometry()       # Vertex coordinates (N_vertices x 3)
    connect = sx.ReadConnect()        # Element connectivity (N_elements x 3)

    print(surfxyz.shape,connect.shape)

    # Convert coordinates from UTM (EPSG:2193, New Zealand Transverse Mercator) to lon/lat
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    myproj = pyproj.Proj(init='epsg:2193', datum='WGS84')

    # Transform from NZTM2000 to WGS84 lon/lat
    surf = pyproj.transform(myproj,lla,surfxyz[:,0],surfxyz[:,1], surfxyz[:,1], radians=False)

    # Compute element centers (centroids of triangular elements)
    centers = (surfxyz[connect[:,0]] + surfxyz[connect[:,1]] + surfxyz[connect[:,2]])/3.

    # Compute epicentral distance (horizontal distance from origin to each element center)
    depi = np.sqrt( (centers[:,0]-origin[0])**2 + (centers[:,1]-origin[1])**2)

    # Create triangulation for plotting in geographic coordinates (lon/lat)
    triang = tri.Triangulation(surf[0],surf[1],connect)
    # Alternative: plot in Cartesian coordinates
    # triang = tri.Triangulation(surfxyz[:,0],surfxyz[:,1],connect)

    # Load ground motion intensity measures from XDMF file
    pga = sx.ReadData('PGA')              # Peak Ground Acceleration
    pgv = sx.ReadData('PGV')              # Peak Ground Velocity
    sa1 = sx.ReadData('SA01.000s')        # Spectral Acceleration at 1.0s
    sa3 = sx.ReadData('SA03.000s')        # Spectral Acceleration at 3.0s
    # sa0_3 = sx.ReadData('SA00.300s')      # Spectral Acceleration at 0.3s

    return pga,pgv,sa1,sa3,triang,depi


import seaborn as sbn
import matplotlib.pyplot as plt


def load_surf_table(xdmfFilename, source_coords, origin=(0,0)):
    """
    Load surface ground motion data and compute distance/azimuth metrics.

    Similar to load_surf_gm but returns a comprehensive pandas DataFrame with
    ground motion intensities, distance metrics (r_rup, r_jb), and azimuth angles.

    Parameters
    ----------
    xdmfFilename : str
        Path to the SeisSol XDMF ground motion output file
    source_coords : array-like
        Fault source coordinates as (x, y, z) in METERS (NZTM2000/EPSG:2193)
        - If 2D array (N x 3): multiple points defining fault surface
        - If 1D array (3,): single point (hypocenter)
        Used for computing r_rup (rupture distance)
        Note: Input in meters, will be converted to km internally
    origin : tuple of float, optional
        (x, y) coordinates of epicenter in METERS (default: (0,0))
        Used for r_epi and r_jb calculations
        Note: Input in meters, will be converted to km internally

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - x (km), y (km), z (km): element center coordinates
        - lon, lat: geographic coordinates (degrees)
        - PGA (%g): Peak Ground Acceleration
        - PGV (cm/s): Peak Ground Velocity
        - SA_1.0 (%g): Spectral Acceleration at 1.0s
        - SA_3.0 (%g): Spectral Acceleration at 3.0s
        - r_epi (km): epicentral distance - horizontal distance from origin
        - r_rup (km): rupture distance - shortest 3D distance to fault surface
        - r_jb (km): Joyner-Boore distance - horizontal distance to fault projection
        - azimuth (deg): azimuth from origin to site (0=North, clockwise)
        - back_azimuth (deg): azimuth from site to origin

    Notes
    -----
    Distance metrics:
    - r_epi: sqrt((x-x0)^2 + (y-y0)^2) - epicentral distance
    - r_rup: min distance from site to any point on fault surface
    - r_jb: horizontal distance from site to surface projection of fault
      (set z=0 for all fault points, compute min horizontal distance)

    Azimuth convention:
    - 0 degrees = North
    - 90 degrees = East
    - 180 degrees = South
    - 270 degrees = West

    Examples
    --------
    >>> # Single hypocenter
    >>> source = np.array([0, 0, -10000])  # x, y, z in meters
    >>> df = load_surf_table('output.xdmf', source, origin=(0, 0))
    >>>
    >>> # Fault surface (multiple points)
    >>> fault_points = np.array([
    ...     [0, 0, -5000],
    ...     [10000, 0, -5000],
    ...     [10000, 0, -15000],
    ...     [0, 0, -15000]
    ... ])
    >>> df = load_surf_table('output.xdmf', fault_points, origin=(5000, 0))
    >>>
    >>> # Analyze results
    >>> print(df[['PGA', 'PGV', 'r_rup', 'r_jb']].describe())
    """
    import pandas as pd
    from obspy.geodetics import gps2dist_azimuth

    # Initialize XDMF reader
    sx = seissolxdmf.seissolxdmf(xdmfFilename)

    # Read mesh metadata
    ndt = sx.ReadNdt()
    surfxyz_m = sx.ReadGeometry()       # Vertex coordinates (N_vertices x 3) in meters
    connect = sx.ReadConnect()          # Element connectivity (N_elements x 3)

    # Convert coordinates from meters to kilometers for distance calculations
    surfxyz = surfxyz_m / 1e3  # Convert m to km

    # Convert coordinates from UTM (EPSG:2193, New Zealand) to lon/lat
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    myproj = pyproj.Proj(init='epsg:2193', datum='WGS84')

    # Compute element centers (centroids of triangular elements) in km
    centers = (surfxyz[connect[:,0]] + surfxyz[connect[:,1]] + surfxyz[connect[:,2]]) / 3.0

    # Transform element centers to lon/lat (convert back to meters for projection)
    centers_lon, centers_lat, _ = pyproj.transform(myproj, lla,
                                                   centers[:,0] * 1e3, centers[:,1] * 1e3, centers[:,2] * 1e3,
                                                   radians=False)

    # Load ground motion intensity measures
    pga = sx.ReadData('PGA')              # Peak Ground Acceleration (m/s²)
    pgv = sx.ReadData('PGV')              # Peak Ground Velocity (m/s)
    sa1 = sx.ReadData('SA01.000s')        # Spectral Acceleration at 1.0s (m/s²)
    sa3 = sx.ReadData('SA03.000s')        # Spectral Acceleration at 3.0s (m/s²)

    # Flatten if needed (handle 2D arrays from XDMF)
    if pga.ndim > 1:
        pga = pga.flatten()
        pgv = pgv.flatten()
        sa1 = sa1.flatten()
        sa3 = sa3.flatten()

    # Ensure source_coords is 2D array (N x 3) and convert to km
    source_coords = np.atleast_2d(source_coords)
    if source_coords.shape[1] != 3:
        raise ValueError("source_coords must have shape (N, 3) with columns [x, y, z]")
    source_coords_km = source_coords / 1e3  # Convert m to km

    # Convert origin to km
    origin_km = (origin[0] / 1e3, origin[1] / 1e3)

    # Compute epicentral distance (r_epi) - horizontal distance from origin (in km)
    r_epi = np.sqrt((centers[:,0] - origin_km[0])**2 + (centers[:,1] - origin_km[1])**2)

    # Compute rupture distance (r_rup) - shortest 3D distance to fault surface (in km)
    # For each element center, find minimum distance to any fault point
    r_rup = np.zeros(len(centers))
    for i, center in enumerate(centers):
        distances = np.sqrt(np.sum((source_coords_km - center)**2, axis=1))
        r_rup[i] = np.min(distances)

    # Compute Joyner-Boore distance (r_jb) - horizontal distance to fault projection (in km)
    # Project fault points to surface (z=0) and compute horizontal distance
    source_xy = source_coords_km[:, :2]  # Take only x, y coordinates (in km)
    centers_xy = centers[:, :2]

    r_jb = np.zeros(len(centers))
    for i, center_xy in enumerate(centers_xy):
        distances = np.sqrt(np.sum((source_xy - center_xy)**2, axis=1))
        r_jb[i] = np.min(distances)

    # Compute azimuth and back-azimuth
    # Convert origin to lon/lat for azimuth calculation (use meters for projection)
    origin_lon, origin_lat, _ = pyproj.transform(myproj, lla,
                                                 origin[0], origin[1], 0.0,
                                                 radians=False)

    azimuth = np.zeros(len(centers))
    back_azimuth = np.zeros(len(centers))

    for i in range(len(centers)):
        # gps2dist_azimuth(lat1, lon1, lat2, lon2)
        # Returns: (distance_m, azimuth, back_azimuth)
        _, az, baz = gps2dist_azimuth(origin_lat, origin_lon,
                                      centers_lat[i], centers_lon[i])
        azimuth[i] = az
        back_azimuth[i] = baz

    # Create DataFrame with all results
    df = pd.DataFrame({
        # 'x (km)': centers[:,0],                # Already in km
        # 'y (km)': centers[:,1],                # Already in km
        # 'z (km)': centers[:,2],                # Already in km
        # 'lon': centers_lon,
        # 'lat': centers_lat,
        'PGA': pga / 9.8,          # Convert m/s² to g
        'PGV': pgv ,               # in m/s
        'pSA_1.0': sa1 / 9.8,       # Convert m/s² to g
        'pSA_3.0': sa3 / 9.8,       # Convert m/s² to g
        'r_epi': r_epi,                   # Already in km
        'r_rup': r_rup,                   # Already in km
        'r_jb': r_jb,                     # Already in km
        'az': azimuth,
        'baz': back_azimuth
    })

    return df


def plot_gm_curve(df_all, df, obvTable, modelname='20240805', gm_flag='pgv'):
    """
    Plot ground motion intensity versus epicentral distance curves.

    Compares simulated ground motion data (SeisSol) with observed data (GeoNet).
    Creates log-log scatter plots showing the distance attenuation relationship.

    Parameters
    ----------
    df_all : pd.DataFrame
        All simulated ground motion data at surface mesh elements
        Must contain columns: "distance (km)" and intensity measure columns
    df : pd.DataFrame
        Simulated ground motion data at specific site locations
        Must contain columns: "distance (km)" and intensity measure columns
    obvTable : pd.DataFrame
        Observed ground motion data from GeoNet stations
        Must contain columns: "Epicentral Distance" and component-specific columns
    modelname : str, optional
        Model identifier for output filename (default: '20240805')
    gm_flag : str, optional
        Ground motion intensity measure to plot: 'pgv', 'sa3', or 'sa1' (default: 'pgv')
        - 'pgv': Peak Ground Velocity
        - 'sa3': Spectral Acceleration at 3.0s period
        - 'sa1': Spectral Acceleration at 1.0s period

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object
    ax : matplotlib.axes.Axes
        Axes object

    Notes
    -----
    - Plot is saved as '{modelname}-{gm_flag}-dist.png' at 300 dpi
    - Both axes use logarithmic scale (base 10)
    - Uses Nature/Science journal colorblind-safe color palette
    - Color scheme:
      * Background mesh elements: light gray (#E0E0E0)
      * Simulated sites: vermillion (#D55E00)
      * Observed vertical: blue (#0173B2)
      * Observed horizontal 1: bluish green (#029E73)
      * Observed horizontal 2: orange (#ECB01E)
    - Markers: circles (o) for simulated, triangles (^) for vertical,
      squares (s) for horizontal 1, diamonds (D) for horizontal 2
    - Journal style: clean spines, subtle grid, professional appearance
    """

    # Create figure with specified size
    fig,ax = plt.subplots(1,1,figsize=(6,4.5))

    if (gm_flag=='pgv'):
        # Plot Peak Ground Velocity (PGV) data
        # Nature/Science style colors: colorblind-safe and print-friendly

        # Simulated: all surface mesh elements (background cloud)
        sbn.scatterplot(data=df_all, x="distance (km)",y="PGV (cm/s)",marker = '.',size=2,
                       color='#E0E0E0', edgecolor='none', alpha=0.6, legend=False,label='model: surface')
        # Simulated: specific receiver sites (highlighted) - Vermillion
        sbn.scatterplot(data=df, x="distance (km)",y="PGV (cm/s)",marker = 'o',
                       color='#D55E00', edgecolor='#7F3700', s=50, linewidth=0.5, legend=False,label='model: sites')

        # Observed: GeoNet station data (3 components) - Nature colorblind-safe palette
        # Vertical component - Blue
        sbn.scatterplot(data=obvTable,x="Epicentral Distance",y="PGV Vertical",marker='^',
                       color='#0173B2', edgecolor='#014A73', s=60, linewidth=0.5, label='Vertical')
        # Horizontal 1 - Bluish Green
        sbn.scatterplot(data=obvTable,x="Epicentral Distance",y="PGV Horizontal 1",marker='s',
                       color='#029E73', edgecolor='#015E45', s=60, linewidth=0.5, label='Horizontal 1')
        # Horizontal 2 - Orange
        sbn.scatterplot(data=obvTable,x="Epicentral Distance",y="PGV Horizontal 2",marker='D',
                       color='#ECB01E', edgecolor='#8D6911', s=60, linewidth=0.5, label='Horizontal 2')

        ax.grid(which='both',linestyle=':',alpha=0.3,color='gray')
        ax.legend(frameon=True, loc='best', fontsize=9)
        ax.set(xlim=(0.1,500),ylim=(0.0001,10))
        # Set axis limits for PGV plot
        

    elif (gm_flag=='sa3'):
        # Plot Spectral Acceleration at 3.0s period
        # Nature/Science style colors: colorblind-safe and print-friendly

        # Simulated: all surface mesh elements
        sbn.scatterplot(data=df_all, x="distance (km)",y="SA 3.0 (%g)",marker = '.',size=2,
                       color='#E0E0E0', edgecolor='none', alpha=0.6, legend=False,label='simulated surface')
        # Simulated: specific receiver sites - Vermillion
        sbn.scatterplot(data=df, x="distance (km)",y="SA 3.0",marker = 'o',
                       color='#D55E00', edgecolor='#7F3700', s=50, linewidth=0.5, legend=False,label='simulated sites')

        # Observed: GeoNet station data (3 components) - Nature colorblind-safe palette
        # Vertical component - Blue
        sbn.scatterplot(data=obvTable,x="Epicentral Distance",y="PSA 3.0 Vertical",marker='^',
                       color='#0173B2', edgecolor='#014A73', s=60, linewidth=0.5, label='Vertical')
        # Horizontal 1 - Bluish Green
        sbn.scatterplot(data=obvTable,x="Epicentral Distance",y="PSA 3.0 Horizontal 1",marker='s',
                       color='#029E73', edgecolor='#015E45', s=60, linewidth=0.5, label='Horizontal 1')
        # Horizontal 2 - Orange
        sbn.scatterplot(data=obvTable,x="Epicentral Distance",y="PSA 3.0 Horizontal 2",marker='D',
                       color='#ECB01E', edgecolor='#8D6911', s=60, linewidth=0.5, label='Horizontal 2')

        ax.grid(which='both',linestyle=':',alpha=0.3,color='gray')
        ax.legend(frameon=True, loc='best', fontsize=9)

        # Set distance range for SA 3.0s plot
        ax.set(xlim=(0.1,500))
        ax.set(xlim=(0.1,500),ylim=(0.00001,0.1))

    else:
        # Plot Spectral Acceleration at 1.0s period (default case)
        # Nature/Science style colors: colorblind-safe and print-friendly

        # Simulated: all surface mesh elements
        sbn.scatterplot(data=df_all, x="distance (km)",y="SA 1.0 (%g)",marker = '.',size=2,
                       color='#E0E0E0', edgecolor='none', alpha=0.6, legend=False,label='simulated surface')
        # Simulated: specific receiver sites - Vermillion
        sbn.scatterplot(data=df, x="distance (km)",y="SA 1.0",marker = 'o',
                       color='#D55E00', edgecolor='#7F3700', s=50, linewidth=0.5, legend=False,label='simulated sites')

        # Observed: GeoNet station data (3 components) - Nature colorblind-safe palette
        # Vertical component - Blue
        sbn.scatterplot(data=obvTable,x="Epicentral Distance",y="PSA 1.0 Vertical",marker='^',
                       color='#0173B2', edgecolor='#014A73', s=60, linewidth=0.5, label='Vertical')
        # Horizontal 1 - Bluish Green
        sbn.scatterplot(data=obvTable,x="Epicentral Distance",y="PSA 1.0 Horizontal 1",marker='s',
                       color='#029E73', edgecolor='#015E45', s=60, linewidth=0.5, label='Horizontal 1')
        # Horizontal 2 - Orange
        sbn.scatterplot(data=obvTable,x="Epicentral Distance",y="PSA 1.0 Horizontal 2",marker='D',
                       color='#ECB01E', edgecolor='#8D6911', s=60, linewidth=0.5, label='Horizontal 2')

        ax.grid(which='both',linestyle=':',alpha=0.3,color='gray')
        ax.legend(frameon=True, loc='best', fontsize=9)
        ax.set(xlim=(0.1,500),ylim=(0.0001,10))

    # Apply log-log scaling for both axes (distance attenuation relationship)
    ax.set_xscale('log',base=10)
    ax.set_yscale('log',base=10)
    
    # Apply Nature/Science style formatting
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.tick_params(width=0.8, direction='out', length=4)

    # Save figure to file with tight layout
    plt.tight_layout()
    plt.savefig(modelname + '-' + gm_flag+ "-dist.png",dpi=300, bbox_inches='tight')

    return fig, ax


