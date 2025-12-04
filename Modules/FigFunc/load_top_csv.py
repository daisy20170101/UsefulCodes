from scipy.interpolate import griddata
import numpy as np
import pandas as pd

def load_elev_csv(csvfile,Nx):
    '''load elevation data and interpolate to my domian;
    and interpolate series (unstructured) data from website onto structured grids
    Usage: mx,my,grd_elev= load_elev_csv(csv_file, Nx(=Ny), Nz) 
    '''
    
    table1 = pd.read_csv(csvfile)
    print(table1.keys())

    grdx = table1['X']
    grdy = table1['Y']
    grdz = table1['Z']   

    grd_x = np.unique(grdx)
    grd_y = np.unique(grdy)
    print(len(grd_x),len(grd_y))

    mx = np.linspace(grd_x[0],grd_x[-1],Nx)
    my = np.linspace(grd_y[0],grd_y[-1],Nx)

    grdx,grdy  = np.meshgrid(mx, my)
    grd_elev = np.zeros((Nx,Nx))

    points = np.array(table1['X'],table1['Y'])
    points = points.transpose()
    # print(points)
    
    values = table1['Z']
    grddata = griddata(points,values, (grdx, grdy), method='nearest')
    grd_elev[:,:] = grddata
        
    return mx,my,grd_elev