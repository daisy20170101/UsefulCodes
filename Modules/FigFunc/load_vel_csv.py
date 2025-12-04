from scipy.interpolate import griddata
import numpy as np
import pandas as pd

def load_vel_csv(csvfile,Nx,Nz):
    '''load velocity data and interpolate to my domian;
    and interpolate series (unstructured) data from website onto structured grids
    Usage: mx,my,grd_z,grdvel,grdrho = load_vel_csv(csv_file, Nx(=Ny), Nz) 
    '''
    
    table1 = pd.read_csv(csvfile)
    print(table1.keys())

    grdx = table1['X']
    grdy = table1['Y']
    grdz = table1['Z']   
    # grdvel = table1['Vs_BLOCK']
    # grdrho = table1['DENSITY']

    grd_z = np.unique(grdz)
    grd_x = np.unique(grdx)
    grd_y = np.unique(grdy)
    print(len(grd_z),len(grd_x),len(grd_y))

    mx = np.linspace(grd_x[0],grd_x[-1],Nx)
    my = np.linspace(grd_y[0],grd_y[-1],Nx)

    grdx,grdy  = np.meshgrid(mx, my)
    # print(grdx,grdy)
    # grdz = np.linspace(grd_z[0],grd_z[-1],Nz)

    grdvel = np.zeros((Nz,Nx,Nx))
    grdrho = np.zeros((Nz,Nx,Nx))

    for k in range(Nz):
        subdata = table1[table1.Z == grd_z[k]]
        points = np.array([subdata['X'][:],subdata['Y'][:]])
        points = points.transpose()
        # print(points)
        
        values = subdata['Vs_BLOCK']
        grddata = griddata(points,values, (grdx, grdy), method='nearest')
        grdvel[k,:,:] = grddata

        values = subdata['DENSITY']
        grddata = griddata(points,values, (grdx, grdy), method='nearest')
        grdrho[k,:,:] = grddata
        
        

    return mx,my,grd_z,grdvel,grdrho