from scipy.interpolate import griddata
import numpy as np
import pandas as pd

def load_bottom_csv(csvfile,Nx,):
    '''load elevation data and interpolate to my domian;
    and interpolate series (unstructured) data from website onto structured grids
    Usage: mx,my,grd_elev= load_elev_csv(csv_file, Nx(=Ny), Nz) 
    '''
    table1 = pd.read_csv(csvfile)
    # elev = pd.read_csv(topofile)

    print(table1.keys())

    grdx = table1['X']
    grdy = table1['Y']
    grdz = table1['Z']   

    grd_z = np.unique(grdz)
    grd_x = np.unique(grdx)
    grd_y = np.unique(grdy)
    
    nz = 129

    mx = np.linspace(grd_x[0],grd_x[-1],Nx)
    my = np.linspace(grd_y[0],grd_y[-1],Nx)

    print('dx,dy =:',str((grd_x[-1]-grd_x[0])/Nx),',',str((grd_y[-1]-grd_y[0])/Nx))

    grdx,grdy  = np.meshgrid(mx,my)  

    gvalues = np.zeros((nz,Nx,Nx))
    gvalues2 = np.zeros((nz,Nx,Nx))

    gfinal = np.zeros((Nx,Nx))

    # points2 = np.array([elev['X'][:],elev[Y][:]])
    # points2 = points2.transpose()
    # values2 = elev['Z']
    # gvalues2 = np.zeros((Nx,Nx))
    # gvalues[:,:]  = griddata(points2,values2, (grdx, grdy), method='nearest')

    for ik, iz in enumerate(grd_z[:]):

        subdata = table1[table1['Z'] == iz]

        points = np.array([subdata['X'][:],subdata['Y'][:]])
        points = points.transpose()
        
        values = subdata['Vs_B']
        gvalues[ik,:,:]  =  griddata(points,values, (grdx, grdy), method='nearest')

    for ix in range(Nx):
        for iy in range(Nx):
            darray2 = gvalues[:,iy,ix]
            # print(darray2)

            ilocate = np.abs(darray2-800.0).argmin()
            gfinal[iy,ix] = grd_z[ilocate]
            # print(darray1[ilocate])
            # print(darray)

    return mx,my,gfinal
