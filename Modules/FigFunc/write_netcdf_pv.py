"""
NetCDF file writing utilities for visualization and SeisSol input.

This module provides functions to write velocity model data and material properties
to NetCDF format for visualization in ParaView and input to SeisSol solver.
"""

import numpy as np
from netCDF4 import Dataset


def write_netcdf(fn, Nx, Ny, Nz, grdvel, grdx=None, grdy=None, grdz=None):
    """
    Write velocity data to NetCDF file for ParaView visualization.

    Parameters:
    -----------
    fn : str
        Output filename
    Nx, Ny, Nz : int
        Grid dimensions
    grdvel : ndarray
        3D velocity array
    grdx, grdy, grdz : ndarray, optional
        Coordinate arrays (must be provided in calling scope if None)
    """
    ncfile = Dataset(fn, 'w', format='NETCDF4')

    ncfile.createDimension('x', Nx)
    ncfile.createDimension('y', Nx)
    ncfile.createDimension('z', Nz)

    ncfile.createVariable('x','f8',('x',))
    ncfile.createVariable('y','f8',('y',))
    ncfile.createVariable('z','f8',('z',))

    if grdx is not None:
        ncfile.variables['x'][:] = np.unique(grdx)
        ncfile.variables['y'][:] = np.unique(grdy)
        ncfile.variables['z'][:] = np.unique(grdz)[0:Nz]*5

    vs_data = ncfile.createVariable('Vs','float64', ('z','y','x'))
    vs_data[:,:,:] = grdvel[:,:,:]
    print(vs_data.dtype)

    print(ncfile)
    ncfile.close()
    print('Dataset is closed!')


def write_netcdf_inp(fname, Nx, Nz, grdvel, grdrho, vpvs, grdx=None, grdy=None, grdz=None):
    """
    Write material properties (rho, mu, lambda) to NetCDF for SeisSol input.

    Parameters:
    -----------
    fname : str
        Output filename
    Nx, Nz : int
        Grid dimensions (assumes Ny = Nx)
    grdvel : ndarray
        3D S-wave velocity array
    grdrho : ndarray
        3D density array
    vpvs : float
        P-wave to S-wave velocity ratio
    grdx, grdy, grdz : ndarray, optional
        Coordinate arrays (must be provided in calling scope if None)
    """
    nu = 0.25

    # Output for ASAGI format
    fout = Dataset(fname, 'w', format='NETCDF4')

    fout.createDimension('x', Nx)
    fout.createDimension('y', Nx)
    fout.createDimension('z', Nz)
    fout.createVariable('x','f8',('x',))
    fout.createVariable('y','f8',('y',))
    fout.createVariable('z','f8',('z',))

    if grdx is not None:
        fout.variables['x'][:] = np.unique(grdx)
        fout.variables['y'][:] = np.unique(grdy)
        fout.variables['z'][:] = np.unique(grdz)

    material = np.dtype([('rho',np.float32),('mu',np.float32),('lambda',np.float32)])
    data = fout.createCompoundType(material,'material')
    vv = fout.createVariable('data',data,dimensions=('z','y','x'))

    for i in np.arange(Nx):
        for j in np.arange(Nx):
            rho_1d = grdrho[:,j,i]
            mu_1d = grdvel[:,j,i] **2 * rho_1d
            lam_1d = (grdvel[:,j,i]*vpvs)**2 * rho_1d - 2*mu_1d

            for k in range(Nz):
                vv[k,j,i] = (rho_1d[k], mu_1d[k], lam_1d[k])

    fout.close()