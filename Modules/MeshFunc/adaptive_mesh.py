"""
Adaptive mesh generation utilities for velocity-aware meshing with Gmsh.

This module provides functions for creating adaptive meshes based on 3D velocity
models, with mesh refinement that considers local velocity structure.
"""

import numpy as np


def create_mesh_tag_function(N):
    """
    Create a helper function to return node tags for mesh generation.

    Parameters:
    -----------
    N : int
        Number of grid points minus 1 (for N+1 x N+1 grid)

    Returns:
    --------
    function
        Tag function that takes (i, j) indices and returns node tag
    """
    def tag(i, j):
        return (N + 1) * i + j + 1
    return tag


def cal_mesh_lc(x0, y0, z0, vs_mesh, z_mesh, x_mesh, y_mesh,
                fre_high=2.0, num_elements_per_wavelength=2.0):
    """
    Calculate adaptive mesh characteristic length based on local velocity.

    Parameters:
    -----------
    x0, y0, z0 : float
        Point coordinates
    vs_mesh : ndarray
        3D velocity mesh
    z_mesh, x_mesh, y_mesh : ndarray
        Coordinate arrays for velocity mesh
    fre_high : float, default=2.0
        Target frequency for meshing
    num_elements_per_wavelength : float, default=2.0
        Number of elements per wavelength

    Returns:
    --------
    float
        Characteristic length for mesh at this point
    """
    lc = 5e3

    if z0 > -15000:
        indz = (np.abs(z0 - z_mesh)).argmin()
        indy = (np.abs(y0 - y_mesh)).argmin()
        indx = (np.abs(x0 - x_mesh)).argmin()

        # Check if point is outside velocity model domain
        if ((np.abs(x0 - x_mesh).min() > 5050.0) |
            (np.abs(y0 - y_mesh).min() > 5050.0) |
            (np.abs(z0 - z_mesh).min() > 5050.0)):
            lc = 3.0 * 1360.0 / fre_high / num_elements_per_wavelength
        else:
            lc = vs_mesh[indz, indy, indx] / fre_high / num_elements_per_wavelength

    return lc


def create_mesh_size_callback(vs_mesh, z_mesh, x_mesh, y_mesh,
                             fre_high=2.0, num_elements_per_wavelength=2.0):
    """
    Create a mesh size callback function for Gmsh.

    Parameters:
    -----------
    vs_mesh : ndarray
        3D velocity mesh
    z_mesh, x_mesh, y_mesh : ndarray
        Coordinate arrays for velocity mesh
    fre_high : float, default=2.0
        Target frequency for meshing
    num_elements_per_wavelength : float, default=2.0
        Number of elements per wavelength

    Returns:
    --------
    function
        Callback function for Gmsh mesh size calculation
    """
    def cal_mesh_size(entity_dim, entity_tag, x, y, z, lc):
        return cal_mesh_lc(x, y, z, vs_mesh, z_mesh, x_mesh, y_mesh,
                          fre_high, num_elements_per_wavelength)

    return cal_mesh_size


def setup_adaptive_mesh_parameters(grdvel, grdz, grdx, grdy, fre_high=2.0):
    """
    Setup parameters for adaptive mesh generation.

    Parameters:
    -----------
    grdvel : ndarray
        3D velocity array from load_vel_csv
    grdz, grdx, grdy : ndarray
        Coordinate arrays
    fre_high : float, default=2.0
        Target frequency

    Returns:
    --------
    dict
        Dictionary containing mesh parameters and callback function
    """
    vs_mesh = grdvel[:, :, :]
    z_mesh = grdz
    x_mesh = grdx
    y_mesh = grdy

    lc = 10e3
    num_elements_per_wavelength = 2.0

    mesh_size_callback = create_mesh_size_callback(
        vs_mesh, z_mesh, x_mesh, y_mesh, fre_high, num_elements_per_wavelength
    )

    return {
        'vs_mesh': vs_mesh,
        'z_mesh': z_mesh,
        'x_mesh': x_mesh,
        'y_mesh': y_mesh,
        'lc': lc,
        'fre_high': fre_high,
        'num_elements_per_wavelength': num_elements_per_wavelength,
        'mesh_size_callback': mesh_size_callback
    }