"""
Basin data loading utilities for Wellington Basin velocity model processing.

This module provides functions to load and process various types of CSV data files
including velocity models, topography, and elevation data.
"""

import numpy as np
import pandas as pd
from scipy.interpolate import griddata


def load_origin_csv(csvfile, yprofile):
    """
    Load and filter velocity model data from CSV for a specific Y profile.

    Parameters:
    -----------
    csvfile : str
        Path to CSV file containing velocity model data
    yprofile : float
        Y coordinate value to filter data around

    Returns:
    --------
    tuple
        (grd_y, filtered_table) - unique Y values and filtered table
    """
    table1 = pd.read_csv(csvfile)
    print(table1.keys())

    grdx = table1['X']
    grdy = table1['Y']
    grdz = table1['Z']
    grdvel = table1['Vs_BLOCK']
    grdrho = table1['DENSITY']

    grd_z = np.unique(grdz)
    grd_x = np.unique(grdx)
    grd_y = np.unique(grdy)

    return grd_y, table1[table1['Y'].apply(np.isclose, b=yprofile, atol=1.0)]


def load_topo_csv(csvfile, yprofile):
    """
    Load topographic data from CSV.

    Parameters:
    -----------
    csvfile : str
        Path to CSV file containing topographic data
    yprofile : float
        Y coordinate value (unused in current implementation)

    Returns:
    --------
    tuple
        (grd_y, table) - unique Y values and full table
    """
    table1 = pd.read_csv(csvfile)
    print(table1.keys())

    grdx = table1['X']
    grdy = table1['Y']
    grdz = table1['Z']

    grd_x = np.unique(grdx)
    grd_y = np.unique(grdy)

    return grd_y, table1


def load_dep_csv(csvfile, yprofile):
    """
    Load depth-related data from CSV.

    Parameters:
    -----------
    csvfile : str
        Path to CSV file containing depth data
    yprofile : float
        Y coordinate value (unused in current implementation)

    Returns:
    --------
    tuple
        (grd_y, table) - unique Y values and full table
    """
    table1 = pd.read_csv(csvfile)
    print(table1.keys())

    grdx = table1['X']
    grdy = table1['Y']
    grdz = table1['Z']

    grd_x = np.unique(grdx)
    grd_y = np.unique(grdy)
    grd_z = np.unique(grdz)

    return grd_y, table1


def load_elev_csv(csvfile, Nx):
    """
    Load elevation data and interpolate to structured grid.

    Parameters:
    -----------
    csvfile : str
        Path to CSV file containing elevation data with X, Y, Z columns
    Nx : int
        Number of grid points in each dimension (creates Nx x Nx grid)

    Returns:
    --------
    tuple
        (mx, my, grd_elev) - X coordinates, Y coordinates, and elevation grid

    Usage:
    ------
    mx, my, grd_elev = load_elev_csv('elevation.csv', 200)
    """
    table1 = pd.read_csv(csvfile)
    print(table1.keys())

    grdx = table1['X']
    grdy = table1['Y']

    grd_x = np.unique(grdx)
    grd_y = np.unique(grdy)
    print(len(grd_x), len(grd_y))

    mx = np.linspace(grd_x[0], grd_x[-1], Nx)
    my = np.linspace(grd_y[0], grd_y[-1], Nx)

    grdx_mesh, grdy_mesh = np.meshgrid(mx, my)
    grd_elev = np.zeros((Nx, Nx))

    points = np.array([table1['X'][:], table1['Y'][:]])
    points = points.transpose()

    values = table1['Z']

    print(points.shape, values.shape)

    grddata = griddata(points, values, (grdx_mesh, grdy_mesh), method='nearest')
    grd_elev[:, :] = grddata

    return mx, my, grd_elev