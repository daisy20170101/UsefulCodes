from palettable.cmocean import sequential as cmapb
from matplotlib.colors import ListedColormap

Cmap2 = ListedColormap(cmapb.Deep_14.mpl_colors)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

import geopandas as gpd


def plot_gme_hik(pgv,sa1,sa3,triang,coast_table, hik_trace, faultTable,modelname='test1', imd='pgv'):

    epi300 = (1740500,5385000)

    figsize_single = (8,6)

    fig,ax0 = plt.subplots(nrows=1,ncols=1,figsize=(figsize_single))

    if imd=='sa1':
        print('plotting SA1.0s')
        sc = ax0.tripcolor(triang,sa1[:]/9.8,cmap=Cmap2, shading='flat',vmin=0,vmax=1.0)
        cl = fig.colorbar(sc,ax=ax0)
        cl.set_label('SA1.0s (g)')

    elif imd=='sa3':
        print('plotting SA3.0s')

        sc = ax0.tripcolor(triang,sa3[:]/9.8,cmap=Cmap2, shading='flat',vmin=0,vmax=0.5)
        cl = fig.colorbar(sc,ax=ax0)
        cl.set_label('SA3.0s (g)')

    else:
        print('plotting PGV')

        sc = ax0.tripcolor(triang,pgv[:],cmap=Cmap2, shading='flat',vmin=0,vmax=1.5)
        cl = fig.colorbar(sc,ax=ax0)
        cl.set_label('PGV (m/s)')
            
    ax0.plot(coast_table['Points:0'],coast_table['Points:1'],'.k',markersize=1)
    ax0.plot(hik_trace['Points:0'],hik_trace['Points:1'],'.',c='blue',alpha=0.7,label='Hikurangi slab',markersize=2.5)

    for ik, ids in enumerate(faultTable['Fault_ID']):
        df = faultTable[faultTable['Fault_ID']==ids]
        fxyz  = df["geometry"].get_coordinates()
        ax0.plot(fxyz['x'],fxyz['y'],'-',c='tomato',alpha=0.7)
        
    ax0.plot(fxyz['x'],fxyz['y'],'-',c='tomato',alpha=0.7, label='crustal fault')
        
    # ax0.set_yticks(lat_ticks)
    ax0.set(xlim=( epi300[0]-130e3, epi300[0]+180e3),ylim=(epi300[1]-110e3,epi300[1]+180e3))

    ax0.set_xlabel('x (m)')
    ax0.set_ylabel('y (m)')
    
    ax0.legend(loc=2)

    outname = imd + '-' + modelname + '.png'
    plt.savefig(outname,dpi=300,transparent=False)
    plt.close()



def plot_gme_crust(pgv,sa1,sa3,triang,coast_table, faultTable,modelname='test1',  imd='pgv'):

    epi300 = (1740500,5385000)

    figsize_single = (8,6)

    fig,ax0 = plt.subplots(nrows=1,ncols=1,figsize=(figsize_single))

    if imd=='sa1':
        print('plotting SA1.0s')
        sc = ax0.tripcolor(triang,sa1[:]/9.8,cmap=Cmap2, shading='flat',vmin=0,vmax=1.)
        cl = fig.colorbar(sc,ax=ax0)
        cl.set_label('SA1.0s (g)')

    elif imd=='sa3':
        print('plotting SA3.0s')

        sc = ax0.tripcolor(triang,sa3[:]/9.8,cmap=Cmap2, shading='flat',vmin=0,vmax=0.5)
        cl = fig.colorbar(sc,ax=ax0)
        cl.set_label('SA3.0s (g)')

    else:
        print('plotting PGV')

        sc = ax0.tripcolor(triang,pgv[:],cmap=Cmap2, shading='flat',vmin=0,vmax=1.5)
        cl = fig.colorbar(sc,ax=ax0)
        cl.set_label('PGV (m/s)')
            
    ax0.plot(coast_table['Points:0'],coast_table['Points:1'],'.k',markersize=1)
    # ax0.plot(hik_trace['Points:0'],hik_trace['Points:1'],'.',c='blue',alpha=0.7,label='Hikurangi slab',markersize=2.5)

    for ik, ids in enumerate(faultTable['Fault_ID']):
        df = faultTable[faultTable['Fault_ID']==ids]
        fxyz  = df["geometry"].get_coordinates()
        ax0.plot(fxyz['x'],fxyz['y'],'-',c='tomato',alpha=0.7)
        
    ax0.plot(fxyz['x'],fxyz['y'],'-',c='tomato',alpha=0.7, label='crustal fault')
        
    # ax0.set_yticks(lat_ticks)
    ax0.set(xlim=( epi300[0]-50e3, epi300[0]+70e3),ylim=(epi300[1]-50e3,epi300[1]+120e3))

    ax0.set_xlabel('x (m)')
    ax0.set_ylabel('y (m)')
    
    ax0.legend(loc=2)

    outname = imd + '-' + modelname + '.png'
    plt.savefig(outname,dpi=300,transparent=False)
    plt.close()