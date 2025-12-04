from palettable.cmocean import sequential as cmapb
from matplotlib.colors import ListedColormap

Cmap2 = ListedColormap(cmapb.Deep_14.mpl_colors)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import geopandas as gpd

import matplotlib
matplotlib.style.use('seaborn-v0_8-talk') # Lovely plotting style



def plot_gme_curve(modelname, gmTable,sa1, rjb, imd='sa1'):
    
    '''plot gme curve with GMPE calculation
    '''

    plt.figure(figsize=(6,4.5))

    # ddep = 9.0
    # Rhyp = np.sqrt((depi/1e3)**2+ddep**2)

    if imd=='sa1':
        plt.plot(rjb/1e3,sa1/9.8,'.',c='royalblue',label='SA1.0s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**gmTable['SA1'],'-',c='tomato',label='Atk22 mean;Vs30=300 m/s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA1']-gmTable['SA1_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA1']+gmTable['SA1_sd']),'--',c='tomato',label='Atk22 +/-std',markersize=0.7)

        plt.ylabel('SA1.0s (g)')
    elif imd=='sa3':
        plt.plot(rjb/1e3,sa1/9.8,'.',c='royalblue',label='SA3.0s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**gmTable['SA3'],'-',c='tomato',label='Atk22 mean;Vs30=300 m/s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA3']-gmTable['SA3_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA3']+gmTable['SA3_sd']),'--',c='tomato',label='Atk22 +/-std',markersize=0.7)
        plt.ylabel('SA3.0s (g)')
    elif imd=='sa0_3':
        plt.plot(rjb/1e3,sa1/9.8,'.',c='royalblue',label='SA0.3s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**gmTable['SA0.3'],'-',c='tomato',label='Atk22 mean;Vs30=300 m/s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA0.3']-gmTable['SA0.3_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA0.3']+gmTable['SA0.3_sd']),'--',c='tomato',label='Atk22 +/-std',markersize=0.7)
        plt.ylabel('SA0.3s (g)')
    else:
        plt.plot(rjb/1e3,pgv,'.',c='royalblue',label='PGV',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**gmTable['SA1'],'-',c='tomato',label='Atk22 mean;Vs30=300 m/s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA1']-gmTable['SA1_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA1']+gmTable['SA1_sd']),'--',c='tomato',label='Atk22 +/-std',markersize=0.7)
        plt.ylabel('PGV (m/s)')


    plt.xlim([1,200])
    plt.ylim([0.001,10])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('hypocentral distance (km)')
    plt.legend(loc=3)

    plt.savefig(modelname +'-'+ imd + '-Rjb.png',dpi=300)
    plt.close()

def plot_gme_hik(pgv,sa1,sa3,triang,coast_table, hik_trace, faultTable,modelname='test1', imd='pgv'):

    epi300 = (1740500,5385000)

    figsize_single = (8,6)

    fig,ax0 = plt.subplots(nrows=1,ncols=1,figsize=(figsize_single))

    if imd=='sa1':

        print('plotting SA1.0s')
        sc = ax0.tripcolor(triang,sa1[:]/9.8,cmap=Cmap2, shading='flat',vmin=0,vmax=1.5)
        cl = fig.colorbar(sc,ax=ax0)
        cl.set_label('SA1.0s (g)')

    elif imd=='sa3':

        print('plotting SA3.0s')

        sc = ax0.tripcolor(triang,sa3[:]/9.8,cmap=Cmap2, shading='flat',vmin=0,vmax=1.0)
        cl = fig.colorbar(sc,ax=ax0)
        cl.set_label('SA3.0s (g)')

    else:
        print('plotting PGV')

        sc = ax0.tripcolor(triang,pgv[:],cmap=Cmap2, shading='flat',vmin=0,vmax=1.5)
        cl = fig.colorbar(sc,ax=ax0)
        cl.set_label('PGV (m/s)')
            
    ax0.plot(coast_table['Points:0'],coast_table['Points:1'],'.k',markersize=1)

    # hikurangi slab edge 
    
    ax0.plot(hik_trace['Points:0'],hik_trace['Points:1'],'.',c='blue',alpha=0.7,label='Hikurangi slab',markersize=2.5)
    ax0.plot(hik_trace['Points:0'][0:1],hik_trace['Points:1'][0:1],'-',c='blue',alpha=0.7,label='Hikurangi slab',markersize=2.5)

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

    epi300 = (1753000,5434000)
    epi230 = (1753000,5434000)


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

    elif imd=='sa0_3':

        print('plotting SA 0.3s')

        sc = ax0.tripcolor(triang,sa0_3[:]/9.8,cmap=Cmap2, shading='flat',vmin=0,vmax=0.5)
        cl = fig.colorbar(sc,ax=ax0)

        cl.set_label('SA 0.3s (g)')

          
    else:
        print('plotting PGV')

        sc = ax0.tripcolor(triang,pgv[:],cmap=Cmap2, shading='flat',vmin=0,vmax=1.5)
        cl = fig.colorbar(sc,ax=ax0)

        cl.set_label('PGV (m/s)')
            
    ax0.plot(coast_table['Points:0'],coast_table['Points:1'],'.k',markersize=1)

    for ik, ids in enumerate(faultTable['Fault_ID']):
        df = faultTable[faultTable['Fault_ID']==ids]
        fxyz  = df["geometry"].get_coordinates()
        ax0.plot(fxyz['x'],fxyz['y'],'-',c='tomato',alpha=0.7)
        
    ax0.plot(fxyz['x'],fxyz['y'],'-',c='tomato',alpha=0.7, label='crustal fault')
    
    # ax0.set_yticks(lat_ticks)
    ax0.set(xlim=( epi230[0]-70e3, epi230[0]+70e3),ylim=(epi230[1]-70e3,epi230[1]+70e3))

    ax0.set_xlabel('x (m)')
    ax0.set_ylabel('y (m)')
    
    ax0.legend(loc=2)

    outname = imd + '-' + modelname + '.png'
    plt.savefig(outname,dpi=300,transparent=False)
    plt.close()