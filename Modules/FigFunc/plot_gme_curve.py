import numpy as np
import matplotlib.tri as tri



def plot_gme_curve(modelname, gmTable,rjb, imd='sa1'):


    plt.figure(figsize=figsize_single)

    # ddep = 9.0
    # Rhyp = np.sqrt((depi/1e3)**2+ddep**2)

    if imd=='sa1':
        plt.plot(rjb/1e3,sa1/9.8,'.',c='royalblue',label='SA1.0s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**gmTable['SA1'],'-',c='tomato',label='Atk22 mean;Vs30=300 m/s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA1']-gmTable['SA1_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA1']+gmTable['SA1_sd']),'--',c='tomato',label='Atk22 +/-std',markersize=0.7)
        plt.ylabel('SA1.0s (g)')
    elif imd=='sa3':
        plt.plot(rjb/1e3,sa3/9.8,'.',c='royalblue',label='SA3.0s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**gmTable['SA3'],'-',c='tomato',label='Atk22 mean;Vs30=300 m/s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA3']-gmTable['SA3_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA3']+gmTable['SA3_sd']),'--',c='tomato',label='Atk22 +/-std',markersize=0.7)
        plt.ylabel('SA3.0s (g)')
    elif imd=='sa0_3':
        plt.plot(rjb/1e3,sa0_3/9.8,'.',c='royalblue',label='SA0.3s',markersize=0.7)
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
    plt.legend(loc=1)

    plt.savefig(modelname +'-'+ imd + '-Rjb.png',dpi=300)
    plt.close()