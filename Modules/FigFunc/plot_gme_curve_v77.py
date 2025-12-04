import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt


def plot_gme_curve(modelname, gmTable, gmTable2,rjb, sa1,sa3,figsize_single,imd='sa1'):


    plt.figure(figsize=figsize_single)

    # ddep = 9.0
    # Rhyp = np.sqrt((depi/1e3)**2+ddep**2)

    if imd=='sa1':

        plt.plot(rjb/1e3,sa1/9.8,'.',c='gray',markersize=1.2)
        plt.plot(gmTable['r_rup'],gmTable['SA1'],'-',c='tomato',label='AtK22 median',markersize=0.7)
        plt.plot(gmTable['r_rup'],(gmTable['SA1_upper'])+ 0.25*(gmTable['SA1_upper']-gmTable['SA1_low']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],(gmTable['SA1_low'])- 0.25*(gmTable['SA1_upper']-gmTable['SA1_low']),'--',c='tomato',label=r'ATK22;+/- 1.5$\sigma$',markersize=0.7)

        plt.plot(gmTable2['r_rup'],gmTable2['SA1'],'-',c='royalblue',label='AG20 median',markersize=0.7)
        plt.plot(gmTable2['r_rup'],(gmTable2['SA1_upper'])+ 0.25*(gmTable2['SA1_upper']-gmTable2['SA1_low']),'--',c='royalblue',markersize=0.7)
        plt.plot(gmTable2['r_rup'],(gmTable2['SA1_low'])- 0.25*(gmTable2['SA1_upper']-gmTable2['SA1_low']),'--',c='royalblue',label=r'ATK22;+/- 1.5$\sigma$',markersize=0.7)
        plt.ylabel('SA1.0s (g)')

    elif imd=='sa3':
        
        plt.plot(rjb/1e3,sa3/9.8,'.', c='gray',markersize=1.2)

        plt.plot(gmTable['r_rup'],gmTable['SA3'],'-',c='tomato',label='ATK22 mean',markersize=0.7)
        plt.plot(gmTable['r_rup'],(gmTable['SA3_upper'])+ 0.25*(gmTable['SA3_upper']-gmTable['SA3_low']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],(gmTable['SA3_low'])-  0.25*(gmTable['SA3_upper']-gmTable['SA3_low']),'--',c='tomato',label=r'ATK22;+/- 1.5$\sigma$',markersize=0.7)

        plt.plot(gmTable2['r_rup'],gmTable2['SA3'],'-',c='royalblue',label='AG20 mean',markersize=0.7)
        plt.plot(gmTable2['r_rup'],gmTable2['SA3_upper']+ 0.25*(gmTable2['SA3_upper']-gmTable2['SA3_low']),'--',c='royalblue',markersize=0.7)
        plt.plot(gmTable2['r_rup'],(gmTable2['SA3_low'])- 0.25*(gmTable2['SA3_upper']-gmTable2['SA3_low']),'--',c='royalblue',label=r'AG20;+/- 1.5$\sigma$',markersize=0.7)
        plt.ylabel('SA3.0s (g)')

    # elif imd=='sa0_3': 
    #     plt.plot(rjb/1e3,sa0_3/9.8,'.',c='royalblue',label='SA0.3s',markersize=0.7)
    #     plt.plot(gmTable['r_rup'],10**gmTable['SA0.3'],'-',c='tomato',label='Atk22 mean;Vs30=300 m/s',markersize=0.7)
    #     plt.plot(gmTable['r_rup'],10**(gmTable['SA0.3']-gmTable['SA0.3_sd']),'--',c='tomato',markersize=0.7)
    #     plt.plot(gmTable['r_rup'],10**(gmTable['SA0.3']+gmTable['SA0.3_sd']),'--',c='tomato',label='Atk22 +/-std',markersize=0.7)
    #     plt.ylabel('SA0.3s (g)')
    else:
        plt.plot(rjb/1e3,sa1,'.',c='royalblue',label='SA1.0s',markersize=0.7)
        plt.plot(gmTable['r_rup'],gmTable['SA1'],'-',c='tomato',label='AG20 mean',markersize=0.7)
        plt.plot(gmTable['r_rup'],(gmTable['SA1']-gmTable['SA1_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],(gmTable['SA1']+gmTable['SA1_sd']),'--',c='tomato',label='AG20 +/-std',markersize=0.7)
        plt.ylabel('PGV (m/s)')


    plt.xlim([1,100])
    plt.ylim([0.0001,10])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('R_jb (km)')
    plt.legend(loc=3,fontsize=9)

    plt.savefig(imd + '-' + modelname  + '-v2.png',dpi=300)
    plt.close()