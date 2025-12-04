import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt


def plot_gme_curve_mean_std(modelname, gmTable,rjb, sa1,sa3,sa0_3,figsize_single,imd='sa1'):


    plt.figure(figsize=figsize_single)

    # ddep = 9.0
    # Rhyp = np.sqrt((depi/1e3)**2+ddep**2)

    # group into bins of Rjb and compute mean and std

    nbins = int(40)
    rjb_bins = np.logspace(0,2,nbins)
    mean_bins = np.zeros(nbins-1)
    std_bins = np.zeros(nbins-1)


    if imd=='sa1':

        for ik, rjb_bin in enumerate(rjb_bins[0:-1]):

            mean_bins[ik] = np.mean(np.log10(sa1[np.where((rjb/1e3 > rjb_bin)&(rjb/1e3 <rjb_bins[ik+1]))]))
            std_bins[ik]  = np.std(np.log10(sa1[np.where((rjb/1e3 > rjb_bin)&(rjb/1e3 <rjb_bins[ik+1]))]))

    elif imd=='sa3':
        
        for ik, rjb_bin in enumerate(rjb_bins[0:-1]):

            mean_bins[ik] = np.mean(np.log10(sa3[np.where((rjb/1e3 > rjb_bin)&(rjb/1e3 <rjb_bins[ik+1]))]))
            std_bins[ik]  = np.std(np.log10(sa3[np.where((rjb/1e3 > rjb_bin)&(rjb/1e3 <rjb_bins[ik+1]))]))

    elif imd=='sa0_3':

        for ik, rjb_bin in enumerate(rjb_bins[0:-1]):

            mean_bins[ik] = np.mean(np.log10(sa0_3[np.where((rjb/1e3 > rjb_bin)&(rjb/1e3 <rjb_bins[ik+1]))]))
            std_bins[ik]  = np.std(np.log10(sa0_3[np.where((rjb/1e3 > rjb_bin)&(rjb/1e3 <rjb_bins[ik+1]))]))
    else:
        print('you have to select one ground motion component!')

    if imd=='sa1':

        plt.plot(rjb/1e3,sa1/9.8,'.',c='skyblue',markersize=0.9)

        plt.plot(gmTable['r_rup'],10**gmTable['SA1'],'-',c='tomato',label='ATK22 mean',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA1'] - gmTable['SA1_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA1'] + gmTable['SA1_sd']),'--',c='tomato',label='ATK22 +/-1.5 std',markersize=0.7)

        plt.plot(rjb_bins[0:-1], 10**mean_bins/9.8,'-',c='royalblue',markersize=0.7,label='simulation mean')
        plt.plot(rjb_bins[0:-1], 10**(mean_bins - 1.28*std_bins)/9.8,'--',c='royalblue',markersize=0.7,label='simulation mean +/-1.5 std')
        plt.plot(rjb_bins[0:-1], 10**(mean_bins + 1.28*std_bins)/9.8,'--',c='royalblue',markersize=0.7)
        
        plt.ylabel('SA 1.0s (g)')
        plt.ylim([0.0001,100])

    elif imd=='sa3':

        plt.plot(rjb/1e3,sa3/9.8,'.',c='skyblue',markersize=0.9)

        plt.plot(gmTable['r_rup'],10**gmTable['SA3'],'-',c='tomato',label='ATK22 mean',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA3']- gmTable['SA3_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA3']+ gmTable['SA3_sd']),'--',c='tomato',label='ATK22 +/- 1.5 std',markersize=0.7)

        plt.plot(rjb_bins[0:-1], 10**mean_bins/9.8,'-',c='royalblue',markersize=0.7,label='simulation mean')
        plt.plot(rjb_bins[0:-1], 10**(mean_bins - 1.28*std_bins)/9.8,'--',c='royalblue',markersize=0.7,label='simulation mean +/-1.5 std')
        plt.plot(rjb_bins[0:-1], 10**(mean_bins + 1.28*std_bins)/9.8,'--',c='royalblue',markersize=0.7)

        plt.ylabel('SA 3.0s (g)')

        plt.ylim([0.000001,100])

    elif imd=='sa0_3':

        plt.plot(rjb/1e3,sa0_3/9.8,'.',c='skyblue',markersize=0.9)

        plt.plot(gmTable['r_rup'],10**gmTable['SA0.3'],'-',c='tomato',label='Atk22 mean;Vs30=300 m/s',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA0.3']-gmTable['SA0.3_sd']),'--',c='tomato',markersize=0.7)
        plt.plot(gmTable['r_rup'],10**(gmTable['SA0.3']+gmTable['SA0.3_sd']),'--',c='tomato',label='ATK22 +/- 1.5 std',markersize=0.7)

        plt.plot(rjb_bins[0:-1], 10**mean_bins/9.8,'-',c='royalblue',markersize=0.7,label='simulation mean')
        plt.plot(rjb_bins[0:-1], 10**(mean_bins-1.28*std_bins)/9.8,'--',c='royalblue',markersize=0.7,label='simulation mean +/- 1.5 std')
        plt.plot(rjb_bins[0:-1], 10**(mean_bins+1.28*std_bins)/9.8,'--',c='royalblue',markersize=0.7)
        
        plt.ylim([0.0001,100])

        plt.ylabel('SA 0.3s (g)')
        
    else:
        print('you have to select one ground motion component!')


    plt.xlim([1,200])
    # plt.ylim([0.0001,100])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('R_jb (km)')
    plt.legend(loc=1,fontsize=9)
    plt.grid(linestyle='-',which='both',linewidth=0.7,c='lightgray')

    plt.savefig(imd + '-' +modelname  + '-rjb-std.png',dpi=300)
    plt.close()