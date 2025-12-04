# plot moment rate and compare

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_mt(foldername = '/Users/DuoL/Documents/SeisSol/Kaikoura/energy/',modelname='nohope2',
            
            ):


    energy_file = foldername +  modelname+'-energy.csv'

    data2 = pd.read_csv(energy_file)
    # print(data2)
    quality = data2.variable.unique()
    # print(quality)

    data_seis = data2[data2['variable']=='seismic_moment']
    seis = data_seis['measurement'].to_numpy()
    time = data_seis['time'].to_numpy()

    seis_rate  = np.diff(seis)/np.diff(time)
    mag  = 2/3*np.log10(seis[-1])-6.07

    fig= plt.figure(figsize=(5,3))

    ax0 = plt.subplot(111)               
    # plt.plot(enf[:,0],enf[:,1]*1e7/1e26,'-',color=colors2[4])   # from MPa*m**2/s to dyne-cm/s
    plt.plot(time[1::],(seis_rate)/1e20,'-',color='tomato',label= modelname +',Mw'+str(round(mag,2)))   # from MPa*m**2/s to dyne-cm/s

    plt.legend(loc=1,fontsize=10)
    plt.xlabel('time (s)')
    plt.ylabel('moment rate (*$10^{20}$ Nm/s)')
    ax0.set(xlim=(0,100))


    return fig, ax0