
# Create discrete colormap with 5 colors
from matplotlib.colors import ListedColormap
import numpy as np
import matplotlib.pyplot as plt

def plot_bias_std_sigma(
        modelname: str='multicrust',
        deltaB_list: list[float]=[1.,1.0],
        deltaS2S_list: list[float]=[0.5,0.5],
        deltaWS_list: list[float]=[0.1,0.1],
        ped_list: list[float]=[1.0,3.0],
        bias_list: list[float]=[1.0,1.0],
):
        
    colors_discrete = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']  # Red, Blue, Green, Purple, Orange
    cmap_discrete = ListedColormap(colors_discrete)
    axeslabels = ['PGA', '0.02', '0.05', '0.1', '0.2', '0.6', '1.0', '2.0', '5.0', '10.0']
    ticklist=[0.01, 0.02, 0.05, 0.1, 0.2, 0.6, 1.0, 2.0, 5.0, 10.0]

    # Calculate total standard deviation
    deltaTotal_list = [np.sqrt(tau**2 + phi_s2s**2 + phi_ws**2) 
                    for tau, phi_s2s, phi_ws in zip(deltaB_list, deltaS2S_list, deltaWS_list)]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    axes[0].plot(ped_list, bias_list, '-^', c='k', linewidth=1.5, zorder=2, 
                label='multi-crust')

    axes[0].set_xscale('log')
    axes[0].set_ylim([-2.0, 2.0])
    axes[0].set_xlabel('Period (s)', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Bias', fontsize=12, fontweight='bold')
    axes[0].axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    axes[0].legend(loc='best', frameon=True, shadow=True, fontsize=10)
    axes[0].grid(which='both', linestyle=':', alpha=0.6)
    axes[0].tick_params(axis='both', which='major', labelsize=10)

    axes[0].set_xticks(ticklist)
    axes[0].set_xticklabels(axeslabels)

    # Second subplot - Standard deviation components + Total
    axes[1].plot(ped_list, deltaB_list, '-o', c='tomato', label=r'$\tau$ (between-event)', 
                linewidth=2.5, markersize=7, alpha=0.9)
    axes[1].plot(ped_list, deltaS2S_list, '-s', c='royalblue', label=r'$\phi_{S2S}$ (site-to-site)', 
                linewidth=2.5, markersize=7, alpha=0.9)
    axes[1].plot(ped_list, deltaWS_list, '-^', c='forestgreen', label=r'$\phi_{0}$ (within-site)', 
                linewidth=2.5, markersize=7, alpha=0.9)
    axes[1].plot(ped_list, deltaTotal_list, '-D', c='darkviolet', label=r'$\sigma$', 
                linewidth=2.5, markersize=7, alpha=0.9, markeredgecolor='black', markeredgewidth=0.5)

    axes[1].set_xscale('log')
    axes[1].set_ylim([0.0, 1.2])
    axes[1].set_xlabel('Period (s)', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Standard Deviation', fontsize=12, fontweight='bold')
    axes[1].legend(loc='best', frameon=True, shadow=True, fontsize=10)
    axes[1].grid(which='both', linestyle=':', alpha=0.6)
    axes[1].tick_params(axis='both', which='major', labelsize=10)

    axes[1].set_xticks(ticklist)
    axes[1].set_xticklabels(axeslabels)

    plt.tight_layout()

    plt.savefig(modelname+'-bias-sigma.png', dpi=300, bbox_inches='tight')

    return fig, axes
