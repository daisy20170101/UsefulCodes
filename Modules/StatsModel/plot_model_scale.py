# Plot Mw vs mean_slip for both datasets
import matplotlib.pyplot as plt
import numpy as np



def plot_model_scaling(df_slipm0_hik,
                       df_slipm0_crt,
                       ):
    fig, ax = plt.subplots(figsize=(6, 4))

    # Hikurangi data (fault-tag 68 from jp3/jp4)
    if len(df_slipm0_hik) > 0:
        ax.scatter(df_slipm0_hik['mean_slip'], df_slipm0_hik['Mw'], 
                s=30, c='tomato', marker='o', alpha=0.7, 
                label=f'Hikurangi ( n={len(df_slipm0_hik)})',
                edgecolors='tomato', linewidth=1.5)

    # Crustal data (other)
    if len(df_slipm0_crt) > 0:
        ax.scatter(df_slipm0_crt['mean_slip'], df_slipm0_crt['Mw'],
                s=30, c='violet', marker='^', alpha=0.7, 
                label=f'Crustal (n={len(df_slipm0_crt)})',
                edgecolors='violet', linewidth=1.5)

    # Add scaling relation line (global average: Mw = 6.93 + 0.82*log10(slip))

    slip_range = np.logspace(-0.3, 2, 100)

    # Mw_global = 6.93 + 0.82 * np.log10(slip_range)
    # ax.plot(slip_range, Mw_global, 'k:', linewidth=2, 
    #         label='Global: Mw=6.93+0.82×log(slip)')

    u_avg = np.arange(0.01,100,10)

    a,a_std,b,b_std = 6.93 ,0.05, 0.82,0.05

    ax.plot((u_avg),a + b*np.log10(u_avg),color='royalblue',linestyle='-',linewidth = 1.5,label=f'WC94')
    ax.fill_between((u_avg),(a-a_std) + (b)*np.log10(u_avg),(a+a_std) + (b)*np.log10(u_avg),alpha=0.1,
                    color='royalblue')

    # Mw={a} + {b}*1og(u_mean) 
    a,a_std,b,b_std = 6.09,0.06, 2,0.01

    ax.plot((u_avg),a + b*np.log10(u_avg),color='forestgreen',linestyle='-',linewidth = 1.5,label=f'DR04')
    ax.fill_between((u_avg),(a-a_std) + (b)*np.log10(u_avg),(a+a_std) + (b)*np.log10(u_avg),alpha=0.1,
                    color='forestgreen')

    
    ax.set_xscale('log')
    ax.set_xlabel('Mean Slip (m)', fontsize=12, weight='bold')
    ax.set_ylabel('Moment Magnitude (Mw)', fontsize=12, weight='bold')
    ax.set_title('Mw-Mean_slip Scaling: Hikurangi vs Crustal Faults', 
                fontsize=14, weight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim([0.1, 100])
    ax.set_ylim([4, 9])

    plt.tight_layout()
    return fig, ax
    