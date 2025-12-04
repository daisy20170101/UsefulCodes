import matplotlib.pyplot as plt
from GMfunc.plot_pga_gmpe import bin_mean_std


def plot_dependency(
    df_event_merge,
    column_x="mag",
    column_y="event_term",
    nbins=10,
    out_name="scaling_depen",
):

    mean_bin, std_bins, rhypo_bins = bin_mean_std(
        nbins, 10 ** df_event_merge[column_y].values, df_event_merge[column_x].values
    )

    fig, ax = plt.subplots(1, 1, figsize=(6, 3))

    ax.plot(
        df_event_merge[column_x],
        df_event_merge[column_y],
        ".",
        alpha=0.1,
        color="lightcoral",
    )

    ax.plot(rhypo_bins, (mean_bin), "-", color="tomato", linewidth=1.7, label="mean")
    ax.fill_between(
        rhypo_bins,
        (mean_bin - std_bins),
        (mean_bin + std_bins),
        color="tomato",
        alpha=0.2,
    )

    ax.set_xlabel(column_x)
    ax.set_ylabel(column_y)
    ax.legend()
    ax.set_ylim([-1.5, 1.5])
    ax.hlines(y=0, xmin=3.8, xmax=7.8, linewidth=2, color="gray")
    ax.grid(which="both")

    plt.tight_layout()
    plt.savefig(out_name + ".png")
    return fig, ax
