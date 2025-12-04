import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def plot_colormap(cmap='viridis', vmin=0, vmax=1, orientation='horizontal', 
                  label='', figsize=(8, 1.2), tick_positions=None, 
                  tick_labels=None, output_file=None, show_numerical=True):
    """
    Plot a standalone colormap/colorbar with customizable parameters.
    
    Parameters:
    -----------
    cmap : str or matplotlib colormap
        Colormap to use (e.g., 'viridis', 'RdYlBu_r', 'plasma', etc.)
    vmin : float
        Minimum value for the colormap normalization
    vmax : float
        Maximum value for the colormap normalization
    orientation : str
        'horizontal' or 'vertical'
    label : str
        Label for the colorbar
    figsize : tuple
        Figure size (width, height)
    tick_positions : list or None
        Custom tick positions. If None, uses default
    tick_labels : list or None
        Custom tick labels. If None, uses tick_positions as labels
    output_file : str or None
        Filename to save the figure. If None, doesn't save
    show_numerical : bool
        Whether to show numerical scale on opposite side
    
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    
    Example:
    --------
    # Basic usage
    fig, ax = plot_colormap(cmap='RdYlBu_r', vmin=-90, vmax=90, 
                           label='Rake Angle (degrees)')
    
    # With custom ticks
    fig, ax = plot_colormap(cmap='plasma', vmin=0, vmax=100,
                           tick_positions=[0, 25, 50, 75, 100],
                           tick_labels=['Low', 'Med-Low', 'Med', 'Med-High', 'High'],
                           label='Intensity', output_file='colorbar.png')
    """

    
    # Adjust figure size based on orientation
    if orientation == 'vertical':
        figsize = (figsize[1], figsize[0])
    
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(bottom=0.5 if orientation == 'horizontal' else 0.2)
    
    # Get colormap
    if isinstance(cmap, str):
        cmap = plt.cm.get_cmap(cmap)
    
    # Create normalization
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    
    # Create the colorbar
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, 
                                    orientation=orientation)
    
    # Set custom ticks if provided
    if tick_positions is not None:
        cb.set_ticks(tick_positions)
        if tick_labels is not None:
            cb.set_ticklabels(tick_labels)
    
    # Set label
    if label:
        cb.set_label(label, fontsize=12, fontweight='bold')
    
    # Increase tick label size
    ax.tick_params(labelsize=11)
    
    # Add numerical scale on opposite side if requested
    if show_numerical and tick_labels is not None and tick_positions is not None:
        if orientation == 'horizontal':
            ax2 = ax.twiny()
            ax2.set_xlim(vmin, vmax)
            # Create reasonable numerical ticks
            num_ticks = 7
            tick_range = np.linspace(vmin, vmax, num_ticks)
            ax2.set_xticks(tick_range)
            ax2.tick_params(labelsize=11)
        else:
            ax2 = ax.twinx()
            ax2.set_ylim(vmin, vmax)
            num_ticks = 7
            tick_range = np.linspace(vmin, vmax, num_ticks)
            ax2.set_yticks(tick_range)
            ax2.tick_params(labelsize=11)
    
    # Save if output file specified
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    plt.show()
    
    return fig, ax