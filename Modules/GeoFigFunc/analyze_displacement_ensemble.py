"""
Analyze ensemble of displacement results from multiple models
Calculate mean and standard deviation across models for each query point
"""

import numpy as np
import matplotlib.pyplot as plt
from extract_xdmf_displacements import extract_displacements_at_positions


def collect_displacement_ensemble(model_list,
                                  query_positions,
                                  xdmf_template,
                                  displacement_fields=['u1', 'u2', 'u3'],
                                  method='linear',
                                  timestep=-1):
    """
    Collect displacements from multiple models at specified query positions

    Parameters:
    -----------
    model_list : list of str
        List of model names
    query_positions : array-like
        Query positions (N, 2) or (N, 3)
    xdmf_template : str
        Template string for XDMF file path with {} placeholder for model name
        Example: '/path/to/data/{}-surface.xdmf'
    displacement_fields : list
        Field names for displacement components
    method : str
        Interpolation method
    timestep : int
        Which timestep to extract (-1 for last)

    Returns:
    --------
    ensemble_data : dict
        Dictionary containing:
        - 'u1': array (n_models, n_points) - u1 for all models
        - 'u2': array (n_models, n_points) - u2 for all models
        - 'u3': array (n_models, n_points) - u3 for all models
        - 'magnitude': array (n_models, n_points) - displacement magnitude
        - 'model_list': list of model names
        - 'query_positions': array of query positions
    """

    n_models = len(model_list)
    query_positions = np.asarray(query_positions)
    n_points = len(query_positions)

    # Initialize storage arrays
    u1_all = np.zeros((n_models, n_points))
    u2_all = np.zeros((n_models, n_points))
    u3_all = np.zeros((n_models, n_points))
    mag_all = np.zeros((n_models, n_points))

    print(f"Collecting displacements from {n_models} models at {n_points} points")
    print("=" * 70)

    successful_models = []

    for id, imodel in enumerate(model_list):
        print(f"\n[{id+1}/{n_models}] Processing model: {imodel}")

        # Generate XDMF file path
        xdmf_file = xdmf_template.format(imodel)
        print(f"  File: {xdmf_file}")

        try:
            # Extract displacements
            results = extract_displacements_at_positions(
                xdmf_file,
                query_positions,
                displacement_fields=displacement_fields,
                method=method,
                timestep=timestep
            )

            # Store results
            u1_all[id, :] = results['u1']
            u2_all[id, :] = results['u2']
            u3_all[id, :] = results['u3']
            mag_all[id, :] = results['displacement_magnitude']

            successful_models.append(imodel)

            print(f"  ✓ Successfully extracted displacements")

        except Exception as e:
            print(f"  ✗ Error: {e}")
            # Set to NaN for failed models
            u1_all[id, :] = np.nan
            u2_all[id, :] = np.nan
            u3_all[id, :] = np.nan
            mag_all[id, :] = np.nan

    print("\n" + "=" * 70)
    print(f"Successfully processed {len(successful_models)}/{n_models} models")

    ensemble_data = {
        'u1': u1_all,
        'u2': u2_all,
        'u3': u3_all,
        'magnitude': mag_all,
        'model_list': model_list,
        'successful_models': successful_models,
        'query_positions': query_positions
    }

    return ensemble_data


def compute_ensemble_statistics(ensemble_data):
    """
    Compute mean and standard deviation across ensemble

    Parameters:
    -----------
    ensemble_data : dict
        Output from collect_displacement_ensemble()

    Returns:
    --------
    stats : dict
        Dictionary containing mean and std for each component:
        - 'u1_mean', 'u1_std': arrays (n_points,)
        - 'u2_mean', 'u2_std': arrays (n_points,)
        - 'u3_mean', 'u3_std': arrays (n_points,)
        - 'mag_mean', 'mag_std': arrays (n_points,)
    """

    print("\nComputing ensemble statistics...")

    # Use nanmean and nanstd to handle any failed models (NaN values)
    stats = {
        'u1_mean': np.nanmean(ensemble_data['u1'], axis=0),
        'u1_std': np.nanstd(ensemble_data['u1'], axis=0),
        'u2_mean': np.nanmean(ensemble_data['u2'], axis=0),
        'u2_std': np.nanstd(ensemble_data['u2'], axis=0),
        'u3_mean': np.nanmean(ensemble_data['u3'], axis=0),
        'u3_std': np.nanstd(ensemble_data['u3'], axis=0),
        'mag_mean': np.nanmean(ensemble_data['magnitude'], axis=0),
        'mag_std': np.nanstd(ensemble_data['magnitude'], axis=0),
        'query_positions': ensemble_data['query_positions']
    }

    print(f"  u1: mean range [{stats['u1_mean'].min():.6e}, {stats['u1_mean'].max():.6e}]")
    print(f"      std  range [{stats['u1_std'].min():.6e}, {stats['u1_std'].max():.6e}]")
    print(f"  u2: mean range [{stats['u2_mean'].min():.6e}, {stats['u2_mean'].max():.6e}]")
    print(f"      std  range [{stats['u2_std'].min():.6e}, {stats['u2_std'].max():.6e}]")
    print(f"  u3: mean range [{stats['u3_mean'].min():.6e}, {stats['u3_mean'].max():.6e}]")
    print(f"      std  range [{stats['u3_std'].min():.6e}, {stats['u3_std'].max():.6e}]")

    return stats


def plot_ensemble_statistics(stats,
                             output_file='ensemble_displacements.png',
                             figsize=(16, 12)):
    """
    Plot mean and standard deviation of displacements across ensemble

    Parameters:
    -----------
    stats : dict
        Output from compute_ensemble_statistics()
    output_file : str
        Output filename
    figsize : tuple
        Figure size
    """

    query_positions = stats['query_positions']
    n_points = len(query_positions)
    point_indices = np.arange(n_points)

    fig, axes = plt.subplots(3, 2, figsize=figsize)

    components = [
        ('u1', 'U1 (East)', axes[0, 0], axes[0, 1]),
        ('u2', 'U2 (North)', axes[1, 0], axes[1, 1]),
        ('u3', 'U3 (Vertical)', axes[2, 0], axes[2, 1])
    ]

    for comp_name, comp_label, ax_mean, ax_std in components:
        mean_key = f'{comp_name}_mean'
        std_key = f'{comp_name}_std'

        mean_vals = stats[mean_key]
        std_vals = stats[std_key]

        # Plot mean with error bars (±1 std)
        ax_mean.errorbar(point_indices, mean_vals, yerr=std_vals,
                        fmt='o-', capsize=5, capthick=2, linewidth=2,
                        markersize=8, label='Mean ± 1σ')
        ax_mean.axhline(0, color='k', linestyle='--', alpha=0.3)
        ax_mean.set_xlabel('Query Point Index', fontsize=12)
        ax_mean.set_ylabel(f'{comp_label} (m)', fontsize=12)
        ax_mean.set_title(f'{comp_label} - Mean', fontsize=14, fontweight='bold')
        ax_mean.grid(True, alpha=0.3)
        ax_mean.legend()

        # Plot standard deviation
        ax_std.plot(point_indices, std_vals, 'o-', linewidth=2, markersize=8,
                   color='orangered', label='Standard Deviation')
        ax_std.set_xlabel('Query Point Index', fontsize=12)
        ax_std.set_ylabel(f'{comp_label} Std Dev (m)', fontsize=12)
        ax_std.set_title(f'{comp_label} - Standard Deviation', fontsize=14, fontweight='bold')
        ax_std.grid(True, alpha=0.3)
        ax_std.legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved figure to: {output_file}")

    return fig


def save_ensemble_data(ensemble_data, stats, output_file='ensemble_data.npz'):
    """
    Save ensemble data and statistics to a compressed NumPy file

    Parameters:
    -----------
    ensemble_data : dict
        Output from collect_displacement_ensemble()
    stats : dict
        Output from compute_ensemble_statistics()
    output_file : str
        Output filename (.npz)
    """

    np.savez_compressed(
        output_file,
        # Raw ensemble data
        u1_all=ensemble_data['u1'],
        u2_all=ensemble_data['u2'],
        u3_all=ensemble_data['u3'],
        magnitude_all=ensemble_data['magnitude'],
        model_list=ensemble_data['model_list'],
        query_positions=ensemble_data['query_positions'],
        # Statistics
        u1_mean=stats['u1_mean'],
        u1_std=stats['u1_std'],
        u2_mean=stats['u2_mean'],
        u2_std=stats['u2_std'],
        u3_mean=stats['u3_mean'],
        u3_std=stats['u3_std'],
        mag_mean=stats['mag_mean'],
        mag_std=stats['mag_std']
    )

    print(f"\n✓ Saved ensemble data to: {output_file}")


def load_ensemble_data(input_file='ensemble_data.npz'):
    """
    Load previously saved ensemble data

    Parameters:
    -----------
    input_file : str
        Input filename (.npz)

    Returns:
    --------
    ensemble_data : dict
    stats : dict
    """

    data = np.load(input_file, allow_pickle=True)

    ensemble_data = {
        'u1': data['u1_all'],
        'u2': data['u2_all'],
        'u3': data['u3_all'],
        'magnitude': data['magnitude_all'],
        'model_list': data['model_list'].tolist(),
        'query_positions': data['query_positions']
    }

    stats = {
        'u1_mean': data['u1_mean'],
        'u1_std': data['u1_std'],
        'u2_mean': data['u2_mean'],
        'u2_std': data['u2_std'],
        'u3_mean': data['u3_mean'],
        'u3_std': data['u3_std'],
        'mag_mean': data['mag_mean'],
        'mag_std': data['mag_std'],
        'query_positions': data['query_positions']
    }

    print(f"✓ Loaded ensemble data from: {input_file}")
    print(f"  Number of models: {len(ensemble_data['model_list'])}")
    print(f"  Number of query points: {len(ensemble_data['query_positions'])}")

    return ensemble_data, stats


# ===========================================================================
# Example usage
# ===========================================================================

if __name__ == "__main__":

    # Define your models
    model_list = ['model1', 'model2', 'model3', 'model4', 'model5']

    # Define query positions
    query_positions = np.array([
        [174.5, -41.0],
        [174.7, -41.1],
        [174.9, -41.2],
        [175.1, -41.3],
        [175.3, -41.4]
    ])

    # XDMF file template
    xdmf_template = '/Volumes/heiterwange/NSHM/Hikurangi/joint3/data-jp3z/{}-surface.xdmf'

    print("=" * 70)
    print("ENSEMBLE DISPLACEMENT ANALYSIS")
    print("=" * 70)

    # Step 1: Collect displacements from all models
    ensemble_data = collect_displacement_ensemble(
        model_list=model_list,
        query_positions=query_positions,
        xdmf_template=xdmf_template,
        displacement_fields=['u1', 'u2', 'u3'],
        method='linear',
        timestep=-1
    )

    # Step 2: Compute statistics
    stats = compute_ensemble_statistics(ensemble_data)

    # Step 3: Save data
    save_ensemble_data(ensemble_data, stats, output_file='ensemble_data.npz')

    # Step 4: Plot results
    fig = plot_ensemble_statistics(stats, output_file='ensemble_displacements.png')

    # Show plots
    plt.show()

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
