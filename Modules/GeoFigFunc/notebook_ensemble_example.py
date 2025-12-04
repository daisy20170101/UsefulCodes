"""
Simple code snippet for Jupyter notebook
Collect ensemble displacements and plot mean ± std
"""

import numpy as np
import matplotlib.pyplot as plt
from extract_xdmf_displacements import extract_displacements_at_positions

# ============================================================================
# YOUR CODE HERE - Setup
# ============================================================================

# Define your model list
model_list = ['model1', 'model2', 'model3', 'model4', 'model5']

# Define query positions
query_positions = np.array([
    [174.5, -41.0],
    [174.7, -41.1],
    [174.9, -41.2],
    # Add more points...
])

# ============================================================================
# Collect displacements from all models
# ============================================================================

n_models = len(model_list)
n_points = len(query_positions)

# Initialize storage arrays
u1_all = []
u2_all = []
u3_all = []

print(f"Processing {n_models} models...")

for id, imodel in enumerate(model_list):
    print(f"[{id+1}/{n_models}] {imodel}")

    # Define your XDMF file
    xdmf_file = f'/Volumes/heiterwange/NSHM/Hikurangi/joint3/data-jp3z/{imodel}-surface.xdmf'

    try:
        # Extract displacements
        results = extract_displacements_at_positions(
            xdmf_file,
            query_positions,
            displacement_fields=['u1', 'u2', 'u3'],
            method='linear'
        )

        # Store results
        u1_all.append(results['u1'])
        u2_all.append(results['u2'])
        u3_all.append(results['u3'])

    except Exception as e:
        print(f"  Error: {e}")
        # Store NaN for failed models
        u1_all.append(np.full(n_points, np.nan))
        u2_all.append(np.full(n_points, np.nan))
        u3_all.append(np.full(n_points, np.nan))

# Convert lists to arrays (shape: n_models × n_points)
u1_all = np.array(u1_all)
u2_all = np.array(u2_all)
u3_all = np.array(u3_all)

print(f"\nData shape: {u1_all.shape}")

# ============================================================================
# Compute mean and standard deviation
# ============================================================================

u1_mean = np.nanmean(u1_all, axis=0)  # Mean across models for each point
u1_std = np.nanstd(u1_all, axis=0)    # Std across models for each point

u2_mean = np.nanmean(u2_all, axis=0)
u2_std = np.nanstd(u2_all, axis=0)

u3_mean = np.nanmean(u3_all, axis=0)
u3_std = np.nanstd(u3_all, axis=0)

print(f"\nU1: mean range [{u1_mean.min():.3e}, {u1_mean.max():.3e}]")
print(f"    std  range [{u1_std.min():.3e}, {u1_std.max():.3e}]")

# ============================================================================
# Plot results
# ============================================================================

fig, axes = plt.subplots(3, 1, figsize=(12, 10))

point_indices = np.arange(n_points)

# U1 component
axes[0].errorbar(point_indices, u1_mean, yerr=u1_std,
                fmt='o-', capsize=5, linewidth=2, markersize=8,
                label='Mean ± 1σ')
axes[0].axhline(0, color='k', linestyle='--', alpha=0.3)
axes[0].set_ylabel('U1 - East (m)', fontsize=12)
axes[0].set_title('U1 Component - Mean and Standard Deviation', fontsize=14, fontweight='bold')
axes[0].grid(True, alpha=0.3)
axes[0].legend()

# U2 component
axes[1].errorbar(point_indices, u2_mean, yerr=u2_std,
                fmt='o-', capsize=5, linewidth=2, markersize=8,
                color='orangered', label='Mean ± 1σ')
axes[1].axhline(0, color='k', linestyle='--', alpha=0.3)
axes[1].set_ylabel('U2 - North (m)', fontsize=12)
axes[1].set_title('U2 Component - Mean and Standard Deviation', fontsize=14, fontweight='bold')
axes[1].grid(True, alpha=0.3)
axes[1].legend()

# U3 component
axes[2].errorbar(point_indices, u3_mean, yerr=u3_std,
                fmt='o-', capsize=5, linewidth=2, markersize=8,
                color='green', label='Mean ± 1σ')
axes[2].axhline(0, color='k', linestyle='--', alpha=0.3)
axes[2].set_xlabel('Query Point Index', fontsize=12)
axes[2].set_ylabel('U3 - Vertical (m)', fontsize=12)
axes[2].set_title('U3 Component - Mean and Standard Deviation', fontsize=14, fontweight='bold')
axes[2].grid(True, alpha=0.3)
axes[2].legend()

plt.tight_layout()
plt.savefig('ensemble_displacements.png', dpi=300, bbox_inches='tight')
plt.show()

# ============================================================================
# Save data (optional)
# ============================================================================

np.savez('ensemble_data.npz',
         u1_all=u1_all,
         u2_all=u2_all,
         u3_all=u3_all,
         u1_mean=u1_mean,
         u1_std=u1_std,
         u2_mean=u2_mean,
         u2_std=u2_std,
         u3_mean=u3_mean,
         u3_std=u3_std,
         model_list=model_list,
         query_positions=query_positions)

print("\n✓ Data saved to ensemble_data.npz")

# To load later:
# data = np.load('ensemble_data.npz', allow_pickle=True)
# u1_mean = data['u1_mean']
# u1_std = data['u1_std']
# etc...
