"""
Diagnostic tools for understanding unusual PCA variance patterns.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def diagnose_pca_variance(pca, variables, threshold=0.05):
    """
    Diagnose unusual PCA variance patterns and suggest explanations.

    Parameters
    ----------
    pca : sklearn PCA object
        Fitted PCA model
    variables : list of str
        Variable names
    threshold : float
        Minimum variance ratio to consider significant (default: 0.05 = 5%)

    Returns
    -------
    dict : Diagnostic information
    """
    var_ratios = pca.explained_variance_ratio_

    diagnostics = {
        'pc1_variance': var_ratios[0],
        'pc2_variance': var_ratios[1],
        'pc1_plus_pc2': var_ratios[0] + var_ratios[1],
        'pc2_greater_than_pc1': var_ratios[1] > var_ratios[0],
        'n_significant_pcs': np.sum(var_ratios > threshold)
    }

    print("="*70)
    print("PCA VARIANCE DIAGNOSTICS")
    print("="*70)
    print(f"\nPC1 explains: {var_ratios[0]:.1%}")
    print(f"PC2 explains: {var_ratios[1]:.1%}")
    print(f"PC1 + PC2:    {diagnostics['pc1_plus_pc2']:.1%}")

    # Check for unusual pattern
    if diagnostics['pc2_greater_than_pc1']:
        print("\n⚠️  WARNING: PC2 > PC1 (unusual!)")
        print("   This suggests:")
        print("   1. Check that data was standardized before PCA")
        print("   2. Multiple independent physical processes")
        print("   3. No single dominant mode of variation")

    # Cumulative variance
    cumvar = np.cumsum(var_ratios)
    print(f"\nCumulative variance:")
    for i in range(min(6, len(variables))):
        print(f"  PC1-PC{i+1}: {cumvar[i]:.1%}")

    # Number of PCs needed
    for target in [0.8, 0.9, 0.95]:
        n_needed = np.argmax(cumvar >= target) + 1
        print(f"\nPCs needed for {target:.0%} variance: {n_needed}")

    # Quality assessment
    print(f"\n{'='*70}")
    pc12_total = diagnostics['pc1_plus_pc2']
    if pc12_total > 0.8:
        quality = "EXCELLENT"
        message = "Two main patterns dominate"
    elif pc12_total > 0.6:
        quality = "GOOD"
        message = "Two main patterns capture most variation"
    elif pc12_total > 0.4:
        quality = "MODERATE"
        message = "Complex system, multiple independent processes"
    else:
        quality = "COMPLEX"
        message = "Very high dimensional, many independent factors"

    print(f"Overall Assessment: {quality}")
    print(f"Interpretation: {message}")
    print(f"{'='*70}\n")

    return diagnostics


def plot_variance_diagnosis(pca, variables, figsize=(12, 4)):
    """
    Create diagnostic plots for PCA variance analysis.

    Parameters
    ----------
    pca : sklearn PCA object
        Fitted PCA model
    variables : list of str
        Variable names
    figsize : tuple
        Figure size

    Returns
    -------
    fig, axes : matplotlib objects
    """
    var_ratios = pca.explained_variance_ratio_
    cumvar = np.cumsum(var_ratios)
    n_components = len(variables)

    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # 1. Individual variance (bar plot)
    colors = ['red' if i < 2 else 'steelblue' for i in range(n_components)]
    axes[0].bar(range(1, n_components+1), var_ratios, color=colors, alpha=0.7)
    axes[0].axhline(0.25, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='25% threshold')
    axes[0].set_xlabel('Principal Component', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Variance Explained', fontsize=11, fontweight='bold')
    axes[0].set_title('Individual PC Variance', fontsize=12, fontweight='bold')
    axes[0].set_ylim(0, max(var_ratios) * 1.1)
    axes[0].grid(True, alpha=0.3, axis='y')
    axes[0].legend()

    # Add percentage labels
    for i, v in enumerate(var_ratios):
        axes[0].text(i+1, v + 0.01, f'{v:.1%}', ha='center', va='bottom', fontsize=9)

    # 2. Cumulative variance
    axes[1].plot(range(1, n_components+1), cumvar, 'o-', linewidth=2, markersize=8, color='darkblue')
    axes[1].axhline(0.8, color='green', linestyle='--', linewidth=1, alpha=0.7, label='80% target')
    axes[1].axhline(0.9, color='orange', linestyle='--', linewidth=1, alpha=0.7, label='90% target')
    axes[1].set_xlabel('Number of Components', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Cumulative Variance Explained', fontsize=11, fontweight='bold')
    axes[1].set_title('Cumulative Variance', fontsize=12, fontweight='bold')
    axes[1].set_ylim(0, 1.05)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    # Add PC1+PC2 annotation
    pc12 = cumvar[1]
    axes[1].plot(2, pc12, 'ro', markersize=12, zorder=10)
    axes[1].annotate(f'PC1+PC2\n{pc12:.1%}', xy=(2, pc12), xytext=(2.5, pc12-0.1),
                    fontsize=10, fontweight='bold', color='red',
                    arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

    # 3. Variance ratio comparison (PC1 vs PC2)
    pc_names = [f'PC{i+1}' for i in range(min(4, n_components))]
    colors_comp = ['red', 'orange', 'steelblue', 'lightblue'][:len(pc_names)]
    axes[2].barh(pc_names, var_ratios[:len(pc_names)], color=colors_comp, alpha=0.7)
    axes[2].set_xlabel('Variance Explained', fontsize=11, fontweight='bold')
    axes[2].set_title('Top PCs Comparison', fontsize=12, fontweight='bold')
    axes[2].grid(True, alpha=0.3, axis='x')

    # Highlight if PC2 > PC1
    if var_ratios[1] > var_ratios[0]:
        axes[2].text(0.5, 0.95, 'PC2 > PC1 ⚠️', transform=axes[2].transAxes,
                    fontsize=11, fontweight='bold', color='red',
                    bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5),
                    verticalalignment='top', horizontalalignment='center')

    plt.tight_layout()
    return fig, axes


def check_data_scaling(X, variable_names):
    """
    Check if data has large scale differences that might affect PCA.

    Parameters
    ----------
    X : np.ndarray
        Data matrix (n_samples x n_features)
    variable_names : list of str
        Variable names

    Returns
    -------
    pd.DataFrame : Summary statistics
    """
    stats = pd.DataFrame({
        'mean': X.mean(axis=0),
        'std': X.std(axis=0),
        'min': X.min(axis=0),
        'max': X.max(axis=0),
        'range': X.max(axis=0) - X.min(axis=0)
    }, index=variable_names)

    print("\n" + "="*70)
    print("DATA SCALING CHECK")
    print("="*70)
    print("\nVariable statistics (before standardization):")
    print(stats.to_string())

    # Check for scale issues
    std_ratio = stats['std'].max() / stats['std'].min()
    range_ratio = stats['range'].max() / stats['range'].min()

    print(f"\nStd deviation ratio (max/min): {std_ratio:.2f}")
    print(f"Range ratio (max/min): {range_ratio:.2f}")

    if std_ratio > 10 or range_ratio > 10:
        print("\n⚠️  WARNING: Large scale differences detected!")
        print("   Standardization is ESSENTIAL before PCA")
    else:
        print("\n✓ Scale differences are reasonable")

    print("="*70 + "\n")

    return stats


def suggest_next_steps(diagnostics, n_features):
    """
    Suggest next analytical steps based on PCA diagnostics.

    Parameters
    ----------
    diagnostics : dict
        Output from diagnose_pca_variance
    n_features : int
        Number of original features
    """
    print("\n" + "="*70)
    print("RECOMMENDED NEXT STEPS")
    print("="*70)

    pc12_total = diagnostics['pc1_plus_pc2']

    if diagnostics['pc2_greater_than_pc1']:
        print("\n1. VERIFY DATA PREPARATION:")
        print("   - Confirm data was standardized (mean=0, std=1)")
        print("   - Check for outliers that might affect scaling")
        print("   - Verify all variables are on comparable scales")

    if pc12_total < 0.6:
        print("\n2. EXAMINE MORE COMPONENTS:")
        print(f"   - Look at PC3, PC4, ... (up to PC{diagnostics['n_significant_pcs']})")
        print("   - Use scree plot to identify 'elbow'")
        print("   - Consider that >2 independent processes exist")

    if pc12_total < 0.8:
        print("\n3. ALTERNATIVE ANALYSES:")
        print("   - Examine correlation matrix directly")
        print("   - Try pair plots for key variable combinations")
        print("   - Consider non-linear dimensionality reduction (t-SNE, UMAP)")

    print("\n4. PHYSICAL INTERPRETATION:")
    print("   - Examine PC loadings to understand what each PC represents")
    print("   - Check if PCs correspond to known physical processes")
    print("   - Consider if complexity reflects real earthquake physics")

    print("\n5. VALIDATION:")
    print("   - Compare PCA results across different subsets of data")
    print("   - Check if pattern is consistent across different models")
    print("   - Test if adding/removing variables changes structure")

    print("="*70 + "\n")


# Example usage function
def full_pca_diagnosis(data_dict, variables):
    """
    Run complete PCA diagnosis workflow.

    Parameters
    ----------
    data_dict : dict
        Dictionary with variable names as keys and numpy arrays as values
    variables : list of str
        Variable names to include

    Returns
    -------
    pca : fitted PCA object
    diagnostics : dict
    """
    # Prepare data
    data_list = [data_dict[var].flatten() for var in variables]
    X = np.column_stack(data_list)

    print("STEP 1: Check data scaling")
    stats = check_data_scaling(X, variables)

    print("\nSTEP 2: Standardize and fit PCA")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca = PCA()
    pca.fit(X_scaled)

    print("STEP 3: Diagnose variance pattern")
    diagnostics = diagnose_pca_variance(pca, variables)

    print("STEP 4: Create diagnostic plots")
    fig, axes = plot_variance_diagnosis(pca, variables)
    plt.show()

    print("STEP 5: Suggest next steps")
    suggest_next_steps(diagnostics, len(variables))

    return pca, diagnostics


if __name__ == "__main__":
    # Example with your specific case
    print("Example: Analyzing PC1=25%, PC2=31.7%")
    print("-" * 70)

    # Simulate your variance pattern
    explained_var = np.array([0.25, 0.317, 0.15, 0.12, 0.08, 0.05])

    class MockPCA:
        def __init__(self, var):
            self.explained_variance_ratio_ = var

    pca = MockPCA(explained_var)
    variables = ['stress', 'slip', 'vr', 'rake_deg', 's2ratio', 'xnuc']

    diag = diagnose_pca_variance(pca, variables)
    fig, axes = plot_variance_diagnosis(pca, variables)
    suggest_next_steps(diag, len(variables))
