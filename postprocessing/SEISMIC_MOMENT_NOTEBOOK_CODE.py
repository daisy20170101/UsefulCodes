"""
CODE TO ADD TO YOUR cov-analys.ipynb NOTEBOOK
==============================================

This shows exactly where to add the seismic moment calculation
in your existing notebook workflow.
"""

# ============================================================================
# ADD THIS TO YOUR NOTEBOOK (after your existing data loading section)
# ============================================================================

# At the top of your notebook, add to imports:
from StatsModel.cov_corr_funcs import (add_coordinates_from_xdmf,
                                       calculate_seismic_moment_by_fault,
                                       plot_seismic_moment_summary)


# Then, after you load ssTable1 and ssTable0, add this section:
# ----------------------------------------------------------------------------
# SEISMIC MOMENT CALCULATION
# ----------------------------------------------------------------------------

# Your existing code:
modelname = 'jp3z'
ssTable1 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/slab2/stress_{modelname}_final.csv')
ssTable0 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/slab2/stress_{modelname}_t0.csv')

print(ssTable1['PSR'].max())

ssTable1['sdrop'] = np.sqrt(ssTable0['Td0']**2+ ssTable0['Ts0']**2) - np.sqrt(ssTable1['Td0']**2+ ssTable1['Ts0']**2)
ssTable1['T0']= np.sqrt(ssTable0['Td0']**2+ ssTable0['Ts0']**2)
ssTable1['Pf'] = np.abs(ssTable0['Pn0'])
ssTable1['r'] = ssTable1['T0']/ssTable1['Pf']

# **NEW CODE STARTS HERE**

# Step 1: Add x, y, z coordinates from XDMF file
xdmf_file = f'/Users/DuoL/Documents/NSHM/Central/Joint3/data-{modelname}/{modelname}-fault.xdmf'

ssTable1 = add_coordinates_from_xdmf(
    fault_data=ssTable1,
    xdmf_file=xdmf_file,
    method='centroid'
)

print(f"\n✓ Added coordinates. Shape: {ssTable1.shape}")
print(f"Coordinate columns: {['x', 'y', 'z']}")
print(ssTable1[['x', 'y', 'z', 'ASl', 'Area']].head())

# Your existing filter (keep this):
ssTable2 = ssTable1[ssTable1['ASl'].between(0.5, 100, inclusive='both')]
print(ssTable1.__len__(), ssTable2.__len__(), ssTable2.keys())

# Step 2: Calculate seismic moment by fault segment
material_file = '/Users/DuoL/Documents/NSHM/Central/cov/mat3d_fault.csv'

df_moment = calculate_seismic_moment_by_fault(
    fault_data=ssTable2,
    material_file=material_file
)

# Step 3: Display results
print("\n" + "="*70)
print(f"SEISMIC MOMENT RESULTS FOR {modelname}")
print("="*70)
print(df_moment)

# Summary
total_M0 = df_moment['total_M0'].sum()
total_Mw = (2/3) * np.log10(total_M0) - 6.07
print(f"\nTotal seismic moment: {total_M0:.3e} N⋅m")
print(f"Equivalent moment magnitude: Mw = {total_Mw:.2f}")

# Step 4: Save results
df_moment.to_csv(f'statsTable/moment-fault-{modelname}.csv', index=False)
print(f"\n✓ Saved to: statsTable/moment-fault-{modelname}.csv")

# Step 5: Plot
fig, axes = plot_seismic_moment_summary(
    df_moment,
    save_path=f'moment-summary-{modelname}.png'
)
plt.show()

# **NEW CODE ENDS HERE**

# Continue with your existing notebook code...


# ============================================================================
# ALTERNATIVE: If you want to do this for multiple models
# ============================================================================

models = ['jp3z', 'jp3xSF2', 'jp3y']  # Your models
moment_results = []

for modelname in models:
    print(f"\n{'='*70}")
    print(f"Processing model: {modelname}")
    print('='*70)

    try:
        # Load data
        ssTable1 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/slab2/stress_{modelname}_final.csv')
        ssTable0 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/slab2/stress_{modelname}_t0.csv')

        # Calculate derived quantities
        ssTable1['sdrop'] = np.sqrt(ssTable0['Td0']**2+ ssTable0['Ts0']**2) - \
                           np.sqrt(ssTable1['Td0']**2+ ssTable1['Ts0']**2)
        ssTable1['T0']= np.sqrt(ssTable0['Td0']**2+ ssTable0['Ts0']**2)
        ssTable1['Pf'] = np.abs(ssTable0['Pn0'])
        ssTable1['r'] = ssTable1['T0']/ssTable1['Pf']

        # Add coordinates
        xdmf_file = f'/Users/DuoL/Documents/NSHM/Central/Joint3/data-{modelname}/{modelname}-fault.xdmf'
        ssTable1 = add_coordinates_from_xdmf(ssTable1, xdmf_file, method='centroid')

        # Filter
        ssTable2 = ssTable1[ssTable1['ASl'].between(0.5, 100)]

        # Calculate moment
        df_moment = calculate_seismic_moment_by_fault(
            ssTable2,
            '/Users/DuoL/Documents/NSHM/Central/cov/mat3d_fault.csv'
        )

        # Add model name
        df_moment['model'] = modelname
        moment_results.append(df_moment)

        # Save individual result
        df_moment.to_csv(f'statsTable/moment-fault-{modelname}.csv', index=False)

    except Exception as e:
        print(f"Error processing {modelname}: {e}")
        continue

# Combine all results
if moment_results:
    df_all_moments = pd.concat(moment_results, ignore_index=True)

    # Save combined
    df_all_moments.to_csv('statsTable/moment-all-models.csv', index=False)

    # Compare by fault-tag
    pivot = df_all_moments.pivot_table(
        index='fault-tag',
        columns='model',
        values=['Mw', 'total_M0', 'mean_slip'],
        aggfunc='first'
    )

    print("\n" + "="*70)
    print("MOMENT MAGNITUDE COMPARISON")
    print("="*70)
    print(pivot['Mw'])

    print("\n" + "="*70)
    print("TOTAL MOMENT COMPARISON (N⋅m)")
    print("="*70)
    print(pivot['total_M0'])


# ============================================================================
# QUICK REFERENCE
# ============================================================================

"""
REQUIRED FILES:
--------------
1. Fault CSV: /Users/DuoL/Documents/NSHM/Central/cov/slab2/stress_{model}_final.csv
2. XDMF file: /Users/DuoL/Documents/NSHM/Central/Joint3/data-{model}/{model}-fault.xdmf
3. Material: /Users/DuoL/Documents/NSHM/Central/cov/mat3d_fault.csv

FUNCTIONS:
---------
1. add_coordinates_from_xdmf(fault_data, xdmf_file, method='centroid')
   → Adds x, y, z columns from XDMF geometry

2. calculate_seismic_moment_by_fault(fault_data, material_file)
   → Calculates M₀ = μ × ASl × Area by fault-tag
   → Returns DataFrame with total_M0, Mw, statistics

3. plot_seismic_moment_summary(df_moment, save_path)
   → Creates 4-panel visualization

OUTPUT:
------
df_moment columns:
  - fault-tag: Segment ID
  - count: Number of elements
  - total_M0: Total moment (N⋅m)
  - Mw: Moment magnitude
  - mean_mu: Average shear modulus (Pa)
  - mean_slip: Average slip (m)
  - total_area: Total area (m²)

PHYSICS:
-------
M₀ = μ × ASl × Area  (seismic moment)
Mw = (2/3) × log₁₀(M₀) - 6.07  (moment magnitude, Hanks & Kanamori 1979)
"""
