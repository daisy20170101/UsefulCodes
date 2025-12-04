# loop for all
from StatsModel.combine_cov_data import group_by_fault_tag
from StatsModel.cov_corr_funcs import add_coordinates_from_xdmf,calculate_seismic_moment_by_fault
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

def plot_all_moment_rates_colored(energy_dir,
                                  pattern="*-energy.csv",
                                  scale_rate=1e20,
                                  cmap_name="plasma",
                                  mw_range=None,          # e.g., (5.0, 9.0); if None, inferred from data
                                  legend_loc=None,        # set to "best" to show legend; None hides it
                                  title="Seismic moment rate vs time",
                                  savepath=None,
                                  show=True):
    """
    Plot moment-rate vs time for all '*-energy.csv' files in energy_dir,
    coloring each line by Mw using the specified colormap.
    """
    energy_dir = Path(energy_dir)
    files = sorted(energy_dir.glob(pattern))
    if not files:
        print(f"[plot_all_moment_rates_colored] No files found in {energy_dir} matching '{pattern}'.")
        return {}

    # First pass: parse, compute (t_mid, rate_scaled, Mw), stash
    traces = []  # list of dicts with keys: model, t_mid, rate_scaled, Mw
    for f in files:
        name = f.stem
        model = name[:-len("-energy")] if name.endswith("-energy") else name

        try:
            df = pd.read_csv(f)
        except Exception as e:
            print(f"[skip] Failed to read {f.name}: {e}")
            continue

        if not {"variable", "measurement", "time"}.issubset(df.columns):
            print(f"[skip] Missing required columns in {f.name}.")
            continue

        data_seis = df[df["variable"] == "seismic_moment"].copy()
        if data_seis.empty:
            print(f"[skip] No 'seismic_moment' rows in {f.name}.")
            continue

        t = pd.to_numeric(data_seis["time"], errors="coerce").to_numpy()
        M0 = pd.to_numeric(data_seis["measurement"], errors="coerce").to_numpy()
        mask = np.isfinite(t) & np.isfinite(M0)
        t, M0 = t[mask], M0[mask]
        if t.size < 3:
            print(f"[skip] Not enough valid points in {f.name}.")
            continue

        order = np.argsort(t)
        t, M0 = t[order], M0[order]

        dt = np.diff(t)
        dM0 = np.diff(M0)
        valid = dt > 0
        if not np.any(valid):
            print(f"[skip] Non-increasing time sequence in {f.name}.")
            continue

        rate = np.full_like(dM0, np.nan, dtype=float)
        rate[valid] = dM0[valid] / dt[valid]
        t_mid = t[1:]
        rate_scaled = rate / scale_rate

        Mw = (2.0/3.0) * np.log10(M0[-1]) - 6.07 if np.isfinite(M0[-1]) and M0[-1] > 0 else np.nan
        print(f"{model}: Mw {Mw:.2f}" if np.isfinite(Mw) else f"{model}: Mw NaN")

        traces.append(dict(model=model, t_mid=t_mid, rate_scaled=rate_scaled, Mw=Mw))

    if not traces:
        print("[plot_all_moment_rates_colored] No valid traces to plot.")
        return {}

    # Set up colormap + normalization
    cmap = cm.get_cmap(cmap_name)
    mw_vals = np.array([tr["Mw"] for tr in traces if np.isfinite(tr["Mw"])])
    if mw_vals.size == 0:
        # fallback range if all Mw are NaN
        norm = Normalize(vmin=6.0, vmax=8.0)
    else:
        if mw_range is None:
            vmin, vmax = float(np.min(mw_vals)), float(np.max(mw_vals))
            if np.isclose(vmin, vmax):
                # expand slightly to avoid zero range
                vmin, vmax = vmin - 0.1, vmax + 0.1
        else:
            vmin, vmax = mw_range
        norm = Normalize(vmin=vmin, vmax=vmax)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    for tr in traces:
        color = cmap(norm(tr["Mw"])) if np.isfinite(tr["Mw"]) else (0.6, 0.6, 0.6, 1.0)  # gray if Mw NaN
        ax.plot(tr["t_mid"][np.isfinite(tr["rate_scaled"])],
                tr["rate_scaled"][np.isfinite(tr["rate_scaled"])],
                '-', lw=1.5,
                label=f"{tr['model']}, Mw {np.round(tr['Mw'], 2)}" if legend_loc else None,
                color=color)

    ax.set_xlabel("Time")
    ax.set_ylabel(f"Seismic moment rate / {scale_rate:g}")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    if legend_loc:
        ax.legend(loc=legend_loc, fontsize=9)

    # Colorbar
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])  # required by older MPL
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Mw")

    fig.tight_layout()

    if savepath is not None:
        savepath = Path(savepath)
        savepath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(savepath, dpi=200, bbox_inches="tight")
        print(f"[saved] {savepath}")

    if show:
        plt.show()

    d = {tr["model"]: tr["Mw"] for tr in traces}
    df = pd.DataFrame.from_dict(d, orient="index", columns=["Mw"])
    df.index.name = "model"
    
    # return model->Mw mapping
    return df,fig,ax



def prepare_fault_moment_data(
        model: str = 'slab2',
        xdmffilefolder: str='/Volumes/heiterwange/NSHM/Hikurangi/crustal/',
        model_list: list[str]=['jp3w','jp3yC'],
        variable_list: list[str]=['rake','ASl','Vr','T0','r','sdrop','Sdir'],
        add_xyz: bool=True,
        outdir: str='slab2test1',
    ):

    for imod, modelname in enumerate(model_list):

        print('start:' +modelname)
        
        # step 1: load fault output
        ssTable1 = pd.read_csv('/Users/DuoL/Documents/NSHM/Central/cov/'+ model+'/stress_'+modelname+'_final.csv')
        ssTable0 = pd.read_csv('/Users/DuoL/Documents/NSHM/Central/cov/'+ model +'/stress_'+modelname+'_t0.csv')

        ssTable1['sdrop'] = np.sqrt(ssTable0['Td0']**2+ ssTable0['Ts0']**2) - np.sqrt(ssTable1['Td0']**2+ ssTable1['Ts0']**2) 

        ssTable1['T0']= np.sqrt(ssTable0['Td0']**2+ ssTable0['Ts0']**2) 
        ssTable1['Pf'] = np.abs(ssTable0['Pn0'])
        ssTable1['r'] = ssTable1['T0']/ssTable1['Pf']
        ssTable1['sdrop'] = np.sqrt(ssTable0['Td0']**2+ ssTable0['Ts0']**2) - np.sqrt(ssTable1['Td0']**2+ ssTable1['Ts0']**2) 


        if add_xyz:
            # importlib.reload(StatsModel.cov_corr_funcs)
            # from StatsModel.cov_corr_funcs import add_coordinates_from_xdmf
            try:
                ssTable1=add_coordinates_from_xdmf(ssTable1,xdmffilefolder+'/'+ modelname +'-fault.xdmf')
            except:
                print('missing:'+modelname)
                pass
        
        # step 2, prepare statistic fault

        # filter data according to final slip
        minslip = 0.01*ssTable1['ASl'].mean()
        maxslip = 10*ssTable1['ASl'].mean()

        ssTable2 = ssTable1[ssTable1['ASl'].between(minslip,maxslip,inclusive='both')]

        # Filter Vr to maximum 8000 m/s
        ssTable2 = ssTable2[ssTable2['Vr'] <= 10000]


        if len(ssTable2)>0:
            # More comprehensive statistics
            df_stats = group_by_fault_tag(
                ssTable2,
                variable_list=variable_list,agg_func='mean'
            )


            # step3 calcualte M0

            # Calculate seismic moment by fault segment
            df_moment = calculate_seismic_moment_by_fault(
                fault_data=ssTable2,
                material_file='/Users/DuoL/Documents/NSHM/Central/cov/resample_mat.csv'
            )

            # Display results
            print(df_moment)

            df_final=pd.merge(df_stats,df_moment,on='fault-tag',how='inner')

            df_final.to_csv('statsTable/'+ outdir + '/'+'statsM0-'+modelname+'.csv')
        

    print('finished!')

from StatsModel.cov_corr_funcs import add_coordinates_from_xdmf

def prepare_fault_hik_crust_data(model_list_crt,
                                 model_list_hik,
                                 ):

    slipm0_hik = []  # Hikurangi: fault-tag=68 in jp3/jp4 models
    slipm0_crt = []  # Crustal: all other cases

    # Material properties file
    material_file = '/Users/DuoL/Documents/NSHM/Central/cov/resample_mat.csv'

    # Process crustal models
    print("Processing crustal models...")
    
    for modelname in model_list_crt:
        try:
            print(f"  Processing {modelname}...", end=' ')
            
            # Load data
            ssTable1 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/crustal/stress_{modelname}_final.csv')
            ssTable0 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/crustal/stress_{modelname}_t0.csv')
            
            # Calculate derived quantities
            ssTable1['T0'] = np.sqrt(ssTable0['Td0']**2 + ssTable0['Ts0']**2)
            ssTable1['Pf'] = np.abs(ssTable0['Pn0'])
            ssTable1['r'] = ssTable1['T0'] / ssTable1['Pf']
            ssTable1['sdrop'] = np.sqrt(ssTable0['Td0']**2+ ssTable0['Ts0']**2) - np.sqrt(ssTable1['Td0']**2+ ssTable1['Ts0']**2) 

            # Add coordinates from XDMF
            xdmf_file = f'/Volumes/heiterwange/NSHM/Hikurangi/crustal/{modelname}-fault.xdmf'
            ssTable1 = add_coordinates_from_xdmf(ssTable1, xdmf_file, method='centroid')
            
            # Filter
            ssTable2 = ssTable1[ssTable1['ASl'].between(0.5, 100)]
            
            # Calculate seismic moment
            df_moment = calculate_seismic_moment_by_fault(ssTable2, material_file)
            
            # Collect mean_slip and Mw for each fault-tag
            for idx, row in df_moment.iterrows():
                fault_tag = row['fault-tag']
                
                data_entry = {
                    'model': modelname,
                    'fault-tag': fault_tag,
                    'mean_slip': row['mean_slip'],
                    'Mw': row['Mw'],
                    'total_M0': row['total_M0'],
                    'count': row['count']
                }
                
                
                slipm0_crt.append(data_entry)
            
            print(f"✓ ({len(df_moment)} segments)")
            
        except FileNotFoundError:
            print(f"✗ Files not found")
        except Exception as e:
            print(f"✗ Error: {e}")

    # Process Hikurangi models
    print("\nProcessing Hikurangi models...")
    for modelname in model_list_hik:
        try:
            print(f"  Processing {modelname}...", end=' ')
            
            # Load data
            ssTable1 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/slab2/stress_{modelname}_final.csv')
            ssTable0 = pd.read_csv(f'/Users/DuoL/Documents/NSHM/Central/cov/slab2/stress_{modelname}_t0.csv')
            
            # Calculate derived quantities
            ssTable1['T0'] = np.sqrt(ssTable0['Td0']**2 + ssTable0['Ts0']**2)
            ssTable1['Pf'] = np.abs(ssTable0['Pn0'])
            ssTable1['r'] = ssTable1['T0'] / ssTable1['Pf']
            
            # Add coordinates from XDMF
            xdmf_file = f'/Volumes/heiterwange/NSHM/Hikurangi/FaultOutput/{modelname}-fault.xdmf'
            ssTable1 = add_coordinates_from_xdmf(ssTable1, xdmf_file, method='centroid')
            
            # Filter
            ssTable2 = ssTable1[ssTable1['ASl'].between(0.5, 100)]
            
            # Calculate seismic moment
            df_moment = calculate_seismic_moment_by_fault(ssTable2, material_file)
            
            # Collect mean_slip and Mw for each fault-tag
            for idx, row in df_moment.iterrows():

                fault_tag = row['fault-tag']
                
                data_entry = {
                    'model': modelname,
                    'fault-tag': fault_tag,
                    'mean_slip': row['mean_slip'],
                    'Mw': row['Mw'],
                    'total_M0': row['total_M0'],
                    'count': row['count']
                }
                
                # Check if fault-tag is 68 AND model starts with jp3 or jp4
                if fault_tag == 68 and (modelname.startswith('jp3') or modelname.startswith('jp4')):
                    slipm0_hik.append(data_entry)
                else:
                    slipm0_crt.append(data_entry)
            
            print(f"✓ ({len(df_moment)} segments)")
            
        except FileNotFoundError:
            print(f"✗ Files not found")
        except Exception as e:
            print(f"✗ Error: {e}")

    # Convert to DataFrames
    df_slipm0_hik = pd.DataFrame(slipm0_hik)
    df_slipm0_crt = pd.DataFrame(slipm0_crt)

    # Summary
    print("\n" + "="*70)
    print("DATA COLLECTION SUMMARY")
    print("="*70)
    print(f"\nHikurangi (fault-tag=68 in jp3/jp4 models): {len(df_slipm0_hik)} records")
    if len(df_slipm0_hik) > 0:
        print(f"  Models: {df_slipm0_hik['model'].unique().tolist()}")
        print(f"  Mw range: [{df_slipm0_hik['Mw'].min():.2f}, {df_slipm0_hik['Mw'].max():.2f}]")
        print(f"  Slip range: [{df_slipm0_hik['mean_slip'].min():.2f}, {df_slipm0_hik['mean_slip'].max():.2f}] m")

    print(f"\nCrustal (other cases): {len(df_slipm0_crt)} records")
    if len(df_slipm0_crt) > 0:
        print(f"  Models: {df_slipm0_crt['model'].unique().tolist()}")
        print(f"  Mw range: [{df_slipm0_crt['Mw'].min():.2f}, {df_slipm0_crt['Mw'].max():.2f}]")
        print(f"  Slip range: [{df_slipm0_crt['mean_slip'].min():.2f}, {df_slipm0_crt['mean_slip'].max():.2f}] m")

    # Display data
    print("\n" + "="*70)
    print("HIKURANGI DATA (fault-tag=68 in jp3/jp4)")
    print("="*70)
    print(df_slipm0_hik)

    print("\n" + "="*70)
    print("CRUSTAL DATA (other)")
    print("="*70)
    print(df_slipm0_crt)

    # Save to CSV
    df_slipm0_hik.to_csv('slipm0_hikurangi.csv', index=False)
    df_slipm0_crt.to_csv('slipm0_crustal.csv', index=False)

    print("\n✓ Saved to: slipm0_hikurangi.csv and slipm0_crustal.csv")

    return slipm0_hik, slipm0_crt,df_slipm0_hik,df_slipm0_crt


def calculate_segment_timing(
    fault_data: pd.DataFrame,
    material_file: str = None,
    rt_col: str = 'RT',
    fault_tag_col: str = 'fault-tag',
    asl_col: str = 'ASl',
    area_col: str = 'Area'
) -> pd.DataFrame:
    """
    Calculate rupture timing and seismic moment for each fault segment.

    For each fault segment (identified by fault-tag), this function:
    1. Filters out non-ruptured elements (RT = 0.0) for timing calculations
    2. Calculates timing statistics from RUPTURED elements only:
       - RT_start: Minimum rupture onset time (first element to rupture)
       - RT_end: Maximum rupture onset time (last element to rupture)
       - duration: RT_end - RT_start (rupture duration)
    3. Calculates seismic moment from elements with 0.5 < ASl < 100 m (if material_file provided):
       - Mw: Moment magnitude
       - mean_slip: Mean slip per segment
       - total_area: Total fault area

    IMPORTANT:
    - Timing: Only elements with RT > 0.0 are used
    - Moment: Only elements with 0.5 < ASl < 100 m are used (filters out very low/high slip)

    Parameters
    ----------
    fault_data : pd.DataFrame
        Fault data with columns including RT (rupture onset time), fault-tag,
        and x, y, z coordinates (if material_file provided)
    material_file : str, optional
        Path to material properties CSV file. If provided, calculates Mw and mean_slip.
        File should have columns: x, y, z, mu
    rt_col : str, optional
        Name of the rupture time column (default: 'RT')
    fault_tag_col : str, optional
        Name of the fault segment ID column (default: 'fault-tag')
    asl_col : str, optional
        Name of the slip column (default: 'ASl')
    area_col : str, optional
        Name of the area column (default: 'Area')

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - fault-tag: Segment identifier
        - count: Number of RUPTURED elements in segment (RT > 0)
        - RT_start: Start time (min RT of ruptured elements)
        - RT_end: End time (max RT of ruptured elements)
        - duration: Rupture duration (RT_end - RT_start)
        - Mw: Moment magnitude (if material_file provided, from 0.5 < ASl < 100 m)
        - mean_slip: Mean slip (if material_file provided, from 0.5 < ASl < 100 m)
        - total_M0: Total seismic moment (if material_file provided, from 0.5 < ASl < 100 m)
        - total_area: Total fault area (if material_file provided, from 0.5 < ASl < 100 m)

    Example
    -------
    >>> ssTable = pd.read_csv('stress_jp3z_final.csv')
    >>> df_timing = calculate_segment_timing(
    ...     ssTable,
    ...     material_file='/path/to/mat3d_fault.csv'
    ... )
    >>> print(df_timing)
       fault-tag  count  RT_start  RT_end  duration    Mw  mean_slip     total_M0
    0         65  23639      2.45   15.23     12.78  7.82       6.78  8.765e+20
    1         68  44737      1.20   18.45     17.25  8.45      15.23  2.456e+21
    2          3  26861      2.34   20.12     17.78  8.12      10.45  1.234e+21

    Note: RT_start values are > 0 since RT=0 elements are filtered out
    """
    print("\n" + "="*70)
    print("CALCULATING SEGMENT TIMING AND MOMENT")
    print("="*70 + "\n")

    # Check required columns
    required_cols = [rt_col, fault_tag_col]
    missing = [col for col in required_cols if col not in fault_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    print(f"Fault data: {len(fault_data)} elements")
    print(f"RT range (all): [{fault_data[rt_col].min():.3f}, {fault_data[rt_col].max():.3f}] s")
    print(f"Number of segments: {fault_data[fault_tag_col].nunique()}")

    # Filter out non-ruptured elements (RT = 0.0)
    fault_data_ruptured = fault_data[fault_data[rt_col] > 0.0].copy()

    n_total = len(fault_data)
    n_ruptured = len(fault_data_ruptured)
    n_not_ruptured = n_total - n_ruptured

    print(f"\nFiltering RT > 0.0 (ruptured elements only):")
    print(f"  Total elements: {n_total}")
    print(f"  Ruptured (RT > 0): {n_ruptured} ({100*n_ruptured/n_total:.1f}%)")
    print(f"  Not ruptured (RT = 0): {n_not_ruptured} ({100*n_not_ruptured/n_total:.1f}%)")

    if n_ruptured == 0:
        print("\n⚠ WARNING: No ruptured elements found (all RT = 0)!")
        # Return empty DataFrame with correct columns
        empty_df = pd.DataFrame(columns=[fault_tag_col, 'count', 'RT_start', 'RT_end', 'duration'])
        return empty_df

    print(f"RT range (ruptured): [{fault_data_ruptured[rt_col].min():.3f}, {fault_data_ruptured[rt_col].max():.3f}] s")

    # Group by fault-tag and calculate timing statistics on RUPTURED elements only
    grouped = fault_data_ruptured.groupby(fault_tag_col)[rt_col]

    df_timing = pd.DataFrame({
        fault_tag_col: fault_data_ruptured.groupby(fault_tag_col)[fault_tag_col].first(),
        'count': grouped.count(),  # Count of ruptured elements only
        'RT_start': grouped.min(),
        'RT_end': grouped.max()
    }).reset_index(drop=True)

    # Calculate duration
    df_timing['duration'] = df_timing['RT_end'] - df_timing['RT_start']

    # Calculate seismic moment if material file provided
    if material_file is not None:
        print("\n" + "-"*70)
        print("Calculating seismic moment...")
        print("-"*70)

        try:
            # Check if coordinates exist
            coord_cols = ['x', 'y', 'z']
            if not all(col in fault_data.columns for col in coord_cols):
                print("⚠ Warning: Missing x, y, z coordinates. Skipping moment calculation.")
                print("  Use add_coordinates_from_xdmf() first to add coordinates.")
            else:
                # Filter by ASl range (0.5 to 100 m) for moment calculation
                fault_data_filtered = fault_data[fault_data[asl_col].between(0.5, 100)].copy()

                n_before = len(fault_data)
                n_after = len(fault_data_filtered)
                print(f"Filtering by ASl (0.5 < ASl < 100 m):")
                print(f"  Elements before filter: {n_before}")
                print(f"  Elements after filter: {n_after} ({100*n_after/n_before:.1f}%)")

                if n_after == 0:
                    print("⚠ Warning: No elements remain after ASl filtering. Skipping moment calculation.")
                else:
                    # Calculate moment using existing function on filtered data
                    df_moment = calculate_seismic_moment_by_fault(
                        fault_data=fault_data_filtered,
                        material_file=material_file,
                        asl_col=asl_col,
                        area_col=area_col,
                        fault_tag_col=fault_tag_col
                    )

                    # Merge timing and moment data
                    df_timing = pd.merge(
                        df_timing,
                        df_moment[[fault_tag_col, 'Mw', 'mean_slip', 'total_M0', 'total_area']],
                        on=fault_tag_col,
                        how='left'
                    )

                    print("✓ Seismic moment calculated and merged")

        except Exception as e:
            print(f"⚠ Warning: Could not calculate seismic moment: {e}")
            print("  Continuing with timing data only.")

    print("\n" + "="*70)
    print("TIMING AND MOMENT RESULTS")
    print("="*70)
    print(df_timing)

    print(f"\nDuration range: [{df_timing['duration'].min():.3f}, {df_timing['duration'].max():.3f}] s")
    print(f"Mean duration: {df_timing['duration'].mean():.3f} s")

    if 'Mw' in df_timing.columns:
        print(f"\nMw range: [{df_timing['Mw'].min():.2f}, {df_timing['Mw'].max():.2f}]")
        print(f"Mean slip range: [{df_timing['mean_slip'].min():.2f}, {df_timing['mean_slip'].max():.2f}] m")

    return df_timing


def process_models_timing(
    model_list: list,
    data_folder: str,
    material_file: str = None,
    xdmf_folder: str = None,
    prefix: str = 'stress_',
    suffix: str = '_final.csv',
    rt_col: str = 'RT',
    fault_tag_col: str = 'fault-tag',
    asl_col: str = 'ASl',
    area_col: str = 'Area',
    add_xyz: bool = True
) -> pd.DataFrame:
    """
    Process multiple models and collect timing and moment data for all segments.

    This function loops through a list of models, loads each stress CSV file,
    optionally adds coordinates from XDMF, calculates timing and seismic moment
    statistics for each fault segment, and combines all results into a single DataFrame.

    NOTE:
    - Timing calculations (RT_start, RT_end, duration): Only RT > 0.0 elements
    - Moment calculations (Mw, mean_slip, total_area): Only 0.5 < ASl < 100 m elements

    Parameters
    ----------
    model_list : list
        List of model names (e.g., ['jp3z', 'jp3xSF2', 'jp3y'])
    data_folder : str
        Path to folder containing stress CSV files
        (e.g., '/Users/DuoL/Documents/NSHM/Central/cov/slab3v')
    material_file : str, optional
        Path to material properties CSV file for moment calculation
        (e.g., '/Users/DuoL/Documents/NSHM/Central/cov/resample_mat.csv')
    xdmf_folder : str, optional
        Path to folder containing XDMF files for coordinate extraction
        (e.g., '/Volumes/heiterwange/NSHM/Hikurangi/FaultOutput/')
        If None, assumes coordinates already exist in CSV
    prefix : str, optional
        File name prefix (default: 'stress_')
    suffix : str, optional
        File name suffix (default: '_final.csv')
    rt_col : str, optional
        Name of rupture time column (default: 'RT')
    fault_tag_col : str, optional
        Name of fault segment ID column (default: 'fault-tag')
    asl_col : str, optional
        Name of slip column (default: 'ASl')
    area_col : str, optional
        Name of area column (default: 'Area')
    add_xyz : bool, optional
        Whether to add coordinates from XDMF (default: True)

    Returns
    -------
    pd.DataFrame
        Combined timing and moment data for all models with columns:
        - model: Model name
        - fault-tag: Segment identifier
        - count: Number of RUPTURED elements (RT > 0)
        - RT_start: Start time (min RT of ruptured elements)
        - RT_end: End time (max RT of ruptured elements)
        - duration: Rupture duration (RT_end - RT_start)
        - Mw: Moment magnitude (if material_file provided, from all elements)
        - mean_slip: Mean slip (if material_file provided, from all elements)
        - total_M0: Total seismic moment (if material_file provided, from all elements)
        - total_area: Total fault area (if material_file provided, from all elements)

    Example
    -------
    >>> model_list_hik = ['jp3z', 'jp3xSF2', 'jp3y']
    >>> df_all = process_models_timing(
    ...     model_list_hik,
    ...     data_folder='/Users/DuoL/Documents/NSHM/Central/cov/slab3v',
    ...     material_file='/Users/DuoL/Documents/NSHM/Central/cov/resample_mat.csv',
    ...     xdmf_folder='/Volumes/heiterwange/NSHM/Hikurangi/FaultOutput/'
    ... )
    >>> print(df_all)
    """
    print("\n" + "="*70)
    print("PROCESSING MULTIPLE MODELS - TIMING AND MOMENT ANALYSIS")
    print("="*70)
    print(f"\nModels to process: {len(model_list)}")
    print(f"Data folder: {data_folder}")
    print(f"File pattern: {prefix}<model>{suffix}")
    if material_file:
        print(f"Material file: {material_file}")
    if xdmf_folder:
        print(f"XDMF folder: {xdmf_folder}")
    print()

    all_timing = []

    for i, modelname in enumerate(model_list, 1):
        print(f"\n[{i}/{len(model_list)}] Processing model: {modelname}")
        print("-" * 70)

        filepath = f"{data_folder}/{prefix}{modelname}{suffix}"

        try:
            # Load fault data
            df = pd.read_csv(filepath)
            print(f"✓ Loaded: {filepath}")
            print(f"  Rows: {len(df)}")

            # Add coordinates from XDMF if requested
            if add_xyz and xdmf_folder is not None:
                xdmf_file = f"{xdmf_folder}/{modelname}-fault.xdmf"
                try:
                    print(f"  Adding coordinates from XDMF...")
                    df = add_coordinates_from_xdmf(df, xdmf_file, method='centroid')
                    print(f"  ✓ Coordinates added")
                except Exception as e:
                    print(f"  ⚠ Could not add coordinates: {e}")

            # Calculate timing (and moment if material_file provided)
            df_timing = calculate_segment_timing(
                df,
                material_file=material_file,
                rt_col=rt_col,
                fault_tag_col=fault_tag_col,
                asl_col=asl_col,
                area_col=area_col
            )

            # Add model name
            df_timing['model'] = modelname

            all_timing.append(df_timing)

            print(f"\n✓ Processed {modelname}: {len(df_timing)} segments")

        except FileNotFoundError:
            print(f"✗ File not found: {filepath}")
            continue
        except Exception as e:
            print(f"✗ Error processing {modelname}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Combine all results
    if not all_timing:
        print("\n⚠ WARNING: No timing data collected!")
        return pd.DataFrame()

    df_all_timing = pd.concat(all_timing, ignore_index=True)

    # Summary
    print("\n" + "="*70)
    print("COMBINED TIMING AND MOMENT RESULTS")
    print("="*70)
    print(f"\nTotal records: {len(df_all_timing)}")
    print(f"Models processed: {df_all_timing['model'].nunique()}")
    print(f"Unique fault segments: {df_all_timing[fault_tag_col].nunique()}")

    print("\nDuration statistics across all models:")
    print(f"  Min: {df_all_timing['duration'].min():.3f} s")
    print(f"  Max: {df_all_timing['duration'].max():.3f} s")
    print(f"  Mean: {df_all_timing['duration'].mean():.3f} s")
    print(f"  Median: {df_all_timing['duration'].median():.3f} s")

    if 'Mw' in df_all_timing.columns:
        print("\nMoment magnitude statistics:")
        print(f"  Min Mw: {df_all_timing['Mw'].min():.2f}")
        print(f"  Max Mw: {df_all_timing['Mw'].max():.2f}")
        print(f"  Mean Mw: {df_all_timing['Mw'].mean():.2f}")

        print("\nMean slip statistics:")
        print(f"  Min: {df_all_timing['mean_slip'].min():.2f} m")
        print(f"  Max: {df_all_timing['mean_slip'].max():.2f} m")
        print(f"  Mean: {df_all_timing['mean_slip'].mean():.2f} m")

    # Reorder columns for clarity
    base_columns = ['model', fault_tag_col, 'count', 'RT_start', 'RT_end', 'duration']

    # Add moment columns if they exist
    if 'Mw' in df_all_timing.columns:
        base_columns.extend(['Mw', 'mean_slip', 'total_M0', 'total_area'])

    # Only include columns that actually exist
    column_order = [col for col in base_columns if col in df_all_timing.columns]
    df_all_timing = df_all_timing[column_order]

    print("\nFirst 10 rows:")
    print(df_all_timing.head(10))

    return df_all_timing


def plot_segment_timing(
    df_timing: pd.DataFrame,
    fault_tag_col: str = 'fault-tag',
    save_path: str = None,
    figsize: tuple = (12, 8),
    dpi: int = 150
):
    """
    Visualize fault segment timing and moment data.

    Creates a 4-panel figure showing:
    1. Timeline of rupture (start/end times by segment)
    2. Duration histogram
    3. Scatter plot: Duration vs Mw (if available) or RT_start vs duration
    4. Summary statistics table (includes Mw and mean_slip if available)

    Parameters
    ----------
    df_timing : pd.DataFrame
        Timing data from calculate_segment_timing() or process_models_timing()
        Must have columns: fault-tag, count, RT_start, RT_end, duration
        Optional: model (for multi-model data), Mw, mean_slip, total_M0
    fault_tag_col : str, optional
        Name of fault segment ID column (default: 'fault-tag')
    save_path : str, optional
        Path to save figure (if None, won't save)
    figsize : tuple, optional
        Figure size in inches (default: (12, 8))
    dpi : int, optional
        Resolution for saved figure (default: 150)

    Returns
    -------
    fig, axes : matplotlib figure and axes objects

    Example
    -------
    >>> # Basic usage (timing only)
    >>> df_timing = calculate_segment_timing(ssTable)
    >>> fig, axes = plot_segment_timing(
    ...     df_timing,
    ...     save_path='segment_timing.png'
    ... )
    >>> plt.show()

    >>> # With seismic moment data
    >>> df_timing = calculate_segment_timing(
    ...     ssTable,
    ...     material_file='/path/to/mat3d_fault.csv'
    ... )
    >>> fig, axes = plot_segment_timing(
    ...     df_timing,
    ...     save_path='timing_moment.png'
    ... )
    >>> plt.show()
    """
    print("\n" + "="*70)
    print("PLOTTING SEGMENT TIMING AND MOMENT")
    print("="*70 + "\n")

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Check if multi-model data and moment data
    has_model = 'model' in df_timing.columns
    has_mw = 'Mw' in df_timing.columns

    # Set title based on data content
    if has_mw:
        fig.suptitle('Fault Segment Rupture Timing and Moment Analysis',
                    fontsize=14, fontweight='bold')
    else:
        fig.suptitle('Fault Segment Rupture Timing Analysis',
                    fontsize=14, fontweight='bold')

    # Panel 1: Timeline (start and end times)
    ax = axes[0, 0]

    if has_model:
        # Group by model for color coding
        for model in df_timing['model'].unique():
            df_model = df_timing[df_timing['model'] == model]

            # Plot horizontal bars for each segment
            for i, (_, row) in enumerate(df_model.iterrows()):
                ax.barh(i, row['duration'], left=row['RT_start'],
                       alpha=0.6, label=model if i == 0 else '')
                ax.text(row['RT_end'], i, f" {row[fault_tag_col]}",
                       va='center', fontsize=8)

        ax.legend(loc='best', fontsize=8)
    else:
        # Single model or no model info
        for i, (_, row) in enumerate(df_timing.iterrows()):
            ax.barh(i, row['duration'], left=row['RT_start'], alpha=0.7)
            ax.text(row['RT_end'], i, f" seg-{row[fault_tag_col]}",
                   va='center', fontsize=8)

    ax.set_xlabel('Time (s)', fontsize=10)
    ax.set_ylabel('Fault Segment', fontsize=10)
    ax.set_title('Rupture Timeline by Segment', fontweight='bold')
    ax.grid(axis='x', alpha=0.3)

    # Panel 2: Duration histogram
    ax = axes[0, 1]
    ax.hist(df_timing['duration'], bins=20, alpha=0.7, edgecolor='black')
    ax.axvline(df_timing['duration'].mean(), color='red', linestyle='--',
              label=f"Mean: {df_timing['duration'].mean():.2f} s")
    ax.axvline(df_timing['duration'].median(), color='orange', linestyle='--',
              label=f"Median: {df_timing['duration'].median():.2f} s")
    ax.set_xlabel('Duration (s)', fontsize=10)
    ax.set_ylabel('Frequency', fontsize=10)
    ax.set_title('Rupture Duration Distribution', fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Panel 3: Scatter plot - RT_start vs duration (or Mw vs duration if available)
    ax = axes[1, 0]

    if has_mw:
        # Plot Mw vs duration if available
        if has_model:
            for model in df_timing['model'].unique():
                df_model = df_timing[df_timing['model'] == model]
                scatter = ax.scatter(df_model['duration'], df_model['Mw'],
                                   s=df_model['mean_slip']*10, alpha=0.6, label=model)
        else:
            scatter = ax.scatter(df_timing['duration'], df_timing['Mw'],
                               s=df_timing['mean_slip']*10, alpha=0.7,
                               c=df_timing[fault_tag_col], cmap='viridis')
            plt.colorbar(scatter, ax=ax, label='fault-tag')

        ax.set_xlabel('Duration (s)', fontsize=10)
        ax.set_ylabel('Moment Magnitude (Mw)', fontsize=10)
        ax.set_title('Duration vs Mw\n(bubble size = mean slip)', fontweight='bold')
    else:
        # Original plot: RT_start vs duration
        if has_model:
            for model in df_timing['model'].unique():
                df_model = df_timing[df_timing['model'] == model]
                scatter = ax.scatter(df_model['RT_start'], df_model['duration'],
                                   s=df_model['count']/10, alpha=0.6, label=model)
        else:
            scatter = ax.scatter(df_timing['RT_start'], df_timing['duration'],
                               s=df_timing['count']/10, alpha=0.7,
                               c=df_timing[fault_tag_col], cmap='viridis')
            plt.colorbar(scatter, ax=ax, label='fault-tag')

        ax.set_xlabel('RT_start (s)', fontsize=10)
        ax.set_ylabel('Duration (s)', fontsize=10)
        ax.set_title('Start Time vs Duration\n(bubble size = element count)', fontweight='bold')

    if has_model:
        ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Panel 4: Summary statistics table
    ax = axes[1, 1]
    ax.axis('off')

    # Compute summary stats
    stats_data = [
        ['Metric', 'Value'],
        ['─' * 20, '─' * 20],
        ['Total segments', f"{len(df_timing)}"],
        ['Total elements', f"{df_timing['count'].sum():,}"],
        ['', ''],
        ['Duration (s):', ''],
        ['  Min', f"{df_timing['duration'].min():.3f}"],
        ['  Max', f"{df_timing['duration'].max():.3f}"],
        ['  Mean', f"{df_timing['duration'].mean():.3f}"],
        ['  Median', f"{df_timing['duration'].median():.3f}"],
        ['', ''],
        ['RT_start (s):', ''],
        ['  Min', f"{df_timing['RT_start'].min():.3f}"],
        ['  Max', f"{df_timing['RT_start'].max():.3f}"],
        ['', ''],
        ['RT_end (s):', ''],
        ['  Min', f"{df_timing['RT_end'].min():.3f}"],
        ['  Max', f"{df_timing['RT_end'].max():.3f}"],
    ]

    # Add moment statistics if available
    if has_mw:
        stats_data.extend([
            ['', ''],
            ['Moment Magnitude:', ''],
            ['  Min Mw', f"{df_timing['Mw'].min():.2f}"],
            ['  Max Mw', f"{df_timing['Mw'].max():.2f}"],
            ['  Mean Mw', f"{df_timing['Mw'].mean():.2f}"],
            ['', ''],
            ['Mean Slip (m):', ''],
            ['  Min', f"{df_timing['mean_slip'].min():.2f}"],
            ['  Max', f"{df_timing['mean_slip'].max():.2f}"],
            ['  Mean', f"{df_timing['mean_slip'].mean():.2f}"],
        ])

    if has_model:
        stats_data.extend([
            ['', ''],
            ['Models', f"{df_timing['model'].nunique()}"],
        ])

    table = ax.table(cellText=stats_data, cellLoc='left', loc='center',
                    colWidths=[0.5, 0.5])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)

    # Style header row
    for i in range(2):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')

    ax.set_title('Summary Statistics', fontweight='bold', pad=20)

    plt.tight_layout()

    # Save if path provided
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"✓ Saved figure to: {save_path}")

    print("✓ Plot created successfully")

    return fig, axes


def plot_moment_histogram_by_model(
    df_timing: pd.DataFrame,
    fault_tag_col: str = 'fault-tag',
    model_col: str = 'model',
    moment_col: str = 'total_M0',
    bins: int = 20,
    log_scale: bool = True,
    save_path: str = None,
    figsize: tuple = (14, 6),
    dpi: int = 150,
    fault_tag_labels: dict = None,
):
    """
    Create histogram plots of seismic moment (total_M0) grouped by model,
    with each fault-tag shown in different colors.

    Parameters
    ----------
    df_timing : pd.DataFrame
        Timing and moment data from process_models_timing()
        Must have columns: model, fault-tag, total_M0
    fault_tag_col : str, optional
        Name of fault segment ID column (default: 'fault-tag')
    model_col : str, optional
        Name of model column (default: 'model')
    moment_col : str, optional
        Name of moment column to plot (default: 'total_M0')
    bins : int, optional
        Number of bins for histogram (default: 20)
    log_scale : bool, optional
        Use log scale for moment axis (default: True, recommended for M0)
    save_path : str, optional
        Path to save figure (if None, won't save)
    figsize : tuple, optional
        Figure size in inches (default: (14, 6))
    dpi : int, optional
        Resolution for saved figure (default: 150)
    fault_tag_labels : dict, optional
        Dictionary mapping fault-tag values to custom labels
        Example: {3: 'Hikurangi Interface', 65: 'Puysegur', 68: 'Alpine Fault'}
        If None, uses default "seg-{fault_tag}" labels

    Returns
    -------
    fig, axes : matplotlib figure and axes objects

    Example
    -------
    >>> df_all = process_models_timing(
    ...     model_list=['jp3z', 'jp3y'],
    ...     data_folder='/path/to/data',
    ...     material_file='mat3d_fault.csv'
    ... )
    >>> # With custom fault-tag labels
    >>> fault_labels = {3: 'Hikurangi Interface', 65: 'Puysegur', 68: 'Alpine Fault'}
    >>> fig, axes = plot_moment_histogram_by_model(
    ...     df_all,
    ...     fault_tag_labels=fault_labels,
    ...     save_path='moment_histogram.png'
    ... )
    >>> plt.show()
    """
    print("\n" + "="*70)
    print("PLOTTING MOMENT HISTOGRAM BY MODEL")
    print("="*70 + "\n")

    # Check required columns
    required_cols = [model_col, fault_tag_col, moment_col]
    missing = [col for col in required_cols if col not in df_timing.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Get unique models and fault tags
    models = sorted(df_timing[model_col].unique())
    fault_tags = sorted(df_timing[fault_tag_col].unique())
    n_models = len(models)

    print(f"Models: {n_models} ({models})")
    print(f"Fault tags: {len(fault_tags)} ({fault_tags})")
    print(f"Total records: {len(df_timing)}")

    # Create color map for fault tags
    colors = plt.cm.tab20(np.linspace(0, 1, len(fault_tags)))
    fault_tag_colors = dict(zip(fault_tags, colors))

    # Create subplots - one for each model
    fig, axes = plt.subplots(1, n_models, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    fig.suptitle(f'Seismic Moment Distribution by Model and Fault Segment',
                fontsize=14, fontweight='bold')

    # Plot histogram for each model
    for i, model in enumerate(models):
        ax = axes[i]
        df_model = df_timing[df_timing[model_col] == model]

        # Get moment values for this model
        moment_values = []
        colors_list = []
        labels_list = []

        # Group by fault-tag
        for fault_tag in fault_tags:
            df_tag = df_model[df_model[fault_tag_col] == fault_tag]
            if not df_tag.empty:
                moment_values.append(df_tag[moment_col].values)
                colors_list.append(fault_tag_colors[fault_tag])

                # Use custom label if provided, otherwise use default
                if fault_tag_labels is not None and fault_tag in fault_tag_labels:
                    label = fault_tag_labels[fault_tag]
                else:
                    label = f"seg-{fault_tag}"
                labels_list.append(label)

        if moment_values:
            # Flatten all moment values to determine bin edges
            all_moments = np.concatenate(moment_values)

            if log_scale and np.all(all_moments > 0):
                # Use log-spaced bins
                bins_edges = np.logspace(np.log10(all_moments.min()),
                                        np.log10(all_moments.max()),
                                        bins + 1)
                ax.set_xscale('log')
            else:
                bins_edges = np.linspace(all_moments.min(), all_moments.max(), bins + 1)

            # Plot stacked histogram with different colors for each fault-tag
            ax.hist(moment_values, bins=bins_edges, alpha=0.7,
                   color=colors_list, label=labels_list, stacked=False, edgecolor='black', linewidth=0.5)

            # Statistics
            total_M0 = df_model[moment_col].sum()
            mean_M0 = df_model[moment_col].mean()

            ax.set_xlabel(f'{moment_col} (N·m)', fontsize=10)
            ax.set_ylabel('Frequency', fontsize=10)
            ax.set_title(f'{model}\nTotal M₀={total_M0:.2e} N·m', fontweight='bold', fontsize=11)
            ax.grid(alpha=0.3, axis='y')

            # Add legend (only for first plot to avoid clutter)
            if i == 0:
                ax.legend(fontsize=8, loc='upper right', ncol=1)

        else:
            ax.text(0.5, 0.5, f'No data for {model}',
                   ha='center', va='center', transform=ax.transAxes)

    # Hide extra subplots if number of models < subplots
    for j in range(n_models, len(axes)):
        axes[j].axis('off')

    plt.tight_layout()

    # Save if path provided
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"✓ Saved figure to: {save_path}")

    print("✓ Histogram plot created successfully\n")

    # Print summary statistics
    print("Summary by model:")
    print("-" * 70)
    for model in models:
        df_model = df_timing[df_timing[model_col] == model]
        total_M0 = df_model[moment_col].sum()
        mean_M0 = df_model[moment_col].mean()
        n_segments = df_model[fault_tag_col].nunique()
        print(f"{model:15s}: Total M₀={total_M0:.3e} N·m, Mean={mean_M0:.3e} N·m, Segments={n_segments}")

    return fig, axes


def plot_variable_histogram(
    df_timing: pd.DataFrame,
    variable: str,
    group_by: str = None,
    bins: int = 30,
    log_scale: bool = False,
    density: bool = False,
    alpha: float = 0.7,
    title: str = None,
    xlabel: str = None,
    save_path: str = None,
    figsize: tuple = (7, 5),
    dpi: int = 300,
    group_labels: dict = None,
):
    """
    Plot histogram for a single variable, optionally grouped by model or fault-tag.

    This is a flexible function that can plot histograms for any numerical column
    (e.g., duration, Mw, mean_slip, RT_start, total_M0, etc.) with optional grouping.

    Parameters
    ----------
    df_timing : pd.DataFrame
        Timing and moment data from process_models_timing()
    variable : str
        Name of the column to plot (e.g., 'duration', 'Mw', 'mean_slip', 'total_M0')
    group_by : str, optional
        Column to group by for color coding (e.g., 'model', 'fault-tag')
        If None, plots a single histogram
    bins : int, optional
        Number of bins for histogram (default: 30)
    log_scale : bool, optional
        Use log scale for x-axis (default: False)
    density : bool, optional
        If True, plot probability density instead of counts (default: False)
    alpha : float, optional
        Transparency of histogram bars (default: 0.7)
    title : str, optional
        Plot title (if None, auto-generated from variable name)
    xlabel : str, optional
        X-axis label (if None, uses variable name)
    save_path : str, optional
        Path to save figure (if None, won't save)
    figsize : tuple, optional
        Figure size in inches (default: (10, 6))
    dpi : int, optional
        Resolution for saved figure (default: 150)
    group_labels : dict, optional
        Dictionary mapping group values to custom labels
        Example: {3: 'Hikurangi Interface', 65: 'Puysegur', 68: 'Alpine Fault'}
        If None, uses default labels

    Returns
    -------
    fig, ax : matplotlib figure and axes objects

    Example
    -------
    >>> # Simple histogram of duration
    >>> fig, ax = plot_variable_histogram(
    ...     df_all_timing,
    ...     variable='duration',
    ...     bins=30
    ... )

    >>> # Histogram of Mw grouped by model
    >>> fig, ax = plot_variable_histogram(
    ...     df_all_timing,
    ...     variable='Mw',
    ...     group_by='model',
    ...     bins=20,
    ...     save_path='Mw_histogram.png'
    ... )

    >>> # Histogram of mean_slip grouped by fault-tag with custom labels
    >>> fault_labels = {3: 'Hikurangi Interface', 65: 'Puysegur', 68: 'Alpine Fault'}
    >>> fig, ax = plot_variable_histogram(
    ...     df_all_timing,
    ...     variable='mean_slip',
    ...     group_by='fault-tag',
    ...     group_labels=fault_labels,
    ...     log_scale=True
    ... )
    """
    print("\n" + "="*70)
    print(f"PLOTTING HISTOGRAM: {variable}")
    print("="*70 + "\n")

    # Check if variable exists
    if variable not in df_timing.columns:
        raise ValueError(f"Column '{variable}' not found in DataFrame. "
                        f"Available columns: {df_timing.columns.tolist()}")

    # Remove NaN values
    df_plot = df_timing[df_timing[variable].notna()].copy()
    n_original = len(df_timing)
    n_valid = len(df_plot)

    if n_valid == 0:
        raise ValueError(f"No valid (non-NaN) values found for '{variable}'")

    if n_valid < n_original:
        print(f"⚠ Removed {n_original - n_valid} NaN values ({100*(n_original-n_valid)/n_original:.1f}%)")

    print(f"Plotting {n_valid} values")
    print(f"Range: [{df_plot[variable].min():.3g}, {df_plot[variable].max():.3g}]")
    print(f"Mean: {df_plot[variable].mean():.3g}")
    print(f"Median: {df_plot[variable].median():.3g}")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Set title and labels
    if title is None:
        if group_by:
            title = f'Distribution of {variable} (grouped by {group_by})'
        else:
            title = f'Distribution of {variable}'

    if xlabel is None:
        xlabel = variable

    # Plot histogram
    if group_by is None:
        # Simple histogram without grouping
        ax.hist(df_plot[variable], bins=bins, alpha=alpha,
               edgecolor='black', density=density, color='steelblue')

        # Add mean and median lines
        ax.axvline(df_plot[variable].mean(), color='red', linestyle='--',
                  linewidth=2, label=f'Mean: {df_plot[variable].mean():.3g}')
        ax.axvline(df_plot[variable].median(), color='orange', linestyle='--',
                  linewidth=2, label=f'Median: {df_plot[variable].median():.3g}')
        ax.legend(fontsize=10)

    else:
        # Grouped histogram
        if group_by not in df_plot.columns:
            raise ValueError(f"Group column '{group_by}' not found in DataFrame")

        groups = sorted(df_plot[group_by].unique())
        print(f"\nGrouping by '{group_by}': {len(groups)} groups")

        # Create color map
        if len(groups) <= 10:
            colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))
        elif len(groups) <= 20:
            colors = plt.cm.tab20(np.linspace(0, 1, len(groups)))
        else:
            colors = plt.cm.viridis(np.linspace(0, 1, len(groups)))

        # Plot each group
        data_by_group = []
        labels_by_group = []
        colors_by_group = []

        for i, group in enumerate(groups):
            df_group = df_plot[df_plot[group_by] == group]
            data_by_group.append(df_group[variable].values)

            # Use custom label if provided, otherwise use default
            if group_labels is not None and group in group_labels:
                label = f"{group_labels[group]} (n={len(df_group)})"
            else:
                label = f"{group} (n={len(df_group)})"

            labels_by_group.append(label)
            colors_by_group.append(colors[i])

            # Print statistics
            display_name = group_labels.get(group, group) if group_labels else group
            print(f"  {display_name}: n={len(df_group)}, "
                  f"mean={df_group[variable].mean():.3g}, "
                  f"median={df_group[variable].median():.3g}")

        # Plot stacked or overlapping histogram
        ax.hist(data_by_group, bins=bins, alpha=alpha, label=labels_by_group,
               color=colors_by_group, edgecolor='black', density=density)

        ax.legend(fontsize=9, loc='best')

    # Apply log scale if requested
    if log_scale:
        ax.set_xscale('log')
        print(f"\nUsing log scale for x-axis")

    # Labels and formatting
    ax.set_xlabel(xlabel, fontsize=11, fontweight='bold')
    if density:
        ax.set_ylabel('Probability Density', fontsize=11, fontweight='bold')
    else:
        ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')

    ax.set_title(title, fontsize=13, fontweight='bold', pad=15)
    ax.grid(alpha=0.3, axis='y')

    plt.tight_layout()

    # Save if path provided
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"\n✓ Saved figure to: {save_path}")

    print("\n✓ Histogram plot created successfully")

    return fig, ax
