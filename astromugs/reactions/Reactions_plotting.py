# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 14:37:45 2026

@author: fxmey
"""

import numpy as np
import pandas as pd
from scipy.interpolate import griddata

def _generate_2d_grid(df, value_column, log_scale=False):
    """Internal helper to interpolate unstructured grid data onto a regular 2D domain (R, Z)."""
    R = df['radius_AU'].values
    Z = df['position_AU'].values
    V = df[value_column].values

    if log_scale:
        V = np.log10(V + 1e-30) # Prevent log10(0) singularities

    # Construct a regular analytical mesh grid for contour mapping
    xi = np.linspace(R.min(), R.max(), 200)
    yi = np.linspace(Z.min(), Z.max(), 200)
    xi, yi = np.meshgrid(xi, yi)

    # Perform linear bivariate interpolation
    zi = griddata((R, Z), V, (xi, yi), method='linear')
    return xi, yi, zi

def print_summary(df_reactions):
    """
    Prints a detailed, non-graphical scientific summary of the unified kinetics database,
    accounting for the R-dependence of the vertical coordinate Z (aspect ratio Z/R).
    """
    if df_reactions.empty:
        print("Warning: Cannot generate summary. The input DataFrame is empty.")
        return

    # Isolate unique spatial nodes to understand the structural grid grid
    df_grid = df_reactions.drop_duplicates(subset=['radius_AU', 'position_AU']).copy()
    
    # Calculate the physical aspect ratio (Z/R) since Z is heavily dependent on R
    df_grid['Z_over_R'] = df_grid['position_AU'] / df_grid['radius_AU']

    print("====================== CHEMICAL NETWORK SUMMARY ======================")
    print(f"Target Species Detected     : {df_reactions['species'].unique().tolist()}")
    print(f"Total Logged Database Rows  : {len(df_reactions):,}")
    print("-" * 70)
    
    print("RADIAL DOMAIN:")
    r_min, r_max = df_grid['radius_AU'].min(), df_grid['radius_AU'].max()
    print(f"  * Radial Range (R)        : {r_min:.1f} to {r_max:.1f} AU")
    print(f"  * Unique Radii Count      : {df_grid['radius_AU'].nunique()} grid points")
    print("-" * 70)
    
    print("VERTICAL GEOMETRY & ASPECT RATIO (Z/R):")
    print(f"  * Absolute Height (Z)     : {df_grid['position_AU'].min():.3f} to {df_grid['position_AU'].max():.2f} AU")
    
    # Extract boundaries at the inner and outer disk to show the geometric flaring
    z_over_r_min = df_grid['Z_over_R'].min()
    z_over_r_max = df_grid['Z_over_R'].max()
    print(f"  * Global Aspect Ratio Z/R : {z_over_r_min:.3f} (Midplane/Inner) to {z_over_r_max:.3f} (Atmosphere/Outer)")
    
    # Local scaling checks for verification
    for r_local in sorted(df_grid['radius_AU'].unique())[:3]:  # Check first 3 radii as samples
        sub_r = df_grid[df_grid['radius_AU'] == r_local]
        print(f"    -> At R = {r_local:3d} AU : Z max = {sub_r['position_AU'].max():.2f} AU (Max Z/R = {sub_r['Z_over_R'].max():.3f})")
    if df_grid['radius_AU'].nunique() > 3:
        print("    ... [Truncated for brevity] ...")
        
    print("-" * 70)
    print("AVERAGE THERMOPHYSICAL CONDITIONS CONTRAST:")
    if 'phys_gas_temp_K' in df_reactions.columns:
        print(f"  * Gas Temperature Range   : {df_reactions['phys_gas_temp_K'].min():.1f} K to {df_reactions['phys_gas_temp_K'].max():.1f} K")
    if 'phys_gas_density' in df_reactions.columns:
        print(f"  * Gas Density Range       : {df_reactions['phys_gas_density'].min():.2e} to {df_reactions['phys_gas_density'].max():.2e} cm^-3")
    if 'phys_visual_extinction_Av' in df_reactions.columns:
        print(f"  * Visual Extinction (Av)  : {df_reactions['phys_visual_extinction_Av'].min():.2f} to {df_reactions['phys_visual_extinction_Av'].max():.2f} mag")
        
    print("-" * 70)
    print("TEMPORAL PROFILE:")
    unique_times = np.sort(df_reactions['time_years'].unique())
    print(f"  * Time Steps Count        : {len(unique_times)} steps")
    print(f"  * Time Domain Boundaries  : {unique_times[0]:.1e} to {unique_times[-1]:.1e} years")
    print("======================================================================")

def get_exact_data(df_reactions, radius_AU, position_AU, time_years, species):
    """Extracts local physical values and explicit reaction profiles from the closest analytical node."""
    df_sub = df_reactions[df_reactions["species"] == species]
    if df_sub.empty:
        return f"Error: Species '{species}' missing from dataset."

    unique_times = df_sub["time_years"].unique()
    closest_time = unique_times[np.abs(unique_times - time_years).argmin()]
    df_time_filtered = df_sub[df_sub["time_years"] == closest_time]

    spatial_distances = np.sqrt(
        (df_time_filtered["radius_AU"] - radius_AU) ** 2
        + (df_time_filtered["position_AU"] - position_AU) ** 2
    )
    closest_idx = spatial_distances.idxmin()
    chosen_point = df_time_filtered.loc[closest_idx]

    r_actual = chosen_point["radius_AU"]
    z_actual = chosen_point["position_AU"]
    sp_actual = chosen_point["spatial_point"]
    t_actual = chosen_point["time_years"]

    reac_sub = df_reactions[
        (df_reactions["radius_AU"] == r_actual)
        & (df_reactions["position_AU"] == z_actual)
        & (df_reactions["time_years"] == t_actual)
        & (df_reactions["species"] == species)
    ]

    reac_grouped = (
        reac_sub.groupby(["reaction_type", "type", "reaction"])[
            ["rate", "percent"]
        ]
        .sum()
        .reset_index()
    )

    prod = reac_grouped[reac_grouped["reaction_type"] == "production"].sort_values(
        by="percent", ascending=False
    )
    dest = reac_grouped[reac_grouped["reaction_type"] == "destruction"].sort_values(
        by="percent", ascending=False
    )

    output = [
        "=" * 60,
        f"    CHEMICAL & PHYSICAL DATA REPORT FOR {species}",
        "=" * 60,
        f"Query Coordinates : R={radius_AU} AU, Z={position_AU} AU, t={time_years} yrs",
        f"Closest Grid Point : R={r_actual:.2f} AU, Z={z_actual:.2f} AU (Index #{int(sp_actual)})",
        f"Actual Time        : t={t_actual:.1f} yrs",
        "-" * 60,
        "PHYSICAL CONDITIONS AT THIS POINT:",
    ]

    for col in chosen_point.index:
        if col.startswith("phys_"):
            output.append(
                f"  * {col.replace('phys_', ''):<22} : {chosen_point[col]}"
            )

    output.extend(["-" * 60, f"PRODUCTION REACTIONS ({len(prod)} pathways) :"])
    for _, row in prod.iterrows():
        output.append(
            f"  [{row['type']:<2}] {row['reaction']:<30} | Rate: {row['rate']:.4e} | {row['percent']:>5}%"
        )

    output.extend(["-" * 60, f"DESTRUCTION REACTIONS ({len(dest)} pathways) :"])
    for _, row in dest.iterrows():
        output.append(
            f"  [{row['type']:<2}] {row['reaction']:<30} | Rate: {row['rate']:.4e} | {row['percent']:>5}%"
        )

    output.append("=" * 60)

    return "\n".join(output)

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_reaction_evolution(df_reactions, species, rmin=None, rmax=None, z_over_r_min=None, z_over_r_max=None, tmin=None, tmax=None, top_n=5, figsize=(16, 6)):
    """
    Plots the temporal evolution of the spatially integrated average reaction rates
    using true grid volume weighting and rigorous temporal dt-weighting to handle
    non-uniform spatial meshes and logarithmic time steps.

    Parameters
    ----------
    df_reactions : pd.DataFrame
        The unified kinetics database containing physical properties, spatial nodes,
        and individual chemical reaction rates.
    species : str
        Target molecular species to analyze (e.g., 'CO', 'H2O').
    rmin : float, optional
        Minimum radial boundary in Astronomical Units (AU). If None, defaults 
        to the minimum radius available in the dataset.
    rmax : float, optional
        Maximum radial boundary in Astronomical Units (AU). If None, defaults 
        to the maximum radius available in the dataset.
    z_over_r_min : float, optional
        Minimum physical aspect ratio (Z/R). If None, defaults to the lowest 
        ratio found within the dataset footprint.
    z_over_r_max : float, optional
        Maximum physical aspect ratio (Z/R). If None, defaults to the highest 
        ratio found within the dataset footprint.
    tmin : float, optional
        Minimum temporal boundary in years. If None, defaults to the earliest
        logged time step.
    tmax : float, optional
        Maximum temporal boundary in years. If None, defaults to the latest
        logged time step.
    top_n : int, optional
        Number of leading pathways to isolate and plot per reaction type. Default is 5.
    figsize : tuple of (float, float), optional
        Dimensions (width, height) of the resulting matplotlib figure. Default is (16, 6).

    Returns
    -------
    None
        Displays a side-by-side subplot tracking production and destruction kinetics over time.
    """
    # Constant: 1 Astronomical Unit in centimeters (cm)
    AU_TO_CM = 1.495978707e13

    # 1. Local copy and aspect ratio calculation (Z/R)
    df = df_reactions[df_reactions['species'] == species].copy()
    if df.empty:
        print(f"Warning: No reaction pathways found for species: {species}")
        return
        
    df['Z_over_R'] = df['position_AU'] / df['radius_AU']

    # 2. Evaluate radial integration boundaries (R)
    r_actual_min, r_actual_max = df['radius_AU'].min(), df['radius_AU'].max()
    r_sel_min = rmin if (rmin is not None and rmin >= r_actual_min) else r_actual_min
    r_sel_max = rmax if (rmax is not None and rmax <= r_actual_max) else r_actual_max
    
    # 3. Evaluate aspect ratio integration boundaries (Z/R)
    z_over_r_actual_min = df['Z_over_R'].min()
    z_over_r_actual_max = df['Z_over_R'].max()
    
    if z_over_r_min == 0 and z_over_r_max == 0:
        z_sel_min, z_sel_max = 0.0, 0.0
    else:
        z_sel_min = z_over_r_min if (z_over_r_min is not None and z_over_r_min >= z_over_r_actual_min) else z_over_r_actual_min
        z_sel_max = z_over_r_max if (z_over_r_max is not None and z_over_r_max <= z_over_r_actual_max) else z_over_r_actual_max

    # 4. Generate spatial masking array
    if z_sel_min == 0.0 and z_sel_max == 0.0:
        spatial_mask = (df['radius_AU'] >= r_sel_min) & (df['radius_AU'] <= r_sel_max) & (df['position_AU'] == 0.0)
        loc_summary = f"Midplane (Z=0), R in [{r_sel_min:.1f}, {r_sel_max:.1f}] AU"
    else:
        spatial_mask = (df['radius_AU'] >= r_sel_min) & (df['radius_AU'] <= r_sel_max) & \
                       (df['Z_over_R'] >= z_sel_min) & (df['Z_over_R'] <= z_sel_max)
        loc_summary = f"R in [{r_sel_min:.1f}, {r_sel_max:.1f}] AU, Z/R in [{z_sel_min:.3f}, {z_sel_max:.3f}]"

    df_filtered = df[spatial_mask].copy()
    if df_filtered.empty:
        print(f"Warning: Spatial integration boundaries returned an empty grid selection.")
        return

    # 5. Spatial Differential Weight Calculation (True Numerical Integration)
    grid_nodes = df_filtered.drop_duplicates(subset=['radius_AU', 'position_AU']).copy()
    
    # Resolve dR elements (mean grid distance to closest neighbors)
    radii = np.sort(grid_nodes['radius_AU'].unique())
    if len(radii) > 1:
        dr_dict = {r: (radii[min(idx+1, len(radii)-1)] - radii[max(idx-1, 0)]) / 2.0 for idx, r in enumerate(radii)}
    else:
        dr_dict = {radii[0]: 1.0}
        
    grid_nodes['dR'] = grid_nodes['radius_AU'].map(dr_dict)
    
    # Resolve dZ vertical elements per localized radial column
    grid_nodes['dZ'] = 1.0
    for r_val in radii:
        sub_z = grid_nodes[grid_nodes['radius_AU'] == r_val].sort_values(by='position_AU')
        z_vals = sub_z['position_AU'].values
        if len(z_vals) > 1:
            dz_vals = (np.roll(z_vals, -1) - np.roll(z_vals, 1)) / 2.0
            dz_vals[0] = z_vals[1] - z_vals[0]
            dz_vals[-1] = z_vals[-1] - z_vals[-2]
            grid_nodes.loc[grid_nodes['radius_AU'] == r_val, 'dZ'] = dz_vals

    # Convert all structural dimensions from AU to cm before computing volume
    r_cm = grid_nodes['radius_AU'].values * AU_TO_CM
    dr_cm = grid_nodes['dR'].values * AU_TO_CM
    dz_cm = grid_nodes['dZ'].values * AU_TO_CM

    # Cylindrical volume element natively resolved in cm^3 (dV = 2 * pi * R * dR * dZ)
    grid_nodes['dV'] = 2.0 * np.pi * r_cm * dr_cm * dz_cm
    total_volume = grid_nodes['dV'].sum()
    
    # Merge volume structures back to weight raw chemical rates
    df_filtered = df_filtered.merge(grid_nodes[['radius_AU', 'position_AU', 'dV']], on=['radius_AU', 'position_AU'], how='left')
    df_filtered['weighted_rate'] = df_filtered['rate'] * df_filtered['dV']

    # 6. Sum weighted rates and normalize across the integration space volume
    df_integrated = df_filtered.groupby(['time_years', 'reaction_type', 'reaction'])['weighted_rate'].sum().reset_index()
    df_integrated['rate'] = df_integrated['weighted_rate'] / total_volume

    # 7. Apply Temporal Boundaries Filtering (tmin, tmax)
    t_actual_min, t_actual_max = df_integrated['time_years'].min(), df_integrated['time_years'].max()
    t_sel_min = tmin if (tmin is not None and tmin >= t_actual_min) else t_actual_min
    t_sel_max = tmax if (tmax is not None and tmax <= t_actual_max) else t_actual_max
    
    df_integrated = df_integrated[(df_integrated['time_years'] >= t_sel_min) & (df_integrated['time_years'] <= t_sel_max)]
    
    if df_integrated.empty:
        print(f"Warning: Temporal boundaries [{t_sel_min:.1e}, {t_sel_max:.1e}] years returned an empty dataset.")
        return

    # 8. Rank leading pathways based on true temporal integration (dt weighting)
    time_steps = np.sort(df_integrated['time_years'].unique())
    if len(time_steps) > 1:
        dt_dict = {}
        for idx, t in enumerate(time_steps):
            if idx == 0:
                dt_dict[t] = time_steps[1] - time_steps[0]
            elif idx == len(time_steps) - 1:
                dt_dict[t] = time_steps[-1] - time_steps[-2]
            else:
                dt_dict[t] = (time_steps[idx+1] - time_steps[idx-1]) / 2.0
    else:
        dt_dict = {time_steps[0]: 1.0}
        
    df_integrated['dt'] = df_integrated['time_years'].map(dt_dict)
    df_integrated['time_weighted_rate'] = df_integrated['rate'] * df_integrated['dt']
    total_duration = sum(dt_dict.values())

    # The ranking is determined by the cumulative area under the curve in the selected time frame
    df_ranking = df_integrated.groupby(['reaction_type', 'reaction'])['time_weighted_rate'].sum().reset_index()
    df_ranking['rate'] = df_ranking['time_weighted_rate'] / total_duration

    # 9. Initialize Subplot Structure
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=True)
    
    for i, rtype in enumerate(['production', 'destruction']):
        ax = axes[i]
        
        rank_sub = df_ranking[df_ranking['reaction_type'] == rtype]
        top_reactions = rank_sub.sort_values(by='rate', ascending=False).head(top_n)['reaction'].tolist()
        
        sub_time = df_integrated[(df_integrated['reaction_type'] == rtype) & (df_integrated['reaction'].isin(top_reactions))]
        if sub_time.empty:
            ax.text(0.5, 0.5, f'No {rtype} pathways logged in this domain', ha='center', va='center', fontsize=12)
            continue
            
        plot_data = sub_time.pivot(index='time_years', columns='reaction', values='rate')
        sorted_columns = rank_sub[rank_sub['reaction'].isin(top_reactions)].sort_values(by='rate', ascending=False)['reaction'].tolist()
        plot_data = plot_data[sorted_columns]
        
        cmap_name = 'Greens_r' if rtype == 'production' else 'Reds_r'
        colors = mpl.colormaps[cmap_name](np.linspace(0.1, 0.6, len(sorted_columns)))

        for col_idx, reaction in enumerate(sorted_columns):
            if reaction in plot_data.columns:
                ax.plot(plot_data.index, plot_data[reaction], label=reaction, 
                        linewidth=2.5, color=colors[col_idx], marker='o', markersize=4, alpha=0.9)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(f'Top {top_n} Leading - {rtype.capitalize()}', fontsize=13, fontweight='bold')
        ax.set_xlabel('Time (years)', fontsize=11, fontweight='bold')
        ax.set_ylabel(r'Volume-Integrated Rate ($cm^{-3} \cdot s^{-1}$)', fontsize=11, fontweight='bold')
        ax.legend(loc='best', fontsize=9, frameon=True, facecolor='white', edgecolor='none')
        ax.grid(True, which="both", ls="--", alpha=0.5)

    plt.suptitle(f"Temporal Evolution of Chemical Pathways for: {species}\nIntegrated Domain: {loc_summary} | Time: [{t_sel_min:.1e}, {t_sel_max:.1e}] years", 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.show()
    
def plot_physical_map(df_reactions, param_names='phys_gas_temp_K', z_over_r_max=None, log_scale=False):
    """
    Generates continuous 2D physical contour maps (R vs Z) from unique spatial nodes.

    Filters the multi-dimensional protoplanetary disk grid using optional vertical 
    aspect ratio thresholds (Z/R) and interpolates the localized numerical 
    attributes to construct smooth physical profiles.

    Parameters
    ----------
    df_reactions : pd.DataFrame
        The unified chemical kinetics database containing columns for radial positions,
        vertical heights, and associated physical parameters.
    param_names : str or list of str, optional
        A single parameter column name or a list of column names to plot.
        Default is 'phys_gas_temp_K'.
    z_over_r_max : float, optional
        Maximum vertical aspect ratio threshold (Z/R). Grid points exceeding 
        this value are pruned to clip upper atmospheric layers. If None, the 
        full simulation domain is displayed.
    log_scale : bool, optional
        If True, applies a base-10 logarithmic transformation to the target 
        physical attributes prior to 2D grid interpolation. Default is False.

    Returns
    -------
    None
        Displays the resulting side-by-side or stacked grid layout of 2D contour maps.
    """
    # 1. Standardize parameter arguments into iterable lists
    if isinstance(param_names, str):
        param_names = [param_names]
    
    # 2. Isolate physical conditions based on unique coordinates, a single time step,
    # and a single species to avoid averaging artifacts caused by structural duplication
    df_grid = df_reactions.drop_duplicates(subset=['radius_AU', 'position_AU', 'time_years', 'species']).copy()
    df_grid = df_grid.drop_duplicates(subset=['radius_AU', 'position_AU'])

    # 3. Apply vertical aspect ratio clipping if a user threshold is specified
    if z_over_r_max is not None:
        df_grid['Z_over_R'] = df_grid['position_AU'] / df_grid['radius_AU']
        df_grid = df_grid[df_grid['Z_over_R'] <= z_over_r_max]
        if df_grid.empty:
            print(f"Error: Vertical boundary filter Z/R <= {z_over_r_max} returned an empty grid selection.")
            return

    # 4. Helper function to translate internal column headers into publication-grade labels
    def clean_name(param):
        if param == "phys_gas_temp_K": return "Gas temperature [K]"
        elif param == 'phys_gas_density': return r"Gas density [$cm^{-3}$]"
        elif param == 'phys_visual_extinction_Av': return r"Visual extinction $A_v$ [mag]"
        elif param == 'phys_x_rate': return 'X-ray ionization rate'
        elif param == 'phys_dust_temp_K': return 'Dust temperature [K]'
        return param
    
    # 5. Sanity check for column key existence within the processed DataFrame
    valid_params = [p for p in param_names if p in df_grid.columns]
    if not valid_params:
        print("Error: No valid physical parameters detected.")
        return

    # 6. Resolve subplot structural layout matrix (2 columns maximum per row)
    num_plots = len(valid_params)
    n_cols = 2 if num_plots > 1 else 1
    n_rows = (num_plots + 1) // 2
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10 if n_cols == 1 else 16, 6 * n_rows), squeeze=False)
    axes_flat = axes.flatten()

    # 7. Core plotting iteration loop across validated parameters
    for i, param_name in enumerate(valid_params):
        ax = axes_flat[i]
        
        # Call the underlying meshgrid generator function to compute continuous 2D matrices
        xi, yi, zi = _generate_2d_grid(df_grid, param_name, log_scale)

        # Dynamically assign colormaps based on physical properties
        # Viridis is used for exponential/logarithmic parameters; Inferno for thermal structures
        cmap = 'viridis' if 'density' in param_name or 'Av' in param_name else 'inferno'
        contour = ax.contourf(xi, yi, zi, levels=50, cmap=cmap)
        
        # Colorbar attachment and scientific labeling
        cbar = fig.colorbar(contour, ax=ax)
        cbar_title = f"Log10({clean_name(param_name)})" if log_scale else clean_name(param_name)
        cbar.set_label(cbar_title, fontsize=11, fontweight='bold')

        # Formatting axes styles and text annotations
        ax.set_xlabel('Radius R (AU)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Height Z (AU)', fontsize=11, fontweight='bold')
        
        filter_str = f" (Z/R <= {z_over_r_max:.2f})" if z_over_r_max is not None else ""
        ax.set_title(f'2D Physical Map: {clean_name(param_name)}{filter_str}', fontsize=12, fontweight='bold')

    # 8. Clean up unused subplot axes panels from the canvas grid
    for j in range(num_plots, len(axes_flat)):
        fig.delaxes(axes_flat[j])

    plt.tight_layout()
    plt.show()
    
    
def plot_chemical_maps(df_reactions, species, time_years=None, z_over_r_max=None, average_time=False, log_scale=True, figsize=None, cmaps=None):
    """
    Renders 2D spatial contour panels mapping localized chemical production, 
    destruction, and net residuals across disk coordinates.

    Processes multi-row panels matching separate rows to different molecular 
    species, evaluates discrete snapshots or cumulative temporal averages, and 
    applies aspect ratio bounds to focus analysis on specific disk zones.

    Parameters
    ----------
    df_reactions : pd.DataFrame
        The unified kinetics database containing columns for radial positions,
        vertical heights, reaction classifications, and raw kinetics rates.
    species : str or list of str
        The molecular target name or list of target names to display. Each 
        species initializes a distinct horizontal panel row.
    time_years : float, optional
        Target time checkpoint in years. If None, the function defaults to 
        the final available time step in the simulation.
    z_over_r_max : float, optional
        Maximum vertical aspect ratio threshold (Z/R). Grid points exceeding 
        this value are trimmed to truncate higher atmospheric disk altitudes.
        If None, the complete vertical column is retained.
    average_time : bool, optional
        If True, aggregates and means all cumulative kinetic states from t=0 
        up to the selected checkpoint. If False, plots an instantaneous 
        snapshot. Default is False.
    log_scale : bool, optional
        If True, applies a base-10 logarithmic transformation to the raw 
        production and destruction rates prior to spatial interpolation. 
        Residual balances remain linear. Default is True.
    figsize : tuple of (float, float), optional
        Dimensions (width, height) of the total canvas output. If None, 
        scales dynamically with the number of rows.
    cmaps : dict, optional
        Custom colormap keys to override visual defaults. Expected structure:
        {'maps': 'colormap_name_1', 'residual': 'colormap_name_2'}.

    Returns
    -------
    None
        Displays a multi-panel grid layout mapping kinetic distributions.
    """
    # 1. Standardize species argument formats into explicit iterable vectors
    species_list = [species] if isinstance(species, str) else list(species)
    num_rows = len(species_list)

    # 2. Assign default sequential and diverging color tables
    default_cmaps = {'maps': 'viridis', 'residual': 'coolwarm'}
    if cmaps is not None:
        default_cmaps.update(cmaps)

    # 3. Handle default dynamic figure sizing allocations
    if figsize is None:
        figsize = (18, 5 * num_rows)

    # 4. Initialize grid axis system
    fig, axes = plt.subplots(num_rows, 3, figsize=figsize, sharex=True, sharey=True)
    if num_rows == 1:
        axes = np.expand_dims(axes, axis=0)

    # 5. Primary vertical loop over the target species entries
    for row_idx, sp in enumerate(species_list):
        sub_sp = df_reactions[df_reactions['species'] == sp].copy()
        if sub_sp.empty:
            print(f"Warning: Species {sp} not found in the dataset.")
            return

        # Apply aspect ratio bounds filtering if requested
        if z_over_r_max is not None:
            sub_sp['Z_over_R'] = sub_sp['position_AU'] / sub_sp['radius_AU']
            sub_sp = sub_sp[sub_sp['Z_over_R'] <= z_over_r_max]
            if sub_sp.empty:
                print(f"Error: Boundary boundary Z/R <= {z_over_r_max} returned an empty matrix row.")
                return

        # 6. Resolve chronological check point boundaries
        unique_times = sub_sp['time_years'].unique()
        selected_time = np.max(unique_times) if time_years is None else unique_times[np.abs(unique_times - time_years).argmin()]
        
        # 7. Isolate networks and condense kinetics based on temporal arguments
        if average_time:
            sub_sp_filtered = sub_sp[sub_sp['time_years'] <= selected_time]
            title_suffix = f"Avg up to {selected_time:.1e} yr"
            
            prod_by_time = sub_sp_filtered[sub_sp_filtered['reaction_type'] == 'production']
            dest_by_time = sub_sp_filtered[sub_sp_filtered['reaction_type'] == 'destruction']
            
            # First sum separate reactions at each individual time node step
            prod_totals = prod_by_time.groupby(['time_years', 'radius_AU', 'position_AU'])['rate'].sum().reset_index()
            dest_totals = dest_by_time.groupby(['time_years', 'radius_AU', 'position_AU'])['rate'].sum().reset_index()
            
            # Mean the condensed snapshots over the cumulative time vector span
            prod_totals = prod_totals.groupby(['radius_AU', 'position_AU'])['rate'].mean().reset_index()
            dest_totals = dest_totals.groupby(['radius_AU', 'position_AU'])['rate'].mean().reset_index()
        else:
            sub_sp_filtered = sub_sp[sub_sp['time_years'] == selected_time]
            title_suffix = f"t={selected_time:.1e} yr"
            
            sub_prod = sub_sp_filtered[sub_sp_filtered['reaction_type'] == 'production']
            sub_dest = sub_sp_filtered[sub_sp_filtered['reaction_type'] == 'destruction']
            
            prod_totals = sub_prod.groupby(['radius_AU', 'position_AU'])['rate'].sum().reset_index()
            dest_totals = sub_dest.groupby(['radius_AU', 'position_AU'])['rate'].sum().reset_index()

        # 8. Merge directional pathways into balanced structural nodes
        merged = pd.merge(prod_totals, dest_totals, on=['radius_AU', 'position_AU'], how='outer', suffixes=('_prod', '_dest')).fillna(0.0)
        merged['residual'] = merged['rate_prod'] - merged['rate_dest']

        # 9. Interpolate point grids to construct smooth 2D mesh arrays
        xi, yi, zi_prod = _generate_2d_grid(merged, 'rate_prod', log_scale)
        xi, yi, zi_dest = _generate_2d_grid(merged, 'rate_dest', log_scale)
        xi, yi, zi_res = _generate_2d_grid(merged, 'residual', log_scale=False)

        # 10. Horizontal rendering iteration over Production, Destruction, and Residuals
        plot_configs = [
            ('Production', zi_prod, default_cmaps['maps']),
            ('Destruction', zi_dest, default_cmaps['maps']),
            ('Net Residual', zi_res, default_cmaps['residual'])
        ]

        for col_idx, (title, grid_data, cmap) in enumerate(plot_configs):
            ax = axes[row_idx, col_idx]
            
            # Avoid crashing if the grid is populated by NaNs or empty blocks
            max_val = np.nanmax(np.abs(grid_data)) if not np.all(np.isnan(grid_data)) else 1.0
            
            # Diverging colormaps require centered levels symmetrically balancing negative/positive regimes
            levels = np.linspace(-max_val, max_val, 51) if col_idx == 2 else 50

            contour = ax.contourf(xi, yi, grid_data, levels=levels, cmap=cmap, extend='both')
            cbar = fig.colorbar(contour, ax=ax, pad=0.02)
            
            # Scientific formatting for colorbars
            if col_idx == 2:
                cbar.set_label(r'Net Rate ($cm^{-3} \cdot s^{-1}$)', fontsize=10, fontweight='bold')
            else:
                unit = r'$\log_{10}(cm^{-3} \cdot s^{-1})$' if log_scale else r'$cm^{-3} \cdot s^{-1}$'
                cbar.set_label(f'Rate ({unit})', fontsize=10, fontweight='bold')

            # Customize tick marks for professional journal standards
            ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=10)
            ax.minorticks_on()

            # Apply labels exclusively to shared perimeter margins
            if row_idx == num_rows - 1:
                ax.set_xlabel('Radius R (AU)', fontsize=11, fontweight='bold')
            if col_idx == 0:
                ax.set_ylabel('Height Z (AU)', fontsize=11, fontweight='bold')
            
            filter_str = f"\n(Z/R <= {z_over_r_max:.2f})" if z_over_r_max is not None else ""
            ax.set_title(f"{sp} - {title}\n({title_suffix}){filter_str}", fontsize=11, fontweight='bold')

    plt.tight_layout()
    plt.show()
    
    
def plot_spatial_trends(df_reactions, species, radius_AU=None, z_over_r_max=None, log_scale=True, figsize=(16, 6)):
    """
    Plots the spatial profiles of chemical rates (production/destruction) and their 
    net residuals as a function of the disk aspect ratio (Z/R) for a specific radius or the global disk.

    Parameters
    ----------
    df_reactions : pd.DataFrame
        The unified kinetics database containing columns for radial positions,
        vertical heights, reaction classifications, and raw kinetics rates.
    species : str or list of str
        Target molecular species to analyze and filter (e.g., 'CO', 'H2O').
    radius_AU : int or float, optional
        Target radial distance in Astronomical Units (AU) to isolate a single 
        vertical column. If None, averages profiles across all available radii.
    z_over_r_max : float, optional
        Maximum vertical aspect ratio threshold (Z/R). Grid points exceeding 
        this value are pruned to clip upper atmospheric layers. If None, the 
        entire vertical domain is parsed.
    log_scale : bool, optional
        If True, applies a base-10 logarithmic scaling to both axes panels.
        If False, uses standard linear axes coordinates. Default is True.
    figsize : tuple of (float, float), optional
        Dimensions (width, height) of the resulting matplotlib figure. Default is (16, 6).

    Returns
    -------
    None
        Displays a side-by-side subplot tracking directional rates and net residuals.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt

    if df_reactions.empty:
        print("Warning: Empty database. Cannot resolve spatial pathways.")
        return

    # 1. Standardize species input into an iterable list
    species_list = [species] if isinstance(species, str) else list(species)

    # 2. Local copy and target species filtering
    df = df_reactions[df_reactions['species'].isin(species_list)].copy()
    if df.empty:
        print(f"Warning: None of the requested species {species_list} were found in the dataset.")
        return

    # 3. Apply radial distance filtering if a specific radius is requested
    if radius_AU is not None:
        # Match the closest available radius integer/float in the database to avoid empty sets
        available_radii = df['radius_AU'].unique()
        closest_r = available_radii[np.abs(available_radii - radius_AU).argmin()]
        df = df[df['radius_AU'] == closest_r]
        radius_summary = f"at R = {closest_r} AU"
        print(f"-> Radial column isolation active. Closest grid radius matched: {closest_r} AU")
    else:
        radius_summary = "Global Disk (Radii-Averaged)"

    # 4. Calculate aspect ratio (Z/R)
    df['Z_over_R'] = df['position_AU'] / df['radius_AU']

    # 5. Apply vertical aspect ratio clipping threshold if requested
    if z_over_r_max is not None:
        df = df[df['Z_over_R'] <= z_over_r_max]
        if df.empty:
            print(f"Error: Boundary boundary Z/R <= {z_over_r_max} returned an empty subset.")
            return

    # 6. Round coordinates to create clean, discrete spatial bins for line smoothing
    df['Z_over_R'] = np.round(df['Z_over_R'], 3)

    # 7. Group to compute directional total rates (Production and Destruction separately)
    # We now group by time and unique spatial signatures to collapse individual reactions first
    df_directional = df.groupby(['species', 'Z_over_R', 'time_years', 'reaction_type'])['rate'].sum().reset_index()
    # Then we average out the temporal variations to extract the core spatial structure
    df_directional = df_directional.groupby(['species', 'Z_over_R', 'reaction_type'])['rate'].mean().reset_index()

    # 8. Compute net residuals (Production - Destruction) by merging subsets
    prod_sub = df_directional[df_directional['reaction_type'] == 'production'].rename(columns={'rate': 'rate_prod'})
    dest_sub = df_directional[df_directional['reaction_type'] == 'destruction'].rename(columns={'rate': 'rate_dest'})
    
    df_residual = pd.merge(prod_sub, dest_sub, on=['species', 'Z_over_R'], how='outer').fillna(0.0)
    df_residual['residual'] = df_residual['rate_prod'] - df_residual['rate_dest']

    # 9. Initialize Side-by-Side Subplots
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    sns.set_theme(style="whitegrid")
    
    filter_suffix = f" (Z/R <= {z_over_r_max:.2f})" if z_over_r_max is not None else ""

    # PANEL 1: Directional Trends (Production vs Destruction)
    sns.lineplot(
        data=df_directional, x='Z_over_R', y='rate', 
        hue='species', style='reaction_type', 
        markers=True, dashes=True, linewidth=2.5, markersize=6, alpha=0.85, ax=axes[0]
    )
    
    if log_scale:
        axes[0].set_xscale('log')
        axes[0].set_yscale('log')
        y_label_prefix = "Log10 Mean"
    else:
        y_label_prefix = "Mean"
        
    axes[0].set_title(f'Directional Rates Profile{filter_suffix}', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Disk Aspect Ratio (Z/R)', fontsize=11, fontweight='bold')
    axes[0].set_ylabel(f'{y_label_prefix} Chemical Rate ($cm^{-3} . s^{-1}$)', fontsize=11, fontweight='bold')
    axes[0].legend(title="Pathways", title_fontproperties={'weight': 'bold'}, loc='best')

    # PANEL 2: Net Residuals (Production - Destruction balance)
    if log_scale:
        df_residual['abs_residual'] = df_residual['residual'].abs()
        df_residual['regime'] = np.where(df_residual['residual'] >= 0, 'Net Production (+)', 'Net Destruction (-)')
        
        df_residual_plot = df_residual[df_residual['abs_residual'] > 0]
        
        sns.lineplot(
            data=df_residual_plot, x='Z_over_R', y='abs_residual', 
            hue='species', style='regime', 
            markers=True, dashes=True, linewidth=2.5, markersize=6, alpha=0.85, ax=axes[1]
        )
        axes[1].set_xscale('log')
        axes[1].set_yscale('log')
        y_res_label = r'Net Residual Value ($|Prod - Dest|$ - $cm^{-3} \cdot s^{-1}$)'
        legend_title = "Net Balance Regime"
    else:
        sns.lineplot(
            data=df_residual, x='Z_over_R', y='residual', 
            hue='species', markers=True, linewidth=2.5, markersize=6, alpha=0.85, ax=axes[1]
        )
        axes[1].axhline(0.0, color='black', linestyle='--', alpha=0.5, linewidth=1.5)
        y_res_label = r'Net Residual Value ($Prod - Dest$ - $cm^{-3} \cdot s^{-1}$)'
        legend_title = "Chemical Species"

    axes[1].set_title(f'Net Residuals Balance{filter_suffix}', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Disk Aspect Ratio (Z/R)', fontsize=11, fontweight='bold')
    axes[1].set_ylabel(y_res_label, fontsize=11, fontweight='bold')
    axes[1].legend(title=legend_title, title_fontproperties={'weight': 'bold'}, loc='best')

    plt.suptitle(f"Spatial Kinetic Profiles Profile: {radius_summary}", fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.show()
    
    
def plot_chemical_equilibrium(df_reactions, species, rmin=None, rmax=None, z_over_r_min=None, z_over_r_max=None, figsize=(10, 6)):
    """
    Plots the absolute total production vs destruction rates over time to diagnose
    whether the selected disk region and specific species have reached chemical equilibrium.
    Ensures missing reaction types are explicitly treated as zero instead of being omitted.

    Parameters
    ----------
    df_reactions : pd.DataFrame
        The unified kinetics database.
    species : str or list of str
        Target molecular species to analyze (e.g., 'CO', or ['CO', 'CN']).
    rmin, rmax : float, optional
        Radial boundaries in AU.
    z_over_r_min, z_over_r_max : float, optional
        Aspect ratio boundaries.
    figsize : tuple, optional
        Figure dimensions. Default is (10, 6).
    """
    import matplotlib.pyplot as plt
    
    # Constant: 1 AU in cm
    AU_TO_CM = 1.495978707e13
    # Absolute physical floor to handle strict zeros on a logarithmic scale smoothly
    LOG_FLOOR = 1e-30

    # 1. Standardize species input format into an explicit iterable list
    species_list = [species] if isinstance(species, str) else list(species)

    # 2. Filter database by the targeted vector of species
    df = df_reactions[df_reactions['species'].isin(species_list)].copy()
    if df.empty:
        print(f"Warning: No data found for requested species: {species_list}")
        return
        
    df['Z_over_R'] = df['position_AU'] / df['radius_AU']

    # 3. Evaluate spatial boundaries
    r_act_min, r_act_max = df['radius_AU'].min(), df['radius_AU'].max()
    r_sel_min = rmin if (rmin is not None and rmin >= r_act_min) else r_act_min
    r_sel_max = rmax if (rmax is not None and rmax <= r_act_max) else r_act_max
    
    if z_over_r_min == 0 and z_over_r_max == 0:
        z_sel_min, z_sel_max = 0.0, 0.0
    else:
        z_sel_min = z_over_r_min if (z_over_r_min is not None) else df['Z_over_R'].min()
        z_sel_max = z_over_r_max if (z_over_r_max is not None) else df['Z_over_R'].max()

    # 4. Apply spatial masking array
    if z_sel_min == 0.0 and z_sel_max == 0.0:
        spatial_mask = (df['radius_AU'] >= r_sel_min) & (df['radius_AU'] <= r_sel_max) & (df['position_AU'] == 0.0)
        loc_summary = f"Midplane (Z=0), R in [{r_sel_min:.1f}, {r_sel_max:.1f}] AU"
    else:
        spatial_mask = (df['radius_AU'] >= r_sel_min) & (df['radius_AU'] <= r_sel_max) & \
                       (df['Z_over_R'] >= z_sel_min) & (df['Z_over_R'] <= z_sel_max)
        loc_summary = f"R in [{r_sel_min:.1f}, {r_sel_max:.1f}] AU, Z/R in [{z_sel_min:.3f}, {z_sel_max:.3f}]"

    df_filtered = df[spatial_mask].copy()
    if df_filtered.empty:
        print("Warning: Selected spatial integration region returned an empty grid subset.")
        return

    # 5. Calculate grid volume elements (dV) for true numerical integration
    grid_nodes = df_filtered.drop_duplicates(subset=['radius_AU', 'position_AU']).copy()
    radii = np.sort(grid_nodes['radius_AU'].unique())
    
    dr_dict = {r: (radii[min(idx+1, len(radii)-1)] - radii[max(idx-1, 0)]) / 2.0 for idx, r in enumerate(radii)} if len(radii) > 1 else {radii[0]: 1.0}
    grid_nodes['dR'] = grid_nodes['radius_AU'].map(dr_dict)
    
    grid_nodes['dZ'] = 1.0
    for r_val in radii:
        sub_z = grid_nodes[grid_nodes['radius_AU'] == r_val].sort_values(by='position_AU')
        z_vals = sub_z['position_AU'].values
        if len(z_vals) > 1:
            dz_vals = (np.roll(z_vals, -1) - np.roll(z_vals, 1)) / 2.0
            dz_vals[0] = z_vals[1] - z_vals[0]
            dz_vals[-1] = z_vals[-1] - z_vals[-2]
            grid_nodes.loc[grid_nodes['radius_AU'] == r_val, 'dZ'] = dz_vals

    grid_nodes['dV'] = 2.0 * np.pi * (grid_nodes['radius_AU'] * AU_TO_CM) * (grid_nodes['dR'] * AU_TO_CM) * (grid_nodes['dZ'] * AU_TO_CM)
    total_volume = grid_nodes['dV'].sum()
    
    df_filtered = df_filtered.merge(grid_nodes[['radius_AU', 'position_AU', 'dV']], on=['radius_AU', 'position_AU'], how='left')
    df_filtered['weighted_rate'] = df_filtered['rate'] * df_filtered['dV']

    # 6. Sum weighted rates across the integration space volume
    df_tot = df_filtered.groupby(['time_years', 'species', 'reaction_type'])['weighted_rate'].sum().reset_index()
    df_tot['rate'] = df_tot['weighted_rate'] / total_volume

    # ==========================================================================
    # CRITICAL FIX: REINDEX MATRIX TO FORCE EXPLICIT ZEROS FOR MISSING DATA
    # ==========================================================================
    # Pivot to align production and destruction on the same row per timestamp
    df_pivot = df_tot.pivot(index=['time_years', 'species'], columns='reaction_type', values='rate').reset_index()
    
    # Force column existence if one type is completely missing from the whole dataset
    if 'production' not in df_pivot.columns: df_pivot['production'] = 0.0
    if 'destruction' not in df_pivot.columns: df_pivot['destruction'] = 0.0
    
    # Fill actual missing timestamps with 0.0
    df_pivot['production'] = df_pivot['production'].fillna(0.0)
    df_pivot['destruction'] = df_pivot['destruction'].fillna(0.0)
    
    # Apply the logarithmic floor to handle zeros safely on the log plot
    df_pivot['production_clip'] = np.where(df_pivot['production'] <= 0.0, LOG_FLOOR, df_pivot['production'])
    df_pivot['destruction_clip'] = np.where(df_pivot['destruction'] <= 0.0, LOG_FLOOR, df_pivot['destruction'])
    # ==========================================================================

    # 7. Render Equilibrium Plots using synchronized coloring
    plt.figure(figsize=figsize)
    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    color_palette = prop_cycle.by_key()['color']

    for idx, sp in enumerate(species_list):
        df_sp = df_pivot[df_pivot['species'] == sp].sort_values(by='time_years')
        if df_sp.empty:
            continue
            
        current_color = color_palette[idx % len(color_palette)]
        
        # Plot production (solid line) and destruction (dashed line)
        plt.plot(df_sp['time_years'], df_sp['production_clip'], color=current_color, 
                 linestyle='-', linewidth=2.5, label=f'{sp} - Production', alpha=0.85)
        plt.plot(df_sp['time_years'], df_sp['destruction_clip'], color=current_color, 
                 linestyle='--', linewidth=2.5, label=f'{sp} - Destruction', alpha=0.85)
    
    plt.xscale('log')
    plt.yscale('log')
    
    # Adjust plot limits slightly above the floor to keep the plot clean
    y_min_data = min(df_pivot['production'].max(), df_pivot['destruction'].max())
    if y_min_data > LOG_FLOOR and not np.all(df_pivot['production'] == 0):
        # Set dynamic bottom limit based on the lowest real non-zero value found
        valid_rates = pd.concat([df_pivot['production'], df_pivot['destruction']])
        real_min = valid_rates[valid_rates > 0.0].min()
        plt.ylim(bottom=max(real_min * 0.1, LOG_FLOOR * 10))

    plt.title(f"Chemical Equilibrium Diagnostic Matrix\nDomain: {loc_summary}", fontsize=13, fontweight='bold')
    plt.xlabel('Time (years)', fontsize=11, fontweight='bold')
    plt.ylabel(r'Volume-Integrated Total Rate ($cm^{-3} \cdot s^{-1}$)', fontsize=11, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Kinetic Tracks", title_fontproperties={'weight': 'bold'})
    plt.grid(True, which="both", ls="--", alpha=0.5)
    
    plt.tight_layout()
    plt.show()