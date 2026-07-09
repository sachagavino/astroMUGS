# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 13:12:33 2026

@author: fxmey
"""

from pathlib import Path
import numpy as np
import pickle
import re
import os
import time
import subprocess
from contextlib import contextmanager
import pandas as pd

import astromugs as mugs
import astromugs.pipeline as pipeline
import astromugs.plotting.plot as mplt

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def parse_nmgc_output(output_text, radius):
    """
    Parses NMGC/Nautilus text output into structured data blocks, including the radius.

    Args:
        output_text (str): Raw text output from the NMGC simulator.
        radius (int|float): The radial distance (AU) to add to each block.

    Returns:
        list[dict]: A list of structured dictionaries with 'radius_AU' included.
    """
    blocks = []
    current_block = None
    current_section = None
    
    header_re = re.compile(
        r"For\s+(?P<species>[A-Za-z0-9+-]+)\s+at\s+output\s+n°(?P<out_num>\d+)\s+"
        r"\((?P<time>[^ ]+)\s+years\)\s+and\s+spatial\s+point\s+n°(?P<pt_num>\d+)\s+"
        r"\((?P<pos>[^ ]+)\s+AU\)"
    )
    
    phys1_re = re.compile(r"Gas density\s*=\s*(?P<dens>[^ ]+).*?Av\s*=\s*(?P<av>[^ ]+).*?X rate\s*=\s*(?P<x_rate>[^ ]+)")
    phys2_re = re.compile(r"Gas temp\s*=\s*(?P<g_temp>[^ ]+).*?Dust temp\s*=\s*(?P<d_temp>[^ ]+)")
    reaction_re = re.compile(r"^\s*(?P<type>\d+)\s+(?P<reaction>.+?)\s+(?P<rate>[0-9E.+-]+)\s+(?P<pct>[0-9.]+)%")

    lines = output_text.splitlines()
    
    for line in lines:
        line_stripped = line.strip()
        
        header_match = header_re.search(line_stripped)
        if header_match:
            if current_block:
                blocks.append(current_block)
                
            current_section = None
            current_block = {
                "radius_AU": radius,
                "species": header_match.group("species"),
                "output_number": int(header_match.group("out_num")),
                "time_years": float(header_match.group("time")),
                "spatial_point": int(header_match.group("pt_num")),
                "position_AU": float(header_match.group("pos")),
                "physical_conditions": {},
                "production_reactions": [],
                "destruction_reactions": []
            }
            continue
            
        if not current_block:
            continue
            
        p1_match = phys1_re.search(line_stripped)
        if p1_match:
            current_block["physical_conditions"].update({
                "gas_density": float(p1_match.group("dens")),
                "visual_extinction_Av": float(p1_match.group("av")),
                "x_rate": float(p1_match.group("x_rate"))
            })
            continue
            
        p2_match = phys2_re.search(line_stripped)
        if p2_match:
            current_block["physical_conditions"].update({
                "gas_temp_K": float(p2_match.group("g_temp")),
                "dust_temp_K": float(p2_match.group("d_temp"))
            })
            continue
            
        if "PRODUCTION" in line_stripped:
            current_section = "PRODUCTION"
            continue
        elif "DESTRUCTION" in line_stripped:
            current_section = "DESTRUCTION"
            continue
        elif "change:" in line_stripped or "Please enter" in line_stripped:
            current_section = None
            continue
            
        if current_section:
            reac_match = reaction_re.search(line_stripped)
            if reac_match:
                rate_str = reac_match.group("rate").strip()
                if 'e' not in rate_str.lower():
                    if '-' in rate_str[1:]:
                        rate_str = rate_str[0] + rate_str[1:].replace('-', 'e-', 1)
                    elif '+' in rate_str[1:]:
                        rate_str = rate_str[0] + rate_str[1:].replace('+', 'e+', 1)

                reac_data = {
                    "type": int(reac_match.group("type")),
                    "reaction": reac_match.group("reaction").strip(),
                    "rate": float(rate_str),
                    "percent": float(reac_match.group("pct"))
                }
                if current_section == "PRODUCTION":
                    current_block["production_reactions"].append(reac_data)
                elif current_section == "DESTRUCTION":
                    current_block["destruction_reactions"].append(reac_data)

    if current_block:
        blocks.append(current_block)
        
    return blocks

def run_nmgc_multi(radius, chemistry_path, output_time=1, reactions="major", molecules="CO", spatial_points=1):
    """
    Executes NMGC for a combination of molecules, spatial points, and times.
    """
    if isinstance(molecules, str):
        molecules = [molecules]
    elif isinstance(molecules, (list, np.ndarray)):
        molecules = list(molecules)
    else:
        raise ValueError("Invalid format for 'molecules' (expected str, list, or array).")
    
    if isinstance(spatial_points, (int, float, list, tuple, np.ndarray)):
        spatial_points = [int(p) for p in np.atleast_1d(spatial_points)]
    else:
        raise ValueError("Invalid format for 'spatial_points' (expected number, list, or array).")

    if isinstance(output_time, (int, float, str)):
        output_times = [output_time]
    else:
        output_times = list(output_time)

    if reactions == "major":
        reac = 0
    else: 
        reac = 1
        if reactions != "all": 
            print(f'Unknown `reactions` parameter {reactions}. Must be "all" or "major". Defaulting to "all".')

    seq = f"{reac}\n" 
    seq += f"{output_times[0]}\n"
    seq += f"{molecules[0]}\n"
    seq += f"{spatial_points[0]}\n"
    
    is_first = True
    for t in output_times:
        for m in molecules:
            for p in spatial_points:
                if is_first:
                    is_first = False
                    continue 
                
                seq += "a\n"
                seq += f"{t}\n"
                seq += f"{m}\n"
                seq += f"{p}\n"

    seq += "q\n"
                
    radius = int(radius)
    target_dir = Path(chemistry_path) / f'{radius}AU' 
    
    start_time = time.time()
    
    process = subprocess.Popen(
        "nmgc major",
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=target_dir
    )

    raw_output, _ = process.communicate(input=seq)
    
    final_results = parse_nmgc_output(raw_output, radius=radius)
    
    check = check_combinations(final_results, radius, molecules, output_times, spatial_points)
    
    if len(check['missing']) > 0 or len(check['surplus']) > 0: 
        print("="*50, "MISSING COMBINATIONS", "="*50)
        print(check['missing'])
        print("="*50, "SURPLUS COMBINATIONS", "="*50)
        print(check['surplus'])
        raise ValueError("Missing or surplus combinations detected between input parameters and simulation output.")

    end_time = time.time()
    print(f" `nmgc major` run completed in {round(end_time - start_time, 2)} s")

    return final_results

def check_combinations(data, radius, molec_list, time_list, height_list):
    """
    Checks for discrepancies between expected and actual (radius, molecule, time, height) combinations.
    """
    expected_combinations = {
        (int(radius), m, float(t), h) 
        for m in molec_list 
        for t in time_list 
        for h in height_list
    }
    
    actual_combinations = set()
    for entry in data:
        r = entry.get('radius_AU', radius)
        specie = entry.get('species')
        time_val = float(entry.get('output_number')) if entry.get('output_number') is not None else None
        height = entry.get('spatial_point') 
        
        actual_combinations.add((int(r), specie, time_val, height))
        
    missing = expected_combinations - actual_combinations
    surplus = actual_combinations - expected_combinations
    
    keys = ['radius_AU', 'species', 'output_time', 'height']
    
    formatted_missing = [dict(zip(keys, c)) for c in missing]
    formatted_surplus = [dict(zip(keys, c)) for c in surplus]
    
    return {
        "surplus": formatted_surplus,
        "missing": formatted_missing
    }

def reactions(chemistry_path, thermal_path, molecules, r_lim=None, t_lim=None, z_arg="all", reac_type="major", write_file=True, overwrite=True):
    """
    Extracts, compresses, and aggregates 2D chemical reaction data across spatial 
    (R, Z) and temporal dimensions, employing incremental disk storage to optimize 
    memory allocation.

    Parameters
    ----------
    chemistry_path : str or Path
        Path to the primary chemical network simulation directory.
    thermal_path : str or Path
        Path to the thermal structure data directory.
    molecules : list of str
        Target molecular species to filter and extract from the raw dataset.
    r_lim : tuple of (float, float), optional
        Radial boundaries (R_min, R_max) in Astronomical Units (AU). If None, 
        all detected radial grid points are processed.
    t_lim : tuple of (float, float), optional
        Temporal boundaries (T_min, T_max) in years. If None, all time steps 
        are computed.
    z_arg : str, optional
        Vertical profile selection profile. Set to 'midplane' to restrict extraction 
        to the dense midplane boundary (Z=0 equivalent), or 'all' for full vertical column.
    reac_type : str, optional
        Classification of reaction filtering (e.g., 'major' or 'all'). Default is 'major'.
    write_file : bool, optional
        If True, serializes the accumulated data chunk-by-chunk into a compressed 
        Parquet file. Default is True.
    overwrite : bool, optional
        If True, clears pre-existing cached Parquet files. If False, bypasses execution 
        and reads the cached binary format.

    Returns
    -------
    pd.DataFrame
        A highly compressed, structured DataFrame containing the accumulated chemical 
        kinetics data across the filtered multi-dimensional space.
    """
    
    def load_z_coordinates_numpy(file_path):
        """Extracts the vertical column heights (Z-grid) using NumPy, discarding comment blocks."""
        return np.loadtxt(file_path, comments='!', usecols=0).tolist()

    path = Path(chemistry_path)
    thermal_path = Path(thermal_path)

    output_file = path / f"{reac_type}_reactions.parquet"

    # --- Cache and File I/O Management ---
    if output_file.is_file():
        if overwrite:
            if write_file:
                output_file.unlink()
                print(f"Existing file removed: {output_file}")
            else:
                print(f"Warning: overwrite=True but write_file=False. {output_file} preserved.")
        else:
            print(f"The file already exists in {output_file} (Loading cached Parquet data...)")
            return pd.read_parquet(output_file)
            
    # 1. Radial Grid Discovery and Lexical Sorting
    r_list = [d.name[:-2] for d in path.iterdir() if d.is_dir() and len(d.name) in (3, 4, 5) and d.name[:-2].isdigit() and d.name[-2:].isalpha()]
    r_list = [int(r) for r in r_list]
    r_list.sort()
    
    # 2. Radial Range Domain Filtering
    if r_lim is not None: 
        if len(r_lim) == 2:
            r_array = np.array(r_list)
            r_list = r_array[(r_array >= r_lim[0]) & (r_array <= r_lim[1])].tolist()
        else: 
            raise ValueError("r_lim must be a tuple (r_min, r_max) or None")

    writer = None
    
    try:
        # 3. Main Computational Loop Over the Radial Space Continuum
        for i, r in enumerate(r_list, 1):
            print("="*20, f"{i}/{len(r_list)}", f"| r = {r} AU", "="*20)

            pipe = pipeline.Interface()
            pipe.add_chemical_path(str(path))
            pipe.add_thermal_path(str(thermal_path))
            pipe.add_chemistry()
            
            time_array = np.array(pipe.chemistry[r]['time'])
            time_index = np.arange(1, len(time_array) + 1)
            
            # 4. Temporal Range Filtering
            if t_lim is not None: 
                if len(t_lim) == 2:
                    t_min_sec = t_lim[0] * 365.25 * 86400
                    t_max_sec = t_lim[1] * 365.25 * 86400
                    mask = (time_array >= t_min_sec) & (time_array <= t_max_sec)
                    time_index = time_index[mask].tolist()
                else: 
                    raise ValueError("t_lim must be a tuple (t_min, t_max) or None")
            else:
                time_index = time_index.tolist()
                    
            # 5. Vertical Coordinate Grid Resolution
            z_points = load_z_coordinates_numpy(path / f"{r}AU/1D_static.dat")
            if z_arg == "midplane": 
                z_list = [1] 
            else: 
                z_list = list(range(1, len(z_points) + 1))
                
            # 6. Execute Subprocess Simulation Engine (NMGC Solver)
            res = run_nmgc_multi(r, str(path), output_time=time_index, reactions=reac_type, 
                                 molecules=molecules, spatial_points=z_list)
            
            if not res:
                continue
            
            # --- FLATTENING LAYER: Unpack nested dictionary blocks on the fly ---
            flattened_reactions = []
            entries = res if isinstance(res, list) else [res]
            
            for entry in entries:
                base = {
                    'radius_AU': entry.get('radius_AU', None),
                    'position_AU': entry.get('position_AU', None),
                    'species': entry.get('species', 'Unknown'),
                    'output_number': entry.get('output_number', None),
                    'time_years': entry.get('time_years', None),
                    'spatial_point': entry.get('spatial_point', None),
                }
                
                if base['radius_AU'] is None or base['position_AU'] is None:
                    continue

                phys_cond = entry.get('physical_conditions', {})
                for k, v in phys_cond.items():
                    base[f'phys_{k}'] = v
                
                # Expand and stamp production pathways
                for rxn in entry.get('production_reactions', []):
                    flattened_reactions.append({
                        **base,
                        'reaction_type': 'production',
                        'type': rxn.get('type', None),
                        'reaction': rxn.get('reaction', 'Unknown'),
                        'rate': rxn.get('rate', 0.0),
                        'percent': rxn.get('percent', 0.0)
                    })
                    
                # Expand and stamp destruction pathways
                for rxn in entry.get('destruction_reactions', []):
                    flattened_reactions.append({
                        **base,
                        'reaction_type': 'destruction',
                        'type': rxn.get('type', None),
                        'reaction': rxn.get('reaction', 'Unknown'),
                        'rate': rxn.get('rate', 0.0),
                        'percent': rxn.get('percent', 0.0)
                    })

            if not flattened_reactions:
                continue

            # Convert newly flattened structure into a pandas chunk DataFrame
            df_chunk = pd.DataFrame(flattened_reactions)
            
            # Enforce a strict, immutable schema for the streaming pipeline.
            # This prevents type fluctuations (e.g., int8 vs int16) between spatial chunks.
            explicit_types = {
                'radius_AU': 'int32',
                'spatial_point': 'int32',
                'output_number': 'int32',
                'type': 'int32',
                'position_AU': 'float64',
                'time_years': 'float64',
                'rate': 'float64',
                'percent': 'float64'
            }
            
            # Apply explicitly defined types only to columns existing in the chunk
            columns_to_cast = {col: dtype for col, dtype in explicit_types.items() if col in df_chunk.columns}
            df_chunk = df_chunk.astype(columns_to_cast)
    
            # 7. Incremental Parquet Stream Storage
            if write_file:
                import pyarrow as pa
                import pyarrow.parquet as pq
                
                table = pa.Table.from_pandas(df_chunk, preserve_index=False)
                if writer is None:
                    writer = pq.ParquetWriter(str(output_file), table.schema, compression='snappy')
                writer.write_table(table)
            
            del res
            del df_chunk

    finally:
        if writer is not None:
            writer.close()

    # --- Data Delivery and Final Sequence ---
    if write_file and output_file.is_file():
        print(f"Chemical data compression completed safely. Global database written to: {output_file}")
        return pd.read_parquet(output_file)
    else:
        print("Warning: Data aggregation executed without continuous disk backup.")
        return pd.DataFrame()