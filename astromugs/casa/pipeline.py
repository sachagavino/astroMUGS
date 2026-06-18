# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:42:35 2026

@author: fxmey
"""

import astromugs.pipeline as pipeline

from casatasks import simobserve, importfits, concat, tclean, exportfits

import radmc3dPy as r3d

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import re
import time
import shutil
import subprocess
import os
import warnings
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.patches import Ellipse

@contextmanager
def silence(verbose=False):
    """Mute all print, errors and warnings if verbose = False"""
    if not verbose:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with open(os.devnull, 'w') as devnull:
                with redirect_stdout(devnull), redirect_stderr(devnull):
                    yield
    else:
        yield
        
def image_name(command, path_lines=""):
    """Takes a radmc3d subprocess command. Returns a corresponding image name"""
    with open(path_lines + "lines.inp", "r") as f:
        mols = [l.split()[0] for l in f.read().strip().split("\n")[2:]]
    idx = re.search(r"imolspec\s+(\d+)", command)
    sub = mols[int(idx.group(1)) - 1] if idx else ""
    clean = re.sub(r"^radmc3d\s+", "", command)
    return "_".join(re.sub(r"imolspec\s+\d+", sub, clean).split())

@contextmanager
def cd(newdir):
    """ Allows to run a cript in a particular directory. 
    Use it like :
    with cd("your/path/"):
        your script
    """
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def writing_input_files(chemistry_path,thermal_path,MOL,datadict,verbose=True):
    """ Writing input files :
     - numberdens_XX.inp
     - gas_velocity.inp
     - lines.inp
     - gas_temperature.inp"""

    chemistry_path = Path(chemistry_path)
    thermal_path = str(Path(thermal_path))
    
    with silence(verbose):
        pipe = pipeline.Interface()
        pipe.add_chemical_path(chemistry_path)  # path to the chemistry directory
        pipe.add_thermal_path(thermal_path)  # path to thermal directory
        pipe.add_chemistry()
        
        # WRITING numberdens_XX.inp FILES
        for molec in MOL:
            pipe.convert_nautilus2radmc(species=molec, numberdens=True)
            
            # WRITING INPUT FILES
        pipe.write_line(control=False, line=True, gasvelocity=True,gastemp='1D_static', species=MOLEC, star_mass=datadict['m_star'])
    if not verbose : print("Done without error")


def find_alma_transitions(selected_bands, therm_path):
    therm_path = Path(therm_path)

    with open(therm_path / "ALMA.dat") as f:
        bands = {int(p[0]): (float(p[4]), float(p[5])) for l in f if not l.startswith('#') and (p := l.split())}
    
    selected_bands = [selected_bands] if isinstance(selected_bands, int) else selected_bands
    results = {}
    
    
    search_pattern = str(therm_path / "molecule_*.inp")
    for file_str in glob.glob(search_pattern):
        file_path = Path(file_str)
        mol = file_path.stem.replace("molecule_", "")
        
        results[mol] = []
        with open(file_path) as f:
            lines = [l.strip() for l in f if l.strip()]
            
        levels = {int(p[0]): p[3] for l in lines if len(p := l.split()) >= 4 and p[0].isdigit() and '_' in p[3]}
        for l in lines:
            p = l.split()
            if len(p) == 6 and p[0].isdigit() and '.' in p[4]: 
                t_id, up, low, _, freq = int(p[0]), int(p[1]), int(p[2]), p[3], float(p[4])
                
                for b in selected_bands:
                    if b in bands and bands[b][0] <= freq <= bands[b][1]:
                        results[mol].append({
                            "bande": b, "id": t_id, "lambdaµm": 1e6*299792458/(freq*1e9), 
                            "trans": f"{levels.get(up, up)} -> {levels.get(low, low)}"
                        })
    return results

def check_and_select_transition(parsed_results, molecule, target_freq_ghz, epsilon_ghz=0.5, target_band=None):
    """
    Checks if a specific molecule has a transition near the target frequency 
    within a given epsilon (tolerance) using the parsed results from find_alma_transitions.
    Optionally filters by ALMA band.
    """
    # Normalize input string to lowercase to avoid naming mismatch
    mol_key = molecule.strip()
    
    print(f"=== Searching for {molecule.upper()} around {target_freq_ghz} GHz (Tol: +/- {epsilon_ghz} GHz) ===")
    
    # Check if the molecule exists in your parsed data
    if mol_key not in parsed_results or not parsed_results[mol_key]:
        print(parsed_results[mol_key])
        print(f"[-] No transitions found for '{molecule}' in your local molecule_*.inp files.")
        return None
        
    matched_transitions = []
    SPEED_OF_LIGHT = 299792458
    
    # Scan all transitions parsed for this molecule
    for sim_trans in parsed_results[mol_key]:
        # Convert wavelength (microns) back to frequency (GHz)
        sim_freq_ghz = (1e6 * SPEED_OF_LIGHT) / (sim_trans["lambdaµm"] * 1e9)
        diff = abs(sim_freq_ghz - target_freq_ghz)
        
        # Check if it falls within our epsilon window
        if diff <= epsilon_ghz:
            # Optional: filter by ALMA band if specified
            if target_band is not None and sim_trans["bande"] != target_band:
                continue
                
            matched_transitions.append({
                "band": sim_trans["bande"],
                "ifreq": sim_trans["id"],
                "freq_ghz": round(sim_freq_ghz, 4),
                "trans_name": sim_trans["trans"],
                "delta_ghz": round(diff, 4)
            })
            
    # Display the results
    if not matched_transitions:
        print(f"[-] No transition found for {molecule.upper()} matching these criteria.")
        return None
        
    print(f"[+] Found {len(matched_transitions)} matching transition(s):")
    for i, match in enumerate(matched_transitions, 1):
        print(f"  {i}. ID {match['ifreq']}: {match['trans_name']}")
        print(f"     -> Frequency: {match['freq_ghz']} GHz (Delta: {match['delta_ghz']} GHz)")
        print(f"     -> ALMA Band: {match['band']}")
        
    # If multiple transitions match, we select the closest one by default
    best_match = min(matched_transitions, key=lambda x: x["delta_ghz"])
    
    if len(matched_transitions) > 1:
        print(f"\n[!] Multiple matches found. Selecting the closest one (ID {best_match['ifreq']}) for your setup.")
        
    return best_match

def run_pipeline(
    available_molecules,
    interactive=True,
    alma_bands=None,
    molecule="CO",
    mode="map",  # "map" or "cube"
    vkms=8.0,
    widthkms=None,
    nlam=6,
    freq=340.0,
    tol=20.0,
    params_dict=None,  # Must contain 'inclination', 'pa' (and optionally 'dpc')
    thermal_path=".",
    index_molec=1,
    sizeau=600,
    npix=200,
    lbd=None,
    lbdrange=None
):
    """Execute the complete RADMC-3D image/cube generation pipeline.

    This function acts as a high-level wrapper orchestrating the workflow 
    between ALMA transition lookups (via `find_alma_transitions`) and RADMC-3D 
    radiative transfer subprocess executions. It natively supports a dual 
    physical framework: Molecular Line Radiative Transfer and Dust Continuum 
    Emission, across both interactive terminal wizards and direct scripting hooks.

    Behavior Matrix
    ---------------
    1. Interactive Mode (`interactive=True`):
       - Prompts the user to isolate the physical domain: Dust Continuum vs. Line Emission.
       - If Continuum is selected: Captures single wavelength (lambda) or a bounding 
         wavelength range (lambdarange), along with canvas dimensions (sizeau, npix).
       - If Line is selected: Interactively captures ALMA receiver bands, targets species, 
         and frequency boundaries.
       - Implements a robust validation loop (max 3 sequential attempts) for type 
         casting and mathematical edge constraints.
       - Executes automated transition queries. If no absolute match is discovered 
         within the tolerance window, a fallback manual catalog table is generated 
         for terminal selection.

    2. Non-Interactive Mode (`interactive=False`):
       - Bypasses terminal interactions, immediately evaluating keyword arguments.
       - If `lbd` or `lbdrange` are populated, immediately routes execution to the 
         Dust Continuum subprocess workflow.
       - If working in Line mode, normalizes molecular strings to uppercase, defaults 
         `alma_bands` to [6, 7] if left unassigned, and forces a deterministic 
         fallback to the first available transition index if automated frequency 
         matching fails.

    Parameters
    ----------
    available_molecules : list of str
        The full inventory of molecular species identifiers supported by the local 
        line database setup (e.g., ['CO', 'H2O', 'HCN']).
    interactive : bool, default True
        If True, initiates the interactive terminal wizard. If False, executes 
        headless automated scripting based entirely on incoming keywords.
    alma_bands : list of int, optional
        Target ALMA receiver bands to search within (integers from 1 to 10). 
        Defaults to [6, 7] in non-interactive line mode if unassigned.
    molecule : str, default "CO"
        Name of the target astronomical molecular species (case-insensitive).
    mode : {"map", "cube"}, default "map"
        Simulation configuration profile used for line imaging:
        - "map": Renders a single velocity channel map at local LSR velocity `vkms`.
        - "cube": Computes a Position-Position-Velocity (PPV) synthetic data cube.
    vkms : float, default 8.0
        The local Local Standard of Rest (LSR) velocity in km/s. Used exclusively 
        when `mode="map"`.
    widthkms : float, optional
        The total symmetric velocity framework window span (km/s). Used 
        exclusively when `mode="cube"`.
    nlam : int, default 6
        Number of discrete wavelength/velocity sampling channels. Used for both 
        the line "cube" mode and the continuum "lambdarange" framework.
    freq : float, default 340.0
        The target baseline frequency in GHz around which the transition line 
        matching algorithm will search.
    tol : float, default 20.0
        The symmetric frequency constraint tolerance window (+/- tol) in GHz.
    params_dict : dict
        Structural dict mapping orientations of the target source. Must contain:
        - params_dict["inclination"] (float/int): Disk inclination angle in degrees.
        - params_dict["pa"] (float/int): Position Angle of the major axis in degrees.
        - params_dict["dpc"] (float/int, optional): Source distance in parsecs. 
          Defaults cleanly to "N/A" if missing.
    thermal_path : str, default "."
        Filesystem path to the working directory housing local molecular data 
        files (`molecule_*.inp`), dust opacity files, and RADMC execution files.
    index_molec : int, default 1
        The tracking integer ID assigned to the molecule inside the `radmc3d.inp` 
        line definition array. Automatically recalculated in interactive line setups.
    sizeau : int, default 600
        The outer spatial scaling bound edge-to-edge canvas size in Astronomical Units (AU).
    npix : int, default 200
        The linear pixel matrix resolution grid (e.g., 200 results in a 200x200 canvas).
    lbd : float, optional
        Target single wavelength in microns (µm) for dust continuum rendering.
    lbdrange : list or tuple of float, optional
        A two-element sequence `[lambda_min, lambda_max]` defining a continuum 
        wavelength range block in microns (µm).

    Returns
    -------
    tuple (str or None, float or None)
        A tuple containing two tracked variables:
        - name (str or None): The final, dynamically updated output file name 
          following the displacement of `image.out`. Returns None if the subprocess fails.
        - widthkms (float or None): Total velocity window span utilized. Returns 
          a float if line `mode="cube"`, and None for maps or continuum runs.

    Raises
    ------
    ValueError
        - Raised if `params_dict` is omitted or resolves to None.
        - Raised if validation inputs fail cumulatively 3 times in interactive mode.
        - Raised if interactive band parsing returns empty collections.
    subprocess.CalledProcessError
        Captured and printed internally if the underlying `radmc3d` binary engine 
        aborts with a non-zero exit runtime status flag.

    Notes
    -----
    - Directory context migrations (`os.chdir`) are entirely isolated inside 
      `try...finally` runtime traps, ensuring the python execution pool safely 
      returns to the base shell working directory if RADMC-3D crashes.
    - Constant factors for micro-to-gigahertz conversions follow the physical 
      speed of light constant: `freq_ghz = 299792.458 / lambda_microns`.
    """
    
    # Check for critical orientation configuration dictionary
    if params_dict is None:
        print("Parameters dictionary is required")
        print("Terminating program.")
        return None, None

    # =========================================================================
    # 0. INTERNAL HELPER FUNCTIONS
    # =========================================================================
    def get_bands(text):
        """Extracts unique ALMA bands between 1 and 10 from text patterns."""
        return sorted(list({int(b) for b in re.findall(r"\d+", text) if 1 <= int(b) <= 10}))

    def ask_value(prompt, cast_type, validation_cond=None, err_msg=""):
        """Prompts the user for a value with up to 3 validation attempts."""
        for i in range(3):
            try:
                user_input = input(prompt).strip()
                # Boolean normalization
                if cast_type == bool:
                    if user_input.lower() in ["y", "yes"]: return True
                    if user_input.lower() in ["n", "no"]: return False
                    raise ValueError()
                # General structural evaluation
                value = cast_type(user_input)
                if validation_cond and not validation_cond(value):
                    raise ValueError()
                return value
            except ValueError:
                print(f"--> {err_msg} (Attempt {i+1}/3)")
        raise ValueError("--> Too many invalid inputs. Terminating program.")

    # =========================================================================
    # 1. PARAMETER MANAGEMENT (INTERACTIVE OR DIRECT)
    # =========================================================================
    if interactive:
        print("\n=== INTERACTIVE CONFIGURATION ===")
        
        # Determine whether execution targets Dust Continuum or Molecular Lines
        continuum_bool = ask_value(
            prompt="Do you want to use radmc3d for disk continuum (i.e. not for a spectral line) ? (y/n): ",
            cast_type=bool,
            err_msg="Please answer by y/n",
        )
        
        if continuum_bool: 
            # Interactive Dust Continuum configuration
            lbd = ask_value(
                prompt="Choose wavelength in µm. If None, you'll need to set a wavelength range \nLambda in µm ? (value/None): ",
                cast_type=lambda x: None if x.lower() == 'none' else float(x),
                err_msg="Please answer by a float or by None",
            )
            if lbd is None:
                lbdrange_min = ask_value(
                    prompt="Choose wavelength range lower bound in µm : ",
                    cast_type=float,
                    err_msg="Please answer by a number",
                )
                lbdrange_max = ask_value(
                    prompt="Choose wavelength range upper bound in µm : ",
                    cast_type=float,
                    err_msg="Please answer by a number",
                )
                
            sizeau = ask_value(
                prompt="Image size (AU): ",
                cast_type=int,
                validation_cond=lambda n: n > 0,
                err_msg="Image size must be a positive integer.",
            )
    
            npix = ask_value(
                prompt="Resolution of the image (number of pixels): ",
                cast_type=int,
                validation_cond=lambda n: n > 0,
                err_msg="npix must be a positive integer.",
            )
            
            name = None  # Prevent UnboundLocalError if subprocess crashes early
            old_cwd = os.getcwd()
            try:
                # Migrate execution to data directory context
                os.chdir(thermal_path)
                t1 = time.time()

                # Build spatial base arguments
                base_cmd = (
                    f"radmc3d image npix {npix} "
                    f"incl {params_dict['inclination']} "
                    f"posang {params_dict['pa']} "
                    f"sizeau {sizeau}"
                )

                # Route command parameters based on single wavelength vs spectral range
                if lbd is not None:
                    cmd = f"{base_cmd} lambda {lbd}"
                else:  
                    nlam = ask_value(
                        prompt="Number of wavelength channels (nlam): ",
                        cast_type=int,
                        validation_cond=lambda n: n > 0,
                        err_msg="Number of channels must be a positive integer.",
                    )
                    cmd = f"{base_cmd} lambdarange {lbdrange_min} {lbdrange_max} nlam {nlam}"

                print(f"\n--> Running terminal command: {cmd}")
                subprocess.run(cmd, shell=True, check=True)

                # Rename the generated file block using pipeline conventions
                name = image_name(cmd, path_lines=thermal_path)
                if os.path.exists("image.out"):
                    os.rename("image.out", name)
                    print(f"--> File successfully renamed to {name}")

                print(f"Process completed in {time.time() - t1:.2f} seconds.")

            except subprocess.CalledProcessError:
                print("--> Error: RADMC-3D subprocess failed.")
            finally:
                # Safely step back out to original directory context
                os.chdir(old_cwd)
                
            print("======================== PARAMETERS ========================")
            print(f'inclination = {params_dict["inclination"]}°')
            print(f'posang = {params_dict["pa"]}°')
            print(f'distance = {params_dict.get("dpc", "N/A")} pc')
            if lbd is not None: 
                print(f'wavelength (µm) = {lbd:.2f}')
            else: 
                print(f'wavelength range : {lbdrange_min:.2f} to {lbdrange_max:.2f} µm')
                print(f"Spectral resolution : {(lbdrange_max-lbdrange_min)/(nlam-1):.2e} µm")
            print("============================================================")
                
            return name, None
        
        else:
            # Interactive Molecular Line configuration
            alma_bands = get_bands(input("Choose ALMA bands (e.g., 6,7): "))
            if not alma_bands:
                print("Please choose one or more bands between 1 and 10.")
                alma_bands = get_bands(input("Choose ALMA bands: "))
                if not alma_bands:
                    raise ValueError("No valid ALMA bands selected.")
            print(f"--> Chosen bands: {alma_bands}")
    
            molecule = ask_value(
                prompt=f"Choose a molecule in {available_molecules}: ",
                cast_type=str,
                validation_cond=lambda m: m.upper() in available_molecules,
                err_msg=f"Please choose a molecule available in {available_molecules}",
            ).upper()
            index_molec = available_molecules.index(molecule) + 1
    
            freq_bool = ask_value(
                prompt=f"Do you want to choose a frequency around which the code will search a transition for {molecule}? (y/n): ",
                cast_type=bool,
                err_msg="Please answer by y/n",
            )
    
            if freq_bool:
                freq = ask_value(
                    prompt="Target frequency in GHz: ",
                    cast_type=float,
                    err_msg="Invalid number for Frequency.",
                )
                tol = ask_value(
                    prompt="Tolerance (+/- tol) in GHz: ",
                    cast_type=float,
                    err_msg="Invalid number for Tolerance.",
                )
                print(f"--> f = {freq:.2f} +/- {tol:.2f} GHz")
            else:
                if isinstance(freq, (float, int)) and isinstance(tol, (float, int)):
                    print(f"--> f = {freq:.2f} +/- {tol:.2f} GHz")
                else:
                    print("Defined values for frequency and tolerance are needed\nTerminating program.")
                    return None, None
    
            # =========================================================================
            # 2. ALMA TRANSITION SEARCH (INTERACTIVE)
            # =========================================================================
            print(f"\n=== Searching for {molecule} in band(s) {alma_bands} ===")
            all_transitions = find_alma_transitions(selected_bands=alma_bands, therm_path=thermal_path)
            selected_setup = check_and_select_transition(
                parsed_results=all_transitions, molecule=molecule, target_freq_ghz=freq, epsilon_ghz=tol
            )
    
            # Manual index matching if automatic selector turns out empty
            if selected_setup is None:
                if not all_transitions or molecule not in all_transitions:
                    print(f"No transitions found for {molecule} in ALMA bands {alma_bands}")
                    return None, None
                else:
                    print("No match within tolerance. Available transitions:")
                    for i, d in enumerate(all_transitions[molecule]):
                        if i == 10 and not ask_value("Show all transitions? (y/n): ", bool, err_msg="y/n"):
                            break
                        print(f"ID: {i} | Band: {d['bande']} | Freq: {299792.458/d['lambdaµm']:.2f} GHz")
    
                    max_id = len(all_transitions[molecule]) - 1
                    chosen_id = ask_value(
                        prompt=f"Choose a transition ID (0 to {max_id}): ",
                        cast_type=int,
                        validation_cond=lambda n: 0 <= n <= max_id,
                        err_msg=f"Invalid ID. Please enter a number between 0 and {max_id}.",
                    )
    
                    selected_dict = all_transitions[molecule][chosen_id]
                    selected_setup = {
                        "band": selected_dict["bande"],
                        "ifreq": selected_dict["id"],
                        "freq_ghz": round(299792.458 / selected_dict["lambdaµm"], 4),
                        "trans_name": selected_dict["trans"],
                    }
    
            # Velocity profiling selection (Cube vs Map)
            mode = ask_value(
                prompt="Choose output mode - Data Cube [cube] or Channel Map [map]: ",
                cast_type=str,
                validation_cond=lambda m: m.lower() in ["cube", "map"],
                err_msg="Invalid mode. Please type 'cube' or 'map'.",
            ).lower()
    
            if mode == "cube":
                widthkms = ask_value(
                    prompt="Specify total velocity range (widthkms) in km/s: ",
                    cast_type=float,
                    validation_cond=lambda w: w > 0,
                    err_msg="Width must be a positive number.",
                )
                nlam = ask_value(
                    prompt="Number of velocity channels (nlam): ",
                    cast_type=int,
                    validation_cond=lambda n: n > 0,
                    err_msg="Number of channels must be a positive integer.",
                )
                vkms = None
            else:  
                vkms = ask_value(
                    prompt="Specify local LSR velocity (vkms) in km/s: ",
                    cast_type=float,
                    err_msg="Invalid number for velocity.",
                )
                widthkms, nlam = None, None
    
            sizeau = ask_value(
                prompt="Image size (AU): ",
                cast_type=int,
                validation_cond=lambda n: n > 0,
                err_msg="Image size must be a positive integer.",
            )
    
            npix = ask_value(
                prompt="Resolution of the image (number of pixels): ",
                cast_type=int,
                validation_cond=lambda n: n > 0,
                err_msg="npix must be a positive integer.",
            )

    else:
        # =========================================================================
        # NON-INTERACTIVE AUTOMATED MODE
        # =========================================================================
        if lbd is not None or lbdrange is not None:
            # Automated Dust Continuum Execution Path
            name = None  # Prevent UnboundLocalError
            old_cwd = os.getcwd()
            try:
                os.chdir(thermal_path)
                t1 = time.time()

                base_cmd = (
                    f"radmc3d image npix {npix} "
                    f"incl {params_dict['inclination']} "
                    f"posang {params_dict['pa']} "
                    f"sizeau {sizeau}"
                )

                if lbd is not None:
                    cmd = f"{base_cmd} lambda {lbd}"
                else:  
                    # Deterministically extract bounds passed in keyword sequence
                    lbdrange_min = lbdrange[0]
                    lbdrange_max = lbdrange[1]
                    cmd = f"{base_cmd} lambdarange {lbdrange_min} {lbdrange_max} nlam {nlam}"

                print(f"\n--> Running terminal command: {cmd}")
                subprocess.run(cmd, shell=True, check=True)

                name = image_name(cmd, path_lines=thermal_path)
                if os.path.exists("image.out"):
                    os.rename("image.out", name)
                    print(f"--> File successfully renamed to {name}")

                print(f"Process completed in {time.time() - t1:.2f} seconds.")

            except subprocess.CalledProcessError:
                print("--> Error: RADMC-3D subprocess failed.")
            finally:
                os.chdir(old_cwd)
                
            print("======================== PARAMETERS ========================")
            print(f'inclination = {params_dict["inclination"]}°')
            print(f'posang = {params_dict["pa"]}°')
            print(f'distance = {params_dict.get("dpc", "N/A")} pc')
            if lbd is not None: print(f'wavelength (µm) = {lbd:.2f}')
            else: 
                print(f'wavelength range : {lbdrange_min:.2f} to {lbdrange_max:.2f} µm')
                print(f"Spectral resolution : {(lbdrange_max-lbdrange_min)/(nlam-1):.2e} µm")
            print("============================================================")
                
            return name, None
        else:
            # Automated Molecular Line Execution Path
            molecule = molecule.upper()
            index_molec = available_molecules.index(molecule) + 1
            if alma_bands is None:
                print("--> ALMA bands set : 6 & 7")
                alma_bands = [6, 7]
            mode = mode.lower()
    
            print(f"\n=== Searching for {molecule} in band(s) {alma_bands} ===")
            all_transitions = find_alma_transitions(selected_bands=alma_bands, therm_path=thermal_path)
            selected_setup = check_and_select_transition(
                parsed_results=all_transitions, molecule=molecule, target_freq_ghz=freq, epsilon_ghz=tol
            )
    
            # Safety fallback: bind index to first available transition record if searching fails
            if selected_setup is None:
                if not all_transitions or molecule not in all_transitions:
                    print(f"No transitions found for {molecule} in ALMA bands {alma_bands}")
                    return None, None
                else:
                    selected_dict = all_transitions[molecule][0]
                    selected_setup = {
                        "band": selected_dict["bande"],
                        "ifreq": selected_dict["id"],
                        "freq_ghz": round(299792.458 / selected_dict["lambdaµm"], 4),
                        "trans_name": selected_dict["trans"],
                    }

    # =========================================================================
    # 3. SUBPROCESS EXECUTION (RADMC-3D) FOR MOLECULAR LINES
    # =========================================================================
    if selected_setup:
        print("[CASA export setup]")
        print(f"  -> ifreq to use: {selected_setup['ifreq']}")
        print(f"  -> nu0 to use: {selected_setup['freq_ghz']} GHz")

        name = None  # Prevent UnboundLocalError
        old_cwd = os.getcwd()
        try:
            os.chdir(thermal_path)
            t1 = time.time()

            # Compile line imaging arguments tracking line index arrays
            base_cmd = (
                f"radmc3d image npix {npix} "
                f"incl {params_dict['inclination']} "
                f"posang {params_dict['pa']} "
                f"sizeau {sizeau} iline {selected_setup['ifreq']} imolspec {index_molec}"
            )

            # Inject appropriate execution keywords depending on mode context
            if mode == "map":
                cmd = f"{base_cmd} vkms {vkms}"
            else:  
                cmd = f"{base_cmd} widthkms {widthkms} linenlam {nlam}"

            print(f"\n--> Running terminal command: {cmd}")
            subprocess.run(cmd, shell=True, check=True)

            # Re-label the final raster array
            name = image_name(cmd, path_lines=thermal_path)
            if os.path.exists("image.out"):
                os.rename("image.out", name)
                print(f"--> File successfully renamed to {name}")

            print(f"Process completed in {time.time() - t1:.2f} seconds.")

        except subprocess.CalledProcessError:
            print("--> Error: RADMC-3D subprocess failed.")
        finally:
            os.chdir(old_cwd)
            
        print("======================== PARAMETERS ========================")
        print(f'inclination = {params_dict["inclination"]}°')
        print(f'posang = {params_dict["pa"]}°')
        print(f'distance = {params_dict.get("dpc", "N/A")} pc')
        print(f'molecule = {molecule}')
        print(f"freq = {selected_setup['freq_ghz']} GHz")
        print(f"lambda = {round(299792.458 / selected_setup['freq_ghz'],2)} µm")
        print(f'iline = {selected_setup["ifreq"]}')
        if widthkms is None: print(f'vkms =  {vkms} km/s') 
        else: 
            print(f'widthkms = {widthkms} km/s')
            print(f'spec res = {(widthkms)/(nlam-1)} km/s  <=> {(selected_setup["freq_ghz"]*(widthkms / 299792.458)/(nlam-1)):e} GHz')
        print("============================================================")
            
        return name, widthkms

    return None, None
    
def casa_pipeline(
    image_folder,
    name,
    main_dict,
    obs_list,
    im,  # Argument sans valeur par défaut remonté pour éviter l'erreur de syntaxe
    widthkms=None,
    verbose=True,
    tc_imsize=[512,512],
    tc_cell="0.025arcsec",
    tc_threshold="2.2mJy",
    tc_interactive=False,
):
    """Execute the complete CASA simulation, measurement set combination, and tclean 

    synthesis imaging pipeline for multi-configuration radio-interferometric data.

    This function acts as an automated workflow inside the CASA environment. It 
    takes a theoretical raw spatial/spectral model (FITS file) generated by a 
    radiative transfer code like RADMC-3D, handles the geometric projection onto 
    equatorial coordinates, and dynamically executes interlinked steps of 
    synthetic observations (`simobserve`), data merging (`concat`), advanced 
    multiscale deconvolution (`tclean`), and standard astronomical FITS export.

    Dynamic Core Steps:
    -------------------
    1. Spectral Analysis: Evaluates the shape of the native RADMC-3D image array 
       to automatically determine the calculation profile ("cube" for multi-channel 
       data, "map" for single multi-frequency synthesis tracking).
    2. Model Alignment: Enforces absolute spatial centering by locking both the 
       telescope pointings (`indirection`) and the model celestial coordinates 
       `direction` to the exact target source structure tracked in `main_dict["coord"]`.
    3. Multi-Array Interferometry: Iterates through an array of diagnostic observation 
       dictionaries, executing distinct base baseline array tracking steps. Generated 
       Measurement Sets (.ms) are cataloged using a dynamic indexing scheme to streamline 
       subsequent multi-configuration concatenation operations.
    4. Advanced Imaging Synthesis: Invokes the `tclean` imaging engine wrapped under 
       a robust multiscale deconvolver model (`scales=[0, 8, 24]`) and a standard Briggs 
       robustness factor (`robust=1.0`) to secure an optimized trade-off between spatial 
       resolution limits and extended gas flux tracking recovery.

    Parameters
    ----------
    image_folder : str or pathlib.Path
        The filesystem base location directory path containing data files and serving as 
        the active workspace. Coerced into an absolute sub-folder (`image_folder / name`).
    name : str
        The primary string identifier model name prefix. Used to search for the initial 
        raw FITS file (`[name].fits`) and to label final deconvolved export files.
    main_dict : dict
        Metadata configuration input tracking dictionary. Must explicitly contain:
        - main_dict["coord"] (str): Target spatial source center equatorial coordinate metrics 
          given in a compliant J2000 standard layout (e.g., '16h28m13.8s -24d31m39s').
    obs_list : list of dict
        A list of observation configuration dictionaries tracking individual baseline setups. 
        Each block dictionary inside the list must map out the following explicit keys:
        - obs["antennalist"] (str): Telescope array layout config tracking path (e.g., 'alma.cycle10.5.cfg').
        - obs["totaltime"] (str): Total integrated observation time span on the target source (e.g., '6985s').
        - obs["refdate"] (str): Astronomical baseline reference evaluation timestamp calendar date (e.g., '2024/06/10').
        - obs["pwv"] (float): Local Precipitable Water Vapor atmospheric weather quality parameter in mm (e.g., 0.897).
        - obs["integration"] (str): Correlator accumulation integration dump interval window time (e.g., '6.048s').
    im : radmc3dPy.image.radmc3dImage object
        The active native RADMC-3D image model instance. Evaluated to dynamically look up 
        the spectral grid array density structure (`im.image.shape[2]`) and extract baseline 
        frequencies (`im.freq`) converted internally from Hertz to Giga-Hertz units.
    widthkms : float, optional
        The total framework velocity coverage span bandwidth expressed in km/s. Required 
        when the script initializes in line spectral data "cube" mode.
    verbose : bool, default True
        If True, displays internal tracking logs, subsystem reports, and processing updates.
    tc_imsize : list of int, optional
        A 2D array matrix size specification list [width, height] defining pixel dimensions 
        for deconvolution synthesis imaging. Defaults to [512, 512] if left unassigned or None.
    tc_cell : str, default "0.025arcsec"
        The precise grid angular dimension sampling pixel mapping scale value assigned to `tclean`.
    tc_threshold : str, default "2.2mJy"
        The flux density stopping threshold intensity limit targeted by the clean loop.
    tc_interactive : bool, default False
        If True, enables the manual interactive deconvolution window interface wizard screens.

    Returns
    -------
    None
        The function yields no native data return type. It processes dataset structures on disk 
        and saves a standard calibrated astronomical FITS data matrix labeled `image_final.fits` 
        inside the image folder workspace.

    Raises
    ------
    SyntaxError
        Raised natively if the function definition positioning order violates the native python 
        parameter constraint architecture rules (non-default keywords following default keyword blocks).
    RuntimeError
        Captured internally from underlying CASA task exceptions if input spatial sampling parameters 
        or rest frequency string labels break interferometric calculation matrices.
    """

    if verbose:
        print("\n=== STARTING CASA SIMULATION PIPELINE ===")
    # Safely determine mode depending on the RADMC-3D grid shape
    nlam = im.image.shape[2]
    mode = "cube" if nlam > 1 else "map"
    
    if verbose and mode == "cube":
        print("Please note that the first and last frequency channels are often set to zero and are therefore lost.")

    tclean_specmode = "cube" if mode == "cube" else "mfs"
    channel_width = (
        f"{widthkms / (nlam - 1):.3f}km/s"
        if (mode == "cube" and widthkms)
        else "0.1km/s"
    )
    num_channels = nlam if mode == "cube" else 1

    # Update path to point to the dedicated object subfolder
    image_folder = image_folder / name

    # Execute the entire workflow inside your 'cd' context manager
    with cd(image_folder):

        # 1. Convert FITS to CASA Image Model
        if verbose:
            print(f"--> Converting {name}.fits to native CASA .image model...")
        importfits(
            fitsimage=f"{name}.fits", imagename=f"{name}.image", overwrite=True
        )

        # 2. Dynamic loop over the cleaned parameters list
        ms_to_concat = []
        for i, obs in enumerate(obs_list):
            # Generate names on the fly using the loop index
            generated_project = f"obs_block_{i}"
            cfg_clean_name = obs["antennalist"].replace(".cfg", "")
            generated_ms = (
                f"{generated_project}/{generated_project}.{cfg_clean_name}.ms"
            )

            if verbose:
                print(
                    f"--> Simulating subarray block {i}: {obs['antennalist']} -> {generated_project}"
                )

            simobserve(
                project=generated_project,
                skymodel=f"{name}.image",
                indirection="J2000 " + str(main_dict["coord"]),
                direction="J2000 " + str(main_dict["coord"]),
                incenter=f"{np.mean(im.freq) * 1e-9:.2f}GHz",
                inwidth="",  
                obsmode="int",
                antennalist=obs["antennalist"],
                totaltime=obs["totaltime"],
                refdate=obs["refdate"],
                thermalnoise="tsys-atm",
                user_pwv=obs["pwv"],
                setpointings=True,
                integration=obs["integration"],
                maptype="ALMA",
                graphics="both" if verbose else "none",
                overwrite=True,
            )

            ms_to_concat.append(generated_ms)

        # 3. Concatenate all generated measurement sets
        if verbose:
            print(f"--> Merging tracking sets: {ms_to_concat}")
        concat(vis=ms_to_concat, concatvis="Combine.ms")

        # 4. Advanced Deconvolution (tclean)
        if verbose:
            print(
                f"--> Executing tclean synthesis imaging ({tclean_specmode.upper()} mode)..."
            )
        tclean(
            vis="Combine.ms",
            imagename=f"{name}_final",
            selectdata=True,
            datacolumn="data",
            imsize=tc_imsize,
            cell=tc_cell,
            specmode=tclean_specmode,
            nchan=num_channels,
            start="",
            width="",
            outframe="LSRK",
            deconvolver="multiscale",
            scales=[0, 8, 24],
            weighting="briggs",
            robust=1.0,
            niter=2000,
            cycleniter=500,
            threshold=tc_threshold,
            interactive=tc_interactive,
        )

        # 5. Export Back to FITS
        if verbose:
            print(f"--> Exporting final deconvolution to {name}_final.fits...")
        exportfits(
            imagename=f"{name}_final.image",
            fitsimage="image_final.fits",
            velocity=True,
            overwrite=True,
        )

    if verbose:
        print("=== CASA PIPELINE SUCCESSFUL ===\n")
        
        

def plot_spectral_cube(
    image_folder,
    filename="flying_saucer_final.fits",
    cmap="jet",
    figsize=(8, 7),
    zoom_window=None,
    channel_range=None,
    show_beam=True,
    beam_xy=(215, 215),
    output_prefix="image_plot",
    save_png=True,
    interactive_show=True,
    verbose=True,
):
    """
    Parse a multi-dimensional spectral FITS cube, dynamically extract WCS mapping,

    and plot sequential channel intensity maps.

    This function opens an astronomical FITS data cube (supporting 3D or 4D sub-arrays) 
    and handles spatial mapping using the World Coordinate System (WCS). It evaluates 
    the header architecture to determine the location and physical unit context of the 
    spectral axis, processes radio-velocity Doppler shifts if necessary, and projects 
    the astronomical synthesized beam configuration as a vector overlay patch over 
    the brightness distribution canvas.

    Dynamic Core Workflows:
    -----------------------
    1. Coordinate Slicing: Bypasses potential non-zero Stokes polarization index profiles 
       by isolating the primary spatial-spectral array framework matrix targets [Channels, Y, X].
    2. Spectral Axis Parsing: Scans standard keyword indices (CTYPE1 to CTYPE4) searching 
       for astronomical frequency or kinematic markers ('FREQ', 'VRAD', 'VOPT'). If coordinates 
       are encoded using radio velocity profiles, it dynamically calculates absolute rest 
       frequencies via the explicit linear radio Doppler tracking convention:
       nu = nu0 * (1.0 - (v_ms / c_ms))
    3. Beam Ellipse Reconstruction: Extracts the major/minor resolution limits (BMAJ, BMIN) 
       and Position Angle (BPA). It translates angular degrees/arcseconds down into absolute 
       pixel grid scale spaces, applying an explicit offset adjustment (bpa + 90) to match 
       the traditional trigonometric orientation mapping convention required by Matplotlib.

    Parameters
    ----------
    image_folder : str or pathlib.Path
        Relative or absolute target directory path hosting the input FITS file and serving 
        as the destination for exported image plots.
    filename : str, default "flying_saucer_final.fits"
        The filename target pointing to the input spectral FITS file.
    cmap : str, default "jet"
        The Matplotlib color map designation profile used to map pixel intensity arrays.
    figsize : tuple of int, default (8, 7)
        Linear dimension bounding constraints (width, height) assigned to each output window.
    zoom_window : tuple of four ints, optional
        A structural sub-array bounding box window given as (xmin, xmax, ymin, ymax). 
        Slices the native coordinate display matrix to isolate a localized central core view. 
        If None, renders the entire un-cropped array canvas.
    channel_range : tuple of two ints, optional
        Index markers indicating the starting and stopping interval array boundaries 
        expressed as (start_channel, end_channel). If None, iterates over the entire cube.
    show_beam : bool, default True
        If True, displays the synthesized elliptical clean beam resolution profile patch overlay.
    beam_xy : tuple of int, default (215, 215)
        Pixel cross-section center anchor coordinates used to plot the clean beam ellipse patch.
    output_prefix : str, default "image_plot"
        The file naming string prefix used to catalog and index exported static graphics files.
    save_png : bool, default True
        If True, dumps high-resolution compressed graphical plots onto the local file system.
    interactive_show : bool, default True
        If True, invokes standard interactive window display screens (plt.show()) inside 
        active terminal/notebook loops. If False, instantly forces canvas structure de-allocation 
        via plt.close() to preserve system volatile memory overhead.
    verbose : bool, default True
        If True, outputs runtime diagnostic updates, auto-detection tracking statuses, 
        and generation logs directly to the active standard terminal output stream.

    Returns
    -------
    None
        The function yields no native data return type. It reads data from files on disk and 
        generates structured static multi-spectral output plots.

    Raises
    ------
    FileNotFoundError
        Raised if the provided directory and filename combination does not map to a physical 
        item on the local storage partition.
    ValueError
        Raised if the target FITS matrix dimensions fail to resolve to a compliant 
        three-dimensional or four-dimensional array structure.
    KeyError
        Handled internally; falls back to an archived default baseline of 0.174 arcseconds 
        if beam parameters cannot be safely recovered from the primary header block or 
        secondary binary table extensions.
    """
    image_folder = Path(image_folder)
    fits_path = image_folder / filename

    if not fits_path.is_file():
        raise FileNotFoundError(f"Target FITS dataset not found at: {fits_path}")

    # 1. Load Data Cube
    hdul = fits.open(fits_path)
    header = hdul[0].header
    raw_data = hdul[0].data

    if raw_data.ndim == 4:
        cube = raw_data[0, :, :, :]  # Drop Stokes axis, keep [Channels, Y, X]
    elif raw_data.ndim == 3:
        cube = raw_data
    else:
        hdul.close()
        raise ValueError("Target FITS file is not a valid 3D or 4D spectral cube.")

    total_channels = cube.shape[0]
    if verbose:
        print(f"Spectral cube detected. Available channels: {total_channels}")

    wcs_2d = WCS(header).celestial

    # 2. Dynamic Spectral Axis Matching Lookup
    axis_idx = None
    for i in range(1, 5):
        ctype_key = f"CTYPE{i}"
        if ctype_key in header:
            ctype_value = str(header[ctype_key]).upper()
            if any(word in ctype_value for word in ["FREQ", "VRAD", "VOPT", "ENER", "WAV"]):
                axis_idx = i
                if verbose:
                    print(f"Spectral axis dynamically localized: Axis {axis_idx} ({ctype_value})")
                break

    if axis_idx is None:
        axis_idx = 4 if "CTYPE4" in header else 3
        if verbose:
            print(f"Non-standard spectral CTYPE profile. Forcing mapping indices to Axis {axis_idx}.")

    crval = header[f"CRVAL{axis_idx}"]
    crpix = header[f"CRPIX{axis_idx}"]
    cdelt = header[f"CDELT{axis_idx}"]

    # Doppler/Velocity mapping tracking fallback verification
    if "RESTFRQ" in header and any(w in str(header[f"CTYPE{axis_idx}"]).upper() for w in ["VRAD", "VOPT"]):
        restfreq = header["RESTFRQ"]
        c_ms = 299792458.0
        frequencies_hz = [
            restfreq * (1.0 - ((crval + (ch - (crpix - 1)) * cdelt) / c_ms))
            for ch in range(total_channels)
        ]
        if verbose:
            print("Kinematic velocity profile detected. Radio Doppler conversion applied.")
    else:
        frequencies_hz = [crval + (ch - (crpix - 1)) * cdelt for ch in range(total_channels)]

    # 3. Dynamic Beam Parameters Parsing Lookup
    if show_beam:
        try:
            if "BMAJ" in header:
                bmaj, bmin, bpa = header["BMAJ"], header["BMIN"], header["BPA"]
                pixel_scale = abs(header["CDELT1"])
                beam_major_pix = bmaj / pixel_scale
                beam_minor_pix = bmin / pixel_scale
            elif len(hdul) > 1:
                beam_table = hdul[1].data
                bmaj, bmin, bpa = beam_table["BMAJ"][0], beam_table["BMIN"][0], beam_table["BPA"][0]
                unit = hdul[1].header.get("TUNIT1", "deg")
                pixel_scale = abs(header["CDELT1"])
                factor = 3600 if "arc" in str(unit).lower() else 1
                beam_major_pix = bmaj / (pixel_scale * factor)
                beam_minor_pix = bmin / (pixel_scale * factor)
            else:
                raise KeyError
        except (KeyError, IndexError, Exception):
            if verbose:
                print("Missing or incomplete BEAMS table data. Applying default synthesized baseline (0.174 arcsec).")
            bpa = 0.0
            pixel_scale_arcsec = abs(header["CDELT1"]) * 3600
            beam_major_pix = beam_minor_pix = 0.174 / pixel_scale_arcsec

    # Establish loop range slicing constraints
    start_ch, end_ch = 0, total_channels
    if channel_range is not None:
        start_ch = max(0, channel_range[0])
        end_ch = min(total_channels, channel_range[1])

    # 4. Processing Loop inside your 'cd' context manager
    with cd(image_folder):
        for ch in range(start_ch, end_ch):
            freq_ghz = frequencies_hz[ch] / 1e9
            lambda_microns = 299792.458 / freq_ghz

            data_channel = cube[ch, :, :]

            # Discard processing frames with completely unmodeled data arrays
            if np.max(data_channel) == 0:
                if verbose:
                    print(f"   -> Channel {ch} skipped because the data slice is empty.")
                continue

            # Initialize Matplotlib Frame
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(1, 1, 1, projection=wcs_2d)

            # Map Data Flux Density intensities
            im_plot = ax.imshow(data_channel, cmap=cmap, origin="lower")
            plt.colorbar(im_plot, label="Intensity (Jy/beam)")

            # Add Clean synthesized Beam Ellipse patch overlay
            if show_beam:
                beam_ellipse = Ellipse(
                    xy=beam_xy,
                    width=beam_major_pix,
                    height=beam_minor_pix,
                    angle=bpa + 90,
                    edgecolor="white",
                    facecolor="grey",
                    linewidth=1.2,
                    alpha=0.85,
                )
                ax.add_patch(beam_ellipse)

            # Labels and titles layout
            ax.set_xlabel("Right Ascension (RA)")
            ax.set_ylabel("Declination (DEC)")
            
            # Enforce raw literal formatting to decouple LaTeX string tokens from Python escape parsing
            ax.set_title(fr"$\lambda$ = {lambda_microns:.2f} $\mu$m (Channel {ch})")

            # Handle dynamic bounding-box spatial zoom slicing windows if provided
            if zoom_window is not None and len(zoom_window) == 4:
                ax.set_xlim((zoom_window[0], zoom_window[1]))
                ax.set_ylim((zoom_window[2], zoom_window[3]))

            # Saving & Screen Display Output control infrastructure pipeline
            if save_png:
                out_name = f"{output_prefix}_ch{ch}_{lambda_microns:.1f}um.png"
                plt.savefig(out_name, dpi=300, bbox_inches="tight")
                if verbose:
                    print(f"   -> Plot frame exported successfully: {out_name}")

            if interactive_show:
                plt.show()
            else:
                plt.close(fig)

    hdul.close()
    if verbose:
        print("Processing finalized. All target visualization profiles generated.\n")
        




def write_fits(
    image_folder,
    name,
    main_dict,
    im,
    interactive=False,
    casa=True,
    verbose=True,
):
    """Manage the output directory filesystem and export the RADMC-3D image structure 

    to a clean, standardized FITS file.

    This utility prepares and structuralizes the target output workspace by resolving 
    potential directory name collisions. Depending on user constraints, it safely purges 
    stale operational data (interactively or via automated script overrides), unlinks 
    pre-existing duplicate tracking files, and invokes the underlying `radmc3dImage` 
    core exporter tool.

    Dynamic Infrastructure Behaviours:
    ----------------------------------
    1. Workspace Sanitation: In interactive mode (`interactive=True`), a terminal prompt 
       guards existing directories. Rejecting the override halts execution immediately 
       with a `return` mechanism to protect historical simulation outputs. In non-interactive 
       mode (`interactive=False`), it automates tree removal to prevent cluster pipeline hangs.
    2. Dynamic Format Optimization: Features an explicit metadata sanitation guard clause. 
       If the third dimension array density of the RADMC-3D payload evaluates to a single 
       channel (`im.image.shape[2] == 1`), the function forces `casa = False`. This bypasses 
       extended Common Astronomy Software Applications coordinate projection headers, 
       ensuring clean compatibility benchmarks for native 2D intensity maps.

    Parameters
    ----------
    image_folder : str or pathlib.Path
        The relative or absolute directory workspace location target path where the resulting 
        FITS product will be compiled. Coerced internally into a `pathlib.Path` instance.
    name : str
        The target output filename string identifier prefix (compiled without the '.fits' extension).
    main_dict : dict
        A structural metadata configuration tracking dictionary. Must explicitly expose 
        the following descriptive operational keys:
        - main_dict["dpc"] (float): Distance to the target astronomical source in parsecs (pc).
        - main_dict["coord"] (str): Center equatorial coordinate string formatted under standard 
          J2000 conventions (e.g., '16h28m13.8s -24d31m39s').
    im : radmc3dPy.image.radmc3dImage object
        The native input image instance vector framework built by the RADMC-3D execution module. 
        Evaluated dynamically to check multi-dimensional array matrix scales (`shape[2]`) 
        and extract baseline frequencies (`im.freq`) converted during export routines.
    interactive : bool, default False
        If True, initializes defensive terminal prompt wizards before undertaking destructive 
        filesystem operations. If False, executes silent automated directory resets.
    casa : bool, default True
        If True, injects specific structural coordinate world keywords into the FITS headers 
        to ensure native, out-of-the-box compatibility with CASA tracking tasks (e.g., `importfits`). 
        Automatically forced to False if the input data represents a 2D single channel continuum map.
    verbose : bool, default True
        If True, pushes active processing records, subsystem reports, directory creation logs, 
        and deletion exceptions directly onto the active standard console output stream.

    Returns
    -------
    None
        The function yields no native data return type. It creates structural directories 
        and outputs a standardized binary FITS file matrix onto the local partition.

    Raises
    ------
    PermissionError
        Captured and handled internally; logs an error message and safely aborts execution 
        if the host Operating System blocks directory tree purging due to insufficient 
        administrative privileges or file-system locking constraints.
    Exception
        Captured and handled internally; catches unexpected runtime OS-level directory 
        creation or file unlinking anomalies, preventing high-level pipeline loop crashes.
    """

    # Ensure image_folder is a Path object
    image_folder = Path(image_folder)

    # Check if the folder already exists
    if image_folder.exists():
        if image_folder.is_dir():
            if interactive:
                response = (
                    input(
                        f"The folder '{image_folder.name}' already exists. \nDo you want to overwrite it? (yes/no): "
                    )
                    .strip()
                    .lower()
                )
                if response in ["yes", "y"]:
                    try:
                        shutil.rmtree(image_folder)
                        image_folder.mkdir(parents=True, exist_ok=True)
                        if verbose:
                            print("Folder successfully reset.")
                    except PermissionError:
                        if verbose:
                            print(
                                "Error: Insufficient permissions to delete this folder."
                            )
                        return
                    except Exception as e:
                        if verbose:
                            print(f"An error occurred during deletion: {e}")
                        return
                else:
                    if verbose:
                        print(
                            "The folder was not deleted. Aborting FITS export."
                        )
                    return  # We stop here to avoid mixing old and new data
            else:
                # Non-interactive mode: auto-overwrite
                try:
                    shutil.rmtree(image_folder)
                    image_folder.mkdir(parents=True, exist_ok=True)
                    if verbose:
                        print("Folder automatically reset.")
                except PermissionError:
                    if verbose:
                        print(
                            "Error: Insufficient permissions to delete this folder."
                        )
                    return
                except Exception as e:
                    if verbose:
                        print(f"An error occurred during deletion: {e}")
                    return
        else:
            if verbose:
                print(
                    "Warning: A file with this name already exists, but it is not a folder."
                )
            return
    else:
        # Create the folder if it doesn't exist
        image_folder.mkdir(parents=True, exist_ok=True)
        if verbose:
            print(
                f"The folder '{image_folder.name}' has been successfully created."
            )

    # Define target FITS path properly using Path operator /
    file_path = image_folder / f"{name}.fits"

    # Remove existing FITS file inside the directory if it's there
    if file_path.is_file():
        file_path.unlink()
        if verbose:
            print("The old FITS file was successfully deleted.")

    # Write the new FITS file using your im1 tool
    if verbose:
        print(f"--> Writing FITS file to: {file_path}")

    if casa and im.image.shape[2]==1: casa = False

    im.writeFits(
        fname=str(file_path),
        dpc=main_dict["dpc"],
        coord=main_dict["coord"],
        casa=casa,
        nu0=np.mean(im.freq),
    )
    


def load_and_plot_image(path, name, main_dict, plot=True, cmap='jet', log=False, arcsec=True, save_png=True):
    """
    Load and plot a RADMC-3D image across all available wavelengths/channels,
    with an option to automatically save each frame to disk.

    Parameters
    ----------
    path : str
        Directory path where the RADMC-3D image file is stored.
    name : str
        Name of the image file (e.g., 'image.out').
    main_dict : dict
        Dictionary containing object/model parameters. Must include the key 
        'dpc' representing the distance to the source in parsecs (float).
    plot : bool, default True
        If True, displays the plots sequentially inside the Jupyter Notebook cells.
        If False, suppresses the display in Jupyter while still allowing download.
    cmap : str, default 'jet'
        The colormap passed to matplotlib for rendering the image intensity.
    log : bool, default False
        If True, plots the image with a logarithmic intensity scale.
    arcsec : bool, default True
        If True, the spatial axes (X and Y) will be displayed in arcseconds. 
        If False, they will be displayed in Astronomical Units (au).
    save_png : bool, default True
        If True, exports and saves each wavelength channel as a high-resolution 
        PNG file in the current working directory.

    Returns
    -------
    im : radmc3dImage
        The raw radmc3dImage object returned by `r3d.image.readImage`.
    """
    # Read the RADMC-3D raw image data
    im = r3d.image.readImage(path + name)
    
    # Get the total number of wavelength channels
    nlam = im.image.shape[2]
    
    if plot or save_png:
        for i in range(nlam):
            
            # 1. Initialize a clean figure
            fig = plt.figure()
            
            # 2. CRITICAL JUPYTER FIX: If plot is False, we use matplotlib's 'ioff' 
            # context manager to block the display ONLY for this specific execution.
            # This leaves the general Jupyter backend completely untouched.
            if not plot:
                plt.ioff()
                
            try:
                # 3. Call the black box plot function
                r3d.image.plotImage(
                    im, 
                    arcsec=arcsec, 
                    au=not arcsec, 
                    log=log,  
                    cmap=cmap, 
                    bunit='inu', 
                    ifreq=i, 
                    dpc=main_dict['dpc']
                )
                
                # 4. Force matplotlib to render the canvas pixels in memory
                plt.draw()
                
                # 5. Handle file export
                if save_png:
                    figure_folder = Path(f"Radmc3d_{name}")
                    os.makedirs(figure_folder, exist_ok=True)
                    filename = f'Radmc3d_{i}_{name}.png'
                    fig.savefig(figure_folder/filename, dpi=300, bbox_inches='tight')
                    print(f"Successfully saved: {filename}")
                
                # 6. Handle Jupyter Notebook display
                if plot:
                    plt.show()  # Explicitly trigger the notebook inline display
                else:
                    plt.close(fig)  # Close to free up memory
                    
            finally:
                # Always restore interactive mode for the rest of the loop/Notebook
                plt.ion()
                
            # 7. Memory cleanup
            plt.close('all')
            
    return im
