# ALMA Observation Simulation Pipeline: RADMC-3D & CASA Integration

This repository provides an automated end-to-end framework to synthesize physical models of protoplanetary and transition disks, compute their radiative transfer signatures, and simulate realistic interferometric observations using the **ALMA (Atacama Large Millimeter/submillimeter Array)** telescope configuration. The pipeline bridges the gap between theoretical astrophysics and observational data by coupling **RADMC-3D** with the **CASA (Common Astronomy Software Applications)** software suite.

■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■

## RADMC-3D and CASA Data Processing Pipeline :  `pipeline.py`

This module provides an integrated wrapper framework to orchestrate radiative transfer simulations via **RADMC-3D**, interferometric synthetic observations via **CASA**, and coordinate-aware astronomical data visualization utilizing **Astropy WCS** and **Matplotlib**.

---

### 1. System and Context Utilities

#### `silence(verbose=False)`
* **Purpose:** Context manager that handles terminal output suppression.
* **Behavior:** Suppresses or allows standard Python outputs (`stdout`, `stderr`) and filters runtime warnings globally based on the evaluation profile parameter.

#### `cd(newdir)`
* **Purpose:** Isolated runtime directory migration supervisor.
* **Behavior:** Safeguards the operating system execution state by migrating to a temporary directory context and automatically guaranteeing a fallback back to the original shell location upon block termination or unexpected crash loops.

#### `image_name(command, path_lines="")`
* **Purpose:** Dynamic file cataloger and parser.
* **Behavior:** Scans local configuration definitions (`lines.inp`) using text pattern matching to translate standard command line interface keywords into descriptive, molecule-specific dataset filenames.

---

### 2. Pipeline Execution Engines

#### `writing_input_files(chemistry_path, thermal_path, MOL, datadict, verbose=True)`
* **Purpose:** Input matrix file compiler for radiative transfer setups.
* **Behavior:** Wraps `astromugs` pipeline interfaces to translate physical chemistry parameters into standard input arrays tracking atomic number densities, spatial velocity vectors, and local 1D thermodynamic temperature fields.

#### `run_pipeline(...)`
* **Purpose:** Master orchestrator for **RADMC-3D** image grid generation.
* **Behavior:** Dual-mode processing engine that supports:
    * *Interactive Mode:* Runs a terminal setup wizard validating numerical types and mathematical constraints.
    * *Automated Scripting Mode:* Immediately processes keywords to trigger dust continuum or spectral lines calculations, dynamically handles subprocess exit codes, and standardizes structural outputs.

#### `casa_pipeline(...)`
* **Purpose:** Interferometric observation and image synthesis emulator.
* **Behavior:** Automated suite executing inside the **CASA** runtime environment. Converts raw spatial models to image models, runs baseline configuration setups via `simobserve`, merges measurement sets via `concat`, runs multiscale deconvolution synthesis via `tclean`, and exports back to FITS structures.

---

### 3. Astronomical Data Parsing and Transition Lookups

#### `find_alma_transitions(selected_bands, therm_path)`
* **Purpose:** Local receiver band line compiler.
* **Behavior:** Scans specific physical input libraries (`ALMA.dat` and `molecule_*.inp`) to map descriptive molecular properties down into matching receiver window target criteria.

#### `check_and_select_transition(parsed_results, molecule, target_freq_ghz, epsilon_ghz=0.5, target_band=None)`
* **Purpose:** Automated spectroscopic target selector.
* **Behavior:** Scans compiled dictionaries for an exact baseline frequency match within an explicit tolerance parameter. Implements a closest-match minimization fallback path when dealing with overlapping spectrum entries.

---

### 4. Visual Analysis and Coordinate Projection

#### `plot_spectral_cube(...)`
* **Purpose:** Multi-dimensional datacube image projection array viewer.
* **Behavior:** Parses 3D/4D FITS dimensions, runs radio Doppler conversions to parse velocity axes profiles, projects coordinate transformations via Astropy WCS, overlays physical clean beam ellipses, and dumps sequential PNG maps.

#### `write_fits(...)`
* **Purpose:** Workspace sanitation manager and exporter tool.
* **Behavior:** Handles structural file system overrides, validates array sizes to drop degenerate dimensional axes, and calls `radmc3dPy` core tools to ensure robust compatibility with downstream CASA tracking packages.

#### `load_and_plot_image(...)`
* **Purpose:** Raw model matrix explorer and memory manager.
* **Behavior:** Reads native code output configurations and steps over channel counts to map pixel variations. Uses dedicated interactive hooks (`plt.ioff`/`plt.ion`) to prevent cluster resource blockages inside notebook cells.

#### `plot_line_with_multiple_continuum_contours(...)`
* **Purpose:** Multi-wavelength publication-ready composite map maker.
* **Behavior:** Loads line emission backdrops alongside continuum maps, reprojects distinct resolution arrays onto a uniform coordinate canvas via `reproject_interp`, sets custom contrast color layers, and projects fractional contour arrays.

#### `moment_map(...)`
* **Purpose:** Mathematical velocity and gas flux integrator.
* **Behavior:** Wraps the CASA `immoments` engine to isolate specific kinematic maps ($M_0$ Integrated Flux, $M_1$ Velocity Vector, $M_2$ Dispersion Limit, or $M_8$ Peak Intensity). Transforms results using coordinate projection layouts, locks scaling boundaries, and appends a telescope beam patch.

■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
## RADMC-3D and CASA Image Generation Pipeline Notebook (`Making_images_radmc3d_CASA_pipeline.ipynb`)

This notebook orchestrates the workflow from astronomical chemistry data to synthetic interférometric observations of disk structures like the "Flying Saucer" edge-on disk.

---

### 1. Environment and Dependencies

#### `Imports`
* Loads `astromugs` for structure configuration and plotting, `radmc3dPy` for radiative transfer, standard Python scientific libraries, and the `astromugs.casa.pipeline` wrapper (`cas`) for CASA imaging tasks.

---

### 2. Global Parameters & Object Setup

#### `Data Paths and Target Setup`
* Defines target molecules and file paths for chemistry and thermodynamics.
* Stores physical parameters for the "Flying Saucer" disk (e.g., coordinates, $M_{\star} = 0.60\ M_{\odot}$, inclination = $87^\circ$).
* Configures ALMA synthetic observing blocks using different antenna configuration lists (`cycle10.5` and `cycle10.3`).

---

### 3. Model Generation & Radiative Transfer

#### `Input Files and Quality Control`
* Calls `cas.writing_input_files(...)` to create RADMC-3D grids tracking density profiles, Keplerian gas velocity fields, and temperatures.
* Contains commented-out plotting hooks for visual structure validation.

#### `RADMC-3D Line Run`
* Automatically scans local database files for the target line (e.g., CO $J=4\rightarrow3$ at 345.796 GHz) and launches the RADMC-3D LTE ray-tracer subprocess to compute raw image intensity data (`image.out`).

---

### 4. CASA Interferometry Engine

#### `Instrument Emulation via cas.casa_pipeline`
Routes raw theoretical physical images through virtual ALMA configurations:

* **Dust Continuum Path:** Generates continuum images at specific wavelengths (e.g., 866 $\mu$m and 1300 $\mu$m) using `cas.run_pipeline`, converts them to FITS format via `cas.write_fits`, and invokes `cas.casa_pipeline` for interferometric visibility sampling (`simobserve`), data merging (`concat`), and multi-frequency cleaning (`tclean`).
* **Molecular Line Path:** Generates a 3D Position-Position-Velocity (PPV) data cube. It calls `cas.casa_pipeline` with the `cube` mode parameter activated to automatically deconvolve and clean each velocity channel under realistic weather noise conditions.

---

### 5. Post-Processing Analysis

#### `Moment Maps & Contours`

> [!WARNING]
> You need to compute continuum for all target wavelengths before running `cas.plot_line_with_multiple_continuum_contours`

* **Kinematics Engine (`cas.moment_map`):** Computes Moment 0 (Integrated Flux intensity maps) and Moment 1 (Velocity field maps) from the cleaned CASA line cube.
* **Multi- $\lambda$ Overlays:** Uses `cas.plot_line_with_multiple_continuum_contours` to reproject the multi-wavelength dust continuum maps onto the celestial coordinate framework (WCS) of the gas raie emission. It overlays discrete continuum fractional contours (10%, 50%, 90% peak) over the moment map to analyze spatial stratification.

■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
## `ALMA.dat` 
This file contains ALMA data, like bands and frequencies.

> [!NOTE]
> It must be in the thermal directory.
