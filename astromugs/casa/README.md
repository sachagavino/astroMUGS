# ALMA Observation Simulation Pipeline: RADMC-3D & CASA Integration

This repository provides an automated end-to-end framework to synthesize physical models of protoplanetary and transition disks, compute their radiative transfer signatures, and simulate realistic interferometric observations using the **ALMA (Atacama Large Millimeter/submillimeter Array)** telescope configuration. The pipeline bridges the gap between theoretical astrophysics and observational data by coupling **RADMC-3D** with the **CASA (Common Astronomy Software Applications)** software suite.

---

## 📁 File Structure & Core Components

### 1. `pipeline.py` (The Computational Engine)
This file acts as the core library module (`astromugs.pipeline`). It encapsulates the low-level data structures, system hooks, and astronomical transformation matrices required to translate physical density structures into synthetic raw visibilities and deconvolved FITS images.

* **Physical Input Generation (`writing_input_files`):** Converts raw outputs from chemical kinetics codes (e.g., Nautilus) into standard RADMC-3D grid files, defining parameters such as molecular number densities (`numberdens_XX.inp`), systematic gas velocities (`gas_velocity.inp`), and 1D/3D thermal structures (`gas_temperature.inp`).
* **Automated Line Transition Parsing (`find_alma_transitions` & `check_and_select_transition`):** Scans localized molecular data files (`molecule_*.inp`) and queries internal ALMA band catalogs. It identifies specific quantum transitions, matching the requested frequency within an adjustable Doppler/frequency tolerance window.
* **Radiative Transfer Execution (`run_pipeline`):** Coordinates background subprocess shells to execute the RADMC-3D core binary. It supports dual profiles: **Dust Continuum Emission** (producing 2D intensity maps across a continuous wavelength range) and **Molecular Line Transfer** (generating 3D Position-Position-Velocity data cubes).
* **Interferometric Simulator (`casa_pipeline`):**
    * Imports the raw, perfect sky-model FITS file and projects it onto equatorial coordinate systems (J2000).
    * Simulates realistic synthesis observations via the `simobserve` task, factoring in specific ALMA antenna configurations (e.g., `alma.cycle10.cfg`), customizable total integration times, reference dates, and atmospheric corruptions via Precipitable Water Vapor (PWV) parameters.
    * Concatenates multi-configuration datasets (`concat`) and feeds them into a multi-scale deconvolution loop (`tclean`) utilizing a Briggs weighting robustness scale to balance resolution against extended flux recovery.
* **Advanced Visualizations (`plot_spectral_cube` & `load_and_plot_image`):** Automatically slices multi-dimensional FITS data cubes along the spectral axis, computes local Doppler velocity shifts, overlays the synthesized clean beam ellipse, and exports static graphical frames to disk.

---

### 2. `Making_images_radmc3d_CASA_pipeline.ipynb` (The Production Interface)
This Jupyter Notebook serves as an interactive implementation guide and high-level wrapper. It provides a reproducible execution sequence demonstrating the pipeline's capabilities on a famous astrophysical benchmark: **The "Flying Saucer" Edge-On Disk**.

The workflow is structured into five distinct operational phases:
1.  **Source Parametrization:** Defines the explicit structural orientation of the target system, including registered J2000 astrometric coordinates, distance (120 pc),  stellar mass ($0.60 M_{\odot}$), and a highly inclined geometry ($87.0^{\circ}$).
2.  **Observational Setup:** Inputs parameters extracted from actual ALMA Archive project tracks (`2023.1.00907.S`), defining a multi-baseline configuration using both extended (`alma.cycle10.5.cfg`) and compact (`alma.cycle10.3.cfg`) antenna arrays.
3.  **Synthetic Data Generation:** Triggers RADMC-3D to calculate the localized non-LTE or LTE population levels and execute ray-tracing for a specific tracer molecule (e.g., Carbon Monoxide, $^{12}\text{CO}$) within the ALMA Band 7 receiver frequency framework.
4.  **CASA Synthesis Imaging:** Conversions to native CASA `.image` formats are performed on disk, followed by visibilities merging, automated Fourier transformations, and multiscale clean loops.
5.  **Data Reduction & Analysis:** Generates final channel maps showing intensity patterns ($\text{Jy/beam}$) as a function of velocity space, instantly exporting individual frames to local directories for comparative studies against real observational data.
