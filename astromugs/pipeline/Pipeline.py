import glob
import os
import sys
import shutil
import numpy as np


from astromugs import radmc3d
from astromugs import nautilus
from astromugs.pipeline.Grid import Grid
from astromugs.modeling.Disk import Disk
from astromugs.utils.struct import StructureParams
from astromugs.utils.thermal import ThermalParams
from astromugs.utils import custom as custom_io

from astromugs.constants.constants import autocm, M_sun, R_sun, c, amu, mu, black_body
from astromugs.modeling.InterstellarRadFields import InterstellarRadFields

import matplotlib.pyplot as plt


class Pipeline:
    """Main user-facing class for setting up and running disk models.

    Orchestrates the full modeling pipeline: dust continuum radiative transfer
    with RADMC-3D, line radiative transfer, local radiation field computation,
    and chemical modeling with Nautilus. Holds the grid, structural parameters,
    thermal parameters, and provides methods to read, write, and run each
    modeling stage.

    Attributes
    ----------
    params : StructureParams
        Structural parameters for the disk model.
    thermalparams : ThermalParams
        Thermal and wavelength parameters for radiative transfer.
    thermalpath : str
        Relative path to the directory containing thermal / RADMC-3D files.
    chempath : str
        Relative path to the directory containing chemistry (Nautilus) files.
    radmc3d_cmd : str
        Command used to invoke RADMC-3D (can be overridden to point to a
        specific binary).
    grid : Grid
        The spatial grid object holding coordinates, dust, and field data.
    nautilus : module
        Reference to the Nautilus chemistry module.
    Ta : ndarray
        Surface-area-weighted dust temperature array, populated after
        ``run_continuum`` or ``write_nautilus``.
    """

    def __init__(self):
        """Initialize the Pipeline with default parameters, grid, and paths."""
        self.params = StructureParams()
        self.thermalparams = ThermalParams()
        self.thermalpath = 'thermal/'
        self.chempath = 'chemistry/'
        self.radmc3d_cmd = 'radmc3d' #in case the user wants to point to a specific radmc3d
        self.nmgc_cmd = 'nmgc'       #in case the user wants to point to a specific nmgc binary
        self.grid = Grid(params=self.params.disk, wave=self.thermalparams.wave)
        self.nautilus = nautilus


    def run_continuum(self, write_control=False, **keywords):
        """Run dust continuum radiative transfer and compute area-weighted temperature.

        Executes the RADMC-3D Monte Carlo dust radiative transfer, reads the
        resulting dust temperature, and computes the surface-area-weighted
        temperature ``self.Ta``. If the output files already exist they are
        read directly; missing files produce a warning but do not halt
        execution (useful when running chemistry without a dust structure).

        Parameters
        ----------
        write_control : bool, optional
            If True, write the RADMC-3D control file before running.
            Default is False.
        **keywords
            Additional keyword arguments forwarded to
            ``run_thermal_radmc3d``.
        """
        # FOLDER WITH THERMAL FILES
        thermpath = self.thermalpath

        # WRITE THERMAL FILES
        if write_control == True:
            self.write_continuum(control=True)

        self.run_thermal_radmc3d(**keywords)




        # THIS SECTION READS THE THERMAL FILES IF THE FILES EXIST.
        # IF THE FILES DO NOT EXIST, A WARNING IS PRINTED BUT THE CODE CONTINUES. ERRORS CAN BE RAISED.
        # WHY DOES THE CODE CONTINUE? BECAUSE SOMETIMES THE USER MAY WANT TO RUN THE CHEMISTRY WITHOUT HAVING A DUST STRUCTURE (E.G. IF THEY WANT TO DERIVE DUST_DENSITY.INP FROM AN EXISTING CHEMISTRY MODEL).
        # ---- READ SECTION ----
        densityfile = thermpath + "dust_density.inp"
        starfile = thermpath + "stars.inp"
        temperaturefile = thermpath + "dust_temperature.dat"

        if os.path.exists(densityfile):

            # READ THERMAL FILES
            nx, ny, nz, x, y, z  = radmc3d.read.grid(thermpath)
            self.grid.set_spherical_grid(r_edge=np.array(x)/autocm, theta_edge=y, phi_edge=z)
            dust_density = radmc3d.read.dust_density(thermpath) # Gives a list of one numpy array. Will be updated when multiple structures
            nbspecies = int(len(dust_density[0])/(nx*ny*nz))
            dust_density = np.reshape(dust_density[0], (nbspecies, nz, ny, nx))
            dust_density[dust_density<=1e-100] = 1e-100 # get a minimum dust density in order to avoid absurd values for chemistry.

            if os.path.exists(temperaturefile):

                self.grid.temperature = radmc3d.read.dust_temperature(thermpath)
                if len(self.grid.temperature) > 0:
                    self.grid.temperature[0] = np.reshape(self.grid.temperature[0], (nbspecies, nz, ny, nx))

                    if nbspecies > 1:

                        a = []
                        dustmass = []

                        for istruc in range(0, len(self.grid.dust)): # create a loop in order to gather all grain sizes into one single array that will be used to compute the area-weighted temperature Ta.

                            dustmodel = self.grid.dust[istruc]
                            #dustmass = dustmodel.grainmass() # gram
                            sizes = dustmodel.sizes()
                            #a = sizes[-1]*1e-4 # cm
                            a.append(sizes[-1]*1e-4) # cm

                            dustmass.append(dustmodel.grainmass()) # gram
                        a = np.hstack(a) # in order to get a numpy list of all sizes no matter how many structures there are.
                        dustmass = np.hstack(dustmass)

                        Ta_num = np.zeros((nz, ny, nx))
                        Ta_denum = np.zeros((nz, ny, nx))

                        for idx in range(len(a)):
                            Ta_num += self.grid.temperature[0][idx, :, :, :]*(dust_density[idx, :, :, :]/dustmass[idx])*a[idx]**2
                            Ta_denum += (dust_density[idx, :, :, :]/dustmass[idx])*a[idx]**2
                        self.Ta = Ta_num/Ta_denum  #Ta is the area-weighted dust temperature.
                    else:
                        self.Ta = self.grid.temperature[0][0,:,:,:]
                else:
                    print(f"\n{temperaturefile} exist but is empty or corrupted.\n")

            else:
                print(f"WARNING: {temperaturefile} does not exist.")

        else:
            print(f"WARNING: {densityfile} does not exist.")

        if os.path.exists(starfile):
            nb_lam, lam, r_star, m_star, Tstar, spectrum = radmc3d.read.stars(thermpath)  #read stars file.
        else:
            print(f"WARNING: {starfile} does not exist.")
        #---- END OF READ SECTION


    def run_line(self, make_image=True, \
                       write_control=False, \
                       lines_format='leiden', \
                       path='chemistry/', #relative or absolute path\
                       incl=90, \
                       npix=800, \
                       iline=None, \
                       lambda_micron=None, \
                       widthkms=10, \
                       linenlam=100, \
                       itime=59, \
                       species='CO', \
                       **keywords):
        """Run line radiative transfer and optionally produce a channel map image.

        Parameters
        ----------
        make_image : bool, optional
            If True, generate a synthetic image via RADMC-3D. Default is True.
        write_control : bool, optional
            If True, write the RADMC-3D control file before running.
            Default is False.
        lines_format : str, optional
            Format of the molecular data file (e.g., 'leiden'). Default is
            'leiden'.
        path : str, optional
            Relative or absolute path to the chemistry model directory.
            Default is 'chemistry/'.
        incl : float, optional
            Disk inclination angle in degrees. Default is 90.
        npix : int, optional
            Number of pixels per side of the output image. Default is 800.
        iline : int or None, optional
            Line transition index for RADMC-3D. Default is None.
        lambda_micron : float or None, optional
            Wavelength in microns at which to produce the image. Default is
            None.
        widthkms : float, optional
            Velocity width of the line channel map in km/s. Default is 10.
        linenlam : int, optional
            Number of wavelength channels across the line. Default is 100.
        itime : int, optional
            Index of the chemistry time output to use. Default is 59.
        species : str, optional
            Chemical species for the line transfer (e.g., 'CO'). Default is
            'CO'.
        **keywords
            Additional keyword arguments forwarded to
            ``run_image_radmc3d``.
        """
        thermpath = self.thermalpath

        if make_image == True:
            self.run_image_radmc3d(incl=incl, npix=npix, iline=iline,  lambda_micron=lambda_micron, **keywords)




    def convert_nautilus2radmc(self, species='CO', dust_density=False, dust_temperature=False, numberdens=False):
        """Convert Nautilus chemistry output to RADMC-3D input files.

        Interpolates chemistry model results onto the radiative transfer grid
        and writes the corresponding RADMC-3D input files. The chemistry model
        must first be attached via ``add_chemistry()``. An ``amr_grid.inp``
        file must already exist (or the grid must be set) because the chemistry
        grid is typically at lower resolution than the RT grid and requires
        interpolation.

        Parameters
        ----------
        species : str or list of str, optional
            Chemical species to convert (e.g., 'CO', 'H2O'). Used for
            creating ``numberdens_XX.inp``. Default is 'CO'.
        dust_density : bool, optional
            If True, write a ``dust_density.inp`` file. Default is False.
        dust_temperature : bool, optional
            If True, write a ``dust_temperature.inp`` file. Default is False.
        numberdens : bool, optional
            If True, write a ``numberdens_XX.inp`` file for the given
            species. Default is False.
        """
        thermpath = self.thermalpath

        if os.path.isfile(thermpath + 'amr_grid.inp'):
            nx, ny, nz, x, y, z  = radmc3d.read.grid(thermpath) # read grid. the amr_grid.inp shows the border of the cells so we convert them.
            self.grid.set_spherical_grid(r_edge=np.array(x)/autocm, theta_edge=y, phi_edge=z)
            x, y, z = self.grid.r, self.grid.theta, self.grid.phi
        elif hasattr(self.grid, 'r_edge') and self.grid.r_edge is not None:
            nx, ny, nz, x, y, z = self.grid.nr, self.grid.ntheta, self.grid.nphi, self.grid.r, self.grid.theta, self.grid.phi
            nx = nx-1
            ny = ny-1
            nz = nz-1
        else:
            print(f"Error: amr_grid not found. Make sure there is one (example: create one using grid.set_spherical_grid())")
            sys.exit(1)  # Exit the program with a non-zero status to indicate an error

        if numberdens==True:
            species_list = [species] if isinstance(species, str) else list(species)
            for sp in species_list:
                numberdens_sph = nautilus.coupling.to_spherical(self.grid.chemmodel[sp], nx, ny, x, y, struct='numberdens_species')
                radmc3d.write.numberdens_mol(numberdens_sph, species=sp, gridstyle="regular", thermpath=thermpath)

        if dust_density==True:
            print('create dust_density.inp')

        if dust_temperature == True:
            print('create dust_temperature.inp')



    def run_localfield(self, nphot_mono=None, write_mcmono=False, run=True, **keywords):
        """Compute the local mean radiation field with RADMC-3D monochromatic Monte Carlo.

        Optionally writes the monochromatic wavelength file, runs the
        RADMC-3D monochromatic Monte Carlo, and reads the resulting
        ``mean_intensity.out`` into ``self.grid.localfield``.

        Parameters
        ----------
        nphot_mono : int or None, optional
            Number of monochromatic photon packages. If None, RADMC-3D uses
            its default. Default is None.
        write_mcmono : bool, optional
            If True, write the ``mcmono_wavelength_micron.inp`` file before
            running. Default is False.
        run : bool, optional
            If True, execute the RADMC-3D monochromatic run. Set to False to
            skip execution and only read existing output. Default is True.
        **keywords
            Additional keyword arguments forwarded to ``write_continuum``.
        """

        self.write_continuum(mcmono=write_mcmono, **keywords)

        if run == True:
            self.run_localfield_radmc3d(nphot_mono=nphot_mono)

        thermpath = self.thermalpath
        nx, ny, nz, x, y, z  = radmc3d.read.grid(thermpath)
        nlam_mono, lam_mono, self.grid.localfield = radmc3d.read.localfield(thermpath)
        lam_mono = 1e6*(c/lam_mono)
        self.grid.localfield = np.reshape(self.grid.localfield, (nlam_mono, nz, ny, nx))


    def run_thermal_radmc3d(self, verbose=True, timelimit=7200, \
            nice=None, **keywords):
        """Run the RADMC-3D thermal Monte Carlo dust radiative transfer.

        Parameters
        ----------
        verbose : bool, optional
            If True, print RADMC-3D output to stdout. Default is True.
        timelimit : int, optional
            Maximum wall-clock time in seconds before the run is killed.
            Default is 7200 (2 hours).
        nice : int or None, optional
            Unix nice priority for the RADMC-3D process. Default is None.
        **keywords
            Additional keyword arguments (currently unused).
        """
        radmc3d.run.thermal(verbose=verbose, timelimit=timelimit, nice=nice, thermpath=self.thermalpath, radmc3d_cmd=self.radmc3d_cmd)


    def run_chemistry(self, verbose=True, timelimit=None):
        """Run the NMGC gas-grain chemistry simulation.

        Calls ``nmgc run`` in ``self.chempath``, which must already contain
        all required NMGC input files. Set ``self.nmgc_cmd`` beforehand if
        the ``nmgc`` binary is not on PATH.

        Parameters
        ----------
        verbose : bool, optional
            If True, NMGC output streams to the terminal. If False, it is
            redirected to ``<chempath>/nmgc.out``. Default is True.
        timelimit : float or None, optional
            Timeout in seconds. None means no timeout. Default is None.
        """
        nautilus.run.run(chempath=self.chempath, nmgc_cmd=self.nmgc_cmd,
                         verbose=verbose, timelimit=timelimit)

    def run_chemistry_custompath(self, custom_path, verbose=True, timelimit=None):
        """Run the NMGC gas-grain chemistry simulation.

        Calls ``nmgc run`` in ``custom_path``, which must already contain
        all required NMGC input files. Set ``self.nmgc_cmd`` beforehand if
        the ``nmgc`` binary is not on PATH.

        Parameters
        ----------
        verbose : bool, optional
            If True, NMGC output streams to the terminal. If False, it is
            redirected to ``<chempath>/nmgc.out``. Default is True.
        timelimit : float or None, optional
            Timeout in seconds. None means no timeout. Default is None.
        """
        nautilus.run.run(chempath=custom_path, nmgc_cmd=self.nmgc_cmd,
                         verbose=verbose, timelimit=timelimit)
    

    def run_localfield_radmc3d(self, nphot_mono=None, verbose=True, timelimit=7200):
        """Run the RADMC-3D monochromatic Monte Carlo for the local radiation field.

        Parameters
        ----------
        nphot_mono : int or None, optional
            Number of monochromatic photon packages. If None, RADMC-3D uses
            its default. Default is None.
        verbose : bool, optional
            If True, print RADMC-3D output to stdout. Default is True.
        timelimit : int, optional
            Maximum wall-clock time in seconds. Default is 7200 (2 hours).
        """
        radmc3d.run.localfield(nphot_mono=nphot_mono, verbose=verbose, timelimit=timelimit, thermpath=self.thermalpath, radmc3d_cmd=self.radmc3d_cmd)


    def run_image_radmc3d(self, npix=300, lambda_micron=None, iline=None, incl=None, verbose=True):
        """Run RADMC-3D to produce a synthetic image.

        Parameters
        ----------
        npix : int, optional
            Number of pixels per side of the image. Default is 300.
        lambda_micron : float or None, optional
            Wavelength of the image in microns. Default is None.
        iline : int or None, optional
            Line transition index. Default is None.
        incl : float or None, optional
            Inclination angle in degrees. Default is None.
        verbose : bool, optional
            If True, print RADMC-3D output to stdout. Default is True.
        """
        radmc3d.run.image(npix=npix, lambda_micron=lambda_micron, iline=iline, incl=incl, verbose=verbose, timelimit=7200, thermpath=self.thermalpath, radmc3d_cmd=self.radmc3d_cmd)


    def write_continuum(self, dens=False, grid=False, opac=False, control=False, stars=False, wave=False, mcmono=False, ext=False):
        """Write RADMC-3D input files for dust continuum radiative transfer.

        Each boolean flag controls whether the corresponding input file is
        written. If no flags are set to True, a warning is printed. The
        thermal directory is created if it does not exist.

        Parameters
        ----------
        dens : bool, optional
            Write ``dust_density.inp``. Default is False.
        grid : bool, optional
            Write ``amr_grid.inp``. Default is False.
        opac : bool, optional
            Write ``dustopac.inp`` (reads existing ``dustkap*`` files in
            the thermal directory). Default is False.
        control : bool, optional
            Write ``radmc3d.inp`` control file. Default is False.
        stars : bool, optional
            Write ``stars.inp`` with stellar properties. Default is False.
        wave : bool, optional
            Write ``wavelength_micron.inp``. Default is False.
        mcmono : bool, optional
            Write ``mcmono_wavelength_micron.inp`` for monochromatic Monte
            Carlo. Default is False.
        ext : bool, optional
            Write ``external_source.inp`` for the interstellar radiation
            field. Default is False.
        """
        #os.system("rm thermal/*.inp")

        thermpath = self.thermalpath

        if not os.path.exists(thermpath):
            os.makedirs(thermpath)

        if control==True:
            print('\nWriting radmc3d.inp:')
            print('----------------------------')
            radmc3d.write.control(self.thermalparams.control, thermpath=thermpath)

        if stars==True:
            print('\nWriting stars.inp:')
            print('----------------------------')
            mstar = []
            rstar = []
            xstar = []
            ystar = []
            zstar = []
            tstar = []
            for i in range(len(self.grid.stars)):
                mstar.append(self.grid.stars[i].mass*M_sun)
                rstar.append(self.grid.stars[i].radius*R_sun)
                xstar.append(self.grid.stars[i].x*autocm)
                ystar.append(self.grid.stars[i].y*autocm)
                zstar.append(self.grid.stars[i].z*autocm)
                tstar.append(self.grid.stars[i].temperature)

            radmc3d.write.stars(rstar, mstar, self.grid.lam, xstar, ystar, zstar, \
                    tstar=tstar, thermpath=thermpath)

        if wave==True:
            print('\nWriting wavelength_micron.inp:')
            print('----------------------------')
            radmc3d.write.wavelength_micron(self.grid.lam, thermpath=thermpath)

        if mcmono==True:
            print('\nWriting mcmono_wavelength_micron.inp:')
            print('----------------------------')
            radmc3d.write.mcmono_wavelength_micron(self.grid.monolam, thermpath=thermpath)

        if ext==True:
            if len(self.grid.isrf) != 0:
                radmc3d.write.external_rad(self.grid.isrf[0], thermpath=thermpath)

        if grid==True:
            print('\nWriting amr_grid.inp:')
            print('----------------------------')
            if self.grid.coordsystem == 'spherical':
                radmc3d.write.amr_grid(self.grid.r_edge*autocm, self.grid.theta_edge, self.grid.phi_edge, gridstyle="regular", coordsystem=self.grid.coordsystem, thermpath=thermpath)

        if dens==True:
            print('\nWriting dust_density.inp:')
            print('----------------------------')
            radmc3d.write.dust_density(self.grid.dustdensity, gridstyle="regular", thermpath=thermpath)

        if opac==True:
            print('\nWriting dustopac.inp:')
            print('----------------------------')
            dustopac = []
            filelist = glob.glob(thermpath + 'dustkap*')
            for files in sorted(filelist):
                dustopac.append(files)
            radmc3d.write.dustopac(dustopac, thermpath=thermpath)

        if len(self.grid.accretionheating) > 0:
            radmc3d.write.accretion_heating(self.grid.w1*autocm, self.grid.w2, self.grid.w3, self.grid.accretionheating[0], gridstyle="regular", thermpath=thermpath)

        if dens == False and grid == False and control == False and opac == False and stars == False and wave == False:
            print('\nWARNING: no RADMC3D file are created.\n')
            print('----------------------------\n')

    def write_line(self, control=False, line=False, gasvelocity=False, gastemp=False, microturb=False, line_format='leiden', species='CO', star_mass=1):
        """Write RADMC-3D input files for line radiative transfer.

        Parameters
        ----------
        control : bool, optional
            Write ``radmc3d.inp`` control file. Default is False.
        line : bool, optional
            Write ``lines.inp`` molecular line data file. Default is False.
        gasvelocity : bool, optional
            Write ``gas_velocity.inp`` assuming Keplerian rotation. Default
            is False.
        gastemp : bool or str, optional
            Controls gas temperature handling.

            - ``False`` (default): no ``gas_temperature.inp`` is written.
              Use ``pipe.thermalparams.control.tgas_eq_tdust = 1`` before
              ``write_line(control=True)`` to let RADMC-3D equate gas and
              dust temperatures internally.
            - ``'1D_static'``: read the gas temperature column ``Tg`` from
              each ``<radius>AU/1D_static.dat`` file in the chemistry
              directory, remap it to the spherical RADMC-3D grid via
              bilinear interpolation (the same method as
              :meth:`convert_nautilus2radmc`), and write
              ``gas_temperature.inp``. This reproduces the gas temperature
              that was actually fed into Nautilus and is independent of
              the chemical timestep.

        microturb : bool, optional
            Write microturbulence file (not yet implemented). Default is
            False.
        line_format : str, optional
            Format of the molecular data file (e.g., 'leiden'). Default is
            'leiden'.
        species : str or list, optional
            Chemical species for the line file. Default is 'CO'.
        star_mass : float, optional
            Stellar mass in solar masses, used for Keplerian velocity
            computation. Default is 1.
        """
        thermpath = self.thermalpath

        if not os.path.exists(thermpath):
            os.makedirs(thermpath)

        if control==True:
            print('\nWriting radmc3d.inp:')
            print('----------------------------')
            radmc3d.write.control(self.thermalparams.control, thermpath=thermpath)

        if line==True:
            print('\nWriting lines.inp:')
            print('----------------------------')
            radmc3d.write.lines(species=species, format=line_format, thermpath=thermpath)

        if gasvelocity == True:
            print('\nWriting gas_velocity.inp:')
            print('----------------------------')
            radmc3d.write.gas_velocity(star_mass=star_mass, r=self.grid.r, theta=self.grid.theta, phi=self.grid.phi, object="disk", thermpath=thermpath)

        if gastemp == '1D_static':
            print('\nWriting gas_temperature.inp from 1D_static.dat:')
            print('----------------------------')
            chempath = self.chempath

            # Discover chemistry radii from folder names
            au_folders = sorted(
                [d for d in os.listdir(chempath)
                 if d.endswith('AU') and os.path.isdir(os.path.join(chempath, d))],
                key=lambda x: int(x.replace('AU', ''))
            )
            if not au_folders:
                print('  [write_line] no AU/ folders found in chempath — skipping gas_temperature.inp')
            else:
                # Build a chemmodel-style dict for Tg, keyed by radius in AU
                # Columns in 1D_static.dat: z(0) nH(1) Tg(2) Av(3) ...
                tg_chemmodel = {}
                for folder in au_folders:
                    r_au = int(folder.replace('AU', ''))
                    static_file = os.path.join(chempath, folder, '1D_static.dat')
                    if not os.path.isfile(static_file):
                        print(f'  [write_line] {folder}/1D_static.dat not found — skipping')
                        continue
                    data = np.loadtxt(static_file, comments='!')
                    tg_chemmodel[r_au] = {
                        'z':  data[:, 0],   # AU, surface → midplane (descending)
                        'Tg': data[:, 2],   # gas temperature [K]
                    }

                if not tg_chemmodel:
                    print('  [write_line] no valid 1D_static.dat files found — skipping')
                else:
                    # Read spherical grid (same pattern as convert_nautilus2radmc)
                    nx, ny, nz, x, y, z = radmc3d.read.grid(thermpath)
                    self.grid.set_spherical_grid(r_edge=np.array(x)/autocm, theta_edge=y, phi_edge=z)
                    x, y = self.grid.r, self.grid.theta

                    tg_sph = nautilus.coupling.to_spherical(
                        tg_chemmodel, nx, ny, x, y, struct='Tg'
                    )
                    # Replace zero/negative cells (inner cavity, outer edge) with a floor of 10 K
                    tg_sph = np.where(tg_sph <= 0, 10.0, tg_sph)
                    radmc3d.write.gas_temperature(tg_sph, thermpath=thermpath)


    def write_nautilus(self, sizes=np.array([[0.1]]),
                       uv_ref=3400,
                       nH_to_AV_conversion=1.600e+21,
                       rsingle=0.1,
                       dtogas=1e-2,
                       ref_radius=100,
                       stop_time=3e6,
                       nb_outputs = 64,
                       tunneling=1,
                       is_h2_formation_rate=0,
                       min_gas_density=1e0,
                       min_av=1e-2,
                       max_uv=None,
                       cap_uv_floor=True,
                       cut_cap=True,
                       max_inv_ab=None,
                       exclude_bins=None,
                       temp_gas='dust',
                       static=True,
                       param=True,
                       element=True,
                       abundances='atomic',
                       network=True,
                       multi_grain=True,
                       tempdecoup=True,
                       coupling_dens=False,
                       coupling_temp=True,
                       coupling_av=True,
                       **keywords):
        """Write Nautilus chemistry input files for each radial point.

        Reads the RADMC-3D thermal output (grid, dust density, dust
        temperature, local radiation field), couples the physical quantities
        onto the chemistry grid, and writes one set of Nautilus input files
        per radial grid point inside the chemistry directory. The existing
        chemistry directory is removed and recreated.

        When ``coupling_dens`` is True the code derives dust and gas number
        densities from the RADMC-3D dust density rather than from an
        analytically added disk/envelope structure.

        Parameters
        ----------
        sizes : ndarray, optional
            Grain size array in microns, shape ``(n_structures, n_bins)``.
            Default is ``np.array([[0.1]])``.
        uv_ref : float, optional
            Reference UV field strength in Habing units. Default is 3400.
        nH_to_AV_conversion : float, optional
            Column density to visual extinction conversion factor in
            cm^-2. Default is 1.6e21.
        rsingle : float, optional
            Single representative grain radius in microns. Default is 0.1.
        dtogas : float, optional
            Dust-to-gas mass ratio. Default is 1e-2.
        ref_radius : float, optional
            Reference radius in AU for the UV scaling. Default is 100.
        stop_time : float, optional
            Chemical evolution stop time in years. Default is 3e6.
        nb_outputs : int, optional
            Number of time outputs written by Nautilus. Default is 64.
        tunneling : int, optional
            Tunneling flag for Nautilus (0 or 1). Default is 1.
        is_h2_formation_rate : int, optional
            H2 formation rate flag for Nautilus. Default is 0.
        min_gas_density : float, optional
            Minimum gas number density in cm^-3 enforced in the chemistry
            input. Default is 1.
        min_av : float, optional
            Minimum visual extinction in mag enforced in the chemistry
            input. Default is 1e-2.
        max_uv : float or None, optional
            Maximum UV field value. If None, no cap is applied. Default is
            None.
        cap_uv_floor : bool, optional
            If True (default), cap the UV factor at 10 Habing for cells at
            the gas density floor. Set to False to keep the full geometric
            UV value in those cells (removes the discontinuity at the
            floor-density boundary visible when plotting the uv column).
        temp_gas : str, optional
            Gas temperature source: 'dust' uses the area-weighted dust
            temperature, 'param' uses a parametrized gas temperature added
            to the grid. Default is 'dust'.
        static : bool, optional
            Write the ``1D_static.dat`` (or multi-grain equivalent) input
            file. Default is True.
        param : bool, optional
            Write the Nautilus ``parameters.in`` file. Default is True.
        element : bool, optional
            Write the ``element.in`` file. Default is True.
        abundances : str, optional
            Initial abundances preset name (e.g., 'atomic'). Default is
            'atomic' for solar-composition. It can also be a filepath.
        network : bool or str, optional
            Controls network file copying. ``False`` skips it entirely.
            ``True`` (default) copies from the built-in
            ``astromugs/nautilus/network/`` directory. A string path copies
            from that directory instead (e.g.
            ``network='/path/to/my_network'``).
        multi_grain : bool, optional
            If True, use multi-grain chemistry mode (NMGC). Default is
            True.
        tempdecoup : bool, optional
            If True, compute a separate area-weighted dust temperature
            from the full size distribution. If False, use the single-bin
            temperature. Default is True.
        coupling_dens : bool, optional
            If True, derive dust and gas densities from the RADMC-3D dust
            density file. Default is False.
        coupling_temp : bool, optional
            If True, interpolate dust temperature from the RADMC-3D output
            onto the chemistry grid. Default is True.
        coupling_av : bool, optional
            If True, compute visual extinction from the local radiation
            field. Default is True.
        **keywords
            Any keyword argument accepted by
            :func:`~astromugs.nautilus.write.parameters_nmgc` can be passed
            here and will be forwarded verbatim to that function. This lets
            you override any entry in the ``parameters.in`` file without
            subclassing. See *Other Parameters* below for the full list.

        Other Parameters
        ----------------
        phase : int, optional
            Chemical phase model. ``0`` = 2-phase, ``1`` = 3-phase
            (default ``1``).
        cr_ionisation_rate : float, optional
            Cosmic-ray ionisation rate [s⁻¹]. Default ``1.3e-17``.
        x_ionisation_rate : float, optional
            X-ray ionisation rate [s⁻¹]. Default ``0.0``.
        is_photodesorb : int, optional
            Photodesorption of ices: ``1`` = on (default), ``0`` = off.
        is_crid : int, optional
            Cosmic-ray induced diffusion (CRID): ``1`` = on, ``0`` = off
            (default).
        is_er_cir : int, optional
            Eley-Rideal and complex-induced reactions: ``1`` = on, ``0`` = off
            (default).
        is_absorption_h2 : int, optional
            H₂ self-shielding (Lee & Herbst 1996): ``1`` = on (default).
        is_absorption_co : int, optional
            CO self-shielding method. ``1`` = Lee & Herbst (1996),
            ``2`` = Visser et al. (2009) (default ``2``).
        is_absorption_n2 : int, optional
            N₂ self-shielding (Li et al. 2013): ``1`` = on (default).
        nb_active_lay : float, optional
            Number of chemically active grain surface layers. Default ``2.0``.
        diff_binding_ratio_surf : float, optional
            Ratio of diffusion barrier to binding energy for surface species.
            Default ``0.4``.
        diff_binding_ratio_mant : float, optional
            Ratio of diffusion barrier to binding energy for mantle species.
            Default ``0.8``.
        cr_peak_grain_temp : float, optional
            Peak grain temperature reached during a cosmic-ray heating event
            [K]. Default ``70.0``.
        cr_peak_duration : float, optional
            Duration of a cosmic-ray heating peak [s]. Default ``1e-5``.
        Fe_ionisation_rate : float, optional
            Cosmic Fe-ion–grain encounter rate [s⁻¹ grain⁻¹]. Default
            ``3e-14``.
        sticking_coeff_neutral : float, optional
            Sticking coefficient for neutral species. Default ``1.0``.
        sticking_coeff_positive : float, optional
            Sticking coefficient for positively charged species. Default
            ``0.0``.
        sticking_coeff_negative : float, optional
            Sticking coefficient for negatively charged species. Default
            ``0.0``.
        surface_site_density : float, optional
            Surface site density on a grain [cm⁻²]. Default ``8e14``.
        chemical_barrier_thickness : float, optional
            Grain reaction activation energy barrier width [cm]. Default
            ``1e-8``.
        relative_tolerance : float, optional
            Relative tolerance of the ODE solver. Default ``1e-4``.
        output_type : str, optional
            Output time spacing: ``'log'`` (default) or ``'linear'``.
        start_time : float, optional
            First output time [yr]. Default ``1.0``.
        modify_rate_flag : int, optional
            Rate modification flag. ``0`` = none (default), ``1`` = H only,
            ``2`` = H + H₂, ``3`` = all species.
        ED_H2 : float, optional
            H₂ binding energy on itself [K]. Default ``23.0``.
        vib_to_dissip_freq_ratio : float, optional
            Ratio of surface-bond vibrational frequency to energy dissipation
            frequency (RRK desorption mechanism). Default ``1e-2``.
        """

        #-----------------------------------------
        # REMOVE IF EXISTS AND CREATE CHEMISTRY FOLDER
        #-----------------------------------------
        thermpath = self.thermalpath
        chempath = self.chempath
        if os.path.exists(chempath):
            shutil.rmtree(chempath)
        os.makedirs(chempath)


        #-----------------------------------------
        # READ THERMAL FILES (just like in function thermal, but sometimes the user does not need to run thermal so this reads the thermal files even though they might already be read by thermal)
        #-----------------------------------------
        nx, ny, nz, x, y, z  = radmc3d.read.grid(thermpath) # read grid

        # x is r_edges in cm, y is theta_edges, z is phi_edges
        self.grid.set_spherical_grid(r_edge=np.array(x)/autocm, theta_edge=y, phi_edge=z)

        dust_density = radmc3d.read.dust_density(thermpath) # read dust_density file
        nbspecies = int(len(dust_density[0])/(nx*ny*nz)) #get number of species
        dust_density = np.reshape(dust_density[0], (nbspecies, nz, ny, nx)) #reshape it
        dust_density[dust_density<=1e-100] = 1e-100 # get a minimum dust density in order to avoid absurd NaN issues.

        nb_lam, lam, r_star, m_star, T_star, spectrum = radmc3d.read.stars(thermpath)  #read stars file.
        external = radmc3d.read.external_source(thermpath)  #read external_source.inp
        # define temperature
        self.grid.temperature = radmc3d.read.dust_temperature(thermpath)
        if len(self.grid.temperature) > 0:
            self.grid.temperature[0] = np.reshape(self.grid.temperature[0], (nbspecies, nz, ny, nx))

            a = []
            dustmass = []
            for istruc in range(0, len(self.grid.dust)): # create a loop in order to gather all grain sizes into one single array that will be used to compute the area-weighted temperature Ta.
                dustmodel = self.grid.dust[istruc]
                sizes = dustmodel.sizes()
                rho_m = dustmodel.rho_m
                #a = sizes[-1]*1e-4 # cm
                a.append(sizes[-1]*1e-4) # convert to cm
                dustmass.append(dustmodel.grainmass()) # gram

            a = np.hstack(a) # in order to get a numpy list of all sizes no matter how many structures there are.
            dustmass = np.hstack(dustmass)

            if tempdecoup == True:
                Ta_num = np.zeros((nz, ny, nx))
                Ta_denum = np.zeros((nz, ny, nx))

                for idx in range(len(a)):
                    Ta_num += self.grid.temperature[0][idx, :, :, :]*(dust_density[idx, :, :, :]/dustmass[idx])*a[idx]**2
                    Ta_denum += (dust_density[idx, :, :, :]/dustmass[idx])*a[idx]**2
                self.Ta = Ta_num/Ta_denum  #Ta is the area-weighted dust temperature.
            else:
                self.Ta = self.grid.temperature[0][0,:,:,:]
        else:
            print('\nNo dust temperature file was found. If coupling_temp is True the chemistry model cannot be created.\n\n')

        nlam_mono, lam_mono, self.grid.localfield = radmc3d.read.localfield(thermpath)
        lam_mono = 1e6*(c/lam_mono) #from Hz to microns. Should be same number as in mcmono_wavelength.inp.

        isrf = InterstellarRadFields(cut=2.e-1, d78=True, vdb82=False) #calculate Draine isrf.


        if len(self.grid.localfield) > 0:
            try:
                self.grid.localfield = np.reshape(self.grid.localfield, (nlam_mono, nz, ny, nx))
            except IOError:
                print('\nPlease, check consistency between the grid size and mean_intensity.out.\n')
                sys.exit(1)
        else:
            print('Warning: No mean_intensity.out file was found. If coupling_av is True then the chemistry model cannot be created. Please check mean_intensity.out file or set coupling_av=False.')


        #-----------------------------------------
        # COUPLING ROUTINES
        #-----------------------------------------
        if coupling_dens == True:
            n_dust, n_gas = nautilus.coupling.dust_density(dtogas, rho_m, a, dust_density, self.grid.rchem*autocm, self.grid.zchem*autocm, self.grid.r, self.grid.theta)

            # If custom sigma, override n_gas with gas density from the custom file, because the n_gas is derived from dtogas in the coupling, which is wrong if sigma_gas is custom.
            if self.params.disk.sigma_compute == 'custom' and len(self.grid.hg_chem) == 0:
                r_custom, _, siggas_table, _ = custom_io.surfacedensities(self.params.disk.sigma_path)
                temp_disk = Disk(params=self.params.disk, dust=self.grid.dust[0])
                hg = temp_disk.scaleheight(self.grid.rchem)  # cm
                for idx, r in enumerate(self.grid.rchem):
                    sig_g = np.interp(r, r_custom, siggas_table)
                    z_cm = self.grid.zchem[idx, :] * autocm
                    n_gas[idx, :] = (sig_g / (np.sqrt(2*np.pi) * hg[idx] * mu * amu)) * np.exp(-z_cm**2 / (2*hg[idx]**2))
            elif self.params.disk.sigma_compute != 'custom':
                print(
                    "\nWARNING: coupling_dens=True but self.params.disk.sigma_compute = "
                    f"'{self.params.disk.sigma_compute}' (default is 'model', not 'custom').\n"
                    "Gas density will be set to (total dust mass) / dtogas everywhere -- a "
                    "fixed dust-to-gas ratio -- NOT from any independent gas surface density.\n"
                    "Did you mean to use a custom gas profile but forgot to re-run your disk "
                    "parameter setup (e.g. pipe1.params.disk.sigma_compute = 'custom' and "
                    ".sigma_path), for instance after restarting the kernel and reusing "
                    "existing RADMC-3D files? Set those before calling write_nautilus() if so.\n"
                )
        if coupling_temp == True:
           if not self.grid.temperature:
               print('coupling_temp==True: The file thermal/dust_temperature.dat is not present or is corrupted. Chemistry model cannot created.')
               sys.exit(1)
           if len(self.grid.hg_chem) > 0:
               T_dust = nautilus.coupling.dust_temperature_disk(self.grid.temperature[0], self.grid.rchem*autocm, self.grid.zchem, self.grid.r, self.grid.theta, self.grid.hg_chem[0]) # dim(a, rchem, zchem)
               T_dust_single = nautilus.coupling.dust_temperature_single_disk(self.Ta, self.grid.rchem*autocm, self.grid.zchem, self.grid.r, self.grid.theta, self.grid.hg_chem[0])
           else:
               T_dust = nautilus.coupling.dust_temperature(self.grid.temperature[0], self.grid.rchem*autocm, self.grid.zchem*autocm, self.grid.r, self.grid.theta) # dim(a, rchem, zchem)
               T_dust_single = nautilus.coupling.dust_temperature_single(self.Ta, self.grid.rchem*autocm, self.grid.zchem*autocm, self.grid.r, self.grid.theta)
        else:
           T_dust = np.expand_dims(self.grid.tgas_chem[0], axis=0) # expend to one extra dimension in order to match the shape of coupled T_dust.
        if coupling_av == True:
            if len(self.grid.hg_chem) > 0:
                av_z = nautilus.coupling.avz_disk(self.grid.localfield, lam_mono, r_star, T_star, self.grid.rchem*autocm, self.grid.zchem, self.grid.r, self.grid.theta, self.grid.hg_chem[0]) # dim(rchem, zchem)
            else:
                av_z = nautilus.coupling.av_z(self.grid.localfield, lam_mono, r_star, T_star, self.grid.rchem*autocm, self.grid.zchem*autocm, self.grid.r, self.grid.theta) # dim(rchem, zchem)
        else:
            av_z = self.grid.avz[0]  #works only if a disk or envelope model is created before.

        #-----------------------------------------
        # WRITE NAUTILUS INPUT FILES FOR EACH RADIUS
        #-----------------------------------------
        print('WRITING NAUTILUS INPUT FILES:')
        print('-----------------------------\n')
        # write input NAUTILUS files
        if multi_grain == True:  #if multi_grain is True, then the code assumes the user will use the NMGC version of Nautilus in multi-grain mode.
             print('Multi-grain mode: Yes\n')
             print('Dust temperature structure: multiple dust temperatures (dust parameters stored in 1D_grain_sizes.in).\n')
        else:
             print('Multi-grain mode: No\n') #if ngmc is False, Nautilus with a single grain bin is assumed to be used. This can be true even if the thermal model has multiple dust bins.
             print('Dust temperature structure: single or area-weigthed (dust parameters stored in 1D_static.dat).\n')


        for idx, r in enumerate(self.grid.rchem):
            path = chempath + '/' + str(int(r)) + 'AU/'
            os.makedirs(path, exist_ok=False)

            #---temporary defining cavity to increase the Av in that area so we don't have convergence issue. This should be removed in the main branch.
            z0 = 200
            phi = (16*np.pi)/180
            zcav = z0*(r/(z0*np.tan(phi/2)))**1.55



            ###!!!! maybe add a factor to account for the excess of UV in the protostar spectrum?
            if len(self.grid.hg_chem) > 0:
                uvfactor = nautilus.write.uv_factordisk(uv_ref, ref_radius, r, self.grid.hg_chem[0][idx]/autocm)
            else:
                uvfactor = nautilus.write.uv_factor(isrf, lam_mono, r_star, T_star, r*autocm, self.grid.zchem[idx,:]*autocm, external)

            avnh_fact = nautilus.write.avnh_factor(nH_to_AV_conversion, dtogas, rsingle, self.grid.nz_chem)

            if temp_gas == 'dust':
                T_gas = T_dust_single[idx,:]
            elif temp_gas == 'param':
                if len(self.grid.tgas_chem) == 0:
                    print("because you set temp_gas='param', you have to add a parametrized Tgas in the model first.")
                    sys.exit(1)
                else:
                    T_gas = self.grid.tgas_chem[0][idx]

            if multi_grain == False:
                nz_actual = self.grid.nz_chem  # default if static not written
                if static == True:
                    if len(self.grid.hg_chem) > 0:
                        dist = self.grid.zchem*self.grid.hg_chem[0][idx]/autocm
                        nH = self.grid.gasdensity_chem[0][idx,:]
                        nd = self.grid.dustdensity_single_chem[0][idx,:]
                    elif coupling_dens == True:
                        dist = self.grid.zchem[idx, :]
                        nH = n_gas[idx, :]
                        md = (4/3)*np.pi*rho_m*(rsingle*1e-4)**3
                        nd = dtogas*n_gas[idx, :]*amu*mu/md
                    nz_actual = nautilus.write.static(path, \
                                    dist, \
                                    nH, \
                                    T_gas, \
                                    av_z[idx, :], \
                                    T_dust_single[idx,:], \
                                    nd, \
                                    rsingle, \
                                    avnh_fact,
                                    uvfactor,
                                    min_gas_density=min_gas_density,
                                    min_av=min_av,
                                    max_uv=max_uv,
                                    cap_uv_floor=cap_uv_floor,
                                    cut_cap=cut_cap,
                                    rho_m=rho_m)
                if param == True:
                    nautilus.write.parameters_nmgc(path, grain_temp='table_1D', nb_outputs=nb_outputs, multi_grain=0, tunneling=tunneling, is_h2_formation_rate=is_h2_formation_rate, resolution=nz_actual, stop_time=stop_time, uv_flux=np.mean(uvfactor), **keywords)
            if multi_grain == True:
                nz_actual = self.grid.nz_chem  # default if static not written
                if static == True:
                    if len(self.grid.hg_chem) > 0:
                        dist = self.grid.zchem*self.grid.hg_chem[0][idx]/autocm
                        nH = self.grid.gasdensity_chem[0][idx,:]
                        nd = self.grid.dustdensity_single_chem[0][idx,:]
                    elif coupling_dens == True:
                        dist = self.grid.zchem[idx, :]
                        nH = n_gas[idx, :]
                        md = (4/3)*np.pi*rho_m*(rsingle*1e-4)**3
                        nd = dtogas*n_gas[idx, :]*amu*mu/md
                    ##!!!!!!! TO BE REMOVED !!!!!!!
                    #av_z[idx, :] = np.where(dist>zcav, av_z[idx, :]*10, av_z[idx, :])
                    ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    nz_actual = nautilus.write.static(path, \
                                    dist, \
                                    nH, \
                                    T_gas, \
                                    av_z[idx, :], \
                                    T_dust_single[idx,:], #if nbspecies > 1, the column is not read.\
                                    nd, #if nbspecies > 1, the column is not read. \
                                    rsingle, #if nbspecies > 1, the column is not read. \
                                    avnh_fact,
                                    uvfactor,
                                    min_gas_density=min_gas_density,
                                    min_av=min_av,
                                    max_uv=max_uv,
                                    cap_uv_floor=cap_uv_floor,
                                    cut_cap=cut_cap,
                                    rho_m=rho_m)
                if param == True:
                    nautilus.write.parameters_nmgc(path, grain_temp='fixed_to_dust_size', nb_outputs=nb_outputs, multi_grain=1, resolution=nz_actual, tunneling=tunneling, is_h2_formation_rate=is_h2_formation_rate, stop_time=stop_time, uv_flux=np.mean(uvfactor), **keywords)

                if nbspecies > 1:
                    if len(self.grid.hg_chem) > 0:
                        nH = self.grid.gasdensity_chem[0][idx,:]
                        nd = self.grid.dustdensity_chem[0][:,idx,:]
                    elif coupling_dens == True:
                        nH = n_gas[idx, :]
                        nd = n_dust[:, idx, :]
                    nautilus.write.grain_sizes(path, sizes, nH, nd, T_dust[:,idx,:], min_gas_density=min_gas_density, cut_cap=cut_cap, max_inv_ab=max_inv_ab, exclude_bins=exclude_bins, rho_m=rho_m)
                else:
                    print('WARNING: multi_grain = True, but the model has only one grain bin. Please, check the number of grain size or switch multi_grain to False.')

            nautilus.write.abundances(path, abundances)
            if network:
                network_path = network if isinstance(network, str) else None
                nautilus.write.network(path, network_path=network_path)
            if element == True:
                nautilus.write.elements(path)
