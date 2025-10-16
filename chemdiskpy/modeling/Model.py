import glob, os, sys, shutil 
import numpy as np

from .. import radmc3d
from .. import nautilus
from .Grid import Grid

from .. constants.constants import autocm, M_sun, R_sun, c, amu, mu, black_body
from ..modeling.InterstellarRadFields import InterstellarRadFields

import matplotlib.pyplot as plt


class Model:

    def __init__(self):
        self.grid = Grid()
        self.nautilus = nautilus

    def run_continuum(self, nphot=1e4, \
                            write_control=False, \
                            **keywords):
        """ 
        Notes:
        run MC dust radiative transfer, open the resulting dust temperature as an array and computes the surface-area weigthed temperature. If run == False, user assumes the RADMC3D output files already exist.
        -----
	    """	
        # FOLDER WITH THERMAL FILES
        thermpath='thermal/' 

        # WRITE THERMAL FILES
        if write_control == True:
            self.write_continuum(nphot_therm=nphot, \
                            control=True, \
                            **keywords)

        self.run_thermal_radmc3d(nphot=nphot, **keywords)




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
            self.grid.set_spherical_grid(x[0], x[-1], nx+1, ny+1, nz+1) #return self.grid.r, self.grid.theta. Necessary in case the user does not create themeselves the radmc3d files.
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
        thermpath='thermal/'

        if make_image == True:
            self.run_image_radmc3d(incl=incl, npix=npix, iline=iline,  lambda_micron=lambda_micron, **keywords)




    def convert_nautilus2radmc(self, species='CO', dust_density=False, dust_temperature=False, numberdens=False):
        """
        -Convert the chemistry grid into radmc3d files. It needs a chemistry model to be run. 
        -The conversion from nautilus to radmc3f files can be done for dust_density.inp, dust_temperature.inp, and numberdens_XX.inp. The user can choice to create all or 1 or 2 files with the flags.
        -The chemistry model has to added to the object using add_chemmodel() function. 
        -note that the conversion to radmc3d REQUIRES the presence of amr_grid.inp file first. Indeed, chemistry models usually have a lower resolution than that for radiative transfer and it is therefore better to define the chosen RT grid. The code will do the rest (interpolation, etc.). 


        Args:
            chempath (str): Relative or absolute path to the chemistry model directory.
            itime (int): Index of the time output to use from the chemistry model. Not yet coded.
            species (str) or list of str: Chemical species of interest (e.g., 'CO', 'H2O') for the creation of numberdens_XX.inp. 
            reader (object, optional): A reader object or module responsible for reading chemistry data. 
                                    Defaults to self.nautilus.read.
        Raises:

        """
        thermpath='thermal/'

        if os.path.isfile(thermpath + 'amr_grid.inp'):
            nx, ny, nz, x, y, z  = radmc3d.read.grid(thermpath) # read grid. the amr_grid.inp shows the border of the cells so we convert them.
            self.grid.set_spherical_grid(x[0]/autocm, x[-1]/autocm, nx+1, ny+1, nz+1, log=True)
            x, y, z = self.grid.r, self.grid.theta, self.grid.phi
            y[0] = 0
            y[-1] = np.pi  # to avoid numerical issues.
        elif hasattr(self.grid, 'r_edge') and self.grid.r_edge is not None:
            nx, ny, nz, x, y, z = self.grid.nr, self.grid.ntheta, self.grid.nphi, self.grid.r, self.grid.theta, self.grid.phi
            nx = nx-1
            ny = ny-1
            nz = nz-1
        else:
            print(f"Error: amr_grid not found. Make sure there is one (example: create one using grid.set_spherical_grid(rin, rout, nr, ntheta, nphi, log=True))")
            sys.exit(1)  # Exit the program with a non-zero status to indicate an error

        if numberdens==True:
            # r_naut, zz_naut, dens_mol_nautilus = nautilus.read.abundance(path=path, itime=itime, species=species) #read abundances of chosen species from Nautilus output files.
            #numberdens_sph = nautilus.coupling.to_spherical(dens_mol_nautilus, nx, ny, nz, x, y, r_naut, zz_naut) #convert abundances into spherical grid.
            numberdens_sph = nautilus.coupling.to_spherical(self.grid.chemmodel[species], nx, ny, x, y, struct='numberdens_species') #convert abundances into spherical grid.
            radmc3d.write.numberdens_mol(numberdens_sph, species=species, gridstyle="regular") #write numberdens_mol.inp file for RADMC-3D.

        if dust_density==True:
            print('create dust_density.inp')

        if dust_temperature == True:
            print('create dust_temperature.inp')


                


    def localfield(self, nphot_mono=1e6, write_mcmono=False, run=True, **keywords):

        self.write_continuum(write=True,\
                           write_dens=False, \
                           write_grid=False, \
                           write_opac=False, \
                           write_control=False, \
                           write_star=False, \
                           write_wave=False, \
                           write_mcmono=write_mcmono, \
                           write_ext=False, **keywords)
        
        if run == True:
            self.run_localfield_radmc3d(nphot_mono=nphot_mono)

        thermpath='thermal/'
        nx, ny, nz, x, y, z  = radmc3d.read.grid(thermpath)
        nlam_mono, lam_mono, self.grid.localfield = radmc3d.read.localfield(thermpath)
        lam_mono = 1e6*(c/lam_mono)
        self.grid.localfield = np.reshape(self.grid.localfield, (nlam_mono, nz, ny, nx))


    def run_thermal_radmc3d(self, nphot=1e6, verbose=True, timelimit=7200, \
            nice=None, **keywords):
        radmc3d.run.thermal(verbose=verbose, timelimit=timelimit, nice=nice)


    def run_localfield_radmc3d(self, nphot_mono=1e6, verbose=True, timelimit=7200):
        radmc3d.run.localfield(nphot_mono=nphot_mono, verbose=verbose, timelimit=timelimit)


    def run_image_radmc3d(self, npix=300, lambda_micron=None, iline=None, incl=None, verbose=True):
        radmc3d.run.image(npix=npix, lambda_micron=lambda_micron, iline=iline, incl=incl, verbose=verbose, timelimit=7200)


    def write_continuum(self, dens=False, grid=False, opac=False, control=False, stars=False, wave=False, mcmono=False, ext=False, **keywords):
        #os.system("rm thermal/*.inp")


        if not os.path.exists('thermal'):
            os.makedirs('thermal')

        if control==True:
            print('\nWriting radmc3d.inp:')
            print('----------------------------')
            radmc3d.write.control(**keywords)

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
                    tstar=tstar)

        if wave==True:
            print('\nWriting wavelength_micron.inp:')
            print('----------------------------')
            radmc3d.write.wavelength_micron(self.grid.lam)

        if mcmono==True:
            print('\nWriting mcmono_wavelength_micron.inp:')
            print('----------------------------')
            radmc3d.write.mcmono_wavelength_micron(self.grid.monolam)

        if ext==True:
            if len(self.grid.isrf) != 0:
                radmc3d.write.external_rad(self.grid.isrf[0])

        if grid==True:
            print('\nWriting amr_grid.inp:')
            print('----------------------------')
            if self.grid.coordsystem == 'spherical':
                radmc3d.write.amr_grid(self.grid.r_edge*autocm, self.grid.theta_edge, self.grid.phi_edge, gridstyle="regular", coordsystem=self.grid.coordsystem)

        if dens==True:
            print('\nWriting dust_density.inp:')
            print('----------------------------')
            radmc3d.write.dust_density(self.grid.dustdensity, gridstyle="regular")

        if opac==True:
            print('\nWriting dustopac.inp:')
            print('----------------------------')
            dustopac = []
            filelist = glob.glob('thermal/dustkap*')
            for files in sorted(filelist):
                dustopac.append(files)
            radmc3d.write.dustopac(dustopac)
            
        if len(self.grid.accretionheating) > 0:
            radmc3d.write.accretion_heating(self.grid.w1*autocm, self.grid.w2, self.grid.w3, self.grid.accretionheating[0], gridstyle="regular")

        if dens == False and grid == False and control == False and opac == False and stars == False and wave == False: 
            print('\nWARNING: no RADMC3D file are created.\n')
            print('----------------------------\n')
                
    def write_line(self, control=False, line=False, gasvelocity=False, gastemp=False, microturb=False, line_format='leiden', species='CO', star_mass=1, **keywords):
        #os.system("rm thermal/*.inp")

        if not os.path.exists('thermal'):
            os.makedirs('thermal')

        if control==True:
            print('\nWriting radmc3d.inp:')
            print('----------------------------')
            radmc3d.write.control(**keywords)

        if line==True:
            print('\nWriting line.inp:')
            print('----------------------------')
            radmc3d.write.lines(species=species, format=line_format)
            
        if gasvelocity == True:
            radmc3d.write.gas_velocity(star_mass=star_mass, r=self.grid.r, theta=self.grid.theta, phi=self.grid.phi, object="disk")   


    # ----WRITE NAUTILUS INPUT FILES----
    # ----------------------------------
    #If coupling_dens == True, the code assumes the user did not add a disk or envelope structure and that there is no calculation of Hg, ng, nd, etc. 
    def write_nautilus(self, sizes=np.array([[0.1]]), uv_ref=3400, nH_to_AV_conversion=1.600e+21, rsingle=0.1, dtogas=1e-2, ref_radius=100,\
                       stop_time=3e6, nb_outputs = 64, tunneling=1, is_h2_formation_rate=0, temp_gas='dust', static=True, param=True, element=True, abundances='atomic', \
                       network=True, multi_grain=True, tempdecoup=True,\
                       coupling_dens=False, coupling_temp=True, coupling_av=True, **keywords):
        
        #-----------------------------------------
        # REMOVE IF EXISTS AND CREATE CHEMISTRY FOLDER
        #-----------------------------------------
        thermpath='thermal/'
        chempath = self.chempath
        if os.path.exists(chempath):
            shutil.rmtree(chempath)
        os.makedirs(chempath)


        #-----------------------------------------
        # READ THERMAL FILES (just like in function thermal, but sometimes the user does not need to run thermal so this reads the thermal files even though they might already be read by thermal)
        #-----------------------------------------
        nx, ny, nz, x, y, z  = radmc3d.read.grid(thermpath) # read grid
        print(type(x), type(y))
        self.grid.set_spherical_grid(x, y, z) #return self.grid.r, self.grid.theta. Necessary in case the user does not create themeselves the radmc3d files.

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
        # isrf = isrf.create_isrf(lam_mono)
        # print(isrf)
        # # #--------------------------------
        # import matplotlib.pyplot as plt
        # from matplotlib.colors import LogNorm
        # fig = plt.figure(figsize=(10, 8.))
        # ax = fig.add_subplot(111)
        # plt.xlabel(r'lamda', fontsize = 17)
        # plt.ylabel(r'isrf', fontsize = 17, labelpad=-7.4)
        # plt.semilogy(isrf[0], isrf[1], c='black')
        # ax.set_ylim(1e-26, 1e-20)
        # plt.show()
        # # #-----------------------------------



        if len(self.grid.localfield) > 0:
            try:
                self.grid.localfield = np.reshape(self.grid.localfield, (nlam_mono, nz, ny, nx)) 
            except IOError:
                print('\nPlease, check consistency between the grid size and external_source size.\n')
                sys.exit(1)
        else:
            print('Warning: No external_source file was found. If coupling_av is True then the chemistry model cannot be created. Please check external_source file or set coupling_av=False.')


        #-----------------------------------------
        # COUPLING ROUTINES
        #-----------------------------------------
        if coupling_dens == True:
            n_dust, n_gas = nautilus.coupling.dust_density(dtogas, rho_m, a, dust_density, self.grid.rchem*autocm, self.grid.zchem*autocm, self.grid.r, self.grid.theta)
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
            path = 'chemistry' + '/' + str(int(r)) + 'AU/'
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
                if param == True:
                    nautilus.write.parameters_nmgc(path, grain_temp='table_1D', nb_outputs=nb_outputs, multi_grain=0, tunneling=tunneling, is_h2_formation_rate=is_h2_formation_rate, resolution=self.grid.nz_chem, stop_time=stop_time, uv_flux=np.mean(uvfactor), **keywords)
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
                    nautilus.write.static(path, \
                                    dist, \
                                    nH, \
                                    T_gas, \
                                    av_z[idx, :], \
                                    T_dust_single[idx,:], \
                                    nd, \
                                    rsingle, \
                                    avnh_fact,
                                    uvfactor)
            if multi_grain == True:
                if param == True:
                    nautilus.write.parameters_nmgc(path, grain_temp='fixed_to_dust_size', nb_outputs=nb_outputs, multi_grain=1, resolution=self.grid.nz_chem, tunneling=tunneling, is_h2_formation_rate=is_h2_formation_rate, stop_time=stop_time, uv_flux=np.mean(uvfactor), **keywords)
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
                    av_z[idx, :] = np.where(dist>zcav, av_z[idx, :]*10, av_z[idx, :])
                    ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    nautilus.write.static(path, \
                                    dist, \
                                    nH, \
                                    T_gas, \
                                    av_z[idx, :], \
                                    T_dust_single[idx,:], #if nbspecies > 1, the column is not read.\
                                    nd, #if nbspecies > 1, the column is not read. \
                                    rsingle, #if nbspecies > 1, the column is not read. \
                                    avnh_fact,
                                    uvfactor)
            
                if nbspecies > 1: 
                    if len(self.grid.hg_chem) > 0:
                        nH = self.grid.gasdensity_chem[0][idx,:]
                        nd = self.grid.dustdensity_chem[0][:,idx,:]
                    elif coupling_dens == True: 
                        nH = n_gas[idx, :]
                        nd = n_dust[:, idx, :]
                    nautilus.write.grain_sizes(path, sizes, nH, nd, T_dust[:,idx,:])
                else:
                    print('WARNING: multi_grain = True, but the model has only one grain bin. Please, check the number of grain size or switch multi_grain to False.')

            nautilus.write.abundances(path, abundances)
            if network == True:
                nautilus.write.network(path)
            if element == True:
                nautilus.write.elements(path)
            # if activ_energies == True:
            #     nautilus.write.activ_energies(path)
            # if surfaces == True:
            #     nautilus.write.surfaces(path)
