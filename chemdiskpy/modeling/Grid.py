"""
_____________________________________________________________________________________________________________
file name: Grid
@author: P.Sheehan. Adapted by S. Gavino for chemistry codes.
last update: Aug 2022
language: PYTHON 3.8
short description:  class Grid for young stellar objects modeling. 
_____________________________________________________________________________________________________________
"""

import numpy as np

class Grid:
    def __init__(self):
        self.density = []
        self.dustdensity = []
        self.gasdensity_chem = []
        self.dustdensity_chem = []
        self.dustdensity_single_chem = []
        self.hg_chem = []
        self.tgas_chem = []
        self.temperature = []
        self.localfield = []
        self.avz = []
        self.stars = []
        self.isrf = []
        self.dust = []
        self.accretionheating = []
        self.chemradii = []
        self.chemparam = []
        self.chemgrid = {}

    #-----
    # ADD/CREATE GRID STRUCTURES THAT EXIST OR THAT HAVE BEEN CREATED.
    #-----

    def add_star(self, star):
        self.stars.append(star)

    def add_isrf(self, isrf):
        self.isrf.append(isrf)

    def add_temperature(self, temperature):
        self.temperature.append(temperature)

    def add_localfield(self, localfield):
        self.localfield.append(localfield)

    def add_density(self, density):
        self.density.append(density)

    def add_dustdensity(self, density):
        self.dustdensity.append(density)

    def add_dustdensity_chem(self, density):
        self.dustdensity_chem.append(density)

    def add_dustdensity_single_chem(self, density):
        self.dustdensity_single_chem.append(density)

    def add_gasdensity_chem(self, density):
        self.gasdensity_chem.append(density)

    def add_gastemperature_chem(self, gas_temperature):
        self.tgas_chem.append(gas_temperature)

    def add_hg_chem(self, hg):
        self.hg_chem.append(hg)

    def add_avz(self, av_z):
        self.avz.append(av_z)

    def add_dust(self, dust):
        self.dust.append(dust)

    def add_accretionheating(self, q_visc):
        self.accretionheating.append(q_visc)

    def add_existingchemradii(self,existingchemradii):
        #args: chemgrid corresponds to an array with the radius and z points.
        self.chemradii.append(existingchemradii) 
    
    def add_existingchemparam(self,existingchemparam):
        #args: chemgrid corresponds to an array with the radius and z points.
        self.chemparam.append(existingchemparam)    

    def add_existingchemgrid(self,existingchemgrid, species):
        #args: chemgrid corresponds to an array with the radius and z points.
        self.chemgrid[species] = existingchemgrid    

        

    #-----
    # SET GRID STRUCTURES.
    #-----

    def set_cartesian_grid(self, xmin, xmax, nx):
        #w1, w2, w3 provide grid with coordinates using the center of each cell.
        self.coordsystem = "cartesian"

        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(xmin, xmax, nx)
        z = np.linspace(xmin, xmax, nx)

        w1 = 0.5*(x[0:x.size-1] + x[1:x.size])
        w2 = 0.5*(y[0:y.size-1] + y[1:y.size])
        w3 = 0.5*(z[0:z.size-1] + z[1:z.size])

        return np.stack((x, y, z)), np.stack((w1, w2, w3))

    def set_spherical_grid(self, rmin, rmax, nr, ntheta, nphi, log=True):
        self.coordsystem = "spherical"
        if log:
            rad = np.logspace(np.log10(rmin), np.log10(rmax), nr, base=10)
        else:
            rad = np.linspace(rmin, rmax, nr)

        theta_angle = np.linspace(0.0, np.pi, ntheta)
        phi_angle = np.linspace(0.0, 2*np.pi, nphi)

        self.r = 0.5*(rad[0:rad.size-1] + rad[1:rad.size])
        self.theta = 0.5*(theta_angle[0:theta_angle.size-1] + theta_angle[1:theta_angle.size])
        self.phi = 0.5*(phi_angle[0:phi_angle.size-1] + phi_angle[1:phi_angle.size])

        self.rad = rad
        self.theta_angle = theta_angle
        self.phi_angle = phi_angle

    def set_chemdisk_grid(self, r, max_H=4, nz_chem=64):
        """
        desc: Spatial grid for disk chemistry model. Defined from a upper limit for the disk atmosphere.
        args:
        -max_H: disk maximum gas scale height below which chemistry grid exists.
        -nz_chem: number of vertical spatial points. Same at all radii.
        -r: radii in au.
        """
        #hg = self.disk.scaleheight(np.array(r))
        pts = np.arange(0, nz_chem, 1)
        #zchem = np.ones((len(rchem), nb_points))

        #hh, ptpt = np.meshgrid(hchem, pts)
        z = (1. - (2.*pts/(2.*nz_chem - 1.)))*max_H#*Hg

        self.rchem = np.array(r)
        self.zchem = z
        self.nz_chem = nz_chem


    def set_chem_grid(self, r, z0=0, zmax=None, msize=None, nbcells=70):
        """
        desc: Custom spatial grid for chemistry model. Can be disk, envelope...
        args:
        -r [au]: radii in au. Must be 1D array of dim = number of radii. Must be inside the radmc3d model.
        -z0 [au]: minimum altitude. Zero by default.
        -zmax [au]: maximum altitude. 
        -msize [au]: modele size in AU if the model is spherical, the user can give msize in order to compute the zmax at each radius.
        -nbcells: number of vertical cells. Same for all radii.  
        """
        self.nz_chem = nbcells
 
        if zmax != None: #if user provides maximum altitude in au, it sets the maximum z at all radii. By default the minimum value z0 is 0. That creates a 2D structure.
            self.rchem = np.array(r)
            z = np.zeros((len(self.rchem), nbcells))
            for idx, rval in enumerate(self.rchem):
                z[idx,:] = np.linspace(z0, zmax, nbcells)
            self.zchem = np.flip(np.array(z), axis=1) #flip because the user gives increasing values and nautilus needs decreasing values.

        if msize != None: #if user wants the altitude max such that the model is a sphere i.e. zmax at each radius follows the spherical structure.
            zmax = np.sqrt(msize**2 - r**2)  # gives max altitude at each radius (polar coordinates).
            zmax = zmax[:-1] # remove the last value because it is zmax = 0 au.

            self.rchem = r[:-1] # because we removed the last zmax.
            z = np.zeros((len(self.rchem), nbcells))

            for idx, zmax_x in enumerate(zmax):
                z[idx,:] = np.linspace(0, zmax_x, nbcells)

            self.zchem = np.flip(z, axis=1) #flip because the user gives increasing values and nautilus needs decreasing values.

            #dz = np.tan(np.pi/nbthetas)*self.rchem #0.0175 is pi/number of thetas.

            # import matplotlib.pyplot as plt
            # fig = plt.figure(figsize=(8, 8.))
            # ax = fig.add_subplot(111)
            # ax.plot(self.rchem, zmax)
            # plt.xlim(0, 5005)
            # plt.ylim(0, 5005)
            # #plt.show()

    def set_wavelength_grid(self, lmin, lmax, nlam, log=False): #microns
        if log:
            self.lam = np.logspace(np.log10(lmin), np.log10(lmax), \
                    nlam)
        else:
            self.lam = np.linspace(lmin, lmax, nlam)


    def set_mcmonowavelength_grid(self, lmin, lmax, nlam, log=False):
        if log:
            self.monolam = np.logspace(np.log10(lmin), np.log10(lmax), \
                    nlam)
        else:
            self.monolam = np.linspace(lmin, lmax, nlam)
