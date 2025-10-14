"""
_____________________________________________________________________________________________________________
file name: Disk
last update: Dec 2021
language: PYTHON 3.8
short description:  Model of a static flared disk adapted for RT and chemistry simulations
_____________________________________________________________________________________________________________
"""
from __future__ import absolute_import
import os, sys, inspect
import numpy as np
#test
from .. constants.constants import mu, autocm, amu, Ggram, kb, M_sun

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import parameters as p
#___________________________________________
#   _______    __    _______    __   ___
#  |   __  \  |  |  |  _____|  |  | /  /
#  |  |  |  | |  |  |  |____   |  |/  /
#  |  |  |  | |  |  |____   |  |     | 
#  |  |__|  | |  |   ____|  |  |  |\  \
#  |_______/  |__|  |_______|  |__| \__\
#___________________________________________
class Disk:
    def __init__(self, ref_radius=p.ref_radius, rin=p.rin, rout=p.rout, star_mass=p.star_mass, disk_mass=p.disk_mass, h0=p.h0, \
                       sigma_gas_ref=p.sigma_gas_ref, Tmidplan_ref=p.Tmidplan_ref, Tatmos_ref=p.Tatmos_ref, sigma_t = p.sigma_t, q_exp=p.q_exp,  \
                       d_exp=p.d_exp, p_exp=p.p_exp, dtogas=p.dtogas, rho_m=p.rho_m, schmidtnumber=p.schmidtnumber, alpha=p.alpha, \
                       settfact=p.settfact, max_H=p.max_H, nz_chem=p.nz_chem, dust_mass=p.disk_dust_mass, q_c = p.q_c, dust=None, \
                       settling=True, isothermal=False, dust_density='g.cm-2', coordsystem='spherical'):                      
        self.ref_radius = ref_radius
        self.rin = rin
        self.rout = rout
        self.star_mass = star_mass
        self.disk_mass = disk_mass
        self.h0 = h0
        self.sigma_gas_ref = sigma_gas_ref
        self.Tmidplan_ref = Tmidplan_ref
        self.Tatmos_ref = Tatmos_ref
        self.sigma_t = sigma_t
        self.q_exp = q_exp
        self.d_exp = d_exp
        self.p_exp = p_exp
        self.dtogas = dtogas
        self.rho_m = rho_m
        self.schmidtnumber = schmidtnumber
        self.alpha = alpha
        self.settfact = settfact
        self.max_H = max_H
        self.nz_chem = nz_chem
        if (dust != None):
            self.dust = dust
        self.dust_mass = dust_mass
        self.settling = settling
        self.isothermal = isothermal
        self.dust_density = dust_density
        self.coordsystem = coordsystem
        self.q_c = q_c

    """
	The following methods give the physical parameters related to the gas phase in the disk.
	The gas is given by Hydrogren nuclei. Equation references are related to the documentation.
	Instance methods:
        A) scaleheight(self, r)
        B) surfacedensity(self, r)
        C) omega2(self, r)
        D) temp_mid(self, r)
        E) temp_atmos(self, r)
        F) temp_altitude(self, r, z)
        G) density_gauss(self, x1, x2, x3=None)
        H) density(self, x1, x2, x3=None)
        I) numberdensity(self, x1, x2, x3=None)
        J) viscous_accretion_heating(self)
    """

    def scaleheight(self, r):
        """ A)
        Returns gas scale height using the power law. Units: [cm].
        Parameters:
            -q_exp:         radial profile exponent. Will set the aspect ratio of the disk. If q_exp = 3, the disk is flat 
        	-r:             distances from the star - can be a list [au]
            -R_ref:         reference radius [au]
        """
        h_exp = (3./2.) - (self.q_exp/2.)
        if self.h0 != None: 
            hgas = self.h0*autocm*(r/(self.ref_radius))**h_exp
        elif self.Tmidplan_ref != None:
            h0 = np.sqrt((kb*self.Tmidplan_ref*(self.ref_radius*autocm)**3)/(mu*amu*Ggram*self.star_mass*M_sun))
            self.h0 = h0/autocm
            #h0 = h0/autocm
            hgas = h0*(r/(self.ref_radius))**h_exp

        return hgas


    def surfacedensity(self, r):
        """ B)
	    Surface density of the gas. Unit: g.cm-2
	    args:
	    	-r:                     distance from the star [au]
        """
        if self.sigma_gas_ref != None:
            sigma_g = 2*self.sigma_gas_ref*(r/(self.ref_radius))**(-self.p_exp)  #2 or not??
        else:
            sigma_g = (1/self.dtogas)*self.sigma_d0*(r/(self.ref_radius))**(-self.p_exp)
        return sigma_g

    def omega2(self, r):
        """ C)
        Returns square of the Keplerian angular velocity.
        Args:
            -r:        distance from the star [au]
        """
        return (Ggram*self.star_mass*M_sun)/(r*autocm)**3


    def temp_mid(self, r):
        """ D)
        Radial temperature profile in the midplane. Unit: Kelvin
        Args: 
	        -r:             distance from the star [au]
        """
        return self.Tmidplan_ref*(r/self.ref_radius)**(-self.q_exp)

    def temp_atmos(self, r):
        """ E)
	    Temperature profile in altitude. Unit: Kelvin   
	    """
        return self.Tatmos_ref*(r/self.ref_radius)**(-self.q_exp)

    def temp_altitude(self, r, z):
        """ F)
        Temperature vertical profile using the definition by Williams and Best (2014). Unit: Kelvin.
        """
        hg = self.scaleheight(r)
        tmid = self.temp_mid(r)
        tatm = self.temp_atmos(r)
        
        hhg, zz = np.meshgrid(hg, z,  indexing='ij')
        ttmid, zz = np.meshgrid(tmid, z, indexing='ij')
        ttatm, zz = np.meshgrid(tatm, z, indexing='ij')
        zz = hhg*zz
        zz0, z = np.meshgrid(zz[:, 0], z,  indexing='ij')
        return ttmid+(ttatm-ttmid)*np.sin((np.pi*zz)/(2*zz0))**(2*self.sigma_t)
	
    def verticaldensity_gauss(self, r):
        """ G)
        Gas number density. Unit: cm-3. ISOTHERMAL ONLY.
        WARNING: We divide here by mass to have the number density instead of mass density.
        """
        sigma = self.surfacedensity()
        sigma[(r >= self.rout) ^ (r <= self.rin)] = 0e0
        H = self.scaleheight(r)
        z = self.altitudes()	
        return (sigma/(mu*amu*H*autocm*np.sqrt(2.*np.pi)))\
               *np.exp(-(z)**2./(2.*(H)**2.))


    def density(self, x1, x2, x3=None):
        """ H)
        Return the gas density. Unit: g/cm^3. The density is computed assumging hydrostatic equilibrium.

        Notes
        -----
        The Gaussian vertical profile is an isothermal approximation. For a non-vertically isothermal 
        profile, the density can be computed iteratively. The profile slightly deviates from a Gaussian profile.
        This is done for cartesian coordinates only for now, will be done for spherical structure in a future update. 
        """
        if self.coordsystem == "spherical":
            rhog = np.ones((len(x1), len(x2), len(x3)))

            rt, tt, pp = np.meshgrid(x1*autocm, x2, x3, indexing='ij')
            rr = rt*np.sin(tt)
            zz = rt*np.cos(tt)
            zzmax = self.rout*autocm*np.cos(tt)
            border = np.greater_equal(abs(zzmax), abs(zz))*1 #used to spherically cut the disk in case the grid is larger than the disk outer radius e.g. if an envelope is present

            sigma_g = self.surfacedensity(rr/autocm)
            hg = self.scaleheight(rr/autocm)

            rhog = (sigma_g/(np.sqrt(2*np.pi)*hg))*np.exp(-(zz[:,:,:]**2)/(2*hg**2))*border


        if self.coordsystem == "nautilus":
            omega2 = self.omega2(x1)
            sigma_g = self.surfacedensity(x1)
            hg = self.scaleheight(x1)
            temp = self.temp_altitude(x1, x2)
            midplane_dens = (sigma_g)/(hg*np.sqrt(2*np.pi))
            rhog = np.ones((len(x1), len(x2)))
            rhog[:, 0] = 0 # init

            for r in range(0, len(x1), 1):
                for z in range(1, len(x2), 1):
                    rhog[r, z] = rhog[r, z-1] - (np.log(temp[r,z]) - np.log(temp[r,z-1])) - ((mu*amu*omega2[r]*hg[r]*x2[z]*(hg[r]*x2[z] - hg[r]*x2[z-1]))/(kb*temp[r,z]))
                rhog[r] = np.exp(rhog[r])
                maximum=np.amax(rhog[r])
                rhog[r] = (rhog[r]*midplane_dens[r])/maximum
        return rhog
        
    def numberdensity(self, x1, x2, x3=None):
        """ I)
        Return the gas number density. Unit: cm^-3. The density is computed assumging hydrostatic equilibrium.

        Notes
        -----
        If vertically isothermal, the profile is Gaussian. If not, the density can be computed iteratively and the profile is not Gaussian.
        """
        ng = self.density(x1, x2, x3)/(mu*amu)            
        return ng

    def viscous_accretion_heating(self, acc_rate, lim_h, x1, x2, x3=None):
        """ I)
        Return the viscous accretion heating. Unit: erg.cm^-3.s-1. The density is computed assumging hydrostatic equilibrium.

        Notes
        -----
        If vertically isothermal, the profile is Gaussian. If not, the density can be computed iteratively and the profile is not Gaussian.
        """
        if self.coordsystem == "spherical":
            yrtosec = 31536000 # number of second in a year. 
            q_visc = np.ones((len(x1), len(x2), len(x3)))
            rt, tt, pp = np.meshgrid(x1*autocm, x2, x3, indexing='ij')
            rr = rt*np.sin(tt)
            zz = rt*np.cos(tt)
            zzmax = self.rout*autocm*np.cos(tt)
            sigma_g = self.surfacedensity(rr/autocm)
            hg = self.scaleheight(rr/autocm)
            omega2 = self.omega2(rr/autocm)
            border_max = np.greater_equal(abs(zzmax), abs(zz))*1  # make border at rout
            border_h = np.greater_equal(abs(lim_h*hg), abs(zz))*1 #make border at lim_h scale height. 
            border = border_max*border_h
            q_visc = ((3*acc_rate*M_sun*omega2)/(4*np.pi*np.sqrt(2*np.pi)*hg*yrtosec))*np.exp(-(zz[:,:,:]**2)/(2*hg**2))*border
            #q_visc[q_visc==0.0] = 1e-30         
        return q_visc

        

    """
    The following methods give the physical quantities related to the grains in the disk.
    Equation references are those of the documentation.
    Instance methods:
	    A) surfacedensity_d(self, r)
        B) scaleheight_d(self, r)
        C) density_d(self, z)
        D) numberdensity_d(self, r)
        E) Qext(self)
        F) av_z(self, lam, nd, r, z)
    """
    def surfacedensity_d(self, r):
        """ A)
        Return dust surface density (g.cm-2)
        """
        sigmad = []
        fraction = self.dust.massfraction()
        #print(r.shape)
        #sig = np.ones((fraction.shape, r[0].shape, r[1].shape, r[2].shape))
        #print(sig)
        if self.sigma_gas_ref != None:
            sigma_di0 = fraction*self.dtogas*self.sigma_gas_ref
            sigma_single0 = self.dtogas*self.sigma_gas_ref
        else:
            self.sigma_d0 = (2-self.p_exp)*self.dust_mass*M_sun/(2*np.pi*(self.ref_radius)**(self.p_exp)) / \
                       (self.rout**(-self.p_exp+2) - self.rin**(-self.p_exp+2)) /autocm**2
            sigma_di0 = fraction*self.sigma_d0
            sigma_single0 = self.sigma_d0
 
        for s in sigma_di0:
            sig = s*(r/self.ref_radius)**(-self.p_exp)
            sig[(r >= self.rout) ^ (r <= self.rin)] = 0e0 # In case sigma is outside disk outer boundaries
            sigmad.append(sig)
        
        sig_single = sigma_single0*(r/self.ref_radius)**(-self.p_exp)
        sig_single[(r >= self.rout) ^ (r <= self.rin)] = 0e0 # In case sigma is outside disk outer boundaries

        return np.array(sigmad), sig_single

    def scaleheight_d(self, r):
        """ B)
	    Return Dust scale height H_d(r, a) dependent on the grain sizes and radii. 

        Notes:
        -----
        2D array (len(r), len(nb_sizes)). Units: [au]
	    """	
        hd = []
        hg = self.scaleheight(r)
        sizes = self.dust.sizes() #grain size in microns
        rsingle = self.dust.rsingle

        if self.settling == True:
            sigma_g = self.surfacedensity(r)
            for a in sizes[-1]:
                stoptime_mid =(np.pi*a*1e-4*self.rho_m)/(2*sigma_g)
                hd.append(hg/(np.sqrt(1 + stoptime_mid*(self.schmidtnumber/self.alpha))))
            stoptime_mid_single =(np.pi*self.dust.rsingle*1e-4*self.rho_m)/(2*sigma_g)
            hd_single = hg/(np.sqrt(1 + stoptime_mid_single*(self.schmidtnumber/self.alpha)))
            return np.array(hd), hd_single

        if self.settling == False:
            for a in sizes[-1]:
                hd.append(hg)
            return np.array(hd), hg
            
    def density_d(self, x1, x2, x3=None):
        """ C)
	    Return dust density rho_d(r, z, a) or rho_d(r, theta, phi, a). 

        Notes:
        -----
        3D array (len(r), len(nb_sizes), len(z)). Units: [g.cm-3]
	    """	
        rhod = []
        rhod_single = []
        if self.coordsystem =='spherical':
            rt, tt, pp = np.meshgrid(x1*autocm, x2, x3, indexing='ij')
            rr = rt*np.sin(tt)
            zz = rt*np.cos(tt)
            zzmax = self.rout*autocm*np.cos(tt)
            sigmad, sigmad_single = self.surfacedensity_d(rr/autocm)
            hd, hd_single = self.scaleheight_d(rr/autocm)
            rhod = np.ones( (len(hd), len(x1), len(x2), len(x3)))
            rhod_single = np.ones((len(x1), len(x2), len(x3)))
            border = np.greater_equal(abs(zzmax), abs(zz))*1 #used to spherically cut the disk in case the grid is larger than the disk outer radius e.g. if an envelope is present. *1 to convert boolean to 0,1 values.
            for i in range(len(hd)): #len(hd) is equal to the number of grain sizes
                rhod[i, :, :, :] = (sigmad[i,:,:,:]/(np.sqrt(2*np.pi)*hd[i,:,:,:]))*np.exp(-(zz[:,:,:]**2)/(2*hd[i,:,:,:]**2))*border


        if self.coordsystem =='nautilus':
            rr, zz = np.meshgrid(x1*autocm, x2, indexing='ij')
            sigmad, sigmad_single = self.surfacedensity_d(rr/autocm)
            hg = self.scaleheight(rr/autocm)
            hd, hd_single = self.scaleheight_d(rr/autocm)
            rhod = np.ones((len(hd), len(x1), len(x2)))
            rhod_single = np.ones((len(x1), len(x2)))
            zz = hg*zz
            for i in range(len(hd)): #len(hd) is equal to the number of grain sizes
                rhod[i, :, :] = (sigmad[i]/(np.sqrt(2*np.pi)*hd[i]))*np.exp(-((zz)**2)/(2*hd[i]**2))

            rhod_single[:, :] = (sigmad_single/(np.sqrt(2*np.pi)*hd_single))*np.exp(-((zz)**2)/(2*hd_single**2))
            self.rhod_single = rhod_single

        if self.dust_density == 'g.cm-2':
            return rhod

    def numberdensity_d(self, x1, x2, x3=None):
        """ D)
	    Return dust number density n_d(r, z, a). Depends on the grain sizes, radii, and altitude. 

        Notes:
        -----
        3D array (len(nb_sizes), len(r), len(z)). Units: [cm-3]
        Example:
        ------- 
            - call n_d[1, 2, :] for densities at third radius and for second grain species.
	    """	
        mass = self.dust.grainmass()
        dens = self.density_d(x1, x2, x3)

        for i in range(len(mass)):
            dens[i] = dens[i]/mass[i]
        
        return dens

    def numberdensity_d_single(self, x1, x2, x3=None):
        """ D)
	    Return dust number density n_d(r, z) in case of single grain for chemistry. Usefull if radmc3d uses multiple sizes and nautilus uses a single one. 

        Notes:
        -----
        2D array (len(r), len(z)). Units: [cm-3]. If nb_species = 1, then numberdensity_d_single[r, z] = numberdensity_d[0, r, z]. We use two seperate functions in case nb_species > 1 and the user needs one size for chemistry.


	    """	
        mass_single = self.dust.grainmass_single()
        dens_single = self.rhod_single

        dens_single = dens_single/mass_single

        return dens_single


    def q_ext(self, lam, A=2, q_c=4):
        """ E)
	    Return extinction efficiency. 

        Notes:
        -----
        3D array (len(nb_sizes), len(r), len(z)). Units: [cm-3]
	    """	
        sizes = self.dust.sizes()
        lambda_c = 2*np.pi*sizes[-1]

        qext = np.ones(( len(sizes[-1]), len(lam)  ))
        for idx_a, a in enumerate(sizes[-1]):
            for idx_wl, wl in enumerate(lam):
                if wl <= np.pi*a:
                    qext[idx_a, idx_wl] = A  #regime A
                elif np.pi*a < wl < 2*np.pi*a:
                    qext[idx_a, idx_wl] = q_c*(wl/lambda_c[idx_a])**(np.log10(q_c)/np.log10(2)) #regime B
                elif wl >= 2*np.pi*a:
                    qext[idx_a, idx_wl] = q_c*(wl/lambda_c[idx_a])**(-2) #regime C
        return qext

    def av_z(self, lam, nd, r, z):
        """ F)
	    Return visual extinction. 
	    """	       
        sizes = self.dust.sizes()
        #extinction efficiency
        lambda_c = 2*np.pi*sizes[-1]
        qext = np.ones(( len(sizes[-1]), len(lam)  ))
        for idx_a, a in enumerate(sizes[-1]):
            for idx_wl, wl in enumerate(lam):
                if wl <= np.pi*a:
                    qext[idx_a, idx_wl] = 1  #regime A
                elif np.pi*a < wl < 2*np.pi*a:
                    qext[idx_a, idx_wl] = self.q_c*(wl/lambda_c[idx_a])**(np.log10(self.q_c)/np.log10(2)) #regime B
                elif wl >= 2*np.pi*a:
                    qext[idx_a, idx_wl] = self.q_c*(wl/lambda_c[idx_a])**(-2) #regime C

        rr, zz = np.meshgrid(r, z, indexing='ij')
        hg = self.scaleheight(rr)    
        zz = hg*zz
        avz = np.zeros( (len(r), len(z)) )
        for idx_r in range(0, len(r), 1):
            for idx_z in range(1, len(z), 1):
                for idx_a in range(0, len(sizes[-1]), 1):
                    avz[idx_r, idx_z] += 1.086*np.pi*nd[idx_a, idx_r, idx_z]*qext[idx_a, 12]*(zz[idx_r, idx_z-1] - zz[idx_r, idx_z])*(sizes[-1][idx_a]*1e-4)**2
            avz[idx_r, :] = np.cumsum(avz[idx_r, :])
            avz[:, 0] = avz[:, 1]
        return avz


    def avz_ism(self, lam, A=2, q_c=4):
        """ E)
	    Return visual extinction using the conversion factor of H column density to Av:. 

        Notes:
        -----
        3D array (len(nb_sizes), len(r), len(z)). Units: [cm-3]
	    """	

    # def av_z2(self ,mass, filelist, ):
    #     #--Av
    #     mass = d.grainmass() #g
    #     kext = np.zeros(len(mass))
    #     filelist = sorted(glob.glob('thermal/dustkap*'))
    #     for ai in range(0, len(mass)):
    #         kappa = pd.read_table(filelist[ai], sep="\s+", comment='#', header=None, skiprows=10)
    #         kext[ai] += (kappa[1].iloc[0] + kappa[2].iloc[0]) # k_extinction at wc ~ 550 nm.

    #     av_mid = np.ones(len(radii))
    #     for ri in range(len(radii)):
    #         z = m.grid.zchem*m.grid.hg_chem[0][ri]
    #         av_z = np.zeros(len(m.grid.zchem))
    #         av_z[0] = 0
    #         for zi in range(1, len(m.grid.zchem)):
    #             for ai in range(0, len(mass)):
    #                 av_z[zi] += 1.086*m.grid.dustdensity_chem[0][ai, ri, zi]*mass[ai]*kext[ai]*(z[zi-1] - z[zi])
    #         av_z = np.cumsum(av_z)
    #         av_mid[ri] = av_z[-1]
    #     print(av_mid)
    #     #--
