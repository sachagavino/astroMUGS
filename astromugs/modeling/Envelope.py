#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_____________________________________________________________________________________________________________
file name: Envelope
last update: may 2021
language: PYTHON 3.8
short description:  geometrical model of a class 0/I envelope. 
                    Every input and output in au, except output for densities returned in g.cm-3. 
_____________________________________________________________________________________________________________
"""


import numpy as np
from scipy.optimize import brenth
from scipy.integrate import trapezoid

from .. constants.constants import mu, autocm, amu, Ggram, kb, M_sun

import parameters as p



#___________________________________________
#   _________   ___    __  ___        ___
#  |   ______| |   \  |  | \  \      /  /
#  |  |______  |    \ |  |  \  \    /  /
#  |   ______| |  |\ \|  |   \  \  /  /
#  |  |______  |  | \    |    \  \/  /
#  |_________| |__|  \___|     \____/
#___________________________________________
class Envelope:
    def __init__(self, rmin=p.rmin, rmax=p.rmax, r_centri=p.r_centri, \
                       acc_rate=p.acc_rate, star_mass=p.star_mass, dust_mass=p.dust_env_mass, dtogas=p.dtogas,  \
                       cavpl=p.cavpl, cav_fact=p.cav_fact, cavz0=p.cavz0, dust=None, dust_density='g.cm-2', coordsystem='spherical'):
        self.rmin=rmin
        self.rmax=rmax
        self.r_centri=r_centri*autocm
        self.acc_rate=acc_rate
        self.star_mass=star_mass
        self.dust_mass=dust_mass
        self.dtogas = dtogas
        self.cavpl=cavpl
        self.cav_fact=cav_fact
        self.cavz0=cavz0
        self.dust_density = dust_density
        self.coordsystem = coordsystem
        if (dust != None):
            self.dust = dust
    """
	The following methods give the physical parameters related to the gas phase and dust in the envelope.
	The gas number density is given relative to Hydrogren nuclei. Equation references are related to the documentation.
	Instance methods:
        A) 
        B) 
        C) 
        D) 

    """

    def density(self, x1, x2, x3=None):
        """ A)
	    Return dust density rho_d(r, z, a) or rho_d(r, theta, phi, a). 

        Notes:
        -----
        3D array (len(r), len(nb_sizes), len(z)). Units: [g.cm-3]
	    """	


        #######
        if self.coordsystem =='spherical':
            rt, tt, pp = np.meshgrid(x1*autocm, x2, x3, indexing='ij')

            mu = np.cos(tt)

            RR = rt * np.sin(tt)
            zz = rt * np.cos(tt)

            mu0 = mu*0
            for ir in range(rt.shape[0]):
                for it in range(rt.shape[1]):
                    mu0[ir,it,0] = brenth(self.solution, -1.0, 1.0, args=(rt[ir,it,0], mu[ir,it,0]))

            rho0 = 1.0
            rho = rho0 * (rt / self.r_centri)**(-1.5) * (1 + mu/mu0)**(-0.5)* (mu/mu0 + 2*mu0**2 * self.r_centri/rt)**(-1)

            mid1 = (np.abs(mu) < 1.0e-10) & (rt < self.r_centri)
            rho[mid1] = rho0 * (rt[mid1] / self.r_centri)**(-0.5) * (1. - rt[mid1] / self.r_centri)**(-1) / 2.

            mid2 = (np.abs(mu) < 1.0e-10) & (rt > self.r_centri)
            rho[mid2] = rho0 * (2.*rt[mid2]/self.r_centri - 1)**(-0.5) * (rt[mid2]/self.r_centri - 1.)**(-1)

            rho[(rt > self.rmax*autocm) ^ (rt < self.rmin*autocm)] = 0e0
                
            if x2.max() > np.pi/2:
                mdot = ((self.dust_mass/self.dtogas)*M_sun)/(2*np.pi*trapz(trapz(rho*rt**2*np.sin(tt),tt,axis=1), \
                        rt[:,0,:],axis=0))[0]
            else:
                mdot = ((self.dust_mass/self.dtogas)*M_sun)/(4*np.pi*trapz(trapz(rho*rt**2*np.sin(tt),tt,axis=1), \
                        rt[:,0,:],axis=0))[0]
            rho *= mdot

            #--- Outflow cavity.
            rho[np.abs(zz)/autocm-self.cavz0/autocm-(RR/autocm)**self.cavpl > 0.0] *= self.cav_fact

            return rho

    def numberdensity(self, x1, x2, x3=None):
        """ B)
        Return the gas number density. Unit: cm^-3.

        Notes
        -----
        If vertically isothermal, the profile is Gaussian. If not, the density can be computed iteratively and the profile is not Gaussian.
        """
        ng = self.density(x1, x2, x3)/(mu*amu)            
        return ng

    def density_d(self, x1, x2, x3=None):
        """ C)
	    Return dust density rho_d(r, z, a) or rho_d(r, theta, phi, a). 

        Notes:
        -----
        3D array (len(r), len(nb_sizes), len(z)). Units: [g.cm-3]
	    """	
        rhog = self.density(x1, x2, x3)
        fraction = self.dust.massfraction()
        rhod = np.ones( (len(fraction), len(x1), len(x2), len(x3)))

        for i in range(len(fraction)): #len(hd) is equal to the number of grain sizes
            rhod[i, :, :, :] = fraction[i]*self.dtogas*rhog[:, :, :]

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
            - call n_d[1][2] for densities at third radius and for second grain species.
	    """	
        mass = self.dust.grainmass()
        dens = self.density_d(x1, x2, x3)
        for i in range(len(mass)):
            dens[i] = dens[i]/mass[i]
            
        return dens

    def solution(self, mu0,r,mu):
        """ E)
	    Return solution

        Notes:
        -----
	    """	
        return mu0**3-mu0*(1-r/self.r_centri)-mu*(r/self.r_centri)






#---old envelope:



#_______________________________________________#
#                     CLASS                     #
#_______________________________________________#

#_____________Outer class: ENVELOPE_____________
#   _________   ___    __  ___        ___
#  |   ______| |   \  |  | \  \      /  /
#  |  |______  |    \ |  |  \  \    /  /
#  |   ______| |  |\ \|  |   \  \  /  /
#  |  |______  |  | \    |    \  \/  /
#  |_________| |__|  \___|     \____/
#_______________________________________________
# class Envelope:
#     def __init__(self, dist, mu, R_centri=p.R_centri, acc_rate=p.acc_rate, star_mass=p.star_mass):
#         self.R_centri=R_centri
#         self.acc_rate=acc_rate
#         self.star_mass=star_mass



#     """
# 	The following methods give the physical parameters related to the gas phase in the envelope.
# 	The gas is given by Hydrogren nuclei. Equation references are related to the documentation.
# 	Instance methods:
#         A) solution (parameter for the analytical solution)
#         B) mu_0()
#         C) env_density(self, acc_rate, star_mass, dist, R_centri, mu, mu_0)
#     """

#     """ A)
#     Equation (1). Initializer. Returns the parameter for the solution of the initial streamline of the infalling material 
#     Arguments:
#         -mu:            cos(theta) where theta is the angle above the midplane 
#         -dist:          distance to the central core in the spherical coordinates (dist, theta, phi) [au]
#         -R_centri:      centrifugal radius [au]
#     """
#     def solution(self, dist, mu, R_centri):
#         sol = (27*mu*dist*R_centri**2+np.sqrt(729*mu**2*dist**2*R_centri**4+108*R_centri**3*(dist-R_centri)**3))**(1/3)/(3*R_centri*2**(1/3))
#         #solution = 729*mu**2*dist**2*R_centri**4 + 108*R_centri**3*(dist-R_centri)**3
#         return sol

#     """ B)
#     Equation (2). Returns the solution of the initial streamline of the infalling material (cos(theta_0))
#     Arguments:
#         -dist:          distance to the central core in the spherical coordinates (dist, theta, phi) [au]
#         -R_centri:      centrifugal radius [au]
#     """	
#     def mu_0(self, dist, R_centri):
#         solution = self.solution()
#         return solution - ((1)/(3*solution))

#     """ C)
#     Equation (2). Returns the density of gas of the envelope [g.cm-3].
#     Arguments:
#         -acc_rate:      accretion rate in solar mass per year (converted in g per second)
#         -dist:          distance to the central core in the spherical coordinates (dist, theta, phi) [au]
#         -R_centri:      centrifugal radius [au]
#         -mu:            cos(theta) where theta is the angle above the midplane
#         -mu0:           initial streamline of the infalling material
#     """	
#     def env_density(self, acc_rate, star_mass, dist, R_centri, mu, mu_0):
#         return (((acc_rate*M_sun)/3.154e7)/(4*np.pi))*(Ggram*star_mass*M_sun*(dist*autocm)**3)**(-1/2)*(1 + mu/mu_0)**(-1/2)*(mu/mu_0 + 2*mu_0**2*(R_centri/dist))**(-1)

#     # #class Dust:
#     #     """
#     #     This Class gives the physical quantities related to the grains in the envelope.
#     #     Equation references are those of the documentation.
#     #     Instance methods:
# 	#         A) initializer ...
#     #     """
        
