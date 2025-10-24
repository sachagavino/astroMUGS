"""
_____________________________________________________________________________________________________________
file name: MRNDistrib
@author: Sacha Gavino
last update: May 2022
language: PYTHON 3.8
short description:  Creates a dust population from MRN.
_____________________________________________________________________________________________________________
"""
import numpy as np
from .. constants.constants import mu, amu



class MRNDistrib:
    
        def __init__(self, rsingle=1.000e-01 , amin=5.000e-03, amax=1.000e+03, \
                           nb_sizes=1, d_exp=3.5, rho_m=2.5, dtogas=0.01, cst_norm=7.41e-26, ext_eff=4):
            self.d_exp = d_exp
            self.rho_m = rho_m
            self.dtogas = dtogas
            self.ext_eff = ext_eff
            self.amin = amin
            self.amax = amax
            self.cst_norm = cst_norm
            self.nb_sizes = nb_sizes
            self.rsingle = rsingle #used for chemistry

        def sizes(self):
            """ A)
            Create an array with min, max, and average value of each interval of grain sizes. 
            Returns
            -------
                Numpy array (amin, amax, average). len(array) = (3, nb_sizes). Units: microns
            """
            if self.nb_sizes == 1:
                sizes_param =  np.array([[self.rsingle]]) #
            if self.nb_sizes > 1:
                a = np.logspace(np.log10(self.amin), np.log10(self.amax), self.nb_sizes+1)
                average = (((1-self.d_exp)/(4-self.d_exp))*((a[1:]**(4-self.d_exp) - a[:-1]**(4-self.d_exp))/(\
	                         a[1:]**(1-self.d_exp) - a[:-1]**(1-self.d_exp))))**(1./3.) 

                sizes_param =  np.array([a[:-1], a[1:], average])  # a[:-1] all but last; a[1:] all but first. Provides the averaged, min and max values of all ranges.

            return sizes_param


        def grainmass(self):
            """ B)
            Create an array with masses of a single grain for each grain population. 
            Returns
            -------
                Numpy array. len(array) = (nb_sizes). Units: gram
            """
            a = self.sizes()
            if self.nb_sizes == 1:
                mass = ((4.*np.pi)/3.)*self.rho_m*(a[-1]*1e-4)**3
                
            if self.nb_sizes > 1:
                mass = ((4.*np.pi)/3.)*self.rho_m*(a[-1]*1e-4)**3

            return mass

        def grainmass_single(self):
            """ B)
            Create mass from single grain (rsingle). 
            Returns
            -------
                Numpy float. Units: gram
            """
            mass_single = ((4.*np.pi)/3.)*self.rho_m*(self.rsingle*1e-4)**3

            return mass_single

        def massfraction(self):
            """ C)
            Calculate the mass fraction of each grain population relative to the total dust mass of the object following a MRN distribution. 
            Returns
            -------
                Numpy array. len(array) = (nb_sizes).
            """
            a = self.sizes()

            if self.nb_sizes == 1:
                fraction = np.array([1])

            if self.nb_sizes > 1:
                cst_norm = (3./(4.*np.pi))*((self.dtogas*mu*amu)/self.rho_m)*(4 - self.d_exp)*\
                           (1/(self.amax**(4-self.d_exp) - self.amin**(4-self.d_exp)))

                mass_density = ((4.*np.pi)/3.*(4.-self.d_exp))*self.rho_m*cst_norm*(a[1]**(4.-self.d_exp) - \
	                             a[0]**(4.-self.d_exp))

                total_density = np.sum(mass_density)
                fraction = mass_density/total_density

            return fraction
