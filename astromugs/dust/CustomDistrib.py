"""
_____________________________________________________________________________________________________________
file name: CustomDistrib
@author: Sacha Gavino
last update: June 2022
language: PYTHON 3.8
short description:  Creates a dust population.
_____________________________________________________________________________________________________________
"""
import numpy as np
from .. constants.constants import mu, amu



class CustomDistrib:

    def __init__(self, rho_m=2.5, units='microns', path='', filename=None):
        """
        The argument 'units' represents the size units given in the file. It can be microns, mm, cm, or m.  
        -------
        """
        self.filename = filename
        self.units = units
        self.path = path
        self.rho_m = rho_m

    def sizes(self):
        """ A)
        Create an array with grain sizes from a given file.
        Returns the list of sizes. Units: microns
        -------
            Numpy array (sizes). len(array) = (nb_sizes).
        """
        if (self.filename == None):
            filename = "dust_sizes.in"
        else:
            filename = self.filename

        sizes = []

        f = open(self.path + filename,"r")


        sizes.append(f.readline().replace("\n",""))

        while len(sizes[-1]) > 0:
            sizes.append(f.readline().replace("\n",""))
        sizes.pop() # remove last empty element

        f.close()

        sizes = np.array([sizes], dtype=float) # we add extra brackets because the array is meant to have more than one dimension, ultimately, with the actual sizes as the last. 

        #convert in microns
        if self.units == 'microns':
            sizes = sizes
        if self.units == 'mm':
            sizes = sizes*1e3
        if self.units == 'cm':
            sizes = sizes*1e4
        if self.units == 'm':
            sizes = sizes*1e6

        return sizes


    def grainmass(self):
        """ B)
        Create an array with masses of a single grain for each grain population. 
        Returns
        -------
            Numpy array. len(array) = (nb_sizes). Units: gram
        """
        a = self.sizes()
        mass = ((4.*np.pi)/3.)*self.rho_m*(a[-1]*1e-4)**3

        return mass