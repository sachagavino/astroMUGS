"""
file name: isrf.
author: Sacha Gavino
date: Nov 2021
language: PYTHON 3.8
_____________________________________________________________________________________________________
short description:  create ISRF as a function of wavelengths.
_____________________________________________________________________________________________________
"""
import numpy as np

class InterstellarRadFields:
    def __init__(self, cut=2.e-1, d78=True, vdb82=True):
        self.cut=cut #micron
        self.d78=d78 #micron
        self.vdb82=vdb82 #micron


    def draine78(self, lam):
        '''
        Desc: Draine (1978) fit from Sterberg & Dalgarno (1995) between 912 and 2000 Angstrom
        Args: lam (microns)
        return: 2D array (wavelength [nm], ISRF [erg.cm-2.s-1.hz-1.sr-1])
        '''
        #lam_d78 = self.lam[self.lam <= 2.e-1]*1e4 #take the right range and convert to angstrom.

        #---- ORIGINAL FLUX ----
        #d78 = (3.2028e13*lam**(-3) - 5.1542e15*lam**(-4) + 2.0546e17*lam**(-5) ) / (4*np.pi) #divide by 4pi to get sr-1
        #d78 *= (2.998e-12)/(2.998e17) # 1photon = 3e-12 ergs and 1nm = 2.998e17 Hz
        #d78 = np.stack([lam, d78])
        #-----------------------


        #---- UPDATED ---- Sterberg & Dalgarno (1995) fit in photons cm-2 s-1 Hz-1 sr-1
        d78 = (1.068e-4*lam**(-1) - 1.719e-2*lam**(-2) + 6.853e-1*lam**(-3))/ (4*np.pi)
        d78 *= 3e-12 #1photon = 3e-12 ergs
        d78 = np.stack([lam, d78])
        #-----------------
        return d78
        

    def van_black82(self, lam):
        '''
        Desc: Extension of van Dishoek & Black (1982) at longer wavelengths (> 2000 Angstrom)
        Args: lam (microns)
        '''
        #lam_vdb82 = self.lam[self.lam > 2.e-1]*1e3
        vdb82 = 3.67e4*lam**(0.7) / (4*np.pi)
        vdb82 *= (2.998e-12)/(2.998e17)
        #vdb82 = 1.38243e-6 * lam**(-0.3)
        ##vdb82 *= (2.998e18)
        vdb82 = np.stack([lam, vdb82])
        return vdb82

    def create_isrf(self, lam):
        lam = lam*1e3 #convert to mm
        if self.d78 is True and self.vdb82 is True:
            lam_d78 = lam[lam <= self.cut*1e3] #take the right range and convert to angstrom.
            lam_vdb82 = lam[lam > self.cut*1e3]
            isrf = np.concatenate((self.draine78(lam_d78), self.van_black82(lam_vdb82)), axis=1)

        elif self.d78 is True and self.vdb82 is False:
            isrf = self.draine78(lam) 

        elif self.d78 is False and self.vdb82 is True:
            isrf = self.van_black82(lam)

        return isrf



