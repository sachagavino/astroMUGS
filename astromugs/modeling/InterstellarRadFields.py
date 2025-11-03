"""
File name: isrf.py
Author: Sacha Gavino
Date: Nov 2023
Language: Python 3.10

Short description:
    Create ISRF as a function of wavelength.
"""

import numpy as np


class InterstellarRadFields:
    def __init__(self, cut=2e-1, d78=True, vdb82=True):
        # cut in microns
        self.cut = cut
        self.d78 = d78
        self.vdb82 = vdb82

    def draine78(self, lam):
        """
        Draine (1978) fit from Sternberg & Dalgarno (1995) between 912 and 2000 Å.

        Parameters
        ----------
        lam : ndarray
            Wavelength in microns.

        Returns
        -------
        ndarray
            2D array: [wavelength (μm), ISRF (erg cm^-2 s^-1 Hz^-1 sr^-1)]
        """
        flux = (
            1.068e-4 * lam ** -1
            - 1.719e-2 * lam ** -2
            + 6.853e-1 * lam ** -3
        ) / (4 * np.pi)

        flux *= 3e-12  # convert photons to ergs

        return np.stack([lam, flux])

    def van_black82(self, lam):
        """
        Extension of van Dishoeck & Black (1982) at wavelengths > 2000 Å.

        Parameters
        ----------
        lam : ndarray
            Wavelength in microns.

        Returns
        -------
        ndarray
            2D array: [wavelength (μm), ISRF (erg cm^-2 s^-1 Hz^-1 sr^-1)]
        """
        flux = 3.67e4 * lam ** 0.7 / (4 * np.pi)
        flux *= 2.998e-12 / 2.998e17  # photon → erg conversion

        return np.stack([lam, flux])

    def create_isrf(self, lam):
        """
        Build the combined ISRF from the selected components.

        Parameters
        ----------
        lam : ndarray
            Wavelength in microns.

        Returns
        -------
        ndarray
            Combined ISRF array.
        """
        lam_mm = lam * 1e3  # convert to nm? (kept original logic)

        if self.d78 and self.vdb82:
            lam_d78 = lam_mm[lam_mm <= self.cut * 1e3]
            lam_vdb82 = lam_mm[lam_mm > self.cut * 1e3]

            isrf_d78 = self.draine78(lam_d78)
            isrf_vdb82 = self.van_black82(lam_vdb82)

            return np.concatenate((isrf_d78, isrf_vdb82), axis=1)

        if self.d78:
            return self.draine78(lam_mm)

        if self.vdb82:
            return self.van_black82(lam_mm)

        return None
