#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
File name: Envelope
Last update: May 2021
Language: Python 3.8

Short description:
    Geometrical model of a Class 0/I envelope.
    Every input and output is in AU, except densities which are
    returned in g cm^-3.
"""

import numpy as np
from scipy.optimize import brenth
from scipy.integrate import trapezoid

from dataclasses import dataclass, fields, is_dataclass
from typing import Optional, Literal, Any

from astromugs.utils.params import EnvelopeParams
from astromugs.constants.constants import (
    mu,
    autocm,
    amu,
    Ggram,
    kb,
    M_sun,
)


# ------------------------------------------------------------------------------
#   _________   ___    __  ___        ___
#  |   ______| |   \  |  | \  \      /  /
#  |  |______  |    \ |  |  \  \    /  /
#  |   ______| |  |\ \|  |   \  \  /  /
#  |  |______  |  | \    |    \  \/  /
#  |_________| |__|  \___|     \____/
# ------------------------------------------------------------------------------
class Envelope:
    def __init__(
        self,
        params: EnvelopeParams,
        dust: Optional[object] = None
    ):
        self.params = params
        self.dust = dust


       # # copy dataclass attributes to instance variables (optional)
        #for field in params.__dataclass_fields__:
            #setattr(self, field, getattr(params, field))

    """
    The following methods give the physical properties related to the gas
    and dust in the envelope.
    The gas number density is relative to hydrogen nuclei.
    """
    def density(self, x1: np.ndarray, x2: np.ndarray, x3: Optional[np.ndarray] = None) -> np.ndarray:
        """
        A)
        Return dust density rho_d(r, z, a) or rho_d(r, theta, phi, a).

        Returns
        -------
        rho : np.ndarray
            3D array with shape (len(r), len(nb_sizes), len(z)), units: g cm^-3
        """
        if self.params.coordsystem == "spherical":
            rt, tt, pp = np.meshgrid(x1 * autocm, x2, x3, indexing="ij")
            mu_val = np.cos(tt)
            RR = rt * np.sin(tt)
            zz = rt * np.cos(tt)
            mu0 = np.zeros_like(mu_val)

            # Find mu0 via root finding
            for ir in range(rt.shape[0]):
                for it in range(rt.shape[1]):
                    mu0[ir, it, 0] = brenth(
                        self.solution,
                        -1.0,
                        1.0,
                        args=(rt[ir, it, 0], mu_val[ir, it, 0]),
                    )

            rho0 = 1.0
            rho = (
                rho0
                * (rt / self.params.r_centri) ** -1.5
                * (1 + mu_val / mu0) ** -0.5
                * (mu_val / mu0 + 2 * mu0**2 * self.params.r_centri / rt) ** -1
            )

            mid1 = (np.abs(mu_val) < 1e-10) & (rt < self.params.r_centri)
            rho[mid1] = (
                rho0
                * (rt[mid1] / self.params.r_centri) ** -0.5
                * (1 - rt[mid1] / self.params.r_centri) ** -1
                / 2
            )

            mid2 = (np.abs(mu_val) < 1e-10) & (rt > self.params.r_centri)
            rho[mid2] = (
                rho0
                * (2 * rt[mid2] / self.params.r_centri - 1) ** -0.5
                * (rt[mid2] / self.params.r_centri - 1) ** -1
            )

            rho[(rt > self.params.rmax * autocm) ^ (rt < self.params.rmin * autocm)] = 0.0

            # Mass normalization
            if x2.max() > np.pi / 2:
                mdot = (
                    (self.params.dust_mass / self.params.dtogas) * M_sun
                ) / (
                    2
                    * np.pi
                    * trapezoid(
                        trapezoid(rho * rt**2 * np.sin(tt), tt, axis=1),
                        rt[:, 0, :],
                        axis=0,
                    )
                )[0]
            else:
                mdot = (
                    (self.params.dust_mass / self.params.dtogas) * M_sun
                ) / (
                    4
                    * np.pi
                    * trapezoid(
                        trapezoid(rho * rt**2 * np.sin(tt), tt, axis=1),
                        rt[:, 0, :],
                        axis=0,
                    )
                )[0]

            rho *= mdot

            # Outflow cavity
            cavity_mask = (
                np.abs(zz) / autocm
                - self.params.cavz0 / autocm
                - (RR / autocm) ** self.params.cavpl
                > 0.0
            )
            rho[cavity_mask] *= self.params.cav_fact

            return rho

    def numberdensity(self, x1, x2, x3=None):
        """
        B)
        Return gas number density in cm^-3.
        """
        ng = self.density(x1, x2, x3) / (mu * amu)
        return ng

    def density_d(self, x1, x2, x3=None):
        """
        C)
        Return dust density rho_d(r, z, a).

        Notes
        -----
        Returns 3D array (len(nb_sizes), len(r), len(z)).
        Units: g cm^-3
        """
        rhog = self.density(x1, x2, x3)
        fraction = self.dust.massfraction()

        rhod = np.ones(
            (len(fraction), len(x1), len(x2), len(x3))
        )

        for i in range(len(fraction)):
            rhod[i] = fraction[i] * self.params.dtogas * rhog

        return rhod

    def numberdensity_d(self, x1, x2, x3=None):
        """
        D)
        Return dust number density n_d(r, z, a).
        Units: cm^-3
        """
        mass = self.dust.grainmass()
        dens = self.density_d(x1, x2, x3)

        for i in range(len(mass)):
            dens[i] /= mass[i]

        return dens

    def solution(self, mu0, r, mu_val):
        """
        E)
        Return the solution function used in the root finding.
        """
        return (
            mu0**3
            - mu0 * (1 - r / self.params.r_centri)
            - mu_val * (r / self.params.r_centri)
        )
