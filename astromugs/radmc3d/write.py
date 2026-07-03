from __future__ import annotations

import numpy as np
from dataclasses import fields as dataclass_fields
from typing import TYPE_CHECKING

from .. constants.constants import M_sun, Ggram, autocm

if TYPE_CHECKING:
    from ..utils.thermal import ControlParams  # only needed for annotation


def control(params: ControlParams, thermpath='thermal/'):
    with open(thermpath + "radmc3d.inp", "w") as f:
        for fld in dataclass_fields(params):
            value = getattr(params, fld.name)
            if value is not None:
                if isinstance(value, float) and value.is_integer():
                    f.write(f"{fld.name} = {int(value)}\n")
                else:
                    f.write(f"{fld.name} = {value}\n")

def stars(rstar, mstar, lam, xstar, ystar, zstar, tstar=None, fstar=None, thermpath='thermal/'):
    nstars = len(rstar)
    nlam = len(lam)

    f = open(thermpath + "stars.inp","w")

    f.write(str(2)+"\n")
    f.write("{0:d}  {1:d}\n".format(nstars, nlam))

    for istar in range (nstars):
        f.write("{0:e}   {1:e}   {2:e}   {3:e}   {4:e}\n".format(rstar[istar], \
                mstar[istar], xstar[istar], ystar[istar], zstar[istar]))

    for ilam in range(nlam):
        f.write("{0:e}\n".format(lam[ilam]))

    for istar in range(nstars):
        if (tstar[istar] != 0):
            f.write("{0:f}\n".format(-tstar[istar]))
        else:
            for i in range(nlam):
                f.write("{0:e}\n".format(fstar[ilam]))

    f.close()

def wavelength_micron(lam, thermpath='thermal/'):
    nlam = len(lam)
    f = open(thermpath + "wavelength_micron.inp","w")
    f.write("{0:d}\n".format(nlam))
    for ilam in range(nlam):
        f.write("{0:f}\n".format(lam[ilam]))
    f.close()

def mcmono_wavelength_micron(lam_mono, thermpath='thermal/'):
    nlam = len(lam_mono)

    f = open(thermpath + "mcmono_wavelength_micron.inp","w")

    f.write("{0:d}\n".format(nlam))
    for ilam in range(nlam):
        f.write("{0:f}\n".format(lam_mono[ilam]))

    f.close()

def amr_grid(x, y, z, gridstyle="regular", coordsystem="cartesian", thermpath='thermal/'):
    nx = x.size-1
    ny = y.size-1
    nz = z.size-1

    incl_x = int(nx > 1)
    incl_y = int(ny > 1)
    incl_z = int(nz > 1)

    f = open(thermpath + "amr_grid.inp","w")

    f.write(str(1)+"\n")

    if (gridstyle == "regular"):
        f.write("0\n")
    elif (gridstyle == "octtree"):
        f.write("1\n")
    elif (gridstyle == "amr"):
        f.write("10\n")

    if (coordsystem == "cartesian"):
        f.write("0\n")
    elif (coordsystem == "spherical"):
        f.write("100\n")
    elif (coordsystem == "cylindrical"):
        f.write("200\n")

    f.write("0\n")
    f.write("{0:d}  {1:d}  {2:d}\n".format(incl_x, incl_y, incl_z))
    f.write("{0:d}  {1:d}  {2:d}\n".format(nx, ny, nz))

    if (gridstyle == "octtree"):
        print("OctTree grids not yet implemented.")
    elif (gridstyle == "amr"):
        print("Layer-style AMR grids not yet implemented.")

    for i in range(nx+1):
        f.write("{0:12.9e}\n".format(x[i]))
    for i in range(ny+1):
        f.write("{0:12.9e}\n".format(y[i]))
    for i in range(nz+1):
        f.write("{0:12.9e}\n".format(z[i]))

    # Insert extra info for octtree and amr grids here...

    f.close()

def dustopac(opacity, thermpath='thermal/'):
    '''
    Desc: write dustopac.inp
    Args: opacity
    '''
    import os
    nspecies = len(opacity)

    f = open(thermpath + "dustopac.inp","w")

    f.write("2\n")
    f.write("{0:d}\n".format(nspecies))
    f.write("==============================================================\n")
    for i in range(nspecies):
        basename = os.path.basename(opacity[i])  # e.g. "dustkappa_silicate.inp"
        filetype = basename.split("_")[0]
        species = basename.split("_")[1].split(".")[0]

        if (filetype == "dustkappa"):
            f.write("1\n")
        elif (filetype == "dustkapscatmat"):
            f.write("10\n")
        elif (filetype == "dustopac"):
            f.write("-1\n")

        f.write("0\n")
        f.write("{0:s}\n".format(species))
        f.write("----------------------------------------------------------\n")

    f.close()


def dust_density(density, gridstyle="regular", thermpath='thermal/'):
    '''
    Desc: write dust_density.inp
    Args: density
    '''
    nstructures = len(density) #disk, envelope...
    print('number of structures (disk, envelope...): ', nstructures, '\n')
    nspecies = 0
    for istruc in range(nstructures):
        nspecies += len(density[istruc]) #nb of species in all structures
        print('number of species in structure {}: '.format(istruc+1), len(density[istruc]))
    print('total number of grain species: ', nspecies)

    if (gridstyle == "regular"):
        nx, ny, nz = density[0][0].shape
        ncells = nx*ny*nz

    f = open(thermpath + "dust_density.inp","w")
    f.write("1\n")
    f.write("{0:d}\n".format(ncells))
    f.write("{0:d}\n".format(nspecies))

    for istruc in range(nstructures): #loop over structures (disk, envelope...)
        for ispec in range(len(density[istruc])): #loop over the dust species within the given structure.
            if (gridstyle == "regular"):
                for iz in range(nz):
                    for iy in range(ny):
                        for ix in range(nx):
                            f.write("{0:e}\n".format(density[istruc][ispec, ix,iy,iz]))
    f.close()


def external_rad(isrf, thermpath='thermal/'):
    '''
    Desc: write external_source.inp
    Args: Interstellar radiation field
    '''
    nlam = len(isrf[0])
    factor = 1e0 #artificially multiply by a factor just to see the impact of different ISRF intensities. Will be removed in future updates
    f = open(thermpath + "external_source.inp","w")

    f.write("2\n")
    f.write("{0:d}\n".format(nlam))
    for ilam in range(nlam):
        f.write("{0:f}\n".format(isrf[0,ilam]*1e-3))
    for ilam in range(nlam):
        f.write("{0:e}\n".format(isrf[1,ilam]*factor))

    f.close()

def accretion_heating(x, y, z, accretionheating, gridstyle="regular", thermpath='thermal/'):
    '''
    Desc: Write viscous accretion heating in the file heatsource.inp
    Args: Viscous accretion heating structure
    '''
    nx = x.size-1
    ny = y.size-1
    nz = z.size-1
    f = open(thermpath + "heatsource.inp","w")
    f.write("1\n")
    f.write("{0:d}\n".format(nx*ny*nz))

    if (gridstyle == "regular"):
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    f.write("{0:e}\n".format(accretionheating[ix,iy,iz]))

    f.close()


def numberdens_mol(numberdens, species='CO', gridstyle="regular", thermpath='thermal/'):
    '''
    Desc: write numberdens_XXX.inp, where XXX is a chemical species
    Args:
    '''
    if (gridstyle == "regular"):
        nx, ny = numberdens.shape
        ncells = nx*ny

    print('writing numberdens_{}.inp...'.format(species))

    f = open(thermpath + "numberdens_{}.inp".format(species),"w")
    f.write("1\n")
    f.write("{0:d}\n".format(ncells))
    if (gridstyle == "regular"):
            for iy in range(ny):
                for ix in range(nx):
                    f.write("{0:e}\n".format(numberdens[ix,iy]))
    f.close()

def gas_temperature(temp, thermpath='thermal/'):
    """Write ``gas_temperature.inp`` for RADMC-3D line radiative transfer.

    Parameters
    ----------
    temp : ndarray, shape (nr, nt)
        Gas temperature in K on the spherical grid, indexed ``[ir, itheta]``.
    thermpath : str, optional
        Output directory. Default is ``'thermal/'``.
    """
    nr, nt = temp.shape
    ncells  = nr * nt
    print(f'writing gas_temperature.inp  ({ncells} cells)...')
    with open(thermpath + 'gas_temperature.inp', 'w') as f:
        f.write("1\n")
        f.write(f"{ncells}\n")
        for iy in range(nt):
            for ix in range(nr):
                f.write(f"{temp[ix, iy]:.6e}\n")


def lines(species='CO', format='leiden', thermpath='thermal/'):
    '''
    Desc: write lines.inp
    Args: species, format
    '''
    f = open(thermpath + "lines.inp","w")
    f.write("2\n")
    f.write("1\n")
    f.write("{} {} 0 0 0".format(species,format))
    f.close()


def gas_velocity(star_mass, r, theta, phi, object="disk", thermpath='thermal/'):
    """Write ``gas_velocity.inp`` for RADMC-3D line transfer (Keplerian disk).

    Computes the Keplerian azimuthal velocity at each grid cell and writes
    it in the spherical-coordinate velocity format expected by RADMC-3D:
    ``(v_r, v_theta, v_phi)`` in cm s :sup:`-1`. For a Keplerian disk
    ``v_r = v_theta = 0`` and ``v_phi = sqrt(G M_star / R_cyl)`` where
    ``R_cyl = r sin(theta)`` is the cylindrical radius.

    Parameters
    ----------
    star_mass : float
        Stellar mass in solar masses.
    r : array_like
        Spherical radial cell centres **in AU**.
    theta : array_like
        Co-latitude cell centres in radians.
    phi : array_like
        Azimuthal cell centres in radians.
    object : str, optional
        Reserved for future use (e.g. envelope kinematics). Default is
        ``'disk'``.
    thermpath : str, optional
        Output directory. Default is ``'thermal/'``.
    """
    nr, ntheta, nphi = len(r), len(theta), len(phi)

    # 3-D spherical grids; rr in cm
    rr, tt, _ = np.meshgrid(r * autocm, theta, phi, indexing='ij')

    # Cylindrical radius R_cyl = r * sin(theta) [cm]
    R_cyl = rr * np.sin(tt)

    # Keplerian azimuthal velocity v_phi = sqrt(G M / R_cyl) [cm/s]
    vphi = np.sqrt(Ggram * star_mass * M_sun / R_cyl)

    # RADMC-3D spherical-grid format: one line per cell with (v_r, v_theta, v_phi)
    # Keplerian rotation: v_r = 0, v_theta = 0, v_phi = v_K
    with open(thermpath + 'gas_velocity.inp', 'w') as f:
        f.write("1\n")                           # iformat
        f.write(f"{nr * ntheta * nphi}\n")
        for k in range(nphi):
            for j in range(ntheta):
                for i in range(nr):
                    f.write(f"0.0e+00 0.0e+00 {vphi[i,j,k]:.6e}\n")
