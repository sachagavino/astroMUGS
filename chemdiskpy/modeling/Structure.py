import sys, os, inspect
import numpy as np

from .Model import Model
from .Star import Star
from .Disk import Disk
from .Envelope import Envelope
from .InterstellarRadFields import InterstellarRadFields
from .. constants.constants import autocm, M_sun, R_sun, L_sun

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import parameters as p

class Structure(Model):

    def add_star(self, mass=0.5, luminosity=1, temperature=4000., x=0., y=0., z=0.):
        self.grid.add_star(Star(mass=mass, luminosity=luminosity, \
                temperature=temperature, x=x, y=y, z=z))

    def add_isrf(self, cut=2.e-1, d78=True, vdb82=True):
        self.isrf = InterstellarRadFields(cut, d78, vdb82)
        self.grid.add_isrf(self.isrf.create_isrf(self.grid.lam))

    def set_spherical_grid(self, rmin, rmax, nr, ntheta, nphi, log=True):
        if log:
            r = np.logspace(np.log10(rmin), np.log10(rmax), nr, base=10)
        else:
            r = np.linspace(rmin, rmax, nr)

        theta = np.linspace(0.0, np.pi, ntheta)
        phi = np.linspace(0.0, 2*np.pi, nphi)

        self.grid.set_spherical_grid(r, theta, phi)


    def add_disk(self, ref_radius=p.ref_radius, rin=p.rin, rout=p.rout, star_mass=p.star_mass, disk_mass=p.disk_mass, h0=p.h0, \
                       sigma_gas_ref=p.sigma_gas_ref, Tmidplan_ref=p.Tmidplan_ref, Tatmos_ref=p.Tatmos_ref, sigma_t = p.sigma_t, q_exp=p.q_exp,  \
                       d_exp=p.d_exp, p_exp=p.p_exp, dtogas=p.dtogas, rho_m=p.rho_m, schmidtnumber=p.schmidtnumber, alpha=p.alpha, \
                       settfact=p.settfact, dust_mass=p.disk_dust_mass, q_c = p.q_c, dust=None, \
                       settling=True, isothermal=False, dust_density='g.cm-2', coordsystem='spherical'):
        self.disk = Disk(ref_radius=ref_radius, rin=rin, rout=rout, star_mass=star_mass, disk_mass=disk_mass, h0=h0, \
                       sigma_gas_ref=sigma_gas_ref, Tmidplan_ref=Tmidplan_ref, Tatmos_ref=Tatmos_ref, sigma_t=sigma_t, q_exp=q_exp,  \
                       d_exp=d_exp, p_exp=p_exp, dtogas=dtogas, rho_m=rho_m, schmidtnumber=schmidtnumber, alpha=alpha, \
                       settfact=settfact, dust_mass=dust_mass, q_c = q_c, dust=dust, \
                       settling=settling, isothermal=isothermal, dust_density=dust_density, coordsystem=coordsystem)

        if (dust != None):
            self.grid.add_dustdensity(self.disk.density_d(self.grid.r, self.grid.theta, self.grid.phi))
        else:
            print('WARNING: no dust model as input.')

    def add_internalheating(self, acc_rate=p.acc_rate, lim_h=p.lim_h):
            self.grid.add_accretionheating(self.disk.viscous_accretion_heating(acc_rate, lim_h, self.grid.r, self.grid.theta, self.grid.phi))


    def add_chemdisk(self, ref_radius=p.ref_radius, rin=p.rin, rout=p.rout, star_mass=p.star_mass, disk_mass=p.disk_mass, h0=p.h0, \
                       sigma_gas_ref=p.sigma_gas_ref, Tmidplan_ref=p.Tmidplan_ref, Tatmos_ref=p.Tatmos_ref, sigma_t = p.sigma_t, q_exp=p.q_exp,  \
                       d_exp=p.d_exp, p_exp=p.p_exp, dtogas=p.dtogas, rho_m=p.rho_m, schmidtnumber=p.schmidtnumber, alpha=p.alpha, \
                       settfact=p.settfact, max_H=p.max_H, nz_chem=p.nz_chem, dust_mass=p.disk_dust_mass, q_c = p.q_c, dust=None, \
                       settling=True, isothermal=False, dust_density='g.cm-2', coordsystem='nautilus'):
        self.chemdisk = Disk(ref_radius=ref_radius, rin=rin, rout=rout, star_mass=star_mass, disk_mass=disk_mass, h0=h0, \
                       sigma_gas_ref=sigma_gas_ref, Tmidplan_ref=Tmidplan_ref, Tatmos_ref=Tatmos_ref, sigma_t=sigma_t, q_exp=q_exp,  \
                       d_exp=d_exp, p_exp=p_exp, dtogas=dtogas, rho_m=rho_m, schmidtnumber=schmidtnumber, alpha=alpha, \
                       settfact=settfact, max_H=max_H, nz_chem=nz_chem, dust_mass=dust_mass, q_c = q_c, dust=dust, \
                       settling=settling, isothermal=isothermal, dust_density=dust_density, coordsystem=coordsystem)

        self.grid.add_dustdensity_chem(self.chemdisk.numberdensity_d(self.grid.rchem, self.grid.zchem))
        self.grid.add_dustdensity_single_chem(self.chemdisk.numberdensity_d_single(self.grid.rchem, self.grid.zchem))
        self.grid.add_gasdensity_chem(self.chemdisk.numberdensity(self.grid.rchem, self.grid.zchem))
        self.grid.add_gastemperature_chem(self.chemdisk.temp_altitude(self.grid.rchem, self.grid.zchem))
        self.grid.add_hg_chem(self.chemdisk.scaleheight(self.grid.rchem))
        self.grid.add_avz(self.chemdisk.av_z(self.grid.lam, self.grid.dustdensity_chem[0], self.grid.rchem, self.grid.zchem))

    def add_envelope(self,rmin=p.rmin, rmax=p.rmax, r_centri=p.r_centri, acc_rate=p.acc_rate, star_mass=p.star_mass, dust_mass=p.dust_env_mass, dtogas=p.dtogas, \
                    cavpl=p.cavpl, cav_fact=p.cav_fact, cavz0=p.cavz0, dust=None, dust_density='g.cm-2', coordsystem='spherical'):
        self.envelope = Envelope(rmin=rmin, rmax=rmax, r_centri=r_centri, \
                       acc_rate=acc_rate, star_mass=star_mass, dust_mass=dust_mass, \
                       dtogas=dtogas, cavpl=cavpl, cav_fact=cav_fact, cavz0=cavz0, dust=dust, dust_density=dust_density, coordsystem=coordsystem)

        if (dust != None):
            self.grid.add_dustdensity(self.envelope.density_d(self.grid.r, self.grid.theta, self.grid.phi))


