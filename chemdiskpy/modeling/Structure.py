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

        if (len(self.grid.dust) > 0):
            self.grid.add_dustdensity(self.disk.density_d(self.grid.r, self.grid.theta, self.grid.phi))
        else:
            print('WARNING: no dust model as input. add your dust model in the grid (grid.add_dust(dust)) before creating the disk structure.')

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


    def add_chemmodel(self, chempath="chemistry/", itime=0, species='CO', reader=None):
        """
        Add an existing chemistry model to the object.

        Args:
            chempath (str): Relative or absolute path to the chemistry model directory.
            structure_type (str): Type of structure ('0D' or '1D'). Currently, only '1D' is supported.
            itime (int): Index of the time output to use from the chemistry model.
            species (str): Chemical species of interest (e.g., 'CO', 'H2O'). Ultimately it will take all species so that arg will disappear.
            reader (object, optional): A reader object or module responsible for reading chemistry data. 
                                    Defaults to self.nautilus.read.

        Raises:
            ValueError: If an unsupported structure type is provided.
            FileNotFoundError: If required files are missing in the specified chempath.
            KeyError: If required parameters are missing in the chemistry data.
            RuntimeError: For any other errors encountered during processing.
        """
        # if structure_type != "1D":
        #     raise ValueError(f"Unsupported structure type: {structure_type}. Only '1D' is supported.")

        if reader is None:
            reader = self.nautilus.read  # Default to the current reader if none is provided

        try:
            radii = reader.radii(chempath=chempath)
            self.grid.add_existingchemradii(radii)

            parameters = reader.parameters(self.grid.chemradii, chempath=chempath)
            if 'nb_grains_1D' not in parameters:
                raise KeyError("The parameter 'nb_grains_1D' is missing in the chemistry parameters.")
            self.grid.add_existingchemparam(parameters)

            grid = reader.grid(radlist=self.grid.chemradii, nb_sizes=int(parameters['nb_grains_1D']), itime=itime, species=species, chempath=chempath)
            self.grid.add_existingchemmodel(grid, species)
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Required file not found in {chempath}: {e}")
        except KeyError as e:
            raise KeyError(f"Missing required parameter in chemistry data: {e}")
        except Exception as e:
            raise RuntimeError(f"An error occurred while adding the chemistry model: {e}")
        



       