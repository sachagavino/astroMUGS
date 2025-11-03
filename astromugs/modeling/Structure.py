import sys 
import os 
import inspect
import numpy as np

from astromugs.utils.params import StructureParams
from astromugs.modeling.Model import Model
from astromugs.modeling.Star import Star
from astromugs.modeling.Disk import Disk
from astromugs.modeling.Envelope import Envelope
from astromugs.modeling.InterstellarRadFields import InterstellarRadFields
from astromugs.constants.constants import autocm, M_sun, R_sun, L_sun

#currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#sys.path.insert(0,parentdir) 


class Structure(Model):
    def __init__(self):
        self.params = StructureParams()
        self.envelope_params = self.params.envelope
        self.disk_params = self.params.disk

        

    def add_chemical_path(self, chemicalpath):
        self.chempath = chemicalpath

    def add_star(self, mass=0.5, luminosity=1, temperature=4000., x=0., y=0., z=0.):
        self.grid.add_star(Star(mass=mass, luminosity=luminosity, \
                temperature=temperature, x=x, y=y, z=z))

    def add_isrf(self, cut=2.e-1, d78=True, vdb82=True):
        self.isrf = InterstellarRadFields(cut, d78, vdb82)
        self.grid.add_isrf(self.isrf.create_isrf(self.grid.lam))

    def add_disk(self,  dust=None, **kwargs):
        # Create a copy of the defaults
        params = DiskParams()

        # Apply overrides given by user
        for key, val in kwargs.items():
            setattr(params, key, val)

        self.disk = Disk(params=params, dust=dust)

        if (len(self.grid.dust) > 0):
            self.grid.add_dustdensity(self.disk.density_d(self.grid.r, self.grid.theta, self.grid.phi))
        else:
            print('WARNING: no dust model as input. add your dust model in the grid (grid.add_dust(dust)) before creating the disk structure.')

    def add_internalheating(self, **kwargs):
        # Create a copy of the defaults
        params = DiskParams()       
        # Apply overrides given by user
        for key, val in kwargs.items():
            setattr(params, key, val)

        self.grid.add_accretionheating(self.disk.viscous_accretion_heating(self.grid.r, self.grid.theta, self.grid.phi))

    def add_chemdisk(self, dust=None, **kwargs):
        # Create a copy of the defaults
        params = DiskParams()

        # Apply overrides given by user
        for key, val in kwargs.items():
            setattr(params, key, val)

        self.chemdisk = Disk(params=params, dust=dust)

        self.grid.add_dustdensity_chem(self.chemdisk.numberdensity_d(self.grid.rchem, self.grid.zchem))
        self.grid.add_dustdensity_single_chem(self.chemdisk.numberdensity_d_single(self.grid.rchem, self.grid.zchem))
        self.grid.add_gasdensity_chem(self.chemdisk.numberdensity(self.grid.rchem, self.grid.zchem))
        self.grid.add_gastemperature_chem(self.chemdisk.temp_altitude(self.grid.rchem, self.grid.zchem))
        self.grid.add_hg_chem(self.chemdisk.scaleheight(self.grid.rchem))
        self.grid.add_avz(self.chemdisk.av_z(self.grid.lam, self.grid.dustdensity_chem[0], self.grid.rchem, self.grid.zchem))

    # def add_envelope(self, dust=None, **kwargs):
    #     # Create a copy of the defaults
    #     params = EnvelopeParams()

    #     # Apply overrides given by user
    #     for key, val in kwargs.items():
    #         setattr(params, key, val)

    #     self.envelope = Envelope(params=params, dust=dust)

    #     if self.grid.dust:
    #         self.grid.add_dustdensity(
    #             self.envelope.density_d(self.grid.r, self.grid.theta, self.grid.phi)
    #         )


    def add_envelope(self, dust=None):
        self.envelope = Envelope(self.envelope_params, dust=dust)
        if dust:
            self.grid.add_dustdensity(
                self.envelope.density_d(self.grid.r, self.grid.theta, self.grid.phi)
            )


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
        



       