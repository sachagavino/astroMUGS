# chemdiskpy: Couple chemistry modeling with thermal modeling of protoplanetary disks with multiple grain sizes


This code performs dust continuum radiative transfer simulations of protoplanetary disks with RADMC3D and allows to use the resulting dust tempeature and local radiation field for gas-grain simulations (NAUTILUS).
The code can convert any relevant RADMC3D input/output files to a NAUTILUS friendly format, and vice-versa. 


## Installation

- latest release from test [https://test.pypi.org/](https://test.pypi.org/):

        pip install -i https://test.pypi.org/simple/ chemdiskpy


## Documentation

        A detailed documentation can be found [here][1]. 

## Quick start


### Set up a directory 
1. Create a working folder where you want to create a model and go to this folder.
2. Add the opacity table in RADM3D format (or use the script provided).
3. Running a model will automatically create a folder **chemistry/** and **thermal/** 

## Set up a model
- All parameters can be written in the file parameters.py in the working directory. See documentation and folder **example_simulation/**. 
- If the user provides all or a few parameters when calling the functions (e.g. if the user calls: **m.add_star(mass=star_mass, luminosity=1., temperature=4500., x=0., y=0., z=0.)**), then the corresponding values in **parameters.py** are 'bypassed'.


## EXAMPLE:

### IMPORT PACKAGE
- Open any notebook or script in the working directory:

```
import chemdiskpy.modeling as modeling
import chemdiskpy.dust as dust
import chemdiskpy.plotting.plot as plot
```

### CREATE AN OBJECT:
```
m = modeling.Structure() 
```

### DISK PARAMETERS
```
nr, ntheta, nphi = 301, 181, 2
rin, rout= 1, 500
dtogas = 1e-2
star_mass = 1. #solar mass
dust_mass = 1e-4 #solar mass
q_c = 12 # extinction efficiency at resonance. Usefull only if you don't compute the local field in RADMC3D.
```

### WAVELENGTH GRID
```
m.grid.set_wavelength_grid(1e-1, 2e3, 100, log=True) # in microns.
m.grid.set_mcmonowavelength_grid(1e-1, 2e3, 100, log=True) # in microns. This grid will be used for the computing of the local radiation field.
```

### RADMC3D GRID
- Typically, the grid is spherical for model with a central object like disks (but can be cartesian if needed):

```
m.set_spherical_grid(rin, rout, nr, ntheta, nphi, log=True)
```

### STAR
- This is a simple black-body spectrum. UV from accreting material will be added in future versions. The star is located at coordinates (0,0,0)
```
m.add_star(mass=star_mass, luminosity=0.8, temperature=3900., x=0., y=0., z=0.)
```

### INTERSTELLAR RADIATION FIELD (optional)
- Draine 1978 between 91.2 and 200 nm, and with the extension of van Dishoeck & Black 1982 at longer wavelengths. Change 'cut' if you want another value.
```
m.add_isrf(cut=2.e-1, d78=True, vdb82=True)
```

### DUST DISK STRUCTURE
- Start by creating a dust population. This library is originally made for multi-grain chemical code but this works fine with a single grain population.
- Set **nb_sizes=1** if you want to use one size. In that case, **amin** and **amax** will not be read and **rsingle** is the grain size.
```
d = dust.DustDistrib(rsingle=0.1, amin=5.000e-03, amax=1.000e+03, nb_sizes=2, rho_m=3.) #microns
sizes = d.sizes() # get the size 
mass = d.grainmass() # get the grain mass for each grain size. 
```

- If **nb_sizes = 1**, then the grain size read is rsingle and amin and amax are not read.
- If **nb_sizes > 1**, then rsingle is not read and the size range is given by amin and amax.

- Then, create a dust disk structure with spherical coordinates for RADMC3D with the values previously set. 
```
m.add_disk(dust=d, rin=rin, rout=rout, star_mass=star_mass, dtogas=dtogas, dust_mass=dust_mass, q_c=q_c, settling=True, coordsystem='spherical')
```

### RUN THERMAL
- The most important parameters are the number of photons **nphot** and **nphot_scat**. 
```
m.thermal(run=True, nphot = 2e7, \
              nphot_scat = 2e7, \
              nphot_spec = 1e5, \
              itempdecoup = 1, \
              istar_sphere = 1, \
              modified_random_walk = 1, \
              rto_style = 1, \
              scattering_mode_max = 2, \
              tgas_eq_tdust = 0)
```

- The parameter **run** is a flag. If True, then **radmc3d mctherm** is run and the dust continuum RT is computed. If False, the user assumes the **dust_temperature.dat** already exists in the model and there is no need to compute the dust temperature. However,  **m.thermal()** is still called because the thermal model will be used to compute the NAUTILUS model. 

### PLOT THE RESULTS (optional)
- Dust temperature:
```
plot.temperature2D() 
```

- Dust density:
```
plot.density2D(mass) # Use the mass previously set up. 
```

- Radial dust temperature profile (in the midplane):
```
plot.midplane_temp() 
```

- Dust opacity (absorption, scattering, angles)
```
plot.opacity() 
```

### RUN LOCAL RADIATION (optional)
```
nphot_mono = 2000000 # choose number of photons sufficiently large.
m.localfield(run=True, nphot_mono = nphot_mono)
```
- Just like thermal(), if **run==True**, then the local radiation field is computed using **radmc3d mcmono**. If it is False, then the user assumes the local radiation field is already computed (mean_intensity.out exists) and does not need to re-compute it again.

- **WARNING**: If you want to compute the local field, make sure the scattering mode is set to a value lower than 2 or greater than 4. To select the mode, set **scattering_mode_max** to the requested value (see RADMC3D documentation) in the **run_thermal** instance above. If you don't set a value, the default value is 1.  


### NAUTILUS GRID
- radii in au where you want to compute chemistry:
```
rchem = [10, 20, 30, 50, 80, 100, 140, 160, 200]
```
- create NAUTILUS grid. 
- The value **max_h** is the maximum scale height (z/H) where chemistry will be computed. The grid is empty above. 
- The value **nz_chem** is the number of vertical points where chemistry is computed for each radius.
```  
m.grid.set_nautilus_grid(rchem, max_h=4, nz_chem=64)
```

### CREATE NAUTILUS DISK STRUCTURE
- Use the same dust model as for RADMC3D. Make sure the disk parameters are the same as provided in the thermal disk model (m.add_disk(...)). 
```
m.add_nautilusdisk(dust=d, rin=rin, rout=rout,  star_mass=star_mass, dust_mass=dust_mass, q_c=q_c, dtogas=dtogas, settling=True)
```

### CREATE NAUTILUS DISK MODEL
- Create a ready-to-use nautilus disk model in **chemistry/**. The dust temperature (and local flux if needed) is extraced from the RADMC3D results.
```
m.write_nautilus(sizes=sizes, uv_ref=3400, dtogas=dtogas, ref_radius=m.disk.ref_radius, cr_ionisation_rate=1.300E-17, network=False, coupling_av=True)
```

- **network=False** implies that the user has to add his/her own chemical network. If True, the kida.14 network will be added in the user's model.
- WARNING: the parameter **coupling_av=True** works only if the local field is computed prior to this.

### PLOT THE RESULTS (optional)
- Vertical dust temperature:
```
plot.vertical_temp(r=100) 
```

### RUN NAUTILUS DISK MODEL
- Next step is to run as you would usually do the nautilus model in **chemistry/**.
- This will be added in a future update.

[1]: https://chemdiskpy.readthedocs.io/en/latest