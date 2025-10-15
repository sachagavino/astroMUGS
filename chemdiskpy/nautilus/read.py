import os
import numpy as np
import pandas as pd



def grid(radlist, nb_sizes, itime, species, chempath):

    columns_static = ['z', 'nH', 'Tg', 'Av', 'diff', 'Td', 'inv_ab', 'conv_factor', 'a']
    grid_dict = {}  # Dictionary to store z arrays

    for radius in radlist[0]:
        static = pd.read_table(f'{chempath}{radius}AU/1D_static.dat', delimiter=r"\s+", comment='!', header=None, engine='python')
        static.columns = columns_static
        z = static['z'].values  # Extract the z column as a NumPy array
        nh = static['nH'].values  # Extract the nH column as a NumPy array
        tg = static['Tg'].values  # Extract the Tg column as a NumPy array

        if nb_sizes == 1:
            td = static['Td'].values
            ab_d = static['inv_ab'].values
            grain_size = static['a'].values
            
            abundance = np.loadtxt(f'{chempath}{radius}AU/ab/{species}.ab', comments='!')
            abundance = np.delete(abundance[itime], 0)
            chem_numdens = nh / abundance

            # Store z, nH, and Tg arrays in a nested dictionary for the current radius
            grid_dict[radius] = {
                'z': z,
                'nH': nh,
                'Tg': tg,
                'nd': ab_d,  # Store the extracted grain sizes
                'Td': td.transpose(grain_temperatures),  # Store the extracted grain sizes
                'numberdens_species' : chem_numdens
            }

        elif nb_sizes > 1:
            # Check if the '1D_grain_sizes.in' file exists
            grain_file_path = f'{chempath}{radius}AU/1D_grain_sizes.in'
            grain_sizes = None  # Default value if the file does not exist
            grain_densities = None
            grain_temperatures = None
            if os.path.exists(grain_file_path):
                # Read the file and extract the first nb_sizes columns
                grain_data = pd.read_table(grain_file_path, delimiter=r"\s+", comment='!', header=None, engine='python')
                grain_abundances = grain_data.iloc[:, nb_sizes:2*nb_sizes].values  # Extract the first nb_sizes columns as a NumPy array
                grain_temperatures = grain_data.iloc[:, 2*nb_sizes:3*nb_sizes].values

                grain_densities = nh / np.transpose(grain_abundances)  # Broadcasting division


            abundance = np.loadtxt(f'{chempath}{radius}AU/ab/{species}.ab', comments='!')
            abundance = np.delete(abundance[itime], 0)
            chem_numdens = nh / abundance

            # Store z, nH, and Tg arrays in a nested dictionary for the current radius
            grid_dict[radius] = {
                'z': z,
                'nH': nh,
                'Tg': tg,
                'grain_numdens': grain_densities,  # Store the extracted grain sizes
                'grain_temperatures': np.transpose(grain_temperatures),  # Store the extracted grain sizes
                'numberdens_species' : chem_numdens
            }

    return grid_dict


def parameters(radlist, chempath):
    chem_parameters = {}  # Dictionary to store the parameters
    with open(f'{chempath}{radlist[0][0]}AU/parameters.in', 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace
            if not line or line.startswith('!'):  # Skip blank lines or comments
                continue

            # Split the line into key and value at the '=' character
            if '=' in line:
                key, value = line.split('=', 1)  # Split only at the first '='
                key = key.strip()  # Remove whitespace around the key
                value = value.split('!')[0].strip()  # Remove comments after the value
                chem_parameters[key] = value  # Add the key-value pair to the dictionary

    return chem_parameters


def radii(chempath):
    radii = os.listdir (chempath) # get all files' and folders' names in the current directory
    radlist = []
    for radius in radii:
        if radius.endswith("AU"):
            rad_d = radius.replace('AU','')
            rad_d = int(rad_d)
            radlist.append(rad_d)

    radlist = sorted(radlist)
    return radlist

# def abundance(path, itime=47, species='CO'):
#     chemmodel = pd.read_table(path + 'disk_t{}.dat'.format(str(itime)), sep=" ", engine='python')
#     chemmodel.dropna(how='all',inplace=True)
#     chemmodel.reset_index(inplace=True)

#     nr = chemmodel['r'].nunique()
#     nz = int(len(chemmodel['r'])/nr)

#     r = chemmodel['r'].unique()
#     z = np.reshape(chemmodel['z'].values, (nr, nz))
#     z = np.fliplr(z)
#     z = np.transpose(z)
#     zz = np.flipud(z)

#     nH = np.reshape(chemmodel['nH'].values, (nr, nz))
#     nH = np.transpose(nH)
#     ab = np.reshape(chemmodel['ab({})'.format(species)].values, (nr, nz))
#     ab = np.transpose(ab)
#     dens_mol = ab*nH #number density cm-3

#     return r, zz, dens_mol