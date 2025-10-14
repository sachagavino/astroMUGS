import numpy as np
import pandas as pd

def abundance(path, itime=47, species='CO'):
    chemmodel = pd.read_table(path + 'disk_t{}.dat'.format(str(itime)), sep=" ", engine='python')
    chemmodel.dropna(how='all',inplace=True)
    chemmodel.reset_index(inplace=True)

    nr = chemmodel['r'].nunique()
    nz = int(len(chemmodel['r'])/nr)

    r = chemmodel['r'].unique()
    z = np.reshape(chemmodel['z'].values, (nr, nz))
    z = np.fliplr(z)
    z = np.transpose(z)
    zz = np.flipud(z)

    nH = np.reshape(chemmodel['nH'].values, (nr, nz))
    nH = np.transpose(nH)
    ab = np.reshape(chemmodel['ab({})'.format(species)].values, (nr, nz))
    ab = np.transpose(ab)
    dens_mol = ab*nH #number density cm-3

    return r, zz, dens_mol