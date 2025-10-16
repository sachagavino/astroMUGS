import numpy as np
from ..constants.constants import c, autocm, amu, mu, black_body

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w

def dust_density(dtogas, rho_m, sizes, dens_radmc3d, rchem, zchem, d, theta): 
    nbspecies, nz, ny, nx = dens_radmc3d.shape
    dens_naut = np.zeros((len(dens_radmc3d), len(rchem), len(zchem[0,:])))
    dens_naut_nd_smooth = np.zeros((len(dens_radmc3d), len(rchem), len(zchem[0,:])))
    dens_naut_nH_smooth = np.zeros((len(rchem), len(zchem[0,:])))
    dens_naut_nd = np.zeros((len(dens_radmc3d), len(rchem), len(zchem[0,:])))
    dens_naut_nH = np.zeros((len(rchem), len(zchem[0,:])))
    for idx_size in range(0, nbspecies, 1):
        for idx, r in enumerate(rchem):
            for idz, z in enumerate(zchem[idx, :]):
                d_pt = np.sqrt(r**2 + z**2)  #convert from cartesian to spherical
                theta_pt = np.arccos(z/d_pt) #convert from cartesian to spherical
                closest_d = min(enumerate(d), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
                closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
                dens_naut[idx_size, idx, idz] = dens_radmc3d[idx_size, 0, closest_t[0], closest_d[0]] 

    #get gas number density:
    for ai in range(0, len(dens_radmc3d), 1):
        dens_naut_nH += dens_naut[ai, :, :]
    dens_naut_nH = (1/dtogas)*dens_naut_nH
    dens_naut_nH = dens_naut_nH/(mu*amu)

    #then convert g.cm-3 to number density for the grains too:
    for ai in range(0, len(dens_radmc3d), 1):
        mass = (4/3)*np.pi*rho_m*sizes[ai]**3
        dens_naut_nd[ai, :, :] = dens_naut[ai, :, :]/mass


    # for idx in range(len(rchem)):
    #     for size_id in range(nbspecies):
    #         dens_naut_nd_smooth[size_id, idx, :] = moving_average(dens_naut_nd[size_id, idx, :], 5) #average the values over a rolling window of 5 points.
    # dens_naut_nd_smooth[:, :, 0:2] = dens_naut_nd[:, :, 0:2]  #clean the boundary effects
    # dens_naut_nd_smooth[:, :, -2:] = dens_naut_nd[:, :, -2:] 


    # for idx in range(len(rchem)):
    #     dens_naut_nH_smooth[idx, :] = moving_average(dens_naut_nH[idx, :], 5) #average the values over a rolling window of 5 points.
    # dens_naut_nH_smooth[:, 0:2] = dens_naut_nH[:, 0:2]  #clean the boundary effects
    # dens_naut_nH_smooth[:, -2:] = dens_naut_nH[:, -2:] 

    # # #--------------------------------
    # import matplotlib.pyplot as plt
    # from matplotlib.colors import LogNorm
    # fig = plt.figure(figsize=(10, 8.))
    # ax = fig.add_subplot(111)
    # plt.xlabel(r'r', fontsize = 17)
    # plt.ylabel(r'z', fontsize = 17, labelpad=-7.4)
    # zz, rr = np.meshgrid(zchem, rchem) #inverse r,z because dim is (len(r), len(z))
    # t = plt.pcolormesh(rr/autocm, zz/autocm, dens_naut_d[0, :,:], cmap='gnuplot2', shading='auto', norm=LogNorm(vmin=1e-35, vmax=1e-15), rasterized=True)
    # #t = plt.contourf(rr/autocm, zz/autocm, dens_naut[0, :,:], levels=[0.1,1,8,10,20,30,40,50,60, 70, 80], cmap='jet', rasterized=True)
    # clr = plt.colorbar(t)
    # plt.show()
    # # #-----------------------------------

    return dens_naut_nd, dens_naut_nH   


def dust_temperature_disk(temp_radmc3d, rchem, zchem, d, theta, hg=None):
    nbspecies, nz, ny, nx = temp_radmc3d.shape

    hhg, zz = np.meshgrid(hg, zchem, indexing='ij')
    zz = hhg*zz
    temp_naut = np.ones((nbspecies, len(rchem), len(zchem)))
    temp_naut_smooth = np.ones((nbspecies, len(rchem), len(zchem)))
    for size_id in range(nbspecies):
        for idx, r in enumerate(rchem):
            for alt in range(len(zchem)):
                d_pt = np.sqrt(r**2 + zz[idx, alt]**2)  #convert from cartesian to spherical
                theta_pt = np.arccos(zz[idx, alt]/d_pt) #convert from cartesian to spherical
                closest_d = min(enumerate(d), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
                closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
                #temp_naut[size_id, idx, alt] = temp_radmc3d[size_id, closest_d[0], closest_t[0], 0] 
                temp_naut[size_id, idx, alt] = temp_radmc3d[size_id, 0, closest_t[0], closest_d[0]] 

    #SMOOTHING TEMPERATURE PROFILE
    for idx in range(len(rchem)):
        for size_id in range(nbspecies):
            temp_naut_smooth[size_id, idx, :] = moving_average(temp_naut[size_id, idx, :], 5) #average the values over a rolling window of 5 points.
    temp_naut_smooth[:, :, 0:2] = temp_naut[:, :, 0:2]  #clean the boundary effects
    temp_naut_smooth[:, :, -2:] = temp_naut[:, :, -2:] 
    
    return temp_naut_smooth 

def dust_temperature_single_disk(temp_radmc3d, rchem, zchem, d, theta, hg=None):
    nz, ny, nx = temp_radmc3d.shape

    hhg, zz = np.meshgrid(hg, zchem, indexing='ij')
    zz = hhg*zz
    temp_naut = np.ones((len(rchem), len(zchem)))
    temp_naut_smooth = np.ones((len(rchem), len(zchem)))

    for idx, r in enumerate(rchem):
        for alt in range(len(zchem)):
            d_pt = np.sqrt(r**2 + zz[idx, alt]**2)  #convert from cartesian to spherical
            theta_pt = np.arccos(zz[idx, alt]/d_pt) #convert from cartesian to spherical
            closest_d = min(enumerate(d), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
            closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
            temp_naut[idx, alt] = temp_radmc3d[0, closest_t[0], closest_d[0]] 

    #SMOOTHING TEMPERATURE PROFILE
    for idx in range(len(rchem)):
        temp_naut_smooth[idx, :] = moving_average(temp_naut[idx, :], 5) #average the values over a rolling window of 5 points.
    temp_naut_smooth[:, 0:2] = temp_naut[:, 0:2]  #clean the boundary effects
    temp_naut_smooth[:, -2:] = temp_naut[:, -2:]    

    return temp_naut_smooth   


def dust_temperature(temp_radmc3d, rchem, zchem, d, theta):
    nbspecies, nz, ny, nx = temp_radmc3d.shape

    temp_naut = np.ones((nbspecies, len(rchem), len(zchem[0, :])))
    temp_naut_smooth = np.ones((nbspecies, len(rchem), len(zchem[0, :])))

    for size_id in range(nbspecies):
        for idx, r in enumerate(rchem):
            for idz, z in enumerate(zchem[idx, :]):
                d_pt = np.sqrt(r**2 + z**2)  #convert from cartesian to spherical
                theta_pt = np.arccos(z/d_pt) #convert from cartesian to spherical
                closest_d = min(enumerate(d), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
                closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
                temp_naut[size_id, idx, idz] = temp_radmc3d[size_id, 0, closest_t[0], closest_d[0]] 

    #SMOOTHING TEMPERATURE PROFILE
    for idx in range(len(rchem)):
        for size_id in range(nbspecies):
            temp_naut_smooth[size_id, idx, :] = moving_average(temp_naut[size_id, idx, :], 5) #average the values over a rolling window of 5 points.
    temp_naut_smooth[:, :, 0:2] = temp_naut[:, :, 0:2]  #clean the boundary effects
    temp_naut_smooth[:, :, -2:] = temp_naut[:, :, -2:] 
    

    # # #--------------------------------
    # import matplotlib.pyplot as plt
    # from matplotlib.colors import LogNorm
    # fig = plt.figure(figsize=(10, 8.))
    # ax = fig.add_subplot(111)
    # plt.xlabel(r'r', fontsize = 17)
    # plt.ylabel(r'z', fontsize = 17, labelpad=-7.4)
    # zz, rr = np.meshgrid(zchem, rchem) #inverse r,z because dim is (len(r), len(z))
    # t = plt.pcolormesh(rr/autocm, zz/autocm, temp_naut_smooth[0,:,:], cmap='gnuplot2', shading='auto', vmin=5, vmax=80, rasterized=True)
    # #t = plt.contourf(rr/autocm, zz/autocm, dens_naut[0, :,:], levels=[0.1,1,8,10,20,30,40,50,60, 70, 80], cmap='jet', rasterized=True)
    # clr = plt.colorbar(t)
    # plt.show()
    # # #-----------------------------------   


    return temp_naut_smooth 

def dust_temperature_single(temp_radmc3d, rchem, zchem, d, theta):
    nz, ny, nx = temp_radmc3d.shape

    temp_naut = np.ones((len(rchem), len(zchem[0,:])))
    temp_naut_smooth = np.ones((len(rchem), len(zchem[0,:])))

    for idx, r in enumerate(rchem):
        for idz, z in enumerate(zchem[idx, :]):
            d_pt = np.sqrt(r**2 + z**2)  #convert from cartesian to spherical
            theta_pt = np.arccos(z/d_pt) #convert from cartesian to spherical
            closest_d = min(enumerate(d), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
            closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
            temp_naut[idx, idz] = temp_radmc3d[0, closest_t[0], closest_d[0]] 

    #SMOOTHING TEMPERATURE PROFILE
    for idx in range(len(rchem)):
        temp_naut_smooth[idx, :] = moving_average(temp_naut[idx, :], 5) #average the values over a rolling window of 5 points.
    temp_naut_smooth[:, 0:2] = temp_naut[:, 0:2]  #clean the boundary effects
    temp_naut_smooth[:, -2:] = temp_naut[:, -2:]    


    # # #--------------------------------
    # import matplotlib.pyplot as plt
    # from matplotlib.colors import LogNorm
    # fig = plt.figure(figsize=(10, 8.))
    # ax = fig.add_subplot(111)
    # plt.xlabel(r'r', fontsize = 17)
    # plt.ylabel(r'z', fontsize = 17, labelpad=-7.4)
    # zz, rr = np.meshgrid(zchem, rchem) #inverse r,z because dim is (len(r), len(z))
    # t = plt.pcolormesh(rr/autocm, zz/autocm, temp_naut_smooth[:,:], cmap='gnuplot2', shading='gouraud', vmin=5, vmax=80, rasterized=True)
    # #t = plt.contourf(rr/autocm, zz/autocm, dens_naut[0, :,:], levels=[0.1,1,8,10,20,30,40,50,60, 70, 80], cmap='jet', rasterized=True)
    # clr = plt.colorbar(t)
    # plt.show()
    # # #-----------------------------------   

    return temp_naut_smooth   


def local_field():
    pass

def avz_disk(field_radmc3d, lam_mono, R_star, T_star, rchem, zchem, d, theta, hg):
    nlam, nph, nt, nr = field_radmc3d.shape
    lamuv = np.where((lam_mono <= 0.2)) # extract the ~ uv
    if len(lamuv[0]) == len(lam_mono):
        extrawave = lam_mono[lamuv[0][-1]] - lam_mono[lamuv[0][-2]]
        extrawave += lam_mono[lamuv[0][-1]]
        lam_mono = np.append(lam_mono, extrawave)
    freq = c/(lam_mono*1e-6)
    fieldint = np.zeros((nt, nr))
    bbint = 0 

    # Integrate over uv frequencies:
    for i in lamuv[0]:
        fieldint += field_radmc3d[i, 0, :, :]*abs(freq[i+1]-freq[i])
        bbint += black_body(T_star, lam_mono[i])*abs(freq[i+1]-freq[i]) #integrate BB spectru over the chosen wavelength range
    fieldint[fieldint==0.0] = 1e-30 #provide arbitrary low values in order to avoid division by zero.

    # Convert from spherical to nautilus grid:
    hhg, zz = np.meshgrid(hg, zchem, indexing='ij')
    zz = hhg*zz
    field_naut, field_naut_smooth = np.ones((len(rchem), len(zchem))), np.ones((len(rchem), len(zchem)))
    avz = np.ones((len(rchem), len(zchem)))
    bbint_map = np.ones((len(rchem), len(zchem)))

    #CREATE UNATTENUATED MAP USING A BB SPECTRUM
    for idx, r in enumerate(rchem):
        for idz in range(len(zchem)):
            bbint_map[idx, idz] = bbint*R_star**2*(1/(r**2+zz[idx, idz]**2))/np.pi

    for idx, r in enumerate(rchem):
        for z in range(len(zchem)):
            d_pt = np.sqrt(r**2 + zz[idx, z]**2)  #convert from cartesian to spherical
            theta_pt = np.arccos(zz[idx, z]/d_pt) #convert from cartesian to spherical
            closest_d = min(enumerate(d), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
            closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
            field_naut[idx, z] = fieldint[closest_t[0], closest_d[0]]
    # Smoothing vertical profiles
    for idx in range(len(rchem)):
        field_naut_smooth[idx, :] = moving_average(field_naut[idx, :], 5) #average the values over a rolling window of 5 points.
    field_naut_smooth[:, 0:2] = field_naut[:, 0:2]  #clean the boundary effects
    field_naut_smooth[:, -2:] = field_naut[:, -2:]   

    #CREATE Av MAP
    for idx in range(len(rchem)):
        for idz in range(len(zchem)):
            avz[idx, idz] = abs(-1.086*np.log(field_naut[idx,idz]/field_naut[idx, idz]))   #field0[idx]))
            avz[idx,1:][avz[idx,1:]==0.0] = np.trim_zeros(avz[idx])[0]
    
    # avz_df = pd.DataFrame(data=avz.transpose())
    # avz_df = avz_df.rolling(window=5, center=True, min_periods=2).mean()
    return avz



def av_z(field_radmc3d, lam_mono, R_star, T_star, rchem, zchem, d, theta):
    nlam, nph, nt, nr = field_radmc3d.shape
    lamuv = np.where((lam_mono <= 0.2)) # extract the ~ uv
    if len(lamuv[0]) == len(lam_mono):
        extrawave = lam_mono[lamuv[0][-1]] - lam_mono[lamuv[0][-2]]
        extrawave += lam_mono[lamuv[0][-1]]
        lam_mono = np.append(lam_mono, extrawave)
    freq = c/(lam_mono*1e-6)
    fieldint = np.zeros((nt, nr))
    bbint = 0 

    # Integrate over uv frequencies:
    for i in lamuv[0]:
          fieldint += field_radmc3d[i, 0, :, :]*abs(freq[i+1]-freq[i])
          bbint += black_body(T_star, lam_mono[i])*abs(freq[i+1]-freq[i]) #integrate BB spectru over the chosen wavelength range
    fieldint[fieldint==0.0] = 1e-30 #provide arbitrary low values in order to avoid division by zero.

    ## Convert from spherical to nautilus grid:
    field_naut, field_naut_smooth = np.ones((len(rchem), len(zchem[0,:]))), np.ones((len(rchem), len(zchem[0,:])))
    avz, avz_smooth = np.ones((len(rchem), len(zchem[0,:]))), np.ones((len(rchem), len(zchem[0,:])))
    bbint_map = np.ones((len(rchem), len(zchem[0,:])))

    #CREATE UNATTENUATED MAP USING A BB SPECTRUM
    for idx, r in enumerate(rchem):
        for idz, z in enumerate(zchem[idx, :]):
            bbint_map[idx, idz] = bbint*R_star**2*(1/(r**2+z**2))/np.pi


    for idx, r in enumerate(rchem):
        for idz, z in enumerate(zchem[idx, :]):
            d_pt = np.sqrt(r**2 + z**2)  #convert from cartesian to spherical
            theta_pt = np.arccos(z/d_pt) #convert from cartesian to spherical
            closest_d = min(enumerate(d), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
            closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
            field_naut[idx, idz] = fieldint[closest_t[0], closest_d[0]]
    # Smoothing vertical profiles
    for idx in range(len(rchem)):
        field_naut_smooth[idx, :] = moving_average(field_naut[idx, :], 5) #average the values over a rolling window of 5 points.
    field_naut_smooth[:, 0:2] = field_naut[:, 0:2]  #clean the boundary effects
    field_naut_smooth[:, -2:] = field_naut[:, -2:] 

    # #--------------------------------
    #import matplotlib.pyplot as plt
    #from matplotlib.colors import LogNorm
    #fig = plt.figure(figsize=(10, 8.))
    #ax = fig.add_subplot(111)
    #plt.xlabel(r'r', fontsize = 17)
    #plt.ylabel(r'z', fontsize = 17, labelpad=-7.4)
    #rr, zz = np.meshgrid(rchem, zchem)
    #t = plt.pcolormesh(rr/autocm, zz/autocm, bbint_map, cmap='gnuplot2', shading='gouraud', norm=LogNorm(vmin=np.min(bbint_map), vmax=np.max(bbint_map)), rasterized=True)
    #clr = plt.colorbar(t)
    #plt.show()
    # #-----------------------------------


    # #--------------------------------
    #import matplotlib.pyplot as plt
    #from matplotlib.colors import LogNorm
    #fig = plt.figure(figsize=(10, 8.))
    #ax = fig.add_subplot(111)
    #plt.xlabel(r'r', fontsize = 17)
    #plt.ylabel(r'z', fontsize = 17, labelpad=-7.4)
    #rr, zz = np.meshgrid(rchem, zchem)
    #t = plt.pcolormesh(rr/autocm, zz/autocm, field_naut, cmap='gnuplot2', shading='gouraud', norm=LogNorm(vmin=np.min(field_naut), vmax=np.max(field_naut)), rasterized=True)
    #clr = plt.colorbar(t)
    #plt.show()
    # #-----------------------------------


    #CREATE Av MAP
    for idx in range(len(rchem)):
        for idz in range(len(zchem[idx,:])):
            avz[idx, idz] = abs(-1.086*np.log(field_naut[idx,idz]/bbint_map[idx, idz]))   #field0[idx]))
            avz[idx, 1:][avz[idx, 1:]==0.0] = np.trim_zeros(avz[idx])[0]

    ##SMOOTHING
    for idx in range(len(rchem)):
        avz_smooth[idx, :] = moving_average(avz[idx, :], 5) #average the values over a rolling window of 5 points.
    for idz in range(len(zchem[0,:])):
        avz_smooth[:, idz] = moving_average(avz_smooth[:, idz], 5) #average the values over a rolling window of 5 points.
    avz_smooth[:, 0:2] = avz[:, 0:2]  #clean the boundary effects
    avz_smooth[:, -2:] = avz[:, -2:]  
    avz_smooth[0:2, :] = avz[0:2, :]  #clean the boundary effects
    avz_smooth[-2:, :] = avz[-2:, :]
    #avz_smooth = np.where(avz_smooth<1, avz_smooth*100, avz_smooth)  


    # # #--------------------------------
    # import matplotlib.pyplot as plt
    # from matplotlib.colors import LogNorm
    # fig = plt.figure(figsize=(10, 8.))
    # ax = fig.add_subplot(111)
    # plt.xlabel(r'r', fontsize = 17)
    # plt.ylabel(r'z', fontsize = 17, labelpad=-7.4)
    # zz, rr = np.meshgrid(zchem, rchem) #inverse r,z because dim is (len(r), len(z))
    # #t = plt.pcolormesh(rr/autocm, zz/autocm, avz, cmap='gnuplot2', shading='gouraud', norm=LogNorm(vmin=1e0, vmax=1e2), rasterized=True)
    # t = plt.contourf(rr/autocm, zz/autocm, avz, levels=[0.1,1,8,10,20,30,40,50,60, 70, 80], cmap='jet')
    # clr = plt.colorbar(t)
    # plt.show()
    # # #-----------------------------------

    return avz_smooth


def to_spherical(chemmodel, nr, nt, dist, theta, struct='numberdens_species'): 
    #self.grid.chemmodel[species], nx, ny, x, y
    spherical_struct = np.zeros((nr, nt))
    r_naut = np.array(list(chemmodel.keys()))
    rcut = r_naut[0]
    for id_thet, thet in enumerate(theta):
        for id_d, d in enumerate(dist):
            r_sph = d*np.sin(thet)
            z_sph = abs(d*np.cos(thet))
            closest_r = min(enumerate(r_naut), key=lambda x: abs(x[1]-r_sph)) #find closest grid point
            closest_z = min(enumerate(chemmodel[closest_r[1]]['z']), key=lambda x: abs(x[1]-z_sph)) #find closest grid point
            if d < rcut:
                spherical_struct[id_d, id_thet] = 0.0
            else:
                spherical_struct[id_d, id_thet] = chemmodel[closest_r[1]][struct][closest_z[0]]

    return spherical_struct

