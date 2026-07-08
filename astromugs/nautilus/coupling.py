import numpy as np
from scipy.interpolate import RegularGridInterpolator
from ..constants.constants import c, autocm, amu, mu, black_body

def moving_average(x, w):
    """Compute a moving average using a uniform convolution kernel.

    Parameters
    ----------
    x : array_like
        Input array to smooth.
    w : int
        Width of the rolling window (number of points).

    Returns
    -------
    numpy.ndarray
        Smoothed array with the same length as `x`.
    """
    return np.convolve(x, np.ones(w), 'same') / w

def dust_density(dtogas, rho_m, sizes, dens_radmc3d, rchem, zchem, d, theta):
    """Interpolate RADMC3D dust density onto the Nautilus chemistry grid.

    Converts dust mass densities from the RADMC3D spherical grid to the
    Nautilus cylindrical grid using bilinear interpolation, preserving
    radial drift gradients and smooth disk boundaries. Returns both dust
    number densities per grain size and gas number density.

    Parameters
    ----------
    dtogas : float
        Dust-to-gas mass ratio.
    rho_m : float
        Material (intrinsic) density of dust grains, in g/cm^3.
    sizes : array_like
        Array of grain radii for each dust species, in cm.
    dens_radmc3d : numpy.ndarray
        Dust mass density from RADMC3D with shape
        ``(nbspecies, nphi, ntheta, nr)``, in g/cm^3.
    rchem : numpy.ndarray
        Radial grid of the Nautilus chemistry model, in cm.
    zchem : numpy.ndarray
        Vertical grid of the Nautilus chemistry model with shape
        ``(nrchem, nzchem)``, in cm. Values can differ per radius.
    d : numpy.ndarray
        Radial (spherical) distance grid of the RADMC3D model, in AU.
    theta : numpy.ndarray
        Co-latitude grid of the RADMC3D model, in radians.

    Returns
    -------
    dens_naut_nd : numpy.ndarray
        Dust number density per grain size, shape
        ``(nbspecies, nrchem, nzchem)``, in cm^-3.
    dens_naut_nH : numpy.ndarray
        Gas hydrogen number density, shape ``(nrchem, nzchem)``, in cm^-3.
    """
    nbspecies, nz, ny, nx = dens_radmc3d.shape
    nrchem = len(rchem)
    nzchem = zchem.shape[1]

    dens_naut = np.zeros((nbspecies, nrchem, nzchem))

    # Convert RADMC3D grid from au to cm to match rchem/zchem units
    d_cm = d * autocm

    # Build query points: convert cylindrical (rchem, zchem) -> spherical (d_pt, theta_pt)
    # zchem has shape (nrchem, nzchem) — z values can differ per radius
    d_pts = np.sqrt(rchem[:, None]**2 + zchem**2)            # (nrchem, nzchem) in cm
    theta_pts = np.arccos(np.clip(zchem / d_pts, -1, 1))     # (nrchem, nzchem)

    # Clamp query points to the RADMC3D grid range to avoid extrapolation
    d_pts = np.clip(d_pts, d_cm[0], d_cm[-1])
    theta_pts = np.clip(theta_pts, theta[0], theta[-1])

    # Stack into (N, 2) array for the interpolator
    query = np.column_stack([d_pts.ravel(), theta_pts.ravel()])

    for idx_size in range(nbspecies):
        # dens_radmc3d shape is (nspecies, nphi, ntheta, nr) — take phi=0
        data_2d = dens_radmc3d[idx_size, 0, :, :]  # (ntheta, nr)
        # RegularGridInterpolator expects data on (d, theta) axes
        # data_2d is indexed as [itheta, ir], so transpose to [ir, itheta]
        interp = RegularGridInterpolator(
            (d_cm, theta), data_2d.T, method='linear', bounds_error=False, fill_value=0.0
        )
        dens_naut[idx_size] = interp(query).reshape(nrchem, nzchem)

    # Gas number density from total dust
    dens_naut_nH = dens_naut.sum(axis=0) / dtogas / (mu * amu)

    # Dust number density per grain size
    dens_naut_nd = np.zeros_like(dens_naut)
    for ai in range(nbspecies):
        mass = (4./3.) * np.pi * rho_m * sizes[ai]**3
        dens_naut_nd[ai] = dens_naut[ai] / mass

    return dens_naut_nd, dens_naut_nH


def dust_temperature_disk(temp_radmc3d, rchem, zchem, d, theta, hg=None):
    """Remap multi-species RADMC3D dust temperatures onto the Nautilus grid using scale heights.

    Converts temperatures from the RADMC3D spherical grid to the Nautilus
    cylindrical grid via bilinear interpolation, then smooths vertical
    profiles with a 5-point moving average.

    Parameters
    ----------
    temp_radmc3d : numpy.ndarray
        Dust temperature from RADMC3D with shape
        ``(nbspecies, nphi, ntheta, nr)``, in K.
    rchem : numpy.ndarray
        Radial grid of the Nautilus chemistry model, in cm.
    zchem : numpy.ndarray
        Normalised vertical grid points (dimensionless, scaled by `hg`).
    d : numpy.ndarray
        Radial (spherical) distance grid of the RADMC3D model, in AU.
    theta : numpy.ndarray
        Co-latitude grid of the RADMC3D model, in radians.
    hg : array_like, optional
        Gas scale height at each radial point, in cm. Used to convert
        normalised `zchem` to physical heights.

    Returns
    -------
    numpy.ndarray
        Smoothed dust temperature on the Nautilus grid with shape
        ``(nbspecies, nrchem, nzchem)``, in K.
    """
    nbspecies, nz, ny, nx = temp_radmc3d.shape
    nrchem = len(rchem)
    nzchem = len(zchem)

    d_cm = d * autocm  # convert RADMC3D grid from au to cm

    # Convert normalised zchem to physical heights using scale heights
    hhg, zz = np.meshgrid(hg, zchem, indexing='ij')
    zz = hhg * zz  # (nrchem, nzchem) in cm

    # Build query points: convert cylindrical (rchem, zz) -> spherical (d_pt, theta_pt)
    d_pts = np.sqrt(rchem[:, None]**2 + zz**2)            # (nrchem, nzchem)
    theta_pts = np.arccos(np.clip(zz / d_pts, -1, 1))     # (nrchem, nzchem)

    # Clamp query points to the RADMC3D grid range
    d_pts = np.clip(d_pts, d_cm[0], d_cm[-1])
    theta_pts = np.clip(theta_pts, theta[0], theta[-1])

    # Stack into (N, 2) array for the interpolator
    query = np.column_stack([d_pts.ravel(), theta_pts.ravel()])

    temp_naut = np.ones((nbspecies, nrchem, nzchem))
    temp_naut_smooth = np.ones((nbspecies, nrchem, nzchem))

    for size_id in range(nbspecies):
        data_2d = temp_radmc3d[size_id, 0, :, :]  # (ntheta, nr)
        interp = RegularGridInterpolator(
            (d_cm, theta), data_2d.T, method='linear', bounds_error=False, fill_value=None
        )
        temp_naut[size_id] = interp(query).reshape(nrchem, nzchem)

    #SMOOTHING TEMPERATURE PROFILE
    for idx in range(nrchem):
        for size_id in range(nbspecies):
            temp_naut_smooth[size_id, idx, :] = moving_average(temp_naut[size_id, idx, :], 5)
    temp_naut_smooth[:, :, 0:2] = temp_naut[:, :, 0:2]  #clean the boundary effects
    temp_naut_smooth[:, :, -2:] = temp_naut[:, :, -2:]

    return temp_naut_smooth

def dust_temperature_single_disk(temp_radmc3d, rchem, zchem, d, theta, hg=None):
    """Remap single-species RADMC3D dust temperature onto the Nautilus grid using scale heights.

    Same as `dust_temperature_disk` but for a single dust species (no
    species dimension). Uses nearest-neighbor lookup followed by 5-point
    moving-average smoothing of vertical profiles.

    Parameters
    ----------
    temp_radmc3d : numpy.ndarray
        Dust temperature from RADMC3D with shape
        ``(nphi, ntheta, nr)``, in K.
    rchem : numpy.ndarray
        Radial grid of the Nautilus chemistry model, in cm.
    zchem : numpy.ndarray
        Normalised vertical grid points (dimensionless, scaled by `hg`).
    d : numpy.ndarray
        Radial (spherical) distance grid of the RADMC3D model, in AU.
    theta : numpy.ndarray
        Co-latitude grid of the RADMC3D model, in radians.
    hg : array_like, optional
        Gas scale height at each radial point, in cm. Used to convert
        normalised `zchem` to physical heights.

    Returns
    -------
    numpy.ndarray
        Smoothed dust temperature on the Nautilus grid with shape
        ``(nrchem, nzchem)``, in K.
    """
    nz, ny, nx = temp_radmc3d.shape

    d_cm = d * autocm  # convert RADMC3D grid from au to cm

    hhg, zz = np.meshgrid(hg, zchem, indexing='ij')
    zz = hhg*zz
    temp_naut = np.ones((len(rchem), len(zchem)))
    temp_naut_smooth = np.ones((len(rchem), len(zchem)))

    for idx, r in enumerate(rchem):
        for alt in range(len(zchem)):
            d_pt = np.sqrt(r**2 + zz[idx, alt]**2)  #convert from cartesian to spherical
            theta_pt = np.arccos(zz[idx, alt]/d_pt) #convert from cartesian to spherical
            closest_d = min(enumerate(d_cm), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
            closest_t = min(enumerate(theta), key=lambda x: abs(x[1]-theta_pt)) #find closest grid point
            temp_naut[idx, alt] = temp_radmc3d[0, closest_t[0], closest_d[0]] 

    #SMOOTHING TEMPERATURE PROFILE
    for idx in range(len(rchem)):
        temp_naut_smooth[idx, :] = moving_average(temp_naut[idx, :], 5) #average the values over a rolling window of 5 points.
    temp_naut_smooth[:, 0:2] = temp_naut[:, 0:2]  #clean the boundary effects
    temp_naut_smooth[:, -2:] = temp_naut[:, -2:]    

    return temp_naut_smooth   


def dust_temperature(temp_radmc3d, rchem, zchem, d, theta):
    """Remap multi-species RADMC3D dust temperatures onto the Nautilus grid.

    Converts temperatures from the RADMC3D spherical grid to the Nautilus
    cylindrical grid via bilinear interpolation, where `zchem` provides
    physical heights per radius. Vertical profiles are smoothed with a
    5-point moving average.

    Parameters
    ----------
    temp_radmc3d : numpy.ndarray
        Dust temperature from RADMC3D with shape
        ``(nbspecies, nphi, ntheta, nr)``, in K.
    rchem : numpy.ndarray
        Radial grid of the Nautilus chemistry model, in cm.
    zchem : numpy.ndarray
        Vertical grid with shape ``(nrchem, nzchem)``, in cm.
        Heights can differ per radial point.
    d : numpy.ndarray
        Radial (spherical) distance grid of the RADMC3D model, in AU.
    theta : numpy.ndarray
        Co-latitude grid of the RADMC3D model, in radians.

    Returns
    -------
    numpy.ndarray
        Smoothed dust temperature on the Nautilus grid with shape
        ``(nbspecies, nrchem, nzchem)``, in K.
    """
    nbspecies, nz, ny, nx = temp_radmc3d.shape
    nrchem = len(rchem)
    nzchem = zchem.shape[1]

    d_cm = d * autocm  # convert RADMC3D grid from au to cm

    # Build query points: convert cylindrical (rchem, zchem) -> spherical (d_pt, theta_pt)
    d_pts = np.sqrt(rchem[:, None]**2 + zchem**2)            # (nrchem, nzchem)
    theta_pts = np.arccos(np.clip(zchem / d_pts, -1, 1))     # (nrchem, nzchem)

    # Clamp query points to the RADMC3D grid range
    d_pts = np.clip(d_pts, d_cm[0], d_cm[-1])
    theta_pts = np.clip(theta_pts, theta[0], theta[-1])

    # Stack into (N, 2) array for the interpolator
    query = np.column_stack([d_pts.ravel(), theta_pts.ravel()])

    temp_naut = np.ones((nbspecies, nrchem, nzchem))
    temp_naut_smooth = np.ones((nbspecies, nrchem, nzchem))

    for size_id in range(nbspecies):
        data_2d = temp_radmc3d[size_id, 0, :, :]  # (ntheta, nr)
        interp = RegularGridInterpolator(
            (d_cm, theta), data_2d.T, method='linear', bounds_error=False, fill_value=None
        )
        temp_naut[size_id] = interp(query).reshape(nrchem, nzchem)

    #SMOOTHING TEMPERATURE PROFILE
    for idx in range(nrchem):
        for size_id in range(nbspecies):
            temp_naut_smooth[size_id, idx, :] = moving_average(temp_naut[size_id, idx, :], 5)
    temp_naut_smooth[:, :, 0:2] = temp_naut[:, :, 0:2]  #clean the boundary effects
    temp_naut_smooth[:, :, -2:] = temp_naut[:, :, -2:]

    return temp_naut_smooth

def dust_temperature_single(temp_radmc3d, rchem, zchem, d, theta):
    """Remap single-species RADMC3D dust temperature onto the Nautilus grid.

    Same as `dust_temperature` but for a single dust species (no species
    dimension). Uses nearest-neighbor lookup followed by 5-point
    moving-average smoothing of vertical profiles.

    Parameters
    ----------
    temp_radmc3d : numpy.ndarray
        Dust temperature from RADMC3D with shape
        ``(nphi, ntheta, nr)``, in K.
    rchem : numpy.ndarray
        Radial grid of the Nautilus chemistry model, in cm.
    zchem : numpy.ndarray
        Vertical grid with shape ``(nrchem, nzchem)``, in cm.
        Heights can differ per radial point.
    d : numpy.ndarray
        Radial (spherical) distance grid of the RADMC3D model, in AU.
    theta : numpy.ndarray
        Co-latitude grid of the RADMC3D model, in radians.

    Returns
    -------
    numpy.ndarray
        Smoothed dust temperature on the Nautilus grid with shape
        ``(nrchem, nzchem)``, in K.
    """
    nz, ny, nx = temp_radmc3d.shape

    d_cm = d * autocm  # convert RADMC3D grid from au to cm

    temp_naut = np.ones((len(rchem), len(zchem[0,:])))
    temp_naut_smooth = np.ones((len(rchem), len(zchem[0,:])))

    for idx, r in enumerate(rchem):
        for idz, z in enumerate(zchem[idx, :]):
            d_pt = np.sqrt(r**2 + z**2)  #convert from cartesian to spherical
            theta_pt = np.arccos(z/d_pt) #convert from cartesian to spherical
            closest_d = min(enumerate(d_cm), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
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
    """Compute the local radiation field (placeholder, not yet implemented)."""
    pass

def avz_disk(field_radmc3d, lam_mono, R_star, T_star, rchem, zchem, d, theta, hg):
    """Compute the vertical visual extinction map using scale-height-based vertical coordinates.

    Integrates the RADMC3D radiation field over UV wavelengths, remaps it
    onto the Nautilus cylindrical grid (with heights derived from gas scale
    heights), and computes Av by comparing the attenuated field to an
    unattenuated blackbody reference.

    Parameters
    ----------
    field_radmc3d : numpy.ndarray
        Monochromatic mean intensity from RADMC3D with shape
        ``(nlam, nphi, ntheta, nr)``, in cgs units.
    lam_mono : numpy.ndarray
        Wavelength grid of the radiation field, in microns.
    R_star : float
        Stellar radius, in cm.
    T_star : float
        Stellar effective temperature, in K.
    rchem : numpy.ndarray
        Radial grid of the Nautilus chemistry model, in cm.
    zchem : numpy.ndarray
        Normalised vertical grid points (dimensionless, scaled by `hg`).
    d : numpy.ndarray
        Radial (spherical) distance grid of the RADMC3D model, in AU.
    theta : numpy.ndarray
        Co-latitude grid of the RADMC3D model, in radians.
    hg : array_like
        Gas scale height at each radial point, in cm.

    Returns
    -------
    numpy.ndarray
        Visual extinction Av on the Nautilus grid with shape
        ``(nrchem, nzchem)``, in magnitudes.
    """
    nlam, nph, nt, nr = field_radmc3d.shape
    d_cm = d * autocm  # convert RADMC3D grid from au to cm
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
            closest_d = min(enumerate(d_cm), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
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
            avz[idx, idz] = abs(-1.086*np.log(field_naut[idx,idz]/bbint_map[idx, idz]))   #field0[idx]))
            avz[idx,1:][avz[idx,1:]==0.0] = np.trim_zeros(avz[idx])[0]
    
    # avz_df = pd.DataFrame(data=avz.transpose())
    # avz_df = avz_df.rolling(window=5, center=True, min_periods=2).mean()
    return avz



def av_z(field_radmc3d, lam_mono, R_star, T_star, rchem, zchem, d, theta):
    """Compute the vertical visual extinction map using physical vertical coordinates.

    Integrates the RADMC3D radiation field over UV wavelengths, remaps it
    onto the Nautilus cylindrical grid (where `zchem` contains physical
    heights per radius), and derives Av by comparing the attenuated field
    to an unattenuated blackbody reference. Both vertical and radial
    profiles are smoothed with a 5-point moving average.

    Parameters
    ----------
    field_radmc3d : numpy.ndarray
        Monochromatic mean intensity from RADMC3D with shape
        ``(nlam, nphi, ntheta, nr)``, in cgs units.
    lam_mono : numpy.ndarray
        Wavelength grid of the radiation field, in microns.
    R_star : float
        Stellar radius, in cm.
    T_star : float
        Stellar effective temperature, in K.
    rchem : numpy.ndarray
        Radial grid of the Nautilus chemistry model, in cm.
    zchem : numpy.ndarray
        Vertical grid with shape ``(nrchem, nzchem)``, in cm.
        Heights can differ per radial point.
    d : numpy.ndarray
        Radial (spherical) distance grid of the RADMC3D model, in AU.
    theta : numpy.ndarray
        Co-latitude grid of the RADMC3D model, in radians.

    Returns
    -------
    numpy.ndarray
        Smoothed visual extinction Av on the Nautilus grid with shape
        ``(nrchem, nzchem)``, in magnitudes.
    """
    nlam, nph, nt, nr = field_radmc3d.shape
    d_cm = d * autocm  # convert RADMC3D grid from au to cm
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
            closest_d = min(enumerate(d_cm), key=lambda x: abs(x[1]-d_pt)) #find closest grid point
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

    return avz_smooth


def to_spherical(chemmodel, nr, nt, dist, theta, struct='numberdens_species'):
    """Convert a Nautilus chemistry structure from cylindrical to spherical coordinates.

    Maps a quantity stored on the Nautilus cylindrical grid (keyed by
    radial position) onto a regular spherical ``(r, theta)`` grid via
    bilinear interpolation. Each chemistry column is first interpolated
    onto a common z grid (built from the union of all column z arrays),
    producing a regular 2-D field in ``(R_cyl, z)`` that is then passed
    to :class:`~scipy.interpolate.RegularGridInterpolator`. Points outside
    the chemistry domain — above the disk surface, beyond the outermost
    radius, or inside the inner radius — receive floor values.

    Parameters
    ----------
    chemmodel : dict
        Dictionary keyed by cylindrical radius **in AU** (int or float).
        Each value is a dict containing:

        - ``'z'``: 1-D array of vertical heights **in AU**, stored
          surface → midplane (descending order).
        - the field named by `struct`: 1-D array of the quantity to remap,
          same length and order as ``'z'``.

    nr : int
        Number of radial cell centres in the output spherical grid.
    nt : int
        Number of co-latitude cell centres in the output spherical grid.
    dist : numpy.ndarray
        Spherical radial distance grid (cell centres) **in AU**, shape
        ``(nr,)``.
    theta : numpy.ndarray
        Co-latitude grid (cell centres) **in radians**, shape ``(nt,)``.
    struct : str, optional
        Key in each ``chemmodel`` entry for the quantity to remap.
        Default is ``'numberdens_species'``.

    Returns
    -------
    numpy.ndarray
        Remapped quantity on the spherical grid, shape ``(nr, nt)``.
        Cells inside the inner radius are set to ``1e-20``; cells outside
        the modelled domain are set to ``0.0``.
    """
    # ── 1. Sorted Nautilus radii (AU) ────────────────────────────────────────
    r_naut = np.array(sorted(chemmodel.keys()), dtype=float)
    r_inner = r_naut[0]

    # ── 2. Build a common z grid from the union of all column z arrays ───────
    # Each column stores z surface→midplane (descending); flip to ascending.
    z_common = np.unique(np.concatenate([chemmodel[r]['z'] for r in r_naut]))
    # z_common is already ascending after np.unique (returns sorted values).

    # ── 3. Interpolate every column onto z_common ─────────────────────────────
    # np.interp requires xp to be ascending, so we flip z (descending→ascending)
    # and the corresponding values.
    #   left  = midplane value  (z < lowest z in column → below midplane)
    #   right = 0.0             (z > highest z in column → above truncation)
    data_2d = np.zeros((len(r_naut), len(z_common)))
    for ir, r in enumerate(r_naut):
        z_asc   = np.asarray(chemmodel[r]['z']).ravel()[::-1]      # ascending
        val_asc = np.asarray(chemmodel[r][struct]).ravel()[::-1]   # corresponding values
        if len(z_asc) != len(val_asc):
            raise ValueError(
                f"to_spherical: length mismatch at r={r} AU — "
                f"z has {len(z_asc)} points but '{struct}' has {len(val_asc)} points. "
                f"Check that 1D_static.dat and abundances.out are consistent."
            )
        data_2d[ir] = np.interp(z_common, z_asc, val_asc,
                                left=val_asc[0], right=0.0)

    # ── 4. Build the 2-D interpolator on (R_cyl, z) ──────────────────────────
    interp = RegularGridInterpolator(
        (r_naut, z_common), data_2d,
        method='linear', bounds_error=False, fill_value=0.0
    )

    # ── 5. Evaluate at all spherical cell centres (vectorised) ───────────────
    dist_g, theta_g = np.meshgrid(dist, theta, indexing='ij')   # (nr, nt)
    R_cyl = dist_g * np.sin(theta_g)                            # (nr, nt) AU
    z_cyl = np.abs(dist_g * np.cos(theta_g))                    # (nr, nt) AU

    query = np.stack([R_cyl.ravel(), z_cyl.ravel()], axis=1)
    result = interp(query).reshape(nr, nt)

    # Inner-radius cavity: set a small floor instead of 0 to avoid log-scale issues
    result[dist_g < r_inner] = 1e-20

    return result

