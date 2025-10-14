################################################################################################################
#                                                     radmc3d                                                  #
################################################################################################################
lmin = 0.1
lmax = 1e5
################################################################################################################
#                                                     nautilus                                                 #
################################################################################################################
rchem = [60,80,100]                         # radial distance where chemistry will be computed
max_H = 4                                   # maximum altitude in fraction of scale height e.g. n=4 means that the largest computed z is 4*Hg
nz_chem = 64                                # number of vertical spatial points for chemistry
################################################################################################################
#                                             disk physical structure                                          #
################################################################################################################
disk_mass = 0.15                            # total mass of the disk in Solar mass
disk_dust_mass = 1.7e-3                     # Total dust mass in the disk in solar mass
ref_radius = 1.000e+02                      # reference radius for parametric laws
cut_radius = 5.000e+02                      # tapered edge radius [au]
rin = 1.000e+00                             # inner radius [au]
rout = 3.000e+02                            # outer radius [au]
Tmidplan_ref = 1.000e+01                    # mid-plan temperature at the reference radius
Tatmos_ref = 5.000e+01                      # atmospheric temperature at the reference radius (where z=n*H)
q_exp = 4.000e-01                           # exponent for the radial variation of temperature and gas scale height
sigma_t = 2.000e+00                         # stiffness of the vertical temperature profile
################################################################################################################
#                                           envelope physical structure                                        #
################################################################################################################
rmin = 3.000e+1                             # outer radius of the envelope [au]
rmax = 1.000e+2                             # inner radius [au]
nr = 4
ntheta = 65
nphi =  2
r_centri = 1.000e+2                         # centrifugal radius (critical radius inside of which the envelope flattens) [au]
acc_rate = 6e-6
dust_env_mass = 1.000e-5                    # total dust mass of the envelope in Solar mass
cavpl = 1.100e0                             # opening angle and shape of the outflow
cav_fact = 2.000e-2                         # factor of decrease of the density in the cavity
cavz0 = 1.000e1                             # [au]

################################################################################################################
#                                                  DISK GAS                                                    #
################################################################################################################
sigma_gas_ref = 3.350e-01             # surface density of the gas at reference radius [g.cm-2]
h0 = 8.21                                   # scale height at reference [au]
p_exp = 1.500e+00                           # surface density exponent
nH_to_AV_conversion = 1.600e+21             # conversion factor of H colmun density to Av (Wagenblast \& Hartquist 1989)
lim_h = 1                                   # value in fraction of scale height until which accretion heating is computed 
################################################################################################################
#                                                    DUST                                                      #
################################################################################################################
d_exp = 3.500e0                             # size distribution exponent. WARNING: don't type 4.
rho_m = 2.500e+00	                        # material density of the grains [g.cm-3]
dtogas = 1.000e-02                          # dust to gas mass ratio in the upper part (usefull only if one grain size!)
schmidtnumber = 1.000e+00	                # Schmidt number. See documentation.  
alpha = 1.000e-02	                        # viscosity coefficient.
settfact = 9.000e-01	                    # settling factor. See documentation. 
ext_eff = 4.000e+00	                        # extinction efficiency at resonance.
nb_sizes = 16                               # number of grain sizes
amin = 5.000e-03                            # minimal grain size (not read if nb_sizes = 1). [microns]   
amax = 1.000e+03                            # max grain size (not read if nb_sizes = 1). [microns]
asingle = 1.000e-05                         # size of small grains if nb_sizes = 1 [microns]
cst_norm = 7.41e-26                         # normalization constant (7.41e-26 in the case of MRN distribution as in MRN1977)
cutoff = 0.000e+00                          # whether yes or no you want a cutoff at cut radius. 1="cutoff", 0="no cutoff"
q_c = 6                                     # Extinction efficiency at resonance
################################################################################################################
#                                                      STAR                                                    #
################################################################################################################
star_mass = 0.579863                        # mass of the central star in Solar mass
uv_ref = 3.400e+03                          # total UV field from the central star at the reference radius

################################################################################################################
#                                           DUST SIZES AND DISTRIBUTION                                        #
################################################################################################################