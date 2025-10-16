import numpy as np
from astropy import units as u
from astropy.modeling.models import BlackBody

#-------------constants-------------   <-- Souldn't be changed
pi = 3.141592653589793
c = 299792458 #m.s-1
h = 6.62607004E-34 #J.s
h_erg = h*1e7 #erg.s
hbar = 1.054571800E-34   #1.054571800E−34 #J.s
kb = 1.38054e-16 #Boltzmann constant in erg.K-1
k = 1.381e-23 #Boltzmann constant in J.K-1
autocm = 1.49597e13 #convert AU to cm
autom = 1.49597e11 #convert AU to m
mu=2.4 #mean molecular of the gas equal to Solar metallicity, can be varied from 2 to 3 depending on the metallicity of the star.
Ggram = 6.6743e-8 #cgs
e = 1.602e-19 #Coulomb

###ASTRONOMICAL CONSTANTS###
M_sun = 1.98847e33 #g
L_sun = 3.9e33    # erg s^-1
R_sun = 6.96e10   # cm
T_sun = 5.78e3    # K
M_earth = 5.9722e27 #g
au = 1.496e11 #m
parsec_in_ly = 3.262 #light-years
parsec_in_m = 3.086e16 #m

###MASS###
#mol & Avogadro
Na = 6.022140857e23 #mol−1
u = 1.660538921e-27 #kg which is 1g/Na mol
amu = 1.66053892e-24   #g (as "u" but in "g")

#nucleous & electron
e_mass = 9.109e-31 #mass of electron in kg
proton_mass_kg = 1.673e-27 #mass of proton in kg
proton_mass_g = 1.673e-24 #mass of proton in g

#atomic mass (kg)
H_mass = 1.6737236e-27
He_mass = 6.6464764e-27
Li_mass = 1.1525801e-26
#Be_mass = 
#B_mass = 
C_mass = 1.9944235e-26
N_mass = 2.3258671e-26
O_mass = 2.6566962e-26
#F_mass = 
#Ne_mass = 
#Na_mass = 
#Mg_mass =
#Al_mass =
Si_mass = 4.6637066e-26
#P_mass = 
#S_mass = 
#Cl_mass = 
#Ar_mass =  

##PLANCK'S LAW (erg.s-1.cm-2,Hz-1.ster-1)
def black_body(T, lam):
    from astropy import units as u
    bb = BlackBody(temperature=T*u.K)
    #bb = ((2*h_erg*freq**3)/c**2)*(1/(np.exp((h_erg*freq)/(kb*T)) - 1))*1e-4

    return bb(lam*u.micron).value
    #return bb