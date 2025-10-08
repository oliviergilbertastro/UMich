import numpy as np
from scipy.integrate import quad

h = 6.626*10**(-27) # [erg.s]
k = 1.380649*10**(-16) # [erg/K]
c = 2.99792458*10**(10) # [cm/s]

def B_E(E:float, kT:float) -> float:
   """Planck function as a function of energy."""
   nu = E/(4.1357E-18) # Convert the energy in keV to frequency in Hertz
   T = kT/(8.617333262E-8) # Get T in kelvin using k in keV/K
   dnu_dE = 1/(4.1357E-18)
   return 2*h*nu**3/c**2 * 1/(np.exp(h*nu/(k*T))-1)*dnu_dE # kT/1000 because it's originally in keV

def get_theoretical_F_03_to_10_keV(kT:float) -> float:
   integrand = lambda E: B_E(E, kT)
   F_3_to_10 = quad(integrand, 0.3, 10)[0]
   return F_3_to_10

def get_infos(kT:float, flux_03_to_10:float, distance:float=3000) -> None:
   print(f"Temperature: {kT/(8.617333262E-8):.3f} K")
   F_3_to_10_theo = get_theoretical_F_03_to_10_keV(kT)
   R_src = distance*(3.086E+13)*np.sqrt(flux_03_to_10/F_3_to_10_theo) # R in km
   print(f"Radius: {R_src:.3f} km")
   print(f"Emitting area: {4*np.pi*R_src**2:.3f} km^2")

src1_kT = 0.760297 # in keV
src1_norm = 1.61602E-07
src1_flux = 1.3491e-14 # ergs/cm^2/s from 0.3-10 keV
get_infos(src1_kT, src1_flux)


src2_kT = 0.741912 # in keV
src2_norm = 1.00526E-06
src2_flux = 8.3923e-14 # ergs/cm^2/s from 0.3-10 keV
get_infos(src2_kT, src2_flux)


src3_kT = 0.905023 # in keV
src3_norm = 5.25623E-07
src3_flux = 4.3773e-14 # ergs/cm^2/s from 0.3-10 keV
get_infos(src3_kT, src3_flux)