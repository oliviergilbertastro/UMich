import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst

def planck(wav, T):
    wav = wav*(u.m)
    return (2*cst.h*cst.c**2*wav**(-5) / (np.exp(cst.h*cst.c/(wav*cst.k_B*T))-1))

def bbody_flux(wav, T):
    return np.pi*planck(wav,T)

f_3000K = lambda wav: bbody_flux(wav,3000*u.K).to(u.erg/u.s/u.m**3).value
f_Vega = lambda wav: bbody_flux(wav,10000*u.K).to(u.erg/u.s/u.m**3).value

def top_hat(wav, wav0, percentage):
    """
    wav0 : central wavelength
    percentage : float, 0 to 1
    actually not used since we directly set the bounds of the integration to wav1 and wav2 instead of integrating from 0 to infinity and adding this to the integral.
    """
    wav1, wav2 = (wav0-wav0*percentage/2).to(u.m), (wav0+wav0*percentage/2).to(u.m)
    return np.where((wav>wav1)&(wav<wav2), 1, 0)

from scipy.integrate import quad

def magnitude_filter(wav0, percentage):
    """
    wav0 : central wavelength of filter
    percentage : float, 0 to 1 for top-hat filter
    """
    wav1, wav2 = (wav0-wav0*percentage/2).to(u.m), (wav0+wav0*percentage/2).to(u.m)
    numerator_flux = quad(lambda wav: f_3000K(wav), wav1.value, wav2.value)[0]
    denominator_flux = quad(lambda wav: f_Vega(wav), wav1.value, wav2.value)[0]
    return -2.5*np.log10(numerator_flux/denominator_flux)

B_01 = magnitude_filter(440*u.nm, percentage=0.01)
V_01 = magnitude_filter(550*u.nm, percentage=0.01)
BV_01 = B_01-V_01

B_20 = magnitude_filter(440*u.nm, percentage=0.20)
V_20 = magnitude_filter(550*u.nm, percentage=0.20)
BV_20 = B_20-V_20

print("(B-V)1%", np.around(BV_01, decimals=3))
print("(B-V)20%", np.around(BV_20, decimals=3))

print("Delta(B-V)", np.around(BV_20-BV_01, decimals=3))