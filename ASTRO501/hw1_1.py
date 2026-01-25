import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst

def planck(wav, T):
    wav = wav.to(u.m)
    return 2*cst.h*cst.c**2*wav**(-5) / (np.exp(cst.h*cst.c/(wav*cst.k_B*T))-1)

def bbody_flux(wav, T):
    return np.pi*planck(wav,T)

f_3000K = lambda wav: bbody_flux(wav,3000)
f_Vega = lambda wav: bbody_flux(wav,10000)

def top_hat(wav, wav0, percentage):
    """
    wav0 : central wavelength
    percentage : float, 0 to 1
    actually not used since we directly set the bounds of the integration to wav1 and wav2 instead of integrating from 0 to infinity and adding this to the integral.
    """
    wav1, wav2 = wav0-wav0*percentage/2, wav0+wav0*percentage/2
    return np.where((wav>wav1)&(wav<wav2), 1, 0)

from scipy.integrate import quad

def magnitude_filter(wav0, percentage):
    """
    wav0 : central wavelength
    percentage : float, 0 to 1
    """
    wav1, wav2 = wav0-wav0*percentage/2, wav0+wav0*percentage/2
    num_flux = quad(lambda wav: f_3000K(wav))

    return -2.5*np.log10()