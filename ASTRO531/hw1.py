import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
from astropy.constants import c

def flux(mag, zp):
    return 10**((mag+48.958+zp)/(-2.5))*(10**(-20)*u.erg*(u.cm)**(-2)*(u.s)**(-1)*(u.Hz)**(-1))

#         U      B      V       R      I     J     H    K
wavs = np.array([0.366, 0.438, 0.545, 0.641, 0.798, 1.22, 1.63, 2.19])*u.um

mags = np.array([14.61, 13.50, 12.11, 11.23, 10.33, 9.15, 8.40, 8.25])
mags_dereddened = np.array([13.86, 12.86, 11.63, 10.82, 10.03, 9.01, 8.31, 8.20])

zps = np.array([-0.152, -0.602, 0.000, 0.555, 1.271, 2.655, 3.760, 4.906])

fluxes = flux(mags_dereddened,zps)
fluxes_uncorr = flux(mags,zps)
freqs = c/wavs
nuFnu = (freqs*fluxes).to(u.erg*(u.cm)**(-2)*(u.s)**(-1))
nuFnu_uncorr = (freqs*fluxes_uncorr).to(u.erg*(u.cm)**(-2)*(u.s)**(-1))

plt.plot(np.log10(wavs.value), np.log10(nuFnu.value), ls="None", marker="o", color="black", label="Corrected")
plt.plot(np.log10(wavs.value), np.log10(nuFnu_uncorr.value), ls="None", marker="o", color="blue", label="Uncorrected")
plt.xlabel(r"$\log (\lambda)$ [$\mathrm{\mu m}$]")
plt.ylabel(r"$\log (\nu F_\nu)$ [$\mathrm{erg/s/cm^2}$]")
plt.legend(fontsize=15)
plt.tight_layout()
plt.show()