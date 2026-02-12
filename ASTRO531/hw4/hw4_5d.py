import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
from tqdm import tqdm
import mesa_web as m

try:
    table = m.read_profile(r"ASTRO531\hw2\MESA-Web_Job_01252661968\profile8.data", as_table=True)
except:
    table = m.read_profile(r"ASTRO531/hw2/MESA-Web_Job_01252661968/profile8.data", as_table=True)
print(table.columns)

Ls = table["luminosity"]*cst.L_sun
Rs = table["radius"]*cst.R_sun
Hs = Ls/(4*np.pi*Rs**2)
Ts = (10**table["logT"])*u.K
surface_flux = cst.sigma_sb*(Ts[-1]**4)

Hs = Hs.to(u.erg/u.cm**2/u.s)
flux_approx = ((Rs[-1]/(Rs))**2*surface_flux).to(u.erg/u.cm**2/u.s)

plt.figure(figsize=(6,6))
plt.plot(Rs.to(u.R_sun), np.log10(Hs.value), ls="-", color="black", marker="None", label="MESA")
plt.plot(Rs.to(u.R_sun), np.log10(flux_approx.value), ls="-", color="red", marker="None", label=r"$H_\mathrm{surface}/r^2$")
plt.xlabel(r"$ R$ $[R_\odot]$")
plt.legend(fontsize=15)
plt.ylabel(r"$\log H$ [$\mathrm{erg/cm^2/s}$]")
plt.tight_layout()
plt.savefig(r"ASTRO531/hw4/hw4_5d.pdf")
plt.show()

