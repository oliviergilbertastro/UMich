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

Rs = table["radius"]*cst.R_sun
Ts = (10**table["logT"])*u.K
logrhos = table["logRho"]
rhos = 10**logrhos*(u.g/u.cm**3)
kappa_Rs = table["opacity"] * u.cm**2/u.g


from scipy.integrate import cumulative_trapezoid
taus = cumulative_trapezoid(((kappa_Rs * rhos).to(1/u.cm)).value, (Rs.to(u.cm)).value, initial=0)
# Taus are from surface to core, and we want the optical depth at the surface to be 2/3, so we adjust:
taus += 2/3
taus = taus[::-1] # Taus are now from core to surface like everything else
print(taus[-1])

plt.figure(figsize=(6,6))
plt.plot(Rs.to(u.R_sun), taus, ls="-", color="black", marker="None")
plt.xlabel(r"$ R$ $[R_\odot]$")
plt.ylabel(r"$\tau$")
plt.tight_layout()
plt.savefig(r"ASTRO531/hw4/hw4_5b.pdf")
plt.show()

