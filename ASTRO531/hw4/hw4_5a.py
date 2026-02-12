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
urs = 4*cst.sigma_sb/cst.c*Ts**4

ratios = (3*Hs/(cst.c*urs)).decompose()

plt.figure(figsize=(6,6))
plt.plot(Rs.to(u.R_sun), np.log10(ratios), ls="-", color="black", marker="None")
plt.xlabel(r"$ R$ $[R_\odot]$")
plt.ylabel(r"$\log( 3H/(u_r c))$")
plt.tight_layout()
plt.savefig(r"ASTRO531/hw4/hw4_5a.pdf")
plt.show()

