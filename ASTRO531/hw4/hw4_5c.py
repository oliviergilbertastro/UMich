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

print(Hs[-1].to(u.erg/u.cm**2/u.s))
T_eff_sun = 5778*u.K
print((cst.sigma_sb*T_eff_sun**4).to(u.erg/u.cm**2/u.s))

