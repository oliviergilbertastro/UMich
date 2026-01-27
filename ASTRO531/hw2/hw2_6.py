import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
import mesa_web as m
from tqdm import tqdm
from scipy.integrate import trapezoid

try:
    hist = m.read_history(r"ASTRO531\hw2\MESA-Web_Job_01252661968\trimmed_history.data", as_table=True)
except:
    hist = m.read_history(r"ASTRO531/hw2/MESA-Web_Job_01252661968/trimmed_history.data", as_table=True)
print(hist)
star_age = list(hist["star_age"])
model_number = list(hist["model_number"])
# Find model closest to 4.5Gyr
star_age_delta = list(np.abs(np.array(star_age)-4.5E9))
idx = star_age_delta.index(min(star_age_delta))
# Model to use along with the stellar age it corresponds to:
print(model_number[idx], star_age[idx])
# The closest model number to 294 in profiles.index is 296, so we use profile 8.

# a) Verify Virial theorem:


try:
    table = m.read_profile(r"ASTRO531\hw2\MESA-Web_Job_01252661968\profile8.data", as_table=True)
except:
    table = m.read_profile(r"ASTRO531/hw2/MESA-Web_Job_01252661968/profile8.data", as_table=True)
print(table.columns)
E_int = table["energy"]*(u.erg/u.g)
M = table["mass"]*(u.M_sun)
R = table["radius"]*u.R_sun
U = trapezoid(E_int, M).to(u.erg) # Integrate the internal energy per unit mass over the mass
E_grav = trapezoid(-cst.G*M/R, M).to(u.erg) # Integrate the gravitational potential energy per unit mass over the mass
print("U:", U)
print("E_grav:", E_grav)
print("U/E_grav:", U/E_grav)
plt.plot(table["mass"])
plt.plot(table["radius"])
plt.show()

plt.plot(table["radius"], table["mass"])
plt.xlabel(r"Radius [$R_\odot$]")
plt.ylabel(r"Mass [$M_\odot$]")
plt.tight_layout()
plt.show()

