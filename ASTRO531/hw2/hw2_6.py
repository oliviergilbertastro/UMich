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



# b)
print((cst.G*u.M_sun**2/u.R_sun).to(u.erg))
alpha = -E_grav/((cst.G*u.M_sun**2/u.R_sun).to(u.erg))
print("alpha:", alpha)


# c)
m_e = cst.m_e
k_B = cst.k_B
h = cst.h
def saha(T,ne,chi=13.6,g_ratio=1/2):
    """returns the fraction n_(i+1)/n_i"""
    return ((2*np.pi*m_e*k_B*T/h**2)**(3/2) * (np.exp(-chi*u.eV/(k_B*T)))/ne*g_ratio).to(u.m**0)
logT = table["logT"]*(u.K)
T = (10**table["logT"])*(u.K)

# Estimate the electron density
n_e = (10**table["logRho"]*(u.g/u.cm**(3)))*table["free_e"]/(1.66E-24*u.g)
#print(n_e.to(u.cm**(-3)))

n_HII_HI = saha(T,n_e)
print(n_HII_HI)
X_HI = 1/(1 + n_HII_HI)
X_HII = n_HII_HI / (1 + n_HII_HI)

n_HeII_HeI = saha(T, n_e, 24.6, 2)
n_HeIII_HeII = saha(T, n_e, 54.4)
print(n_HeII_HeI)
print(n_HeIII_HeII)


He_tot = (1+n_HeII_HeI + n_HeII_HeI * n_HeIII_HeII)
X_HeI = 1 / He_tot
X_HeII = n_HeII_HeI / He_tot
X_HeIII = (n_HeII_HeI * n_HeIII_HeII) / He_tot

HI_mass_frac = X_HI * table["h1"]
HII_mass_frac = X_HII * table["h1"]
HeI_mass_frac = X_HeI * table["he4"]
HeII_mass_frac = X_HeII * table["he4"]
HeIII_mass_frac = X_HeIII * table["he4"]

ax = plt.subplot(111)
ax.plot(logT, HI_mass_frac, ls="-", lw=2, color="red", label="HI")
ax.plot(logT, HII_mass_frac, ls="--", lw=2, color="red", label="HII")
ax.plot(logT, HeI_mass_frac, ls="-", lw=2, color="blue", label="HeI")
ax.plot(logT, HeII_mass_frac, ls="--", lw=2, color="blue", label="HeII")
ax.plot(logT, HeII_mass_frac, ls="-.", lw=2, color="blue", label="HeIII")
ax.legend()
ax.invert_xaxis()
ax.set_xlabel(r"$\log T$ [$\mathrm{K}$]")
ax.set_ylabel(r"Mass Fraction")
plt.tight_layout()
plt.savefig(r"ASTRO531\hw2\hw2_6c.pdf")
plt.show()

ax = plt.subplot(111)
ax.plot(table["radius"], HI_mass_frac, ls="-", lw=2, color="red", label="HI")
ax.plot(table["radius"], HII_mass_frac, ls="--", lw=2, color="red", label="HII")
ax.plot(table["radius"], HeI_mass_frac, ls="-", lw=2, color="blue", label="HeI")
ax.plot(table["radius"], HeII_mass_frac, ls="--", lw=2, color="blue", label="HeII")
ax.plot(table["radius"], HeII_mass_frac, ls="-.", lw=2, color="blue", label="HeIII")
ax.legend()
ax.set_xlabel(r"Radius [$R_\odot$]")
ax.set_ylabel(r"Mass Fraction")
plt.tight_layout()
plt.savefig(r"ASTRO531\hw2\hw2_6c2.pdf")
plt.show()