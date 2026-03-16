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


rmo = table["opacity"]*(u.cm**2/u.g)
logT = (table["logT"])

plt.figure(figsize=(10,8))
ax = plt.subplot(211)
ax.plot(logT, rmo)
ax.invert_xaxis()
ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"RMO [$\mathrm{cm^2/g}$]")




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

HI_mass_frac = X_HI #* table["h1"]
HII_mass_frac = X_HII #* table["h1"]
HeI_mass_frac = X_HeI #* table["he4"]
HeII_mass_frac = X_HeII #* table["he4"]
HeIII_mass_frac = X_HeIII #* table["he4"]

ax2 = plt.subplot(212, sharex=ax)
ax2.plot(logT, HI_mass_frac, ls="-", lw=2, color="red", label="HI")
ax2.plot(logT, HII_mass_frac, ls="--", lw=2, color="red", label="HII")
ax2.plot(logT, HeI_mass_frac, ls="-", lw=2, color="blue", label="HeI")
ax2.plot(logT, HeII_mass_frac, ls="--", lw=2, color="blue", label="HeII")
ax2.plot(logT, HeII_mass_frac, ls="-.", lw=2, color="blue", label="HeIII")
ax2.legend()
ax2.set_xlabel(r"$\log T$ [$\mathrm{K}$]")
ax2.set_ylabel(r"Number fraction")








plt.tight_layout()
plt.savefig(r"ASTRO531\hw5\plot2a.pdf")
plt.show()