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
import pandas as pd
try:
    fermi_ints = pd.read_csv(r"ASTRO531\hw3\FermiDiracIntegrals.txt", delim_whitespace=True)
except:
    fermi_ints = pd.read_csv(r"ASTRO531/hw3/FermiDiracIntegrals.txt", delim_whitespace=True)


fermi_ints["2/3 F_3/2"] = fermi_ints["2/3"]
fermi_ints["F_1/2"] = fermi_ints["F_3/2"]
del fermi_ints["F_3/2"]
del fermi_ints["2/3"]

def F_12(alpha):
    return (fermi_ints[fermi_ints["alpha"] == alpha])["F_1/2"].values
def F_32(alpha):
    return (fermi_ints[fermi_ints["alpha"] == alpha])["2/3 F_3/2"].values

from scipy.interpolate import CubicSpline
spline = CubicSpline(fermi_ints["F_1/2"], fermi_ints["alpha"])
def get_alpha_from_F_12(F_12):
    """Finds the closest alpha value in the table"""
    return spline(F_12)
try:
    hist = m.read_history(r"ASTRO531\hw2\MESA-Web_Job_01252661968\trimmed_history.data", as_table=True)
except:
    hist = m.read_history(r"ASTRO531/hw2/MESA-Web_Job_01252661968/trimmed_history.data", as_table=True)

star_age = list(hist["star_age"])
model_number = list(hist["model_number"])
# Find model closest to 4.5Gyr
star_age_delta = list(np.abs(np.array(star_age)-4.5E9))
idx = star_age_delta.index(min(star_age_delta))
# Model to use along with the stellar age it corresponds to:

# The closest model number to 294 in profiles.index is 296, so we use profile 8.


try:
    table = m.read_profile(r"ASTRO531\hw2\MESA-Web_Job_01252661968\profile8.data", as_table=True)
except:
    table = m.read_profile(r"ASTRO531/hw2/MESA-Web_Job_01252661968/profile8.data", as_table=True)
print(table.columns)

def get_F12_from_output(rho, T, n_e):
    """rho: density, T: temperature, n_e: # of free electrons/nucleon"""
    mu_e = rho/(n_e*u.u)
    return (cst.h**3/(4*np.pi)*(rho/(mu_e*u.u))*(2*cst.m_e*cst.k_B*T)**(-3/2)).decompose()

rhos = (10**table["logRho"])*(u.g/u.cm**3)
temperatures = (10**table["logT"])*(u.K)
n_free = (table["free_e"])/(u.cm**3)
mus = table["mu"]
nes = rhos/(mus*u.u * (1+1/n_free.value)) # Get n_e from n_free
mu_es = (rhos/(nes*u.u)).decompose()
print(cst.m_e)
print(cst.k_B)
print(cst.h)
import matplotlib.colors as colors
fig = plt.figure()
ax1 = plt.subplot(211)
ax2 = plt.subplot(212, sharex=ax1)
alphas = fermi_ints["alpha"].values
etas = table["eta"]
table_alpha = -etas
ax1.plot(alphas, F_32(alphas)/F_12(alphas), color="black", ls="-", lw=3, marker="None")
ymin, ymax = ax1.get_ylim()
xmin, xmax = ax1.get_xlim()
ax1.fill_betweenx([ymin, ymax], xmin, 0, color="blue", alpha=0.4, label="Degeneracy")
ax1.set_ylim(ymin, ymax)
ax1.set_xlim(xmin, xmax)
ax1.legend(fontsize=15)
print(len(rhos),len(mu_es),len(temperatures))
ax2_line = ((rhos/(mu_es)*temperatures**(-3/2)).decompose()).to(u.g/u.K**(3/2)/u.cm**3)


print(len(ax2_line), len(alphas))
ax2.plot(fermi_ints["alpha"],np.log10(fermi_ints["F_1/2"])-8.044, label="Clayton 2-6", color="black", ls="-", lw=3, marker="None")
ax2.plot(table_alpha, np.log10(ax2_line.value), color="red", ls="-", lw=3, marker="None", label="Solar model", zorder=2)
ax2.set_xlabel(r"$\alpha$ [$-$]")
ax1.set_ylabel(r"$2/3 F_{3/2}/F_{1/2}$ [$-$]")
ax2.set_ylabel(r"$\log(\rho/\mu_e T^{-3/2})$ [$\mathrm{g/cm^3/K^{3/2}}$]")
ax2.legend(fontsize=15)
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig(r"ASTRO531\hw3\3b.pdf")
plt.show()