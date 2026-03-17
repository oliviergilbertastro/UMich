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
grad_a = table["grada"]
grad_R = table["gradr"]
grad_T = table["gradT"]
L = table["luminosity"]*u.L_sun
rho = 10**table["logRho"]*(u.g/u.cm**3)
R = table["radius"]*u.R_sun
M = table["mass"]*u.M_sun
P = table["pressure"]*(u.dyn/u.cm**2)
logT = (table["logT"])
T = 10**logT * u.K

g = (cst.G * M / R**2).to(u.cm / u.s**2) # gravity
H_P = P / (rho * g) # Pressure scale height
H_P = H_P.to(u.R_sun)
temp_excess = T*(grad_T-grad_a)
temp_excess = temp_excess.to(u.K)
from scipy.interpolate import interp1d
x_primary = np.array(logT)
x_top_vals = np.array(R)
# your data


plt.figure(figsize=(8,6))

ax2 = plt.subplot(111)
ax2.plot(logT, temp_excess, ls="-", lw=2, color="red")
ax2.legend(fontsize=15)
ax2.invert_xaxis()
ax2.set_xlabel(r"$\log T$ [$\mathrm{K}$]")
ax2.set_ylabel(r"$\Delta \nabla T \ell$ [$\mathrm{K}$]")
#ax2.set_yscale("log")
# create top axis
secax = ax2.twiny()
secax.set_xlim(ax2.get_xlim())

# choose clean points
xticks = ax2.get_xticks()
logT_arr = np.array(logT)
inds = [np.argmin(np.abs(logT_arr - xt)) for xt in xticks]
inds = np.unique(inds)  # avoid duplicates
inds = inds[1:-1]
print(len(x_primary))
print(inds)
print(x_primary[inds])
secax.set_xticks(x_primary[inds])
secax.set_xticklabels([f"{x_top_vals[i]:.2f}" for i in inds])
secax.set_xlabel(r"Radius [$R_\odot$]")
ax2.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig(r"ASTRO531\hw5\plot3h.pdf")
plt.show()