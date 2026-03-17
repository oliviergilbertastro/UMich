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
logT = (table["logT"])
R = table["radius"]

from scipy.interpolate import interp1d
x_primary = np.array(logT)
x_top_vals = np.array(R)
# your data


plt.figure(figsize=(8,6))

ax2 = plt.subplot(111)
ax2.plot(logT, grad_a, ls="-", lw=2, color="red", label=r"$\nabla_a$")
ax2.plot(logT, grad_R, ls="-", lw=2, color="blue", label=r"$\nabla_R$")
ax2.plot(logT, grad_T, ls="-", lw=2, color="green", label=r"$\nabla_T$")
ax2.legend(fontsize=15)
ax2.invert_xaxis()
ax2.set_xlabel(r"$\log T$ [$\mathrm{K}$]")
ax2.set_ylabel(r"Gradient")
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
plt.ylim(0, 1.2)
plt.tight_layout()
plt.savefig(r"ASTRO531\hw5\plot3a.pdf")
plt.show()