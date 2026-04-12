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




profile_tip_rgb = m.read_profile(r"ASTRO531/hw2/MESA-Web_Job_01252661968/profile189.data", as_table=True)

def get_attributes(table):
    rmo = table["opacity"]*(u.cm**2/u.g)
    grad_a = table["grada"]
    grad_R = table["gradr"]
    grad_T = table["gradT"]
    rho = 10**table["logRho"]*(u.g/u.cm**3)
    R = table["radius"]*u.R_sun
    logT = (table["logT"])
    Pgas = table['pgas']
    Prad = table['prad']
    Ptot = Pgas + Prad
    beta = Pgas / Ptot
    gamma2 = (32 - 24*beta-3*beta**2) / (24-18*beta-3*beta**2)
    return logT, R, rho, grad_a, grad_R, grad_T, rmo, gamma2

atr_tip_rgb = get_attributes(profile_tip_rgb)

plt.figure(figsize=(6,5))
for i, atr in enumerate([atr_tip_rgb]):
    ax = plt.subplot(1,1,i+1)
    ax.plot(atr[0], atr[3], ls=":", marker="None", color="red", lw=3, label=r"$\nabla_a$")
    ax.plot(atr[0], atr[4], ls="--", marker="None", color="blue", lw=3, label=r"$\nabla_R$")
    ax.plot(atr[0], atr[5], ls="-", marker="None", color="black", lw=3, label=r"$\nabla$")
    conv_mask = atr[4] > atr[3]
    conv_inds = np.where(conv_mask)[0]
    splits = np.split(conv_inds, np.where(np.diff(conv_inds) != 1)[0] + 1)
    conv_region = (atr[0][conv_inds[0]], atr[0][conv_inds[-1]])
    x_lims = ax.get_xlim()
    y_lims = ax.get_ylim()
    ax.fill_betweenx([*y_lims],*conv_region,alpha=0.3,color="black",label="Convective")
    x_lims = ax.get_xlim()
    #ax.text((x_lims[1]-x_lims[0])/4+x_lims[0], 0.8, ["1.54Myr","8.43Myr","20.3Myr"][i],fontsize=17)
    ax.invert_xaxis()
    ax.set_ylim(0,1)
    ax.legend(fontsize=15)
ax.set_ylabel(r"")
ax.set_xlabel(r"$\log T$ [K]")
plt.tight_layout()
plt.savefig("ASTRO531/hw8/4a.pdf")
plt.show()

plt.figure(figsize=(6,5))
ax = plt.subplot(111)
for i, atr in enumerate([atr_tip_rgb]):
    ax.plot(atr[0], np.log10(atr[6].value), ls="-", lw=3, marker="None")
    ax.invert_xaxis()
ax.set_ylabel(r"$\log$ RMO [$\mathrm{cm^2/g}$]")
ax.legend(fontsize=15)
ax.set_xlabel(r"$\log T$ [K]")
plt.tight_layout()
plt.savefig("ASTRO531/hw8/4b.pdf")
plt.show()

plt.figure(figsize=(6,5))
ax = plt.subplot(111)
for i, atr in enumerate([atr_tip_rgb]):
    ax.plot(atr[0], atr[7].value, ls="-", marker="None")
    ax.invert_xaxis()
ax.axhline(5/3, label="Monoatomic gas", color="black", ls=":", lw=3)
ax.axhline(4/3, label="Onset of instability", color="red", ls=":", lw=3)
ax.legend(fontsize=15, loc="center right")
ax.set_ylabel(r"$\Gamma_2$")
ax.set_xlabel(r"$\log T$ [K]")
plt.tight_layout()
plt.savefig("ASTRO531/hw8/4c.pdf")
plt.show()

plt.figure(figsize=(6,5))
ax = plt.subplot(111)
for i, atr in enumerate([atr_tip_rgb]):
    ax.plot(atr[1].value,atr[0], ls="-", lw=3, marker="None")
    #ax.invert_xaxis()
ax.legend(fontsize=15)
ax.set_xlabel(r"$R$ [$\mathrm{R_\odot}$]")
ax.set_ylabel(r"$\log T$ [K]")
plt.tight_layout()
plt.savefig("ASTRO531/hw8/4d.pdf")
plt.show()

