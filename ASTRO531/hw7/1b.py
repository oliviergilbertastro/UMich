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
    hist = m.read_history(r"ASTRO531\hw2\MESA-Web_Job_01252661968\trimmed_history.data", as_table=True)
except:
    hist = m.read_history(r"ASTRO531/hw2/MESA-Web_Job_01252661968/trimmed_history.data", as_table=True)
print(hist)
star_age = list(hist["star_age"])
model_number = list(hist["model_number"])

hist_pre_main = hist[hist["star_age"]<4.5E9]



try:
    table = m.read_profile(r"ASTRO531\hw7\MESA-Web_Job_03292665778\profile8.data", as_table=True)
except:
    table = m.read_profile(r"ASTRO531/hw7/MESA-Web_Job_03292665778/profile8.data", as_table=True)

# Find models closest to 1.5, 8.6, 20.3 Myr:

try:
    hist = m.read_history(r"ASTRO531\hw7\MESA-Web_Job_03292665778\trimmed_history.data", as_table=True)
except:
    hist = m.read_history(r"ASTRO531/hw7/MESA-Web_Job_03292665778/trimmed_history.data", as_table=True)
print(hist)
star_age = list(hist["star_age"])
model_number = list(hist["model_number"])
# Find model closest to 4.5Gyr
for age in (1.5E6, 8.6E6, 20.3E6, 4.5E9):
    star_age_delta = list(np.abs(np.array(star_age)-age))
    idx = star_age_delta.index(min(star_age_delta))
    # Model to use along with the stellar age it corresponds to:
    print(age, model_number[idx], star_age[idx])

model_186 = m.read_profile(r"ASTRO531/hw7/MESA-Web_Job_03292665778/profile94.data", as_table=True)
model_204 = m.read_profile(r"ASTRO531/hw7/MESA-Web_Job_03292665778/profile103.data", as_table=True)
model_219 = m.read_profile(r"ASTRO531/hw7/MESA-Web_Job_03292665778/profile111.data", as_table=True)
idxs = [186,204,219]

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

atr_186, atr_204, atr_219 = get_attributes(model_186), get_attributes(model_204), get_attributes(model_219)

plt.figure(figsize=(8,9))
for i, atr in enumerate([atr_186, atr_204, atr_219]):
    ax = plt.subplot(3,1,i+1)
    ax.plot(atr[0], atr[3], ls=":", marker="None", color="red", label=r"$\nabla_a$")
    ax.plot(atr[0], atr[4], ls="--", marker="None", color="blue", label=r"$\nabla_R$")
    ax.plot(atr[0], atr[5], ls="-", marker="None", color="black", label=r"$\nabla$")
    x_lims = ax.get_xlim()
    ax.text((x_lims[1]-x_lims[0])/4+x_lims[0], 0.8, ["1.54Myr","8.43Myr","20.3Myr"][i],fontsize=17)
    ax.invert_xaxis()
    ax.set_ylim(0,1)
    ax.legend(fontsize=15)
ax.set_ylabel(r"")
ax.set_xlabel(r"$\log T$ [K]")
plt.tight_layout()
plt.savefig("ASTRO531/hw7/1bi.pdf")
plt.show()

plt.figure(figsize=(7,6))
ax = plt.subplot(111)
for i, atr in enumerate([atr_186, atr_204, atr_219]):
    ax.plot(atr[0], np.log10(atr[6].value), ls="-", marker="None", label=["1.54Myr","8.43Myr","20.3Myr"][i])
    ax.invert_xaxis()
ax.set_ylabel(r"$\log$ RMO [$\mathrm{cm^2/g}$]")
ax.legend(fontsize=15)
ax.set_xlabel(r"$\log T$ [K]")
plt.tight_layout()
plt.savefig("ASTRO531/hw7/1bii.pdf")
plt.show()

plt.figure(figsize=(7,6))
ax = plt.subplot(111)
for i, atr in enumerate([atr_186, atr_204, atr_219]):
    ax.plot(atr[0], atr[7].value, ls="-", marker="None", label=["1.54Myr","8.43Myr","20.3Myr"][i])
    ax.invert_xaxis()
ax.axhline(5/3, label="Monoatomic gas", color="black", ls=":")
ax.legend(fontsize=15)
ax.set_ylabel(r"$\Gamma_2$")
ax.set_xlabel(r"$\log T$ [K]")
plt.tight_layout()
plt.savefig("ASTRO531/hw7/1biii.pdf")
plt.show()

plt.figure(figsize=(7,6))
ax = plt.subplot(111)
for i, atr in enumerate([atr_186, atr_204, atr_219]):
    ax.plot(atr[0], atr[1].value, ls="-", marker="None", label=["1.54Myr","8.43Myr","20.3Myr"][i])
    ax.invert_xaxis()
ax.legend(fontsize=15)
ax.set_ylabel(r"$R$ [$\mathrm{R_\odot}$]")
ax.set_xlabel(r"$\log T$ [K]")
plt.tight_layout()
plt.savefig("ASTRO531/hw7/1biv.pdf")
plt.show()


plt.figure(figsize=(7,6))
ax = plt.subplot(111)
for i, atr in enumerate([atr_186, atr_204, atr_219]):
    ax.plot(atr[0], atr[2].value, ls="-", marker="None", label=["1.54Myr","8.43Myr","20.3Myr"][i])
    ax.invert_xaxis()
ax.legend(fontsize=15)
ax.set_ylabel(r"$\rho$ [$g/cm^3$]")
ax.set_yscale("log")
ax.set_xlabel(r"$\log T$ [K]")
plt.tight_layout()
plt.savefig("ASTRO531/hw7/1bv.pdf")
plt.show()