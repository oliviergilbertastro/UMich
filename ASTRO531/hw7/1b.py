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

print(model_186.columns)
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
    gamma2 = (32 - 24*beta-3*beta**2) / (24-14*beta-3*beta**2)
    return logT, R, rho, grad_a, grad_R, grad_T, rmo, gamma2

atr_186, atr_204, atr_219 = get_attributes(model_186), get_attributes(model_204), get_attributes(model_219)

plt.figure(figsize=(6,6))
ax = plt.subplot(111)
ax.plot(atr_186[0], atr_186["log_L"], ls="None", marker="o", color="red", label=("Selected PMS models" if i==0 else None))
ax.invert_xaxis()
ax.set_xlabel(r"$\log T_\mathrm{eff}$ [K]")
ax.set_ylabel(r"$\log L$ [$L_\odot$]")
ax.legend(fontsize=15)
plt.tight_layout()
plt.savefig("ASTRO531/hw7/1b.pdf")
plt.show()