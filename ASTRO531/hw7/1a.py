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


# Find models closest to 1.5, 8.6, 20.3 Myr:

try:
    hist = m.read_history(r"ASTRO531\hw2\MESA-Web_Job_01252661968\trimmed_history.data", as_table=True)
except:
    hist = m.read_history(r"ASTRO531/hw2/MESA-Web_Job_01252661968/trimmed_history.data", as_table=True)
print(hist)
star_age = list(hist["star_age"])
model_number = list(hist["model_number"])
# Find model closest to 4.5Gyr
for age in (1.5E6, 8.6E6, 20.3E6, 4.5E9):
    star_age_delta = list(np.abs(np.array(star_age)-age))
    idx = star_age_delta.index(min(star_age_delta))
    # Model to use along with the stellar age it corresponds to:
    print(age, model_number[idx], star_age[idx])

rmo = table["opacity"]*(u.cm**2/u.g)
grad_a = table["grada"]
grad_R = table["gradr"]
grad_T = table["gradT"]
logT = (table["logT"])

hist_pre_main = hist[hist["star_age"]<4.5E9]
plt.figure(figsize=(6,6))
ax = plt.subplot(111)
ax.plot(hist_pre_main["log_Teff"], hist_pre_main["log_L"], color="black", label="Evolutionary track")
ax.invert_xaxis()
ax.set_xlabel(r"$\log T_\mathrm{eff}$ [K]")
ax.set_ylabel(r"$\log L$ [$L_\odot$]")
ax.legend(fontsize=15)
plt.show()
