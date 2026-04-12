


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

hist_sun = m.read_history(r"ASTRO531/hw2/MESA-Web_Job_01252661968/trimmed_history.data", as_table=True)
hist_Z0001 = m.read_history(r"ASTRO531/hw8/MESA-Web_Job_04052666824/trimmed_history.data", as_table=True)

def cut_pre_main_sequence(hist):
    """Only keep the history starting from the main sequence"""
    Lh = 10**hist["log_LH"]
    L  = 10**hist["log_L"]
    # ZAMS = when H burning dominates total luminosity
    zams_index = np.argmax(Lh / L > 0.99)
    return hist[zams_index:]

ms_sun = cut_pre_main_sequence(hist_sun)
ms_Z0001 = cut_pre_main_sequence(hist_Z0001)








# log Rho log T diagram
plt.figure(figsize=(18,6))
ax = plt.subplot(111)
ax.plot(ms_sun["log_center_T"], ms_sun["log_center_Rho"], color="black", ls="-", lw=3, alpha=1, label=r"$Z=0.02$")
ax.plot(ms_Z0001["log_center_T"], ms_Z0001["log_center_Rho"], color="blue", ls="-", lw=3, alpha=1, label=r"$Z=0.001$")
ax.set_xlabel(r"$\log T_\mathrm{core}$ [K]")
ax.set_ylabel(r"$\log \rho_\mathrm{core}$ [$\mathrm{g/cm^3}$]")
ax.legend(fontsize=15)#, loc=("lower left" if i!=2 else "center left"))
#ax.set_ylim(-2,5.25)
#ax.set_xlim(3.4,5.1)
logrho = np.linspace(*ax.get_ylim(), 200)
logT_deg = 5.3571 + 0.4857 * (logrho)

ax.plot(logT_deg, logrho, 'b--')
ax.fill_betweenx(
    logrho,
    logT_deg,
    ax.get_xlim()[0],
    color="blue",
    alpha=0.2
)
ax.text(6.2,4,"Non-relativistic\ndegeneracy",color="blue",weight="bold",fontsize=16)
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/3b.pdf")
plt.show()