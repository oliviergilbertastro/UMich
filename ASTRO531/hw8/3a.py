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

plt.figure(figsize=(6,6))
ax = plt.subplot(111)
ax.plot(ms_sun["log_Teff"], ms_sun["log_L"], color="black", ls="-", lw=5, alpha=0.6, label=r"$Z=0.02$")
ax.plot(ms_Z0001["log_Teff"], ms_Z0001["log_L"], color="blue", ls="-", lw=5, alpha=0.6, label=r"$Z=0.001$")
ax.set_xlabel(r"$\log T_\mathrm{eff}$ [K]")
ax.set_ylabel(r"$\log L$ [$L_\odot$]")
ax.legend(fontsize=15, loc="lower left")
#ax.set_ylim(-2,5.25)
#ax.set_xlim(3.4,5.1)
ax.invert_xaxis()
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/3a.pdf")
plt.show()






def get_evolution_tracks(hist):
    Xc = hist["center_h1"]  # hydrogen
    Yc = hist["center_he4"] # helium
    Cc = hist["center_c12"] # carbon
    Oc = hist["center_o16"] # oxygen
    main_seq = Xc > 1e-3
    rgb = (Xc < 1e-3) & (Yc > 1e-2)
    hb = (Yc > 1e-3) & (Cc < 1e-3)
    early_agb = (Yc < 1e-3) & ((Cc + Oc) > 1e-3)