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

def cut_pre_main_sequence(hist):
    """Only keep the history starting from the main sequence"""
    Lh = 10**hist["log_LH"]
    L  = 10**hist["log_L"]
    # ZAMS = when H burning dominates total luminosity
    zams_index = np.argmax(Lh / L > 0.99)
    return hist[zams_index:]

ms_sun = cut_pre_main_sequence(hist_sun)

plt.figure(figsize=(6,6))
ax = plt.subplot(111)
#ax.plot(ms_sun["log_Teff"], ms_sun["log_L"], color="black", ls="-", lw=5, alpha=0.6)
line = ax.scatter(ms_sun["log_Teff"], ms_sun["log_L"], c=ms_sun["model_number"], alpha=1, label=r"$1M_\odot$")
for i in range(len(ms_sun)//1000):
    plt.text(ms_sun["log_Teff"][i*1000], ms_sun["log_L"][i*1000], ms_sun["model_number"][i*1000])
plt.colorbar(line, label="Model number")
ax.set_xlabel(r"$\log T_\mathrm{eff}$ [K]")
ax.set_ylabel(r"$\log L$ [$L_\odot$]")
ax.legend(fontsize=15, loc="lower left")
ax.set_ylim(-2,5.25)
ax.set_xlim(3.4,5.1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()