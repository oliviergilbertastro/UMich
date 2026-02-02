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

def occupation_ind(E_kT, alpha):
    return 1/(np.exp(alpha+E_kT)+1)

E_kTs = np.linspace(1E-10, 1E4, 1000)
for alph in [10,0,-4,-30,-100]:
    plt.plot(np.log10(E_kTs), occupation_ind(E_kTs, alph), ls="-", marker="None", label=rf"$\alpha=${alph}")

plt.xlabel(r"$\log(E/kT)$ [$-$]")
plt.ylabel(r"P(p) [$-$]")
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig(r"ASTRO531\hw3\2b.pdf")
plt.show()