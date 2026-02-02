import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
import pandas as pd

fermi_ints = pd.read_csv(r"ASTRO531\hw3\FermiDiracIntegrals.txt", delim_whitespace=True)
print(fermi_ints.columns)
fermi_ints["2/3 F_3/2"] = fermi_ints["2/3"]
fermi_ints["F_1/2"] = fermi_ints["F_3/2"]
del fermi_ints["F_3/2"]
del fermi_ints["2/3"]
print(fermi_ints)
def F_12(alpha):
    return (fermi_ints[fermi_ints["alpha"] == alpha])["F_1/2"].values
def F_32(alpha):
    return (fermi_ints[fermi_ints["alpha"] == alpha])["2/3 F_3/2"].values

def n_e(alpha, T):
    return 4*np.pi/cst.h**3 *(2*cst.m_e*cst.k_B*T)**(2/3) * F_12(alpha)

def P_e(alpha, T):
    return n_e(alpha,T)*cst.k_B*T*F_32(alpha)/F_12(alpha)
print(cst.m_e)
print(cst.k_B)
print(cst.h)

alphas = fermi_ints["alpha"].values
for T in [1E7, 1E8, 1E9]:
    n_es = n_e(alphas, T*u.K).to(u.cm**-3)
    Pressures = P_e(alphas, T*u.K)
    print(n_es)
    plt.plot(n_es, Pressures, ls="-", marker="None", label=r"$T=10^{"+f"{np.around(np.log10(T))}"+r"}K")
    pass

plt.xlabel(r"$n_e$ [$-$]")
plt.ylabel(r"P_e [$-$]")
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig(r"ASTRO531\hw3\2c.pdf")
plt.show()