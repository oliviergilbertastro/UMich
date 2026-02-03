import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
import pandas as pd

try:
    fermi_ints = pd.read_csv(r"ASTRO531\hw3\FermiDiracIntegrals.txt", delim_whitespace=True)
except:
    fermi_ints = pd.read_csv(r"ASTRO531/hw3/FermiDiracIntegrals.txt", delim_whitespace=True)

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
    return 4*np.pi/cst.h**3 *(2*cst.m_e*cst.k_B*T)**(3/2) * F_12(alpha)

def P_e(alpha, T):
    return (n_e(alpha,T)*cst.k_B*T*F_32(alpha)/F_12(alpha)).to(u.dyne/u.cm**2)

mu_e = 1 # for Hydrogen, which we assume
def P_e_non_degen(n_e, T):
    return (n_e*cst.k_B*T).to(u.dyne/u.cm**2)
def P_e_full_degen(n_e, T):
    return (cst.h**2/(20*cst.m_e)*(3/np.pi)**(2/3)*n_e**(5/3)).to(u.dyne/u.cm**2)


print(cst.m_e)
print(cst.k_B)
print(cst.h)
import matplotlib.colors as colors
fig = plt.figure()
alphas = fermi_ints["alpha"].values
norm = colors.Normalize(vmin=np.quantile(alphas, 0.5), vmax=np.max(alphas))
plot_nes = np.linspace(1E14,1E29, 1000)
for T in [1E7, 1E8, 1E9]:
    n_es = n_e(alphas, T*u.K).to(u.cm**-3)
    Pressures = P_e(alphas, T*u.K)
    print(n_es)
    plt.plot(np.log10(n_es.value), np.log10(Pressures.value), ls="-", lw=3, marker="None", label=r"$T=10^{"+f"{np.around(np.log10(T)):.0f}"+r"}$K")
    col = plt.gca().lines[-1].get_color()
    plt.plot(np.log10(plot_nes), np.log10(P_e_non_degen(plot_nes*u.cm**-3, T*u.K).value), color=col, ls="--", lw=2)
    plt.plot(np.log10(plot_nes), np.log10(P_e_full_degen(plot_nes*u.cm**-3, T*u.K).value), color=col, ls=":", lw=2)
    alpha_color = plt.scatter(np.log10(n_es.value),  np.log10(Pressures.value), c=alphas, cmap="inferno", s=20, norm=norm, zorder=2, edgecolors="black")
    pass

plt.plot(np.log10([1E-16,1E-15]), np.log10(P_e_non_degen(np.array([1E-16,1E-15])*u.cm**-3, T*u.K).value), color="grey", ls="--", lw=3, label="Non-degenerate")
plt.plot(np.log10([1E-16,1E-15]), np.log10(P_e_full_degen(np.array([1E-16,1E-15])*u.cm**-3, T*u.K).value), color="grey", ls=":", lw=3, label="Fully-degenerate")
fig.colorbar(alpha_color, label=r"$\alpha$")
plt.xlabel(r"$\log n_e$ [$\mathrm{cm^{-3}}$]")
plt.ylabel(r"$\log P_e$ [$\mathrm{dyne/cm^{2}}$]")
plt.xlim(14, 29)
plt.ylim(4, 23)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig(r"ASTRO531\hw3\2d.pdf")
plt.show()