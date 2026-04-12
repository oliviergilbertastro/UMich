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


profile_Z0001_zams = m.read_profile(r"ASTRO531/hw8/MESA-Web_Job_04052666824/profile6.data", as_table=True)
profile_Z002_zams = m.read_profile(r"ASTRO531/hw2/MESA-Web_Job_01252661968/profile8.data", as_table=True)

print(profile_Z0001_zams.columns)

print(profile_Z0001_zams["radius"][-1])
print(profile_Z002_zams["radius"][-1])

table = profile_Z0001_zams

rmo = table["opacity"]*(u.cm**2/u.g)
logT = (table["logT"])

plt.figure(figsize=(10,8))
ax = plt.subplot(211)
ax.plot(logT, np.log10(rmo.value), lw=2, color="blue", label=r"MESA $Z=0.001$ star")
ax.invert_xaxis()
ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$\log$ RMO [$\mathrm{cm^2/g}$]")




# c)
m_e = cst.m_e
k_B = cst.k_B
h = cst.h
def saha(T,ne,chi=13.6,g_ratio=1/2):
    """returns the fraction n_(i+1)/n_i"""
    return ((2*np.pi*m_e*k_B*T/h**2)**(3/2) * (np.exp(-chi*u.eV/(k_B*T)))/ne*g_ratio).to(u.m**0)
logT = table["logT"]*(u.K)
T = (10**table["logT"])*(u.K)

# Estimate the electron density
n_e = (10**table["logRho"]*(u.g/u.cm**(3)))*table["free_e"]/(1.66E-24*u.g)
#print(n_e.to(u.cm**(-3)))

n_HII_HI = saha(T,n_e)
print(n_HII_HI)
X_HI = 1/(1 + n_HII_HI)
X_HII = n_HII_HI / (1 + n_HII_HI)

n_HeII_HeI = saha(T, n_e, 24.6, 2)
n_HeIII_HeII = saha(T, n_e, 54.4)
print(n_HeII_HeI)
print(n_HeIII_HeII)


He_tot = (1+n_HeII_HeI + n_HeII_HeI * n_HeIII_HeII)
X_HeI = 1 / He_tot
X_HeII = n_HeII_HeI / He_tot
X_HeIII = (n_HeII_HeI * n_HeIII_HeII) / He_tot

HI_mass_frac = X_HI #* table["h1"]
HII_mass_frac = X_HII #* table["h1"]
HeI_mass_frac = X_HeI #* table["he4"]
HeII_mass_frac = X_HeII #* table["he4"]
HeIII_mass_frac = X_HeIII #* table["he4"]

ax2 = plt.subplot(212, sharex=ax)
ax2.plot(logT, HI_mass_frac, ls="-", lw=2, color="red", label="HI")
ax2.plot(logT, HII_mass_frac, ls="--", lw=2, color="red", label="HII")
ax2.plot(logT, HeI_mass_frac, ls="-", lw=2, color="blue", label="HeI")
ax2.plot(logT, HeII_mass_frac, ls="--", lw=2, color="blue", label="HeII")
ax2.plot(logT, HeIII_mass_frac, ls="-.", lw=2, color="blue", label="HeIII")
ax2.legend()
ax2.set_xlabel(r"$\log T$ [$\mathrm{K}$]")
ax2.set_ylabel(r"Number fraction")



sigma_T = 6.65E-25*u.cm**2
mu_e = 1.1435
kappa_s = sigma_T/(mu_e*u.u)
kappa_s = kappa_s.to(u.cm**2/u.g)
print(kappa_s)
ax.axhline(np.log10(kappa_s.value), ls="--", color="black", label=r"$\kappa_T$")
ax.legend(fontsize=12)

plt.tight_layout()
plt.savefig(r"ASTRO531/hw8/figures/3c_rmoZ0001.pdf")
plt.show()





def plot_gradients(table):
    rmo = table["opacity"]*(u.cm**2/u.g)
    grad_a = table["grada"]
    grad_R = table["gradr"]
    grad_T = table["gradT"]
    logT = (table["logT"])
    R = table["radius"]/table["radius"][-1]

    conv_mask = grad_R > grad_a
    conv_inds = np.where(conv_mask)[0]
    if len(conv_inds) == 0:
        print("No convective region found.")
    else:
        # contiguous regions (in case there are multiple zones)
        splits = np.split(conv_inds, np.where(np.diff(conv_inds) != 1)[0] + 1)

        conv_regions = []
        conv_regions_T = []
        for region in splits:
            r_min = R[region[0]]
            r_max = R[region[-1]]
            logT_min = logT[region[0]]
            logT_max = logT[region[-1]]
            conv_regions.append((r_min, r_max))
            conv_regions_T.append((logT_min, logT_max))

        print("Convective regions (R/R_*):")
        for rmin, rmax in conv_regions:
            print(f"{rmin:.4f} to {rmax:.4f}  (ΔR = {rmax - rmin:.4f})")

    from scipy.interpolate import interp1d
    x_primary = np.array(logT)
    x_top_vals = np.array(R)

    plt.figure(figsize=(8,6))

    ax2 = plt.subplot(111)
    ax2.plot(logT, grad_a, ls="-", lw=2, color="red", label=r"$\nabla_a$")
    ax2.plot(logT, grad_R, ls="-", lw=2, color="blue", label=r"$\nabla_R$")
    ax2.plot(logT, grad_T, ls="-", lw=2, color="green", label=r"$\nabla_T$")
    ax2.invert_xaxis()
    ax2.set_xlabel(r"$\log T$ [$\mathrm{K}$]")
    ax2.set_ylabel(r"Gradient")
    xlims, ylims = ax2.get_xlim(), ax2.get_ylim()
    for T_min, T_max in conv_regions_T:
        ax2.fill_betweenx([*ylims],T_min,T_max,alpha=0.3,color="black",label="Convective")
    ax2.legend(fontsize=15)
    ax2.set_xlim(*xlims)
    ax2.set_ylim(*ylims)
    # create top axis
    secax = ax2.twiny()
    secax.set_xlim(ax2.get_xlim())

    # choose clean points
    xticks = ax2.get_xticks()
    logT_arr = np.array(logT)
    inds = [np.argmin(np.abs(logT_arr - xt)) for xt in xticks]
    inds = np.unique(inds)  # avoid duplicates
    inds = inds[1:-1]
    print(len(x_primary))
    print(inds)
    print(x_primary[inds])
    secax.set_xticks(x_primary[inds])
    secax.set_xticklabels([f"{x_top_vals[i]:.2f}" for i in inds])
    secax.set_xlabel(r"$R/R_\star$")
    plt.ylim(0, 1.2)
    plt.tight_layout()


plot_gradients(profile_Z002_zams)
plt.title(r"$Z=0.02$",fontsize=16)
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/3c_Z002_gradients.pdf")
plt.show()
plot_gradients(profile_Z0001_zams)
plt.title(r"$Z=0.001$",fontsize=16)
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/3c_Z0001_gradients.pdf")
plt.show()












import pandas as pd
# ALPHA DEGEN
try:
    fermi_ints = pd.read_csv(r"ASTRO531\hw3\FermiDiracIntegrals.txt", delim_whitespace=True)
except:
    fermi_ints = pd.read_csv(r"ASTRO531/hw3/FermiDiracIntegrals.txt", delim_whitespace=True)


fermi_ints["2/3 F_3/2"] = fermi_ints["2/3"]
fermi_ints["F_1/2"] = fermi_ints["F_3/2"]
del fermi_ints["F_3/2"]
del fermi_ints["2/3"]

def F_12(alpha):
    return (fermi_ints[fermi_ints["alpha"] == alpha])["F_1/2"].values
def F_32(alpha):
    return (fermi_ints[fermi_ints["alpha"] == alpha])["2/3 F_3/2"].values

from scipy.interpolate import CubicSpline
spline = CubicSpline(fermi_ints["F_1/2"], fermi_ints["alpha"])
def get_alpha_from_F_12(F_12):
    """Finds the closest alpha value in the table"""
    return spline(F_12)

def get_F12_from_output(rho, T, n_e):
    """rho: density, T: temperature, n_e: # of free electrons/nucleon"""
    mu_e = rho/(n_e*u.u)
    return (cst.h**3/(4*np.pi)*(rho/(mu_e*u.u))*(2*cst.m_e*cst.k_B*T)**(-3/2)).decompose()

def plot_alpha_degen(table, label, ax):
    rhos = (10**table["logRho"])*(u.g/u.cm**3)
    temperatures = (10**table["logT"])*(u.K)
    n_free = (table["free_e"])/(u.cm**3)
    mus = table["mu"]
    nes = rhos/(mus*u.u * (1+1/n_free.value)) # Get n_e from n_free
    etas = table["eta"]
    table_alpha = -etas
    F_12s = get_F12_from_output(rhos, temperatures, nes)
    alphas = get_alpha_from_F_12(F_12s)
    print(alphas)
    print(table_alpha)
    ax.plot(alphas, table_alpha, lw=2, ls="-", label=label)

ax2 = plt.subplot(111)
plot_alpha_degen(profile_Z002_zams, label=r"$Z=0.02$", ax=ax2)
plot_alpha_degen(profile_Z0001_zams, label=r"$Z=0.001$", ax=ax2)
xlims, ylims = ax2.get_xlim(), ax2.get_ylim()
ax2.plot([*xlims],[*xlims],ls="--", color="black", lw=2)
ax2.legend(fontsize=15)
ax2.set_xlim(*xlims)
ax2.set_ylim(*ylims)
plt.xlabel(r"$\alpha$ [$-$]")
plt.ylabel(r"$\alpha=-\eta$ from MESA [$-$]")
plt.tight_layout()
plt.savefig(r"ASTRO531/hw8/figures/3c_alpha.pdf")
plt.show()