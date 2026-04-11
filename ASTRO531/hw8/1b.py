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



hist_05_Msun = m.read_history(r"ASTRO531/hw8/MESA-Web_Job_04052666820/trimmed_history.data", as_table=True)
hist_5_Msun = m.read_history(r"ASTRO531/hw8/MESA-Web_Job_04052666822/trimmed_history.data", as_table=True)
#plt.plot(hist_05_Msun["model_number"],hist_05_Msun["center_he4"]) # profile 7 is when the model starts to change the He mass frac
#plt.plot(hist_5_Msun["model_number"],hist_5_Msun["center_he4"]) # profile 8 is when the model starts to change the He mass frac
#plt.show()

profile_05_Msun = m.read_profile(r"ASTRO531/hw8/MESA-Web_Job_04052666820/profile7.data", as_table=True)
profile_5_Msun = m.read_profile(r"ASTRO531/hw8/MESA-Web_Job_04052666822/profile8.data", as_table=True)
print(profile_05_Msun.columns)

enclosed_M_05 = profile_05_Msun["mass"]*(u.M_sun)
total_M_05 = enclosed_M_05[-1]
enclosed_M_norm_05 = enclosed_M_05/total_M_05

enclosed_M_5 = profile_5_Msun["mass"]*(u.M_sun)
total_M_5 = enclosed_M_5[-1]
enclosed_M_norm_5 = enclosed_M_5/total_M_5

figsize=(6,5)
radius_norm_05 = profile_05_Msun["radius"]/(profile_05_Msun["radius"][-1])
radius_norm_5 = profile_5_Msun["radius"]/(profile_5_Msun["radius"][-1])
plt.figure(figsize=figsize)
ax = plt.subplot(111)
ax.plot(profile_05_Msun["logT"],radius_norm_05,marker="None",ls="-",label=r"0.5$M_\odot$")
ax.plot(profile_5_Msun["logT"],radius_norm_5,marker="None",ls="-",label=r"5$M_\odot$")
ax.set_ylabel(r"$r/R_\star$")
ax.set_xlabel(r"$\log T$ [$\mathrm{K}$]")
ax.legend()
ax.invert_xaxis()
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/1b_rad_logT.pdf")
plt.show()







m_e = cst.m_e
k_B = cst.k_B
h = cst.h
def plot_ionization_states_table(table):
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

    ax2 = plt.subplot(111)
    ax2.plot(logT, HI_mass_frac, ls="-", lw=2, color="red", label="HI")
    ax2.plot(logT, HII_mass_frac, ls="--", lw=2, color="red", label="HII")
    ax2.plot(logT, HeI_mass_frac, ls="-", lw=2, color="blue", label="HeI")
    ax2.plot(logT, HeII_mass_frac, ls="--", lw=2, color="blue", label="HeII")
    ax2.plot(logT, HeIII_mass_frac, ls="-.", lw=2, color="blue", label="HeIII")
    ax2.legend()
    ax2.set_xlabel(r"$\log T$ [$\mathrm{K}$]")
    ax2.set_ylabel(r"Number fraction")
    ax2.invert_xaxis()


plot_ionization_states_table(profile_05_Msun)
plt.title(r"0.5$M_\odot$",fontsize=16)
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/1b_05_ionization.pdf")
plt.close()
plot_ionization_states_table(profile_5_Msun)
plt.title(r"5$M_\odot$",fontsize=16)
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/1b_5_ionization.pdf")
plt.close()