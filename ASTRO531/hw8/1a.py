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
# Plot the 
plt.figure(figsize=figsize)
ax = plt.subplot(111)
ax.plot(enclosed_M_norm_05,(profile_05_Msun["logRho"]),marker="None",ls="-",label=r"0.5M_\odot")
ax.plot(enclosed_M_norm_5,(profile_5_Msun["logRho"]),marker="None",ls="-",label=r"5M_\odot")
ax.set_xlabel(r"$M(r)/M_\star$")
ax.set_ylabel(r"$\log \rho$ [$\mathrm{g/cm^3}$]")
ax.legend()
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/1a_logrho.pdf")
plt.show()

plt.figure(figsize=figsize)
ax = plt.subplot(111)
ax.plot(enclosed_M_norm_05,(profile_05_Msun["logT"]),marker="None",ls="-",label=r"0.5M_\odot")
ax.plot(enclosed_M_norm_5,(profile_5_Msun["logT"]),marker="None",ls="-",label=r"5M_\odot")
ax.set_xlabel(r"$M(r)/M_\star$")
ax.set_ylabel(r"$\log T$ [$\mathrm{K}$]")
ax.legend()
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/1a_logT.pdf")
plt.show()


plt.figure(figsize=figsize)
ax = plt.subplot(111)
ax.plot(enclosed_M_norm_05,np.log10(profile_05_Msun["eps_nuc"]),marker="None",ls="-",label=r"0.5M_\odot")
ax.plot(enclosed_M_norm_5,np.log10(profile_5_Msun["eps_nuc"]),marker="None",ls="-",label=r"5M_\odot")
ax.set_xlabel(r"$M(r)/M_\star$")
ax.set_ylabel(r"$\log \epsilon_\mathrm{nuc}$ [$\mathrm{erg/g/s}$]")
ax.legend()
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/1a_epsnuc.pdf")
plt.show()

luminosity_norm_05 = profile_05_Msun["luminosity"]/(profile_05_Msun["luminosity"][-1])
luminosity_norm_5 = profile_5_Msun["luminosity"]/(profile_5_Msun["luminosity"][-1])
plt.figure(figsize=figsize)
ax = plt.subplot(111)
ax.plot(enclosed_M_norm_05,luminosity_norm_05,marker="None",ls="-",label=r"0.5M_\odot")
ax.plot(enclosed_M_norm_5,luminosity_norm_5,marker="None",ls="-",label=r"5M_\odot")
ax.set_xlabel(r"$M(r)/M_\star$")
ax.set_ylabel(r"$L(r)/L_\star$")
ax.legend()
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/1a_luminosity.pdf")
plt.show()

radius_norm_05 = profile_05_Msun["radius"]/(profile_05_Msun["radius"][-1])
radius_norm_5 = profile_5_Msun["radius"]/(profile_5_Msun["radius"][-1])
plt.figure(figsize=figsize)
ax = plt.subplot(111)
ax.plot(radius_norm_05,enclosed_M_norm_05,marker="None",ls="-",label=r"0.5M_\odot")
ax.plot(radius_norm_5,enclosed_M_norm_5,marker="None",ls="-",label=r"5M_\odot")
ax.set_xlabel(r"$r/R_\star$")
ax.set_ylabel(r"$M(r)/M_\star$")
ax.legend()
plt.tight_layout()
plt.savefig("ASTRO531/hw8/figures/1a_rad_mass.pdf")
plt.show()