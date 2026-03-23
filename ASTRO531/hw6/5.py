import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
from tqdm import tqdm


L = cst.L_sun
M = cst.M_sun
R_0 = 4*cst.R_sun
G = cst.G
alpha = 6/7

D_n = 5.991
B_n = 0.206

mu = 0.61
m_H = cst.m_p
k = cst.k_B

def R(t):
    """Radius as function of time"""
    return R_0/(1+(2*L*t*R_0)/(alpha*G*M**2))

def T_c(t):
    """Central temperature as function of time"""
    return 3**(1/3)*mu*m_H/k * B_n*D_n**(1/3)*G*M/R(t)


times = np.linspace(0,2*1E7,1000)*(u.year)
radii = R(times).to(u.R_sun)
temps = T_c(times).to(u.K)

plt.figure(figsize=(6,6))
ax2 = plt.subplot(111)
ax2.plot(times/1E6, radii, ls="-", lw=2, color="black")
ax2.legend(fontsize=17)
ax2.set_xlabel(r"$t$ [$\mathrm{Myr}$]")
ax2.set_ylabel(r"$R(t)$ [$R_\odot$]")
plt.tight_layout()
plt.savefig(r"ASTRO531\hw6\plot_radius.pdf")
plt.show()

plt.figure(figsize=(6,6))
ax2 = plt.subplot(111)
ax2.plot(times/1E6, temps, ls="-", lw=2, color="black")
ax2.axhline(1E7, ls="--",color="red", lw=2,label="H-burning start")
ax2.legend(fontsize=17)
ax2.set_xlabel(r"$t$ [$\mathrm{Myr}$]")
ax2.set_ylabel(r"$T_c(t)$ [$\mathrm{K}$]")
plt.tight_layout()
plt.savefig(r"ASTRO531\hw6\plot_temp.pdf")
plt.show()

t_KH = G*M**2/(2*R_0*L)
print(t_KH.to(u.Myr))