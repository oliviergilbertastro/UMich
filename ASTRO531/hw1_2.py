import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
from astropy.constants import c, L_sun, G, M_sun, R_sun, sigma_sb, m_e, k_B, h

N_e = 1E13 *(u.cm)**(-3)
def saha(T):
    return (2*np.pi*m_e*k_B*T/h**2)**(3/2) * (np.exp(-13.6*u.eV/(k_B*T))/N_e)
def fraction(T):
    return (4*np.exp(-10.2*u.eV/(k_B*T))) / ( 1 + saha(T) )
temps = np.linspace(3000,30000,1000)*u.K
f2 = fraction(temps)

plt.plot(temps, f2, ls="-", marker="None", color="black")
plt.xlabel("Temperature [K]")
plt.ylabel(r"$\log f_2$")
plt.tight_layout()
plt.savefig("ASTRO531/hw1_2.pdf")
plt.show()
f_nus = np.array([1.790, 4.063, 3.636, 3.064, 2.416, 1.589, 1.021, 0.640 ])