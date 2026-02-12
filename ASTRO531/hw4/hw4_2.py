import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
from astropy.table import Table
import astropy.constants as cst
import astropy.units as u

t = Table.read(
    r"ASTRO531/hw4/table7.txt",
    format="ascii.fixed_width",
    header_start=3,
    data_start=5,
    names= ["Sol SDSS","Type","MJD","Plate","Fiber","Teff","e_Teff","logg","e_logg","Mwd","e_Mwd","dwd","e_dwd","MType","dsec","e_dsec","Flag","Notes"],
    data_end=-1
)
t = t[t['Mwd'] != 0.00]

Ms = t["Mwd"]*u.M_sun
print(Ms)
Rs = np.sqrt(cst.G*Ms/(10**t["logg"]*(u.cm/u.s**2)))
Rs = Rs.to(u.R_sun)
x,y = np.log10(Rs.value), np.log10(Ms.value)
mask = ~np.isnan(x) & ~np.isnan(y) & (x<-2)
x_fit = x[mask]
y_fit = y[mask]
from scipy.optimize import curve_fit
res = curve_fit(lambda x,a: -x/3+a, x_fit,y_fit)[0]
x_sim = np.linspace(min(x),max(x),100)
y_sim = -x_sim/3+res[0]
plt.figure()
plt.plot(x,y, color="black", lw=2, ls="None", marker=".", label="Data from Rebassa-Mansergas et al.")
plt.plot(x_sim,y_sim, color="blue", lw=2, ls="-", marker="None", label="Fit with -1/3 slope")
plt.axhline(np.log10(1.44), lw=2, ls="--", marker="None", color="red", label="Chandrasekhar Mass")
plt.xlabel(r"$\log R$ [$R_\odot$]")
plt.ylabel(r"$\log M$ [$M_\odot$]")
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig(r"ASTRO531/hw4/hw4_2.pdf")
plt.show()