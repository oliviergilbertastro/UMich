import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
import pandas as pd
import mesa_web as m

rho_avg = (cst.M_sun/(4*np.pi*cst.R_sun**3/3)).to(u.g/u.cm**3)
print(rho_avg)

csts = 0.0184561759676
rho_c = rho_avg/csts
print(rho_c)

xi1 = 6.89685
a = cst.R_sun/xi1
print(a.to(u.cm))
print((a*xi1)/cst.R_sun)

K = a**2*np.pi *cst.G * rho_c**(2/3)
print(K.decompose())

polytrope = pd.read_csv(r"ASTRO531/hw4/Polytrope_n3.txt", header=1, delim_whitespace=True)
xi = polytrope["Xi"].values
theta = polytrope["Theta"].values

n=3
r = (xi*a).to(u.cm)
rho = (rho_c*theta**n).to(u.g/u.cm**3)
pressure = K*rho**(4/3)
pressure = pressure.to(u.dyn/u.cm**2)


try:
    table = m.read_profile(r"ASTRO531\hw2\MESA-Web_Job_01252661968\profile8.data", as_table=True)
except:
    table = m.read_profile(r"ASTRO531/hw2/MESA-Web_Job_01252661968/profile8.data", as_table=True)
print(table.columns)

plt.figure()
plt.plot(r, rho, color="black", ls="-", marker="None", label="Polytrope")
plt.plot((table["radius"]*u.R_sun).to(u.cm), 10**(table["logRho"]), color="red", ls="--", marker="None", label="MESA")
plt.xlabel(r"$r$ [$cm$]")
plt.ylabel(r"$\rho$ [$\mathrm{g/cm^3}$]")
plt.legend(fontsize=15)
#plt.xscale("log")
plt.tight_layout()
plt.savefig(r"ASTRO531/hw4/hw4_3b_rho.pdf")
plt.show()
plt.figure()
plt.plot(r, pressure, color="black", ls="-", marker="None", label="Polytrope")
plt.plot((table["radius"]*u.R_sun).to(u.cm), table["pressure"], color="red", ls="--", marker="None", label="MESA")
#plt.xscale("log")
plt.xlabel(r"$r$ [$cm$]")
plt.ylabel(r"$P$ [$\mathrm{dyn/cm^2}$]")
plt.legend(fontsize=15)
plt.tight_layout()
plt.savefig(r"ASTRO531/hw4/hw4_3b_P.pdf")
plt.show()