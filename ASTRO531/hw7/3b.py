import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 

file = "ASTRO531/hw7/tableb2.dat"

df = pd.read_csv(file,delim_whitespace=True,comment='#',header=None)
mass = df[4] # already in solar masses
luminosity = df[10] # already in solar luminosities
print(df.head())

from scipy.optimize import curve_fit


mask = (mass > 0) & (luminosity > 0) & np.isfinite(mass) & np.isfinite(luminosity) & (np.log10(mass) > -0.2) 

mass_clean = mass[mask]
luminosity_clean = luminosity[mask]

lms = df[(np.log10(mass)<0.3) & (mask)]
ums = df[(np.log10(mass)>0.3) & (mask)]

def line_fit(x,a,b):
    return a*x+b
ums_fit = curve_fit(lambda x,b:line_fit(x,3,b),np.log10(ums[4]),(ums[10]))[0]
lms_fit = curve_fit(lambda x,b:line_fit(x,5.46,b),np.log10(lms[4]),(lms[10]))[0]

plt.figure(figsize=(7,7))
ax = plt.subplot(111)
x_fit = np.linspace(np.min(np.log10(mass_clean)),np.max(np.log10(mass_clean)),1000)
ax.plot(np.log10(mass_clean), (luminosity_clean), ls="None", marker=".", color="black", label="Malkov+2007")
ax.plot(x_fit, line_fit(x_fit,3,ums_fit[0]), ls="--", marker="None", color="red", label=r"$L\propto M^{3}$")
ax.plot(x_fit, line_fit(x_fit,5.46,lms_fit[0]), ls="--", marker="None", color="blue", label=r"$L\propto M^{5.46}$")
ax.set_ylabel(r"$\log L/L_\odot$")
ax.set_xlabel(r"$\log M/M_\odot$")
ax.legend(fontsize=16)
plt.tight_layout()
plt.savefig("ASTRO531/hw7/3b.pdf")
plt.show()