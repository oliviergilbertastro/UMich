import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import pandas as pd
import astropy.units as u
import astropy.constants as cst

constants = 4*np.pi* ( 1/8*(3/np.pi)**(1/3)*cst.h*cst.c/(np.pi*cst.G*u.u**(4/3)) )**(3/2)
xi2_dxi_dtheta_term = 2.01824
chandrashekar_lim = constants.to(u.M_sun)*xi2_dxi_dtheta_term
print(chandrashekar_lim)
print(chandrashekar_lim/4)
print(chandrashekar_lim/8)


polytrope = pd.read_csv(r"ASTRO531/hw4/Polytrope_n3.txt", header=1, delim_whitespace=True)
xi1 = polytrope["Xi"].values[-1]
dtheta_dxi = np.gradient(polytrope["Theta"], polytrope["Xi"])
dxi_dtheta = np.gradient(polytrope["Xi"], polytrope["Theta"])