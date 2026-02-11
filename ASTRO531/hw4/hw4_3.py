import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
import pandas as pd

rho_avg = (cst.M_sun/(4*np.pi*cst.R_sun**3/3)).to(u.g/u.cm**3)
print(rho_avg)

csts = 0.0184561759676
rho_c = rho_avg/csts
print(rho_c)

xi1 = 6.89685
a = xi1/cst.R_sun
print(a.to(u.cm**-1))

K = a**2*np.pi *cst.G * rho_c**(2/3)
print(K.decompose())