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