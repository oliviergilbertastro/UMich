import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst

def P(r):
    """returns P in units of P_c as long as r is in units of R"""
    return (1-(r)**2)

Rs = np.linspace(0,1,100)
plt.plot(Rs, P(Rs), ls="-", lw=3, color="black")
plt.xlabel(r"Radius [$r/R$]")
plt.ylabel(r"Pressure [$P/P_c$]")
plt.tight_layout()
plt.savefig(r"ASTRO531\hw2\hw2_1d.pdf")
plt.show()