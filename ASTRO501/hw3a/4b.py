import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
from scipy.optimize import brentq
from scipy.special import j1
import astropy.units as u

def f(x,V):
    """Func to find root """
    return j1(x)/x-V/2
x = brentq(f,0.01,10,args=0.3)
lam = 1.65*(u.um)
B = 300*(u.m)
theta = (x*lam/(np.pi*B)).decompose()

print(theta)