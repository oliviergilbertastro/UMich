import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
from tqdm import tqdm

Chi_H = 13.6*u.eV
def saha(rho, T):
    """returns the fraction y^2/(1-y)"""
    return (cst.m_p/rho * (2*np.pi*cst.m_e*cst.k_B*T/(cst.h**2))**(3/2)*np.exp(-Chi_H/(cst.k_B*T))).decompose()

def y(C):
    return (-C+np.sqrt(C**2+4*C))/2
def D(y):
    return y*(1-y)/((2-y)*(1+y))
def Gamma3(D,T):
    return 1 + (2+2*D*(3/2+Chi_H/(cst.k_B*T)))/(3+2*D*(3/2+Chi_H/(cst.k_B*T))**2)

resolution = 100
log_rhos = np.linspace(-15,1,resolution)
rhos = (10**log_rhos)*(u.g/u.cm**3)
log_Ts = np.linspace(3.4,5,resolution)
T_s = (10**log_Ts)*(u.K)

Gammas = np.empty((resolution,resolution))
Ys = np.empty((resolution, resolution))
D_for_01 = D(0.01)
for x2,T in tqdm(enumerate(T_s)):
    for x1,p in enumerate(rhos):
        C = saha(p,T)
        y_ = y(C)
        D_ = D(y_)
        Gamma3_ = Gamma3(D_, T)
        Gammas[x2,x1] = Gamma3_
        Ys[x2,x1] = y_

plt.figure(figsize=(8,8))
img = plt.imshow(Gammas, cmap="inferno",origin="lower", extent=(-15,1,3.4,5), aspect="auto")
cb = plt.colorbar(img)
cb.set_label(r"$\Gamma_3$", fontsize=15)
plt.contour(Gammas,levels=[4/3],colors='white',linewidths=2,extent=(-15,1,3.4,5))
plt.contour(Ys,levels=[0.01, 0.99],colors=['blue', 'red'],linewidths=2,origin="lower",extent=(-15,1,3.4,5))
plt.xlabel(r"$\log \rho$ [$\mathrm{g/cm^3}$]")
plt.ylabel(r"$\log T$ [$\mathrm{K}$]")
plt.tight_layout()
plt.savefig(r"ASTRO531/hw4/hw4_4c.pdf")
plt.show()

