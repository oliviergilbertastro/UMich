import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
from scipy.optimize import brentq
from scipy.special import j1
import astropy.units as u

def get_theta(V):
    def f(x):
        """Func to find root """
        return j1(x)/x-V/2
    x = brentq(f,0.01,10)
    lam = 1.65*(u.um)
    B = 300*(u.m)
    theta = (x*lam/(np.pi*B)).decompose()
    return (theta*u.rad).to(u.mas)

V = 0.30
theta = get_theta(V)
theta_min = get_theta(V+0.03)
theta_max = get_theta(V-0.03)
print(theta, theta_min, theta_max,[theta-theta_min],[theta_max-theta])

for V in [0.33,0.30,0.27]:
    theta = get_theta(V).value
    theta_min = get_theta(V+0.03).value
    theta_max = get_theta(V-0.03).value
    plt.errorbar([theta], [V], xerr=([theta-theta_min],[theta_max-theta]), yerr=[0.03], marker="o",ls="None", label=rf"$V=${V}$\pm0.03$")
vs = np.linspace(0.24,0.36,100)
thetas = []
for v in vs:
    thetas.append(get_theta(v).value)
thetas = np.array(thetas)
print(thetas.shape)
plt.plot(thetas, vs, ls="-", color="black")
plt.xlabel(r"$\theta$ [$\mathrm{mas}$]")
plt.ylabel(r"Visibility Amplitude $V$")
plt.legend(fontsize=13)
plt.tight_layout()
plt.savefig(r"ASTRO501/hw3a/4b.pdf")
plt.show()
