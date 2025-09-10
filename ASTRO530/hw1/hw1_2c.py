import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expn
from scipy.integrate import quad
from matplotlib import cm

sigma = 5.6704E-5 # [erg cm^-2 Hz K^-4]
sigma = 1
R_star2_T_eff_over_d2 = 1

def flux(cos_alpha):
    if cos_alpha >= 0:
        return R_star2_T_eff_over_d2*sigma*cos_alpha
    return 0

from mpl_toolkits.mplot3d import Axes3D

import numpy as np
theta = np.linspace(-np.pi/2, np.pi/2, 1000)
phi = np.linspace(0, 2*np.pi, 1000)
theta, phi = np.meshgrid(theta, phi)
r = 1

x = np.cos(theta)*np.cos(phi)*r
y = np.cos(theta)*np.sin(phi)*r
z = np.sin(theta)

cos_alpha = np.cos(phi)*np.cos(theta)
flux_vectorized = np.vectorize(flux)
print(cos_alpha.shape)


scalar_field = flux_vectorized(cos_alpha)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
cb = ax.plot_surface(x, y, z, facecolors=cm.inferno((scalar_field - np.min(scalar_field)) /
                          (np.max(scalar_field) - np.min(scalar_field))), cmap='inferno', alpha=0.8, ccount=100, rcount=100)

norm = plt.Normalize(vmin=scalar_field.min(), vmax=scalar_field.max())
colors = plt.cm.inferno(norm(scalar_field))

mappable = plt.cm.ScalarMappable(cmap="inferno", norm=norm)
mappable.set_array([])  # required for ScalarMappable without explicit data
cbar = fig.colorbar(mappable, ax=ax, shrink=0.8, label=r"$F(\alpha)$ [$R^2 T_\mathrm{eff}^4 \sigma / d^2$]")
cbar.set_label(
    r"$F(\alpha)$ [$R^2 T_\mathrm{eff}^4 \sigma / d^2$]", 
    fontsize=18
)

ax.set_xlabel(r"$x$ $[R_\mathrm{P}]$", fontsize=17)
ax.set_ylabel(r"$y$ $[R_\mathrm{P}]$", fontsize=17)
ax.set_zlabel(r"$z$ $[R_\mathrm{P}]$", fontsize=17)
ax.set_box_aspect((1, 1, 1))
plt.tight_layout()
plt.savefig(f"ASTRO530/hw1_2c1.pdf")
plt.show()


def temperature(cos_alpha): # Returns temperature in units
    if cos_alpha >= 0:
        return cos_alpha**(1/4)
    return 0

temperature_vectorized = np.vectorize(temperature)

scalar_field = temperature_vectorized(cos_alpha)
scalar_field = flux_vectorized(cos_alpha)**(1/4)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
cb = ax.plot_surface(x, y, z, facecolors=cm.inferno((scalar_field - np.min(scalar_field)) /
                          (np.max(scalar_field) - np.min(scalar_field))), cmap='inferno', alpha=1, ccount=20, rcount=20)

norm = plt.Normalize(vmin=scalar_field.min(), vmax=scalar_field.max())
colors = plt.cm.inferno(norm(scalar_field))

mappable = plt.cm.ScalarMappable(cmap="inferno", norm=norm)
mappable.set_array([])  # required for ScalarMappable without explicit data
cbar = fig.colorbar(mappable, ax=ax, shrink=0.8, label=r"$T_\mathrm{P,eff}(\alpha)$  [$\sqrt{\frac{R}{d}} T_\mathrm{eff}$]")
cbar.set_label(
    r"$T_\mathrm{P,eff}(\alpha)$  $\left[\sqrt{\frac{R}{d}} T_\mathrm{eff}\right]$", 
    fontsize=18
)

ax.set_xlabel(r"$x$ $[R_\mathrm{P}]$", fontsize=17)
ax.set_ylabel(r"$y$ $[R_\mathrm{P}]$", fontsize=17)
ax.set_zlabel(r"$z$ $[R_\mathrm{P}]$", fontsize=17)
ax.set_box_aspect((1, 1, 1))
plt.tight_layout()
plt.savefig(f"ASTRO530/hw1_2c.pdf")
plt.show()

fig = plt.figure(figsize=(7, 5))
alphas = np.linspace(0,np.pi,100)
plt.plot(alphas, [temperature(np.cos(alphas[i])) for i in range(len(alphas))])
plt.hlines(1/np.sqrt(2), 0, np.pi, ls="--", color="red", label=r"$T_\mathrm{P,eff}=\sqrt{\frac{R}{2d}} T_\mathrm{eff}$")
plt.ylabel(r"$T_\mathrm{P,eff}(\alpha)$  $\left[\sqrt{\frac{R}{d}} T_\mathrm{eff}\right]$", fontsize=17, color="black")
plt.xlabel(r"$\alpha$  $[\mathrm{rad}]$", fontsize=17)
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig(f"ASTRO530/hw1_2c3.pdf")
plt.show()