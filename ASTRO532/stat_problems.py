import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def factorial(x:int) -> int:
    if x > 1:
        return x*factorial(x-1)
    return x

def prob_poisson(k:float, lam:float) -> float:
    return lam**k * np.exp(-lam) / (factorial(k))

def make_2d(x, y, z, per_pixel_limit=None, cmap=None, title=None, x_scale="linear", y_scale="linear"):

    fig, ax = plt.subplots(figsize=(6,5))
    im = plt.imshow(count_stat, cmap=("rainbow" if cmap is None else cmap), vmin=per_pixel_limit, vmax=np.nanmax(count_stat), origin="lower", aspect="auto", extent=[x_limits[0], x_limits[1], y_limits[0], y_limits[1]], alpha=None)
    
    cbar = plt.colorbar(im)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.set_label(r"P", fontsize=14, rotation=270)
    plt.xscale(x_scale)
    plt.yscale(y_scale)
    plt.tight_layout()
    plt.show()

N = 25

ks = np.arange(N)+1
lams = np.arange(N)+1

Ps = np.ndarray((N,N))
for k in ks:
    for lam in lams:
        Ps[lam-1,k-1] = prob_poisson(lam,k)

plt.imshow(Ps, origin="lower", extent=[1, N, 1, N], cmap="rainbow")
plt.xlabel(r"$k$", fontsize=16)
plt.ylabel(r"$\lambda$", fontsize=16)
plt.show()