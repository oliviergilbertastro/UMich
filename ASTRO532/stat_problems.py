import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def factorial(x:int) -> int:
    if x > 1:
        return x*factorial(x-1)
    return x

def prob_poisson(k:float, lam:float) -> float:
    if lam > 15:
        print(np.exp(-lam))
    return lam**k * np.exp(-lam) / (factorial(k))
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
lams = list(lams)
for lam in lams:
    plt.plot(ks, Ps[(lams.index(lam)),:], label=f"$\lambda$={lam}")
    plt.legend()
    plt.show()