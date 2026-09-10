import numpy as np
import astropy.constants as cst
from astropy import units as u

R = 15*u.pc
V = (4*np.pi*R**3/3).to(u.cm**3)
print(V)

n_H2 = 300*(u.cm**-3)
N_H2 = n_H2*V
print(N_H2)