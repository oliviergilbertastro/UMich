import astropy.units as u
import astropy.constants as cst
import numpy as np

T_B = 500*u.m**2 * 1*u.Jy / cst.k_B

print(T_B.to(u.K))

T_B = np.pi*u.m**2 * 1*u.Jy / cst.k_B
print(T_B.to(u.K))