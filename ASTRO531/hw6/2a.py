import astropy.constants as cst
import astropy.units as u
import numpy as np
m_p = cst.m_p

S_0 = 3.78E-22

prefactor = 2.979E32

r_pp_const = prefactor * S_0 #/ m_p**2

print(r_pp_const/1E11)

rho_H = 50*(u.g/u.cm**3)
X_H = 0.5
n_H = rho_H*X_H/m_p
n_H = n_H.to(1/u.cm**3)
print(n_H.to(1/u.cm**3))
T_6 = 15.
r_pp = 1.1261E11*rho_H**2* X_H**2*(T_6)**(-2/3)*np.exp(-33.81*T_6**(-1/3)) *(u.cm**-3*u.s**-1)
print(r_pp.value)

t_pp = n_H.value/(2*r_pp.value)*u.s
print(t_pp.to(u.Gyr))