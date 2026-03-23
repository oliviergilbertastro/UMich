import astropy.constants as cst
import astropy.units as u

m_p = cst.m_p

S_0 = 3.78E-22

prefactor = 2.979E32

r_pp_const = prefactor * S_0 #/ m_p**2

print(r_pp_const/1E11)