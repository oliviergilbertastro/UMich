import astropy.units as u

factor = 500*u.m**2 / (1.38E-23*(u.J/u.K))
print(factor.to(u.K/u.Jy))