import astropy.units as u

factor = (1.38E-23*(u.J/u.K))/(500*u.m**2)

print(factor.to(u.Jy/u.K))