import astropy.units as u

factor = (1.38E-23*(u.J/u.K))/(500*u.m**2)

print(factor.to(u.Jy/u.K))

print(2.76/109544.51, 4.564E-4*2.76)
print(3*2.76/109544.51, 3*4.564E-4*2.76)