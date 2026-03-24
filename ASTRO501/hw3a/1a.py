import astropy.units as u

factor = 500*u.m**2 / (1.38E-23*(u.J/u.K))

print(factor.to(u.K/u.Jy))

F_nu = 1.0*u.mJy
T_A = (factor*F_nu).to(u.K)
print(T_A)

# A_eff*F_nu = k*T_B*A_eff*DeltaOmega/lambda**2
# T_B = F_nu*lambda**2/(k*DeltaOmega)
T_B = 1.0*(u.mJy)*(1*u.cm)**2/(((3*u.deg**2).to(u.sr))*(1.38E-23*(u.J/u.K)))
print(T_B.to(u.K/u.sr))

theta_beam = (1*u.cm)**2/(500*u.m**2)
print(theta_beam.decompose())

import numpy as np

img_fov = np.array([650,500])*5*u.arcsec
img_fov= img_fov.to(u.degree)
print(img_fov)