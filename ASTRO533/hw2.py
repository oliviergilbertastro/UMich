import numpy as np
import astropy.constants as cst
from astropy import units as u


print((cst.h).to(u.erg * u.s))
print((cst.c).to(u.um/u.s))
print((cst.k_B).to(u.erg/u.K))

lambda_peak = (cst.h).to(u.erg * u.s)*(cst.c).to(u.um/u.s)/(cst.k_B).to(u.erg/u.K)/2.8
print(lambda_peak)


lam = 1216.67*u.AA
A_ul = 6E8/u.s

FWHM = A_ul/(2*np.pi)
FWHM_lam = FWHM*lam**2/cst.c
print(FWHM_lam.to(u.AA))
FWHM_kms = FWHM*lam
print(FWHM_kms.to(u.km/u.s))

b = 1.28*u.km/u.s * (6000/100)**(1/2)
print(b)
FWHM_lam = b*lam/cst.c
print(FWHM_lam.to(u.AA))