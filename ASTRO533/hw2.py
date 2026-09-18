import numpy as np
import astropy.constants as cst
from astropy import units as u


print((cst.h).to(u.erg * u.s))
print((cst.c).to(u.um/u.s))
print((cst.k_B).to(u.erg/u.K))

lambda_peak = (cst.h).to(u.erg * u.s)*(cst.c).to(u.um/u.s)/(cst.k_B).to(u.erg/u.K)/2.8
print(lambda_peak)
