from astropy import units as u
from astropy import constants as cst
import numpy as np

M_sun = 2E33*u.g
M_particle = 8E-25*u.g
R_sun = cst.R_sun.to(u.cm)
print("R_sun:", R_sun)
V_sun = 4*np.pi*R_sun**3/3
print("V_sun:", V_sun)
N_particles = M_sun/M_particle
n = N_particles/V_sun
print("N_particles:", N_particles)
print("Number density:", n)
T = ((1/18)**2*n)**(1/3)
T = T.value *u.K
print("T:", T)