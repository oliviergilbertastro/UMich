import numpy as np
import astropy.constants as cst
from astropy import units as u

R = 15*u.pc
V = (4*np.pi*R**3/3).to(u.cm**3)
print(V)

n_H2 = 300*(u.cm**-3)
N_H2 = n_H2*V
print(N_H2)

m_H2 = (2*u.u).to(u.g)
print(m_H2)

m_cloud = N_H2*m_H2
print(m_cloud)

m_all_clouds_in_MW = (1.5E9*(u.M_sun)).to(u.g)

N_clouds_in_MW = m_all_clouds_in_MW/m_cloud
print("N_clouds_in_MW:", N_clouds_in_MW)


dust_radius = 0.1*u.um
dust_radius = dust_radius.to(u.cm)
print(dust_radius)

dust_bulk_density = 2.2*(u.g/u.cm**3)

mass_dust_grain = 4*np.pi*dust_radius**3/3 * dust_bulk_density
print(mass_dust_grain)

mass_dust_cloud = 0.01*m_cloud

N_dust_grains = mass_dust_cloud/mass_dust_grain
print(N_dust_grains)