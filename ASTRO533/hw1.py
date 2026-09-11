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

n_dust_grains = N_dust_grains/V
print("n_dust_grains:",n_dust_grains)

sigma_dust = np.pi*dust_radius**2
print("Cross-section:", sigma_dust)

mfp = 1/(sigma_dust*n_dust_grains)
print("Mean free path:",mfp)

mass_metals_cloud = m_cloud*0.02
print(m_cloud)
print(mass_metals_cloud)


V_galaxy = np.pi*(15*u.kpc)**2*150*u.pc
V_galaxy = V_galaxy.to(u.cm**3)
print(V_galaxy)
print(V_galaxy.to(u.kpc**3))

n_clouds = N_clouds_in_MW/V_galaxy
n_clouds = n_clouds.to(u.kpc**-3)
print(n_clouds)

LOS_distance = 8*u.kpc
sigma_cloud = np.pi*(15*u.pc)**2
sigma_cloud = sigma_cloud.to(u.kpc**2)
print("sigma_cloud:", sigma_cloud)

mfp_cloud = 1/(sigma_cloud*n_clouds)
print(mfp_cloud)

expected_nb_of_clouds_in_LOS = LOS_distance/mfp_cloud
print(expected_nb_of_clouds_in_LOS)

print(np.exp(-expected_nb_of_clouds_in_LOS))


P = 4E-13*(u.dyn/u.cm**2)
print(P*3/2)

print((1*u.dyn).to(u.eV/u.cm))

print((P*3/2).to(u.eV/u.cm**3))