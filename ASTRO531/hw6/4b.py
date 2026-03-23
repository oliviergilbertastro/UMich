import astropy.constants as cst
import astropy.units as u

m_p = cst.m_p
Delta_E = 26.73121*u.MeV
epsi_H = (Delta_E/(4*m_p)).to(u.erg/u.g)
print(epsi_H)

m_He = 4.002603254*u.u
Delta_E = 7.27425*u.MeV
epsi_He = (Delta_E/(3*m_He)).to(u.erg/u.g)
print(epsi_He)


time_H = epsi_H*0.2*cst.M_sun/cst.L_sun
print(time_H.to(u.Myr))

time_He = epsi_He*0.2*cst.M_sun/cst.L_sun
print(time_He.to(u.Myr))