import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
from tqdm import tqdm
import mesa_web as m

try:
    table = m.read_profile(r"ASTRO531\hw2\MESA-Web_Job_01252661968\profile8.data", as_table=True)
except:
    table = m.read_profile(r"ASTRO531/hw2/MESA-Web_Job_01252661968/profile8.data", as_table=True)
print(table.columns)


rmo = table["opacity"]*(u.cm**2/u.g)
grad_a = table["grada"]
grad_R = table["gradr"]
grad_T = table["gradT"]
L = table["luminosity"]*u.L_sun
rho = 10**table["logRho"]*(u.g/u.cm**3)
R = table["radius"]*u.R_sun
M = table["mass"]*u.M_sun
P = table["pressure"]*(u.dyn/u.cm**2)
logT = (table["logT"])
v_conv = table["conv_vel"]*(u.cm/u.s)
gamma1 = table["gamma1"]
cs = (np.sqrt(gamma1*P/rho)).to(u.cm/u.s) # sound speed
T = 10**logT * u.K

g = (cst.G * M / R**2).to(u.cm / u.s**2) # gravity
H_P = P / (rho * g) # Pressure scale height
H_P = H_P.to(u.R_sun)


conv = grad_R > grad_a

l_ = np.mean(H_P[conv])
v_conv_ = np.mean(v_conv[conv])

timescale = (l_/v_conv_).to(u.s)
print(timescale)
timescale =timescale.to(u.hour)
print(timescale)