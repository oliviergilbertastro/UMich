import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
from astropy.constants import c, L_sun
import mesa_web as m
from tqdm import tqdm

hist = m.read_history("", as_table=True)
table = m.read_history(r"ASTRO531\hw2\MESA-Web_Job_01252661968\profile1.data", as_table=True)
print(table.columns)


# a) Verify Virial theorem:
Us = []
for i in tqdm(range(1,290)):
    try:
        table_i = m.read_history(r"ASTRO531\hw2\MESA-Web_Job_01252661968\profile"+f"{i}"+".data", as_table=True)
        E_int = table_i["energy"] # as a cumulative function of radius
        U = E_int[-1]
        Us.append(U)
    except:
        pass

plt.plot(Us)
plt.show()
