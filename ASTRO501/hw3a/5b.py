import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 

phi = np.deg2rad(40)
HA = np.linspace(-6, 6, 200) * 15 * np.pi/180  # radians

# example baseline (replace with real ones)
baselines = [()] # east, north
Be = 100   # east component (m)
Bn = 50    # north component (m)

u = Be*np.cos(HA) - Bn*np.sin(phi)*np.sin(HA)
v = Be*np.sin(HA) + Bn*np.sin(phi)*np.cos(HA)

plt.plot(u, v)

# repeat for all baselines...

plt.xlabel('u')
plt.ylabel('v')
plt.axis('equal')
plt.show()