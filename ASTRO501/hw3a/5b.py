import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 

phi = np.deg2rad(40)
HA = np.linspace(-6, 6, 200) * 15 * np.pi/180  # radians

baselines = [(1,0.5),(5,2),(4,1.5),(5.5,1),(0.5,-1),(1.5,-0.5)] # east, north
for i,B in enumerate(baselines):
    Be = B[0]   # east component
    Bn = B[1]    # north component

    u = Be*np.cos(HA) - Bn*np.sin(phi)*np.sin(HA)
    v = Be*np.sin(HA) + Bn*np.sin(phi)*np.cos(HA)

    plt.plot(u, v, label=f"Baseline {i}")

plt.xlabel('u')
plt.ylabel('v')
plt.axis('equal')
plt.legend(fontsize=12)
plt.show()