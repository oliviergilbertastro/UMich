import numpy as np
import matplotlib.pyplot as plt

periods = [5,10,20] # days
t = np.linspace(0,30,1000) # days
offset=10 # days

stars = [np.cos(2*np.pi/periods[i]*t-2*np.pi*offset) for i in range(len(periods))]

for i in range(len(periods)):
    plt.plot(t, stars[i], label=f"Period={periods[i]}")
plt.axvline(offset, color="red", ls="--", lw=3)
plt.axvline(offset+10, color="red", ls="--", lw=3)
plt.xlabel("Time in days", fontsize=15)
plt.ylabel("Brightness")
plt.legend(fontsize=14)
plt.savefig("ASTRO102/lab_cepheids.pdf")
plt.show()