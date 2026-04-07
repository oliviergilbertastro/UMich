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



hist_05_Msun = m.read_history(r"ASTRO531/hw8/MESA-Web_Job_04052666820/trimmed_history.data", as_table=True)
hist_5_Msun = m.read_history(r"ASTRO531/hw8/MESA-Web_Job_04052666822/trimmed_history.data", as_table=True)
plt.plot(hist_05_Msun["model_number"],hist_05_Msun["center_he4"]) # profile 7 is when the model starts to change the He mass frac
plt.plot(hist_5_Msun["model_number"],hist_5_Msun["center_he4"]) # profile 8 is when the model starts to change the He mass frac
plt.show()