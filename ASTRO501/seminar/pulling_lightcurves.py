import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
import numpy as np


ra, dec = 166.1138,	38.2088

import pandas as pd

url = f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE {ra} {dec} 0.001&FORMAT=csv"
df = pd.read_csv(url)

print(df["objectid"].unique())