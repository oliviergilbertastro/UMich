import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
import numpy as np
import pandas as pd
from urllib.parse import urlencode

ra, dec = 166.1138,	38.2088

params = {
    "POS": f"CIRCLE {ra} {dec} 0.001",
    "FORMAT": "csv"
}

url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?" + urlencode(params)
print(url)
df = pd.read_csv(url)
df.to_csv("ASTRO501/seminar/ztf_lc.csv")
print(df.head())