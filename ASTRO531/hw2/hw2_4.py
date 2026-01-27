import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 
import astropy.units as u
import astropy.constants as cst
import pandas as pd

df = pd.read_fwf(r"ASTRO531\hw2\Table_abundances.txt", colspecs=[
        (0, 12),   # Element
        (22, 25),  # Symbol
        (32, 35),  # Atomic Number
        (39, 45),  # Atomic Weight
        (48, 55),  # A09
        (57, 65),  # G98
    ], skiprows=4, skipfooter=4, names=["Element","Symbol","Atomic Number","Atomic Weight","A09","G98"], engine="python")
print(df)

# a) + b)
for comp in ["A09","G98"]:
    H_abundance = 10**((df[df["Symbol"] == "H"])[comp].values[0]) * (df[df["Symbol"] == "H"])["Atomic Weight"].values[0]
    He_abundance = 10**((df[df["Symbol"] == "He"])[comp].values[0]) * (df[df["Symbol"] == "He"])["Atomic Weight"].values[0]
    metal_abundance = np.sum(10**((df[(df["Symbol"] != "H") & (df["Symbol"] != "He")])[comp].values)* (df[(df["Symbol"] != "H") & (df["Symbol"] != "He")])["Atomic Weight"].values)
    total_abundance = H_abundance+He_abundance+metal_abundance

    X = H_abundance/total_abundance
    Y = He_abundance/total_abundance
    Z = metal_abundance/total_abundance
    Z_X = Z/X
    print(comp, X, Y, Z, Z_X)

    X_i = (10**(df[comp].values) * df["Atomic Weight"].values)/total_abundance
    inverse_mu = np.sum([X_i/df["Atomic Weight"].values * (1+df["Atomic Number"].values)])
    mu = 1/inverse_mu
    print(comp, mu)

    plt.plot(df["Atomic Number"].values, df[comp].values/(df[comp].values[0]), label=comp, ls="-", lw=3)
plt.xlabel("Atomic number")
plt.ylabel("Abundance/(Hydrogen abundance)")
plt.legend()
plt.tight_layout()
plt.savefig(r"ASTRO531\hw2\hw2_4c1.pdf")
plt.show()
for comp in ["A09","G98"]:
    H_abundance = 10**((df[df["Symbol"] == "H"])[comp].values[0]) * (df[df["Symbol"] == "H"])["Atomic Weight"].values[0]
    He_abundance = 10**((df[df["Symbol"] == "He"])[comp].values[0]) * (df[df["Symbol"] == "He"])["Atomic Weight"].values[0]
    metal_abundance = np.sum(10**((df[(df["Symbol"] != "H") & (df["Symbol"] != "He")])[comp].values)* (df[(df["Symbol"] != "H") & (df["Symbol"] != "He")])["Atomic Weight"].values)
    total_abundance = H_abundance+He_abundance+metal_abundance
    plt.plot(df["Atomic Weight"].values, df[comp].values/(df[comp].values[0]), label=comp, ls="-", lw=3)
plt.xlabel("Atomic weight")
plt.ylabel("Abundance/(Hydrogen abundance)")
plt.legend()
plt.tight_layout()
plt.savefig(r"ASTRO531\hw2\hw2_4c2.pdf")
plt.show()