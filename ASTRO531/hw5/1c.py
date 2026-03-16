import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
plt.rc("axes", labelsize=16) 

def read_opacity_table(table_num):
    """
    Read a specific OP Rosseland opacity table.

    Parameters
    ----------
    filename : str
        Path to opacity table file
    table_num : int
        Table number (1-126)

    Returns
    -------
    logT : np.ndarray
    logR : np.ndarray
    kappa : np.ndarray
        2D array (logT, logR) of log opacity
    """

    start_line = 240 + (table_num - 1) * 77
    end_line = start_line + 77

    with open(r"ASTRO531/hw5/GS98.txt") as f:
        lines = f.readlines()[start_line:end_line]

    # logR header
    logr_line = lines[4]
    logR = np.array([float(x) for x in logr_line.split()[1:]])

    logT = []
    kappa = []

    for line in lines[6:]:
        parts = line.split()
        if len(parts) < 2:
            continue

        logT.append(float(parts[0]))
        row = []

        for v in parts[1:]:
            val = float(v)
            if val == 9.999:
                row.append(np.nan)
            else:
                row.append(val)
        for i in range(19-len(parts[1:])):
            row.append(np.nan)
        kappa.append(row)

    logT = np.array(logT)
    kappa = np.array(kappa)

    return (logT, logR, kappa)

table_1 = read_opacity_table(table_num=1)
table_73 = read_opacity_table(table_num=73)
table_13 = read_opacity_table(table_num=13)

R_indices = np.arange(0, len(table_73[1]), 2)

plt.figure(figsize=(10,7))
for i in R_indices:
    plt.plot(table_13[0], table_13[2][:,i], ls="--", label=rf"He + high Z $\log R$={table_13[1][i]:.1f}")
    last_col = plt.gca().lines[-1].get_color()
    plt.plot(table_73[0], table_73[2][:,i], color=last_col, label=f"Solar $\log R$={table_73[1][i]:.1f}")
    plt.plot(table_1[0], table_1[2][:,i], ls=":", color=last_col, label=f"He $\log R$={table_1[1][i]:.1f}")

plt.xlabel(r"$\log(T)$ [$\mathrm{10^6 K}$]")
plt.ylabel(r"$\log$ RMO [$\mathrm{cm^2/g}$]")
plt.legend(ncol=3, fontsize=10)
plt.tight_layout()
plt.savefig(r"ASTRO531/hw5/plot1c.pdf")
plt.show()