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

table_66 = read_opacity_table(table_num=66)
table_73 = read_opacity_table(table_num=73)
table_78 = read_opacity_table(table_num=78)

R_indices = np.arange(0, len(table_66[1]), 2)

plt.figure(figsize=(10,7))
for i in R_indices:
    plt.plot(table_66[0], table_66[2][:,i], ls="--", label=rf"Z=0 $\log R$={table_66[1][i]:.1f}")
    last_col = plt.gca().lines[-1].get_color()
    plt.plot(table_73[0], table_73[2][:,i], color=last_col, label=f"Z=0.02 $\log R$={table_73[1][i]:.1f}")
    plt.plot(table_78[0], table_78[2][:,i], ls=":", color=last_col, label=f"Z=0.1 $\log R$={table_78[1][i]:.1f}")

plt.xlabel(r"$\log(T)$ [$\mathrm{10^6 K}$]")
plt.ylabel(r"$\log$ RMO [$\mathrm{cm^2/g}$]")
plt.legend(ncol=3, fontsize=10)
plt.tight_layout()
plt.savefig(r"ASTRO531/hw5/plot1b.pdf")
plt.show()