import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

t = Table.read(
    "table7.dat",
    format="ascii.fixed_width",
    header_start=0,
    data_start=4
)
