import numpy as np
import matplotlib.pyplot as plt

def flux(mag, zp):
    return 10**((mag+48.958+zp)/(-2.5))

#         U      B      V       R      I     J     H    K
wavs = [0.366, 0.438, 0.545, 0.641, 0.798, 1.22, 1.63, 2.19] # microns