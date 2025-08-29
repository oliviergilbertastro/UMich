import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as spi


def simpson(func, bounds, N):
    '''
    Function which solves the integral within the bounds using the Simpson method
    bounds: tuple of the bounds [e.g. (a,b)]
    N: int representing the number of slices we want
    '''
    a, b = bounds
    h = (b-a)/N
    res = func(a)+func(b)
    for k in range(2,N,2):
      #Even
      res += 2*func(a+k*h)
    for k in range(1,N,2):
      #Odd
      res += 4*func(a+k*h)
    return h/3*res

def E2(x:float) -> float:
    return simpson(lambda t : np.exp(-x*t)/t**2, bounds=(1,np.inf), N=100)
print(E2(1))
def f_nu(tau_nu:float):
    pass