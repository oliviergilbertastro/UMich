import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import expn

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
    return expn(2, x)

def S_nu(nu:float, tau_nu:float) -> float:
   """We can use Planck's function as the source function because of local thermodynamic equilibrium."""
   T = ((3/4)*T_eff**4*(tau_nu+2/3))**(-4)
   return 2*h*nu**3/c**2 * 1/(np.exp(h*nu/(k*T))-1)

def func_inside_int(t_nu, tau_nu):
   return S_nu(t_nu)*E2(t_nu-tau_nu)

def f_nu(nu:float, tau_nu:float) -> float:
    return 2*np.pi*simpson(lambda t_nu: func_inside_int(t_nu, tau_nu), bounds=(tau_nu,np.inf), N=1000) - \
           2*np.pi*simpson(lambda t_nu: func_inside_int(t_nu, tau_nu), bounds=(tau_nu,np.inf), N=1000)

# CONSTANTS in gaussian units:
h = 6.626*10**(-27) # [erg.s]
k = 1.380649*10**(-16) # [erg/K]
c = 2.99792458*10**(10) # [cm/s]

taus = [3,1,0.3]
T_eff = 5500 # [K]

wavs = np.logspace(1e-5, 1e-3, 100) # in cm so that when converted to microns it goes from -1 to 1 in logspace
freqs = c/wavs
print(wavs)
print(freqs)
plt.plot(np.log(wavs*1e4), freqs, marker=".", ls="None")
plt.show()
plt.plot(freqs, marker=".", ls="None")
plt.show()