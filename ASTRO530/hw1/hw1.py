import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expn
from scipy.integrate import quad

def E2(x:float) -> float:
    return expn(2, x)

def S_nu(tau_nu:float, nu:float) -> float:
   """We can use Planck's function as the source function because of local thermodynamic equilibrium."""
   T = ((3/4)*T_eff**4*(tau_nu+2/3))**(1/4)
   #print(T, T_eff, tau_nu)
   return 2*h*nu**3/c**2 * 1/(np.exp(h*nu/(k*T))-1)

def f_nu(nu:float, tau_nu:float) -> float:
    #print(quad(lambda t_nu: S_nu(nu, t_nu)*E2(t_nu-tau_nu), a=tau_nu, b=np.inf)[0], quad(lambda t_nu: S_nu(nu, t_nu)*E2(tau_nu-t_nu), a=0, b=tau_nu)[0])
    return 2*np.pi*quad(lambda t_nu: S_nu(t_nu, nu)*E2(t_nu-tau_nu), a=tau_nu, b=np.inf)[0] - 2*np.pi*quad(lambda t_nu: S_nu(t_nu,nu)*E2(tau_nu-t_nu), a=0, b=tau_nu)[0]

# CONSTANTS in gaussian units:
h = 6.626*10**(-27) # [erg.s]
k = 1.380649*10**(-16) # [erg/K]
c = 2.99792458*10**(10) # [cm/s]

taus = [3,1,0.3, .1] # Grey atmosphere, so kappa_nu = kappa, so tau is independent of nu
T_eff = 5500 # [K]

wavs = np.logspace(-5, -3, 100) # in cm so that when converted to microns it goes from -1 to 1 in logspace
freqs = c/wavs

for tau in taus:
  Fs = []
  for f in freqs:
    Fs.append(f_nu(f, tau))
  plt.plot(np.log10(wavs*1e4), np.log10(np.array(Fs)*freqs), marker="None", ls="-", label=r"$\tau$"+f"$ = ${tau:.1f}")
plt.xlabel(r"$\log \lambda$ [$\mathrm{\mu m}$]", fontsize=15)
plt.ylabel(r"$\log \nu F_\nu$ [$\mathrm{erg \, cm^{-2} \, s^{-1} \, Hz^{-1}}$]", fontsize=15)
plt.legend()
plt.savefig("ASTRO530/hw1_a.pdf")
plt.show()

# Integrate F over frequency at each tau:

N_samples = np.linspace(0,1,100)*200
error_list = []
for N in N_samples:
  wavs = np.logspace(-5, -3, int(N)) # in cm so that when converted to microns it goes from -1 to 1 in logspace
  freqs = c/wavs # redefine frequencies

  sigma = 5.6704E-5 # [erg cm^-2 Hz K^-4]
  integrated_fluxes = []
  for tau in taus:
    Flux_nu = lambda f: f_nu(f, tau)
    Fs = [] # Sample the flux at different freqs
    for f in freqs:
      Fs.append(f_nu(f, tau))
    Flux = np.trapezoid(Fs[::-1], freqs[::-1])
    integrated_fluxes.append(Flux)
  error_array = np.abs(np.array(integrated_fluxes)-sigma*T_eff**4)/(sigma*T_eff**4)
  error = np.mean(error_array)
  error_list.append(error*100)
  print(r"accuracy (average % error)", error)
  if len(freqs) == 100: # Check to only make this plot for N=100
    plt.plot(taus, integrated_fluxes, ls="None", marker="o", color="black")
    plt.hlines(sigma*T_eff**4, 0, 3.1, linestyles="--", color="red", label=r"$\sigma T_\mathrm{eff}^4$")
    plt.fill_between([0, 3.1], y1=(sigma*T_eff**4)*0.99, y2=(sigma*T_eff**4)*1.01, color="red", alpha=0.5)
    plt.fill_between([0, 3.1], y1=(sigma*T_eff**4)*0.98, y2=(sigma*T_eff**4)*1.02, color="red", alpha=0.3)
    plt.fill_between([0, 3.1], y1=(sigma*T_eff**4)*0.97, y2=(sigma*T_eff**4)*1.03, color="red", alpha=0.1)
    plt.xlabel(r"$\tau$ [-]", fontsize=15)
    plt.ylabel(r"$F$ [$\mathrm{erg \, cm^{-2} \, s^{-1}}$]", fontsize=15)
    plt.legend(loc="upper right", fontsize=16)
    plt.savefig("ASTRO530/hw1_b.pdf")
    plt.show()

print(error_list)
plt.plot(N_samples, error_list, ls="-", marker=".", color="black")
plt.xlabel(r"$N$ number of samples", fontsize=15)
plt.ylabel(r"Error [%]", fontsize=15)
plt.savefig("ASTRO530/hw1_c.pdf")
plt.show()