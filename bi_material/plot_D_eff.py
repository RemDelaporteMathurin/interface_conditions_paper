import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.stats import linregress

try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass


def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


size = 4e-3
c_0 = 1e20
Ts = np.linspace(300, 700, num=10)
diffusion_coeffs = []
for T in Ts:
    data = np.genfromtxt('flux_at_several_temperatures/derived_quantities_{:.0f}.csv'.format(T), delimiter=',', names=True)
    val_flux = -data['Flux_surface_2_solute']
    diffusion_coeffs.append(val_flux*size/c_0)
diffusion_coeffs = np.array(diffusion_coeffs)
plt.scatter(1/Ts, diffusion_coeffs)
plt.yscale("log")
slope, intercept, r_value, p_value, std_err = linregress(1/Ts, np.log10(diffusion_coeffs))
k_B = 8.6e-5
D_0 = 10**intercept
E_D = -slope*np.log(10)*k_B
plt.plot(1/Ts, D_0*np.exp(-E_D/k_B/Ts), linestyle="dashed", label=r"${} \exp(-{:.2f}/(k_B T))$".format(latex_float(D_0), E_D))

D_0_W = 2.9e-7*0.8165
E_D_W = 0.39
plt.plot(1/Ts, D_0_W*np.exp(-E_D_W/k_B/Ts), color="tab:grey", label="$D_\mathrm{W}$")

D_0_Cu = 6.6e-7
E_D_Cu = 0.387
plt.plot(1/Ts, D_0_Cu*np.exp(-E_D_Cu/k_B/Ts), color="tab:orange", label="$D_\mathrm{Cu}$")

plt.legend()
plt.xlabel(r"1/T (K$^{-1})$")
plt.ylabel(r"$D\mathrm{eff}$ (m$^2$ s$^{-1}$)")
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)



plt.show()
