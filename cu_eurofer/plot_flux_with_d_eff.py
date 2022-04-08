import numpy as np
import matplotlib.pyplot as plt
try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass

data_mu = np.genfromtxt('chemical_pot_conservation/derived_quantities.csv', delimiter=',', names=True)
data_D_eff = np.genfromtxt('equivalent_D_eff/derived_quantities.csv', delimiter=',', names=True)
data_c = np.genfromtxt('concentration_continuity/derived_quantities.csv', delimiter=',', names=True)

# plot 3 curves
plt.figure()
plt.plot(data_c['ts'], -data_c['Flux_surface_2_solute'], label=r"continuity of $c_\mathrm{m}$")
plt.plot(data_mu['ts'], -data_mu['Flux_surface_2_solute'], label="continuity of $\mu$", linestyle="dashed")
plt.plot(data_D_eff['ts'], -data_D_eff['Flux_surface_2_solute'], label=r"with effective D")
plt.xscale("log")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Outgassing flux (H m$^{-2}$ s$^{-1}$)")
# plt.ylim(bottom=-0.1e11)
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
plt.show()
