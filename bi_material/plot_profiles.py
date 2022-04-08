import numpy as np
import matplotlib.pyplot as plt
try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass

data = np.genfromtxt('solute_profiles.csv', delimiter=',', names=True)
plt.plot(data["arc_length"], data['solute_c'], label=r"continuity of $c_\mathrm{m}$", linestyle="-")
plt.plot(data["arc_length"], data['solute_mu'], label="continuity of $\mu$", linestyle="dashed")
plt.plot([data["arc_length"][0], data["arc_length"][-1]], [data['solute_c'][0], data['solute_c'][-1]], label="with effective $D$")
# plt.ylim(bottom=y_lims[0], top=y_lims[1])
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
plt.xlabel("depth (m)")
plt.ylabel("Solute concentration (m$^{-3}$)")
plt.legend()
plt.show()