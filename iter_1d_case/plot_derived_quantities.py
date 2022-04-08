import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass

data_mu = np.genfromtxt('chemical_pot_conservation/derived_quantities.csv', delimiter=',', names=True)
data_c = np.genfromtxt('concentration_continuity/derived_quantities.csv', delimiter=',', names=True)

# plot inventory
plt.figure()
ret_mu = data_mu['Total_retention_volume_8'] + data_mu['Total_retention_volume_7'] + data_mu['Total_retention_volume_6']
ret_c = data_c['Total_retention_volume_8'] + data_c['Total_retention_volume_7'] + data_c['Total_retention_volume_6']
plt.plot(data_c['ts'], ret_c, label=r"continuity of $c_\mathrm{m}$")
plt.plot(data_mu['ts'], ret_mu, label="continuity of $\mu$", linestyle="dashed")

plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Inventory (H m$^{-2}$)")
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
plt.show()

# plot relative diff
plt.figure()

ret_c_interp = interp1d(data_c["ts"], ret_c)
ret_mu_interp = interp1d(data_mu["ts"], ret_mu)

diff = abs(ret_c_interp(data_c["ts"]) - ret_mu_interp(data_c["ts"])) /  ret_mu_interp(data_c["ts"])
plt.plot(data_c["ts"], diff)
plt.yscale("log")
plt.xscale("log")
plt.ylabel("Absolute relative difference")
plt.xlabel("Time (s)")
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
plt.show()

# plot inventories and diff

fig, axs = plt.subplots(2, 1, sharex="col")

max_t = 2.4e7

axs[0].plot(data_c['ts'][data_c["ts"] < max_t], ret_c[data_c["ts"] < max_t], label=r"continuity of $c_\mathrm{m}$")
axs[0].plot(data_mu['ts'][data_mu["ts"] < max_t], ret_mu[data_mu["ts"] < max_t], label="continuity of $\mu$", linestyle="dashed")
axs[0].set_ylabel("Inventory (H m$^{-2}$)")
axs[0].legend()
axs[1].plot(data_c["ts"][data_c["ts"] < max_t], diff[data_c["ts"] < max_t])
axs[1].set_ylabel("Absolute relative difference")
for ax in axs:
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.minorticks_on()
    ax.grid(which='minor', alpha=0.3)
    ax.grid(which='major', alpha=0.7)
axs[1].set_xlabel("Time (s)")
plt.show()

# plot flux
plt.figure()
plt.plot(data_c['ts'], -data_c['Flux_surface_2_solute'], label=r"continuity of $c_\mathrm{m}$")
# plt.plot(data_mu['ts'], -data_mu['Flux_surface_2_solute'], label="continuity of $\mu$", linestyle="dashed")
# plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Outgassing flux at coolant (H m$^{-2}$ s$^{-1}$)")
plt.ylim(bottom=-0.1e11)
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
plt.show()
