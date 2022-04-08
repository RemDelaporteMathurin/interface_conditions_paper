import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass

field = "retention"
chemical_pot_conservation = False
if chemical_pot_conservation:
    folder = "profiles_"+field+"_mu_continuity"
    folder_t = "chemical_pot_conservation"
else:
    folder = "profiles_"+field+"_c_continuity"
    folder_t = "concentration_continuity"


def create_2d_data(folder, num_steps, header="0"):
    data_2d = []
    for i in range(num_steps):
        data_1d = np.genfromtxt(folder + '/{:.0f}.csv'.format(i), delimiter=',', names=True)
        data_2d.append(data_1d[header])
    return np.array(data_2d)


x = np.genfromtxt(folder + '/{:.0f}.csv'.format(0), delimiter=',', names=True)['arc_length']
t = np.genfromtxt(folder_t+'/derived_quantities.csv', delimiter=',', names=True)["ts"][:-1]
num_steps = len(t)

if field == "solute":
    header = "0"
    minimum_value = 1e20
    maximum_value = 2e24
elif field == "retention":
    header = "retention"
    minimum_value = 1e22
    maximum_value = 4e26

data_2d = create_2d_data(folder, num_steps, header)

data_2d[data_2d < minimum_value] = minimum_value
levels = np.logspace(
        np.log10(minimum_value),
        np.log10(maximum_value),
        1000)
locator = ticker.LogLocator(base=10)

xx, tt = np.meshgrid(x, t)
cs = plt.contourf(tt, xx, data_2d, levels=levels, locator=locator, vmin=minimum_value, extend="min")
cb = plt.colorbar(cs, ticks=locator)
if field == "0":
    cb.set_label('$c_\mathrm{m}$ (m$^{-3}$)')
elif field == "retention":
    cb.set_label('retention (m$^{-3}$)')

plt.gca().invert_yaxis()
plt.xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('x (m)')
for c in cs.collections:  # for avoiding white lines in pdf
    c.set_edgecolor("face")
plt.show()
