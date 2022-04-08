import numpy as np
import matplotlib.pyplot as plt
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


def plot_2_profiles(field, csvlabels, y_label, y_lims):
    fig, axs = plt.subplots(ncols=2, nrows=1, sharey=True, figsize=[6.4, 4.8])

    for t, ax in zip(times, axs):
        data = np.genfromtxt(field + "_profiles_t={:.1e}s.csv".format(t), delimiter=',', names=True)
        ax.plot(data["arc_length"], data[csvlabels[0]], label=r"continuity of $c_\mathrm{m}$", linestyle="-")
        ax.plot(data["arc_length"], data[csvlabels[1]], label="continuity of $\mu$", linestyle="dashed")
        ax.set_yscale("log")
        ax.set_title("$t = {}$ s".format(latex_float(t)))
        ax.set_ylim(bottom=y_lims[0], top=y_lims[1])
        ax.minorticks_on()
        ax.grid(which='minor', alpha=0.3)
        ax.grid(which='major', alpha=0.7)
        ax.set_xlabel("depth (m)")
    axs[0].set_ylabel(y_label)
    plt.legend()




times = [2.4e6, 2.4e7]
fields = ['solute', 'retention']
y_labels = ["Solute concentration (m$^{-3}$)", "Retention (m$^{-3}$)"]
y_lims = [[1e21, 1e24], [1e22, 1e26]]

# Plot one field at a time
for field, y_label, y_lim in zip(fields, y_labels, y_lims):
    plot_2_profiles(
        field, csvlabels=[field + '_c', field + '_mu'],
        y_label=y_label, y_lims=y_lim)

# plot 2 fields in one plot
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, sharey='row', sharex='col', figsize=[6.4, 4.8*1.5])
axs_top = (ax1, ax2)
axs_bottom = (ax3, ax4)
for i, t in enumerate(times):
    for field, ax in zip(fields, [axs_top[i], axs_bottom[i]]):
        data = np.genfromtxt(field + "_profiles_t={:.1e}s.csv".format(t), delimiter=',', names=True)
        csvlabels = [field + '_c', field + '_mu']
        ax.plot(data["arc_length"], data[csvlabels[0]], label=r"continuity of $c_\mathrm{m}$", linestyle="-")
        ax.plot(data["arc_length"], data[csvlabels[1]], label="continuity of $\mu$", linestyle="dashed")
        ax.set_yscale("log")
        ax.minorticks_on()
        ax.grid(which='minor', alpha=0.3)
        ax.grid(which='major', alpha=0.7)

for axs, y_lim, y_label in zip([axs_top, axs_bottom], y_lims, y_labels):
    axs[0].set_ylim(bottom=y_lim[0], top=y_lim[1])
    axs[0].set_ylabel(y_label)

for ax, t in zip(axs_top, times):
    ax.set_title("$t = {}$ s".format(latex_float(t)))
for ax in axs_bottom:
    ax.set_xlabel("depth (m)")
axs_top[0].legend()
plt.show()
