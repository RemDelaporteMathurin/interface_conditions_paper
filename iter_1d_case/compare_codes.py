import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass

color_cyle = plt.rcParams['axes.prop_cycle'].by_key()['color']

HEADERS_CONCENTRATIONS_TMAP = \
    ['Mobile_atm3', 'Trap_1_atm3', 'Trap_2_atm3',
     'Trap_3_atm3', 'Trap_4_atm3', 'Trap_5_atm3']


def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


def read_tmap_results(t):
    data = np.genfromtxt('results_tmap/{:.1e}/profiles.txt'.format(t), delimiter="\t", names=True)
    return data


def read_retention_tmap(t):
    data = read_tmap_results(t)
    x = data['x_m'] * 1e3  # in mm
    retention = 0
    for header in HEADERS_CONCENTRATIONS_TMAP:
        retention += data[header]
    return x, retention


def read_solute_abaqus(t):
    data = np.genfromtxt('results_abaqus/abaqus{:.1e}.csv'.format(t), delimiter=",", names=True)
    x = data['xm'] * 1e3  # in mm
    solute = data["solutem3"]
    return x, solute


def read_retention_abaqus(t):
    data = np.genfromtxt('results_abaqus/abaqus{:.1e}.csv'.format(t), delimiter=",", names=True)
    x = data['xm'] * 1e3  # in mm
    retention = data["retentionm3"]
    return x, retention


def read_solute_tmap(t):
    data = read_tmap_results(t)
    x = data['x_m'] * 1e3  # in mm
    solute = data['Mobile_atm3']
    return x, solute


def read_retention_festim(t):
    data = np.genfromtxt('chemical_pot_conservation/retention_{:.1f}s.txt'.format(t), delimiter=" ")
    # remove bunny ears
    indexes = np.where(data[:, 1] > 0)
    x = data[:, 0][indexes] * 1e3  # in mm
    retention = data[:, 1][indexes]
    return x, retention


def read_solute_festim(t):
    data = np.genfromtxt('chemical_pot_conservation/0_{:.1f}s.txt'.format(t), delimiter=" ")
    x = data[:, 0] * 1e3  # in mm
    solute = data[:, 1]
    return x, solute


if __name__ == "__main__":
    piyears = np.pi*365.25*24*3600
    # times = [2.4e6, 2.4e7, piyears]
    times = [2.4e6, 2.4e7]
    y_labels = ["Solute concentration (m$^{-3}$)", "Retention (m$^{-3}$)"]
    y_lims = [[1e21, 3e24], [1e22, 4e26]]
    fig, (axs_top, axs_bottom) = plt.subplots(
        ncols=len(times), nrows=2, sharey='row', sharex='col', figsize=[6.4, 4.8])

    # # plot Abaqus
    # for axs, read_profile in zip([axs_top, axs_bottom], [read_solute_abaqus, read_retention_abaqus]):
    #     for ax, t in zip(axs, times):
    #         x, res = read_profile(t)
    #         ax.scatter(x, res, label="ABAQUS", marker='+', zorder=3, color=color_cyle[1])

    # plot TMAP
    for axs, read_profile in zip([axs_top, axs_bottom], [read_solute_tmap, read_retention_tmap]):
        for ax, t in zip(axs, times):
            x, res = read_profile(t)
            ax.scatter(x, res, label="TMAP7", marker='+', zorder=3, color=color_cyle[0])

    # plot FESTIM
    for axs, read_profile in zip([axs_top, axs_bottom], [read_solute_festim, read_retention_festim]):
        for ax, t in zip(axs, times):
            x, res = read_profile(t)
            ax.plot(x, res, label="FESTIM", zorder=3, color=color_cyle[1])

    # layout
    for axs, y_lim, y_label in zip([axs_top, axs_bottom], y_lims, y_labels):
        axs[0].set_ylim(bottom=y_lim[0], top=y_lim[1])
        axs[0].set_ylabel(y_label)
        axs[0].set_yscale('log')
        # for ax in axs:
        #     ax.minorticks_on()
        #     ax.grid(which='minor', alpha=0.3)
        #     ax.grid(which='major', alpha=0.7)

        h = y_lim[1] - y_lim[0]
        for ax in axs:
            rect1 = Rectangle((0, y_lim[0]), 6, h, facecolor="tab:grey", alpha=0.5, edgecolor="none")
            ax.add_patch(rect1)

            rect2 = Rectangle((6, y_lim[0]), 1, h, facecolor=(228/255, 146/255, 64/255), alpha=0.5, edgecolor="none")
            ax.add_patch(rect2)

            rect3 = Rectangle((7, y_lim[0]), 1.55, h, facecolor=(180/255, 95/255, 6/255), alpha=0.5, edgecolor="none")
            ax.add_patch(rect3)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

    for ax, t in zip(axs_top, times):
        ax.set_title("$t = {}$ s".format(latex_float(t)))
    for ax in axs_bottom:
        ax.set_xlabel("depth (mm)")
    plt.legend(bbox_to_anchor=(1.05, 1.3))
    plt.tight_layout()
    plt.show()