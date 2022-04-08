import matplotlib.pyplot as plt
import numpy as np

# try:
#     plt.rc('text', usetex=True)
#     plt.rc('font', family='serif', size=12)
# except:
#     pass
k_B = 8.6e-5

D_0_W = 2.9e-7*0.8165
E_D_W = 0.39


def D_W(T, main=False):
    return D_0_W*np.exp(-E_D_W/k_B/T)


D_0_Cu = 6.6e-7
E_D_Cu = 0.387


def D_Cu(T, main=False):
    return D_0_Cu*np.exp(-E_D_Cu/k_B/T)


D_0_CuCrZr = 3.92e-7
E_D_CuCrZr = 0.418


def D_CuCrZr(T, main=False):
    return D_0_CuCrZr*np.exp(-E_D_CuCrZr/k_B/T)


def D_Eurofer(T, main=False):
    return 1.5e-7*np.exp(-0.15/k_B/T)


S_0_W = 1.87e24
E_S_W = 1.04


def S_W(T, main=False):
    return S_0_W*np.exp(-E_S_W/k_B/T)


S_0_Cu = 3.14e24
E_S_Cu = 0.572


def S_Cu(T, main=False):
    return S_0_Cu*np.exp(-E_S_Cu/k_B/T)


S_0_CuCrZr = 4.28e23
E_S_CuCrZr = 0.387


def S_CuCrZr(T, main=False):
    return S_0_CuCrZr*np.exp(-E_S_CuCrZr/k_B/T)


def S_Eurofer(T, main=False):
    return 6.14e20*np.exp(-0.25/k_B/T)


if __name__ == "__main__":
    def tick_function(X):
        V = 1000/(X)
        return ["%.0f" % z for z in V]

    T = np.arange(300, 1300, step=1)
    T_Cu = np.arange(470, 1200, step=1)
    T_CuCrZr = np.arange(561, 769, step=1)
    T_Eurofer = np.arange(300, 1300, step=1)
    # T_W = np.arange()
    grey_W = "tab:grey"
    orange_Cu = (228/255, 146/255, 64/255)
    yellow_CuCrZr = (180/255, 95/255, 6/255)
    blue_eurofer = "tab:blue"

    fig, axs = plt.subplots(2, 1, sharex="col")
    ax1 = axs[0]
    ax1_1 = ax1.twiny()
    ax1.plot(1000/T, D_W(T, main=True), label="W", color=grey_W)
    ax1.plot(1000/T_Cu, D_Cu(T_Cu, main=True), label="Cu", color=orange_Cu)
    ax1.plot(1000/T_CuCrZr, D_CuCrZr(T_CuCrZr, main=True), label="CuCrZr", color=yellow_CuCrZr)
    ax1.plot(1000/T_Eurofer, D_Eurofer(T_Eurofer, main=True), label="Eurofer", color=blue_eurofer)

    new_tick_locations = np.array([1, 1.5, 2, 2.5, 3])
    ax1_1.set_xlim(ax1.get_xlim())
    ax1_1.set_xticks(new_tick_locations)
    ax1_1.set_xticklabels(tick_function(new_tick_locations))
    ax1.set_yscale("log")
    ax1.minorticks_on()
    ax1.grid(which='minor', alpha=0.3)
    ax1.grid(which='major', alpha=0.7)
    # ax1.legend()
    ax1.set_ylabel(r"$D$ (m$^{2}$ s$^{-1}$)", fontsize=15)
    ax1_1.set_xlabel(r"T (K)", fontsize=15)

    ax2 = axs[1]
    ax2.plot(1000/T, S_W(T, main=True), label="W", color=grey_W)
    ax2.plot(1000/T_Cu, S_Cu(T_Cu, main=True), label="Cu", color=orange_Cu)
    ax2.plot(1000/T_CuCrZr, S_CuCrZr(T_CuCrZr, main=True), label="CuCrZr", color=yellow_CuCrZr)
    ax2.plot(1000/T_Eurofer, S_Eurofer(T_Eurofer, main=True), label="Eurofer", color=blue_eurofer)

    ax2.annotate("W", (1000/T[0], S_W(T[0], main=True)), color=grey_W)
    ax2.annotate("Cu", (1000/T_Cu[0], S_Cu(T_Cu[0], main=True)), color=orange_Cu)
    ax2.annotate("CuCrZr", (1000/T_CuCrZr[0], S_CuCrZr(T_CuCrZr[0], main=True)), color=yellow_CuCrZr)
    ax2.annotate("EUROFER", (1000/T_Eurofer[0], S_Eurofer(T_Eurofer[0], main=True)), color=blue_eurofer)

    ax1.annotate("W", (1000/T[0], D_W(T[0], main=True)), color=grey_W)
    ax1.annotate("Cu", (1000/T_Cu[0], D_Cu(T_Cu[0], main=True)), color=orange_Cu)
    ax1.annotate("CuCrZr", (1000/T_CuCrZr[0] - 0.5, 0.2*D_CuCrZr(T_CuCrZr[0], main=True)), color=yellow_CuCrZr)
    ax1.annotate("EUROFER", (1000/T_Eurofer[0], D_Eurofer(T_Eurofer[0], main=True)), color=blue_eurofer)

    ax2.set_yscale("log")
    ax2.minorticks_on()
    ax2.grid(which='minor', alpha=0.3)
    ax2.grid(which='major', alpha=0.7)
    # ax2.legend()
    ax2.set_xlim(right=3.75)
    ax2.set_ylabel(r"$S$ (m$^{-3}$ Pa$^{-0.5}$)", fontsize=15)
    ax2.set_xlabel(r"1000/T (K$^{-1}$)", fontsize=15)
    plt.tight_layout()
    plt.show()