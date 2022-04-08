from solve_analytical_chemical_pot import D1, D2, S1, S2, a1, res, u_0, u_L, x_L, x_int, f

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker


flux = -D1*res[a1]
L = 1
flux = flux.subs(x_L, L)
c_0 = 1
flux = flux.subs(u_0, c_0).subs(u_L, 0)
flux = flux.subs(f, 0)

D_1, S_1 = 1, 1

flux = flux.subs(S1, S_1).subs(D1, D_1)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)
fig.suptitle("Absolute relative difference on outgassing flux")
axs = [ax1, ax2, ax3, ax4]
for x_int_val, ax in zip([0.2, 0.4, 0.6, 0.8], axs):

    diffs = []
    alphas = np.logspace(-5, 5, num=30)
    betas = np.logspace(-5, 5, num=100)
    i = 0
    for alpha in alphas:
        diffs.append([])
        D_2 = alpha*D_1
        for beta in betas:
            flux_mu = flux.subs(S2, beta*S_1).subs(D2, D_2).subs(x_int, x_int_val*L)
            flux_c = flux.subs(S2, S_1).subs(D2, D_2).subs(x_int, x_int_val*L)
            diff_val = abs(flux_mu - flux_c)/flux_mu
            diffs[-1].append(diff_val)
            print(i, end="\r")
            i +=1

    ax.set_title(r"$x_\mathrm{int}/L = $" + "{:.1f}".format(x_int_val))

    diffs = np.array(diffs)

    levels = np.logspace(
        -6,
        5,
        100)
    locator = ticker.LogLocator()
    X, Y = np.meshgrid(betas, alphas)
    CS = ax.contourf(X, Y, diffs, levels=levels, locator=locator)
    cbar = fig.colorbar(CS, ticks=locator, ax=ax)
    formatter = ticker.LogFormatter(10, labelOnlyBase=False, minor_thresholds=(np.inf, np.inf))
    cbar.formatter = formatter
    ax.set_xscale("log")
    ax.set_yscale("log")
    for c in CS.collections:
        c.set_edgecolor("face")
ax3.set_xlabel(r"$S_2/S_1$")
ax4.set_xlabel(r"$S_2/S_1$")
ax1.set_ylabel(r"$D_2/D_1$")
ax3.set_ylabel(r"$D_2/D_1$")



plt.show()
