from solve_analytical_chemical_pot import D1, D2, S1, S2, a1, res, u_0, u_L, x_L, x_int, f

import numpy as np
import matplotlib.pyplot as plt

k_B = 8.6e-5


def arhenius(D_0, E_D, T):
    return D_0*np.exp(-E_D/k_B/T)


flux = -D1*res[a1]
L = 1
flux = flux.subs(x_L, L)
c_0 = 1e20
flux = flux.subs(u_0, c_0).subs(u_L, 0)
flux = flux.subs(f, 0)

D_1_0 = 2.9e-7*0.8165
E_D_1 = 0.39

D_2_0 = 6.6e-7
E_D_2 = 0.397

S_1_0 = 1.87e24
E_S_1 = 1.04

S_2_0 = 3.14e24
E_S_2 = 0.572
fluxes = []
Ts = np.linspace(300, 900, num=50)
interfaces = np.linspace(0.1, 0.9, num=5)
for x_interface in interfaces:
    fluxes.append([])
    flux_ = flux.subs(x_int, x_interface*L)
    for T in Ts:
        D_1, D_2 = arhenius(D_1_0, E_D_1, T), arhenius(D_2_0, E_D_2, T)
        S_1, S_2 = arhenius(S_1_0, E_S_1, T), arhenius(S_2_0, E_S_2, T)
        flux_val = flux_.subs(S1, S_1).subs(S2, S_2).subs(D1, D_1).subs(D2, D_2)
        flux_val *= L/arhenius(D_1_0, E_D_1, 300)/c_0
        fluxes[-1].append(flux_val)

    plt.plot(Ts, fluxes[-1], label="$x_\mathrm{int}/L=" + "{:.1f}$".format(x_interface))

plt.legend()
plt.yscale('log')
plt.xlabel("$T$ (K)")
plt.ylabel(r"$\varphi_\mathrm{outgassing}\cdot L / (c_0 \cdot D(300 \mathrm{K})) $")


plt.figure()
fluxes = []
Ts = np.linspace(300, 900, num=5)
interfaces = np.linspace(0.1, 0.9, num=50)
for T in Ts:
    fluxes.append([])
    D_1, D_2 = arhenius(D_1_0, E_D_1, T), arhenius(D_2_0, E_D_2, T)
    S_1, S_2 = arhenius(S_1_0, E_S_1, T), arhenius(S_2_0, E_S_2, T)
    flux_ = flux.subs(S1, S_1).subs(S2, S_2).subs(D1, D_1).subs(D2, D_2)
    for x_interface in interfaces:

        flux_val = flux_.subs(x_int, x_interface*L)
        flux_val *= L/arhenius(D_1_0, E_D_1, 300)/c_0
        fluxes[-1].append(flux_val)

    plt.plot(interfaces, fluxes[-1], label="$T=" + "{:.0f}$ K".format(T))
plt.legend()
plt.yscale('log')
plt.xlabel("$x_\mathrm{int}/L$")
plt.ylabel(r"$\varphi_\mathrm{outgassing}\cdot L / (c_0 \cdot D(300 \mathrm{K})) $")
plt.show()
