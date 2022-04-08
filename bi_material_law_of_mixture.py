from bi_material_flux import parameters
import numpy as np
import FESTIM
import matplotlib.pyplot as plt

atom_density_W = 6.28e28  # 6.3222e28  # atomic density m^-3
atom_density_Cu = 8.43e28  # 8.4912e28  # atomic density m^-3
NL_W = 6*atom_density_W
NL_Cu = 6*atom_density_Cu

T = parameters["temperature"]["value"]
W, Cu = parameters["materials"][0], parameters["materials"][1]
S_W = W["S_0"] * np.exp(-W["E_S"]/FESTIM.k_B/T)
S_Cu = Cu["S_0"] * np.exp(-Cu["E_S"]/FESTIM.k_B/T)

alphas = np.linspace(0, 1)
p0 = 1e5
S_Ws, S_Cus, ratios = [], [], []
for alpha in alphas:
    # S_W_mixture = (alpha/atom_density_W + (1-alpha)/(S_W*p0**0.5))**(-1)
    # S_Cu_mixture = (alpha/atom_density_Cu + (1-alpha)/(S_Cu*p0**0.5))**(-1)
    S_W_mixture = (1/NL_W*(alpha + (1-alpha)/(np.exp(-W["E_S"]/FESTIM.k_B/T))))**(-1)
    S_Cu_mixture = (1/NL_W*(alpha + (1-alpha)/(np.exp(-W["E_S"]/FESTIM.k_B/T))))**(-1)
    S_Ws.append(S_W_mixture)
    S_Cus.append(S_Cu_mixture)
    ratios.append(S_W_mixture/S_Cu_mixture)
# plt.plot(alphas, S_Ws)
# plt.plot(alphas, S_Cus)
plt.plot(alphas, ratios)
# plt.yscale('log')
# plt.xscale('log')
plt.savefig('out.png')
# print(S_W, S_W_mixture)
# print(S_Cu, S_Cu_mixture)

W["S_0"] = S_W_mixture
W["E_S"] = 0
Cu["S_0"] = S_Cu_mixture
Cu["E_S"] = 0

parameters["exports"]["derived_quantities"]["folder"] = \
    "bi_material_law_of_mixtures/{:.1f}".format(alpha)
del parameters["exports"]["xdmf"]

# FESTIM.generic_simulation.run(parameters)
