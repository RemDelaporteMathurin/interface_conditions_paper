import cu_eurofer
import FESTIM
import numpy as np

from solve_analytical_chemical_pot import D1, D2, S1, S2, a1, res, u_0, u_L, x_L, x_int, f

folder = "cu_eurofer/equivalent_D_eff"
parameters = cu_eurofer.parameters
del parameters["exports"]["xdmf"]

parameters["exports"]["derived_quantities"]["folder"] = folder
T = parameters["temperature"]["value"]
size = parameters["mesh_parameters"]["size"]
D_eff = D1*res[a1]*x_L/(u_L - u_0)

D1_simu = parameters["materials"][0]["D_0"]*np.exp(-parameters["materials"][0]["E_D"]/FESTIM.k_B/T)
D2_simu = parameters["materials"][1]["D_0"]*np.exp(-parameters["materials"][1]["E_D"]/FESTIM.k_B/T)
S1_simu = parameters["materials"][0]["S_0"]*np.exp(-parameters["materials"][0]["E_S"]/FESTIM.k_B/T)
S2_simu = parameters["materials"][1]["S_0"]*np.exp(-parameters["materials"][1]["E_S"]/FESTIM.k_B/T)
u_0_simu = parameters["boundary_conditions"][0]["value"]
u_L_simu = parameters["boundary_conditions"][1]["value"]
print(S1_simu, S2_simu)
x_interface = parameters["materials"][0]["borders"][1]
D_eff = D_eff.subs(x_L, size).subs(D1, D1_simu).subs(D2, D2_simu)
D_eff = D_eff.subs(x_int, x_interface).subs(S1, S1_simu).subs(S2, S2_simu)
D_eff = D_eff.subs(u_0, u_0_simu).subs(u_L, u_L_simu)
D_eff = D_eff.subs(f, 0)

parameters["materials"] = [
    {
        "D_0": float(D_eff),
        "E_D": 0,
        "borders": [0, size],
        "id": 1
    }
]
FESTIM.generic_simulation.run(parameters)
