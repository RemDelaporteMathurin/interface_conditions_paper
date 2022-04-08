import bi_material_flux
import FESTIM
import numpy as np

print('check that chemcial pot is True!')
folder = "bi_material/flux_at_several_temperatures"
parameters = bi_material_flux.parameters
del parameters["exports"]["xdmf"]
parameters["solving_parameters"]["type"] = "solve_stationary"
parameters["exports"]["derived_quantities"]["folder"] = folder

for T in np.linspace(300, 700, num=10):
    parameters["temperature"]["value"] = T
    parameters["exports"]["derived_quantities"]["file"] = "derived_quantities_{:.0f}.csv".format(T)

    FESTIM.generic_simulation.run(parameters)
