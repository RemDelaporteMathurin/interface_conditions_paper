import FESTIM

chemical_pot_conservation = False
if chemical_pot_conservation:
    folder = "cu_eurofer/chemical_pot_conservation"
else:
    folder = "cu_eurofer/concentration_continuity"

size = 4e-3
parameters = {
    "mesh_parameters": {
        "size": size,
        "initial_number_of_cells": 1000,
        "refinements": [
        ]

        },
    "materials": [
        {
            # Cu
            "D_0": 6.6e-7,  # OK TMAP eq(9)
            "E_D": 0.387,  # OK TMAP eq(9)
            "S_0": 3.14e24,  # OK TMAP eq(8) at/m3.Pa0.5
            "E_S": 0.572,  # OK TMAP eq(8) eV
            "borders": [0, size/2],  # OK TMAP
            "id": 1,
        },
        {
            # Eurofer
            "D_0": 1.5e-7,  # Aiello
            "E_D": 0.149,  # Aiello
            "S_0": 1.02e-3*6.022e23,  # Aiello at/m3.Pa0.5
            "E_S": 0.246,  # Aiello eV
            "borders": [size/2, size],  # OK TMAP
            "id": 2,
        },
        ],
    "traps": [
        ],
    "boundary_conditions": [
        {
           "type": "dc",
           "surfaces": 1,
           "value": 1e20
        },

        {
            "type": "dc",
            "surfaces": 2,
            "value": 0,
        },
        ],
    "temperature": {
        "type": "expression",
        "value": 500
        },
    "solving_parameters": {
        "final_time": 1e6,
        "initial_stepsize": 100,
        "adaptive_stepsize": {
            "stepsize_change_ratio": 1.2,
            "t_stop": 1e8,
            "stepsize_stop_max": 1e7,
            "dt_min": 1e-8,
            },
        "newton_solver": {
            "absolute_tolerance": 1e5,
            "relative_tolerance": 1e-10,
            "maximum_iterations": 20,
        },
        "traps_element_type": "DG"
        },
    "exports": {
        "xdmf": {
            "functions": ["solute"],
            "labels": ["solute"],
            "folder": folder,
            "checkpoint": False,
        },

        "derived_quantities": {

            "surface_flux": [
                {
                    "surfaces": [1, 2],
                    "field": "solute"
                }
            ],
            "file": "derived_quantities.csv",
            "folder": folder
        }
    }
}

if not chemical_pot_conservation:
    for mat in parameters["materials"]:
        del mat['S_0']
        del mat['E_S']
if __name__ == "__main__":
    FESTIM.generic_simulation.run(parameters)