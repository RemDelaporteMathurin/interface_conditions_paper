import FESTIM
import fenics
import sympy as sp
import matplotlib.pyplot as plt
import numpy as np

# Parameters
## Case 1
interface = 0.6
size = 1
S_01, S_02 = 1, 1.5
D_01, D_02 = 1, 2
E_D1, E_D2 = 0, 0
E_S1, E_S2 = 0, 0
c_0, c_L = 0.1, 0.2
f = 1/(size**2/(D_01*D_02)**0.5/c_0)

## Case 2
# interface = 0.5
# size = 1.25
# S_01, S_02 = 1, 0.5
# D_01, D_02 = 3, 2
# E_D1, E_D2 = 0, 0
# E_S1, E_S2 = 0, 0
# c_0, c_L = 0.4, 0.1

# f = 2/(size**2/(D_01*D_02)**0.5/c_0)
print(f*(size**2/(D_01*D_02)**0.5/c_0))

# # case 3
# interface = 0.6
# size = 1
# S_01, S_02 = 1.87e24, 3.14e24
# D_01, D_02 = 2.9e-7*0.8165, 6.6e-7
# E_D1, E_D2 = 0.39, 0.387
# E_S1, E_S2 = 1.04, 0.572
# f = 0
# c_0, c_L = 1e20, 0
T = 300

parameters = {
    "materials": [
        {
            "D_0": D_01,
            "E_D": E_D1,
            "S_0": S_01,
            "E_S": E_S1,
            "id": 1,
            "borders": [0, interface],
        },
        {
            "D_0": D_02,
            "E_D": E_D2,
            "S_0": S_02,
            "E_S": E_S2,
            "id": 2,
            "borders": [interface, size],
        },
    ],
    "source_term": [
        {
            "value": f,
            "volumes": [1, 2]
        },
    ],
    "temperature": {
        "type": "expression",
        "value": T
    },
    "traps": [],
    "boundary_conditions": [
        {
            "type": "dc",
            "value": c_0,
            "surfaces": [1]
        },
        {
            "type": "dc",
            "value": c_L,
            "surfaces": [2]
        },
    ],
    "mesh_parameters": {
        "initial_number_of_cells": 50,
        "size": size
    },
    "exports": {
        "xdmf": {
                "functions": ["solute"],
                "labels": ["solute"],
                "folder": "MES_chemical_pot",
                "checkpoint": False,
            },
        "error": {}
    },
    "solving_parameters": {
        "final_time": 1,
        "initial_stepsize": 0.05,
        "type": "solve_stationary",
        "newton_solver": {
            "absolute_tolerance": 1e-10,
            "relative_tolerance": 1e-10,
            "maximum_iterations": 15
        }
    }
}

out = FESTIM.generic_simulation.run(parameters)
u_sim = out["u"]

def arhenius(D_0, E_D, T):
    k_B = 8.6e-5
    return D_0*np.exp(-E_D/k_B/T)

D1, D2 = arhenius(D_01, E_D1, T), arhenius(D_02, E_D2, T)

S1, S2 = arhenius(S_01, E_S1, T), arhenius(S_02, E_S2, T)

x_int = interface
u0, uL = c_0, c_L
a1, a2, b1, b2 = sp.symbols('a1, a2, b1, b2')
x = sp.Symbol("x")
u1 = -1/2*f/D1*x**2 + a1*x + b1
u2 = -1/2*f/D2*x**2 + a2*x + b2
xL = size
x0 = 0

list_of_equations = [
    u1.subs(x, x_int)/S1 - u2.subs(x, x_int)/S2,
    D1*sp.diff(u1, x) - D2*sp.diff(u2, x),
    u1.subs(x, x0) - u0,
    u2.subs(x, xL) - uL
]
res = sp.solve(list_of_equations, a1, a2, b1, b2)

u1 = u1.subs(a1, res[a1]).subs(b1, res[b1])
u2 = u2.subs(a2, res[a2]).subs(b2, res[b2])

def u_exact(x):
    if x <= x_int:
        return u1.subs('x', x)
    elif x > x_int:
        return u2.subs('x', x)

error = 0
dx_ = size/parameters["mesh_parameters"]["initial_number_of_cells"]
diff = []
for x_val in np.arange(0, size, dx_):
    error += abs(u_exact(x_val) - u_sim(x_val))
    diff.append(float(u_exact(x_val) - u_sim(x_val)))

error = error/parameters["mesh_parameters"]["initial_number_of_cells"]
# print(np.linalg.norm(diff)*dx_)
print(error, dx_)
V = u_sim.function_space()

u_exact_expr = u1*(x < x_int) + u2*(x >= x_int)
u_exact_expr = u_exact_expr.subs('x', FESTIM.x)
u_exact_expr = fenics.Expression(sp.printing.ccode(u_exact_expr), degree=1)
u_e = fenics.interpolate(u_exact_expr, V)
error_L2 = fenics.errornorm(u_e, u_sim, 'L2')
print("L2 error ", error_L2)
vertex_values_u_D = u_e.compute_vertex_values(V.mesh())
vertex_values_u = u_sim.compute_vertex_values(V.mesh())
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
print("Max error ", error_max)

# plot
X1 = np.linspace(0, interface, num=50, endpoint=True)
X2 = np.linspace(interface + 1e-5, size)

C1 = np.array([u_exact(val) for val in X1])
C2 = np.array([u_exact(val) for val in X2])
plt.plot(X1/xL, C1/u0, color="tab:blue", label="Analytical")
plt.plot(X2/xL, C2/u0, color="tab:blue")

# rescale solution
u_scaled = fenics.project(u_sim/u0, u_sim.function_space())
u_scaled.function_space().mesh().coordinates()[:, 0] *= 1/xL
fenics.plot(u_scaled, linestyle="None", marker="+", label="FESTIM")

plt.xlabel(r"x/L")
plt.ylabel(r"$c_\mathrm{m}/c_0$")
plt.legend()
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
plt.savefig("out_MES_case1.pdf")
