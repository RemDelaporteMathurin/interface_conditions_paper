import FESTIM
import fenics
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# # # Meshing
def mesh_and_mark(dim=1, N=100, size=1, interface=0.5):
    nx = ny = nz = N

    coords = [FESTIM.x]
    if dim >= 2:
        coords.append(FESTIM.y)
    if dim == 3:
        coords.append(FESTIM.z)

    if dim == 1:
        mesh = fenics.UnitIntervalMesh(nx)
    if dim == 2:
        mesh = fenics.UnitSquareMesh(nx, ny)
    if dim == 3:
        mesh = fenics.UnitCubeMesh(nx, ny, nz)
    x = mesh.coordinates()
    x[:, :] *= size

    tol = 1E-14
    subdomain_1 = fenics.CompiledSubDomain('x[0] <= interface + tol',
                                           tol=tol, size=size,
                                           interface=interface)
    subdomain_2 = fenics.CompiledSubDomain('x[0] >= interface - tol', tol=tol,
                                           size=size, interface=interface)
    top_l = fenics.CompiledSubDomain('on_boundary && near(x[1], size, tol) && x[0] <= interface + tol',
                                     tol=tol, size=size, interface=interface)
    top_r = fenics.CompiledSubDomain('on_boundary && near(x[1], size, tol) && x[0] >= interface - tol',
                                     tol=tol, size=size, interface=interface)
    bottom_l = fenics.CompiledSubDomain('on_boundary && near(x[1], 0, tol) && x[0] <= interface + tol',
                                        tol=tol, size=size, interface=interface)
    bottom_r = fenics.CompiledSubDomain('on_boundary && near(x[1], 0, tol) && x[0] >= interface - tol',
                                        tol=tol, size=size, interface=interface)
    left = fenics.CompiledSubDomain('on_boundary && near(x[0], 0, tol)',
                                    tol=tol, size=size, interface=interface)
    right = fenics.CompiledSubDomain('on_boundary && near(x[0], size, tol)',
                                     tol=tol, size=size, interface=interface)
    front_l = fenics.CompiledSubDomain('on_boundary && near(x[2], 0, tol) && x[0] <= interface + tol',
                                       tol=tol, size=size, interface=interface)
    front_r = fenics.CompiledSubDomain('on_boundary && near(x[2], 0, tol) && x[0] >= interface - tol',
                                       tol=tol, size=size, interface=interface)
    back_l = fenics.CompiledSubDomain('on_boundary && near(x[2], size, tol) && x[0] <= interface + tol',
                                      tol=tol, size=size, interface=interface)
    back_r = fenics.CompiledSubDomain('on_boundary && near(x[2], size, tol) && x[0] >= interface - tol',
                                      tol=tol, size=size, interface=interface)
    vm = fenics.MeshFunction("size_t", mesh, mesh.topology().dim())
    sm = fenics.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    vm.set_all(1)
    subdomain_1.mark(vm, 1)
    subdomain_2.mark(vm, 2)

    sm.set_all(0)
    left.mark(sm, 1)
    right.mark(sm, 2)
    if dim >= 2:
        top_l.mark(sm, 3)
        top_r.mark(sm, 4)
        bottom_l.mark(sm, 5)
        bottom_r.mark(sm, 6)
    if dim == 3:
        front_l.mark(sm, 7)
        front_r.mark(sm, 8)
        back_l.mark(sm, 9)
        back_r.mark(sm, 10)
    return mesh, vm, sm, coords


interface = 0.5
size = 1
mesh, vm, sm, coords = \
    mesh_and_mark(dim=2, N=100, size=size, interface=interface)

# Parameters
T = 500 + 30*sp.cos(2*fenics.pi*FESTIM.x)*sp.cos(2*fenics.pi*FESTIM.y)*sp.cos(2*fenics.pi*FESTIM.t)
S_01, S_02 = 1, 1.05
D_01, D_02 = 1, 2
E_D1, E_D2 = 0.1, 0.2
E_S1, E_S2 = 0.1, 0.09
D1 = D_01*sp.exp(-E_D1/FESTIM.k_B/T)
D2 = D_02*sp.exp(-E_D2/FESTIM.k_B/T)
S1 = S_01*sp.exp(-E_S1/FESTIM.k_B/T)
S2 = S_02*sp.exp(-E_S2/FESTIM.k_B/T)
c1 = 2 + FESTIM.t + sp.cos(2*fenics.pi*FESTIM.x)*sp.cos(2*fenics.pi*FESTIM.y)
c2 = S2/S1*c1

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

    "temperature": {
        "type": "expression",
        "value": T
    },
    "initial_conditions": [
        {
            "value": c1 * (FESTIM.x <= interface) + c2 * (FESTIM.x > interface),
            "component": 0
        },
    ],
    "traps": [],
    "boundary_conditions": [
        {
            "type": "dc",
            "value": c1,
            "surfaces": [1, 3, 5, 7, 9]
        },
        {
            "type": "dc",
            "value": c2,
            "surfaces": [2, 4, 6, 8, 10]
        },
    ],
    "mesh_parameters": {
        "mesh": mesh,
        "meshfunction_cells": vm,
        "meshfunction_facets": sm,
    },
    "exports": {
        "xdmf": {
                "functions": ["solute"],
                "labels": ["solute"],
                "folder": "MMS_chemical_pot",
                "checkpoint": False,
            },
        "error": {}
    },
    "solving_parameters": {
        "final_time": 0.01,
        "initial_stepsize": 0.01,
        # "type": "solve_stationary",
        "newton_solver": {
            "absolute_tolerance": 1e-10,
            "relative_tolerance": 1e-10,
            "maximum_iterations": 15
        }
    }
}

# create source terms
D = [D1, D2]
cs = [c1, c2]
f = []
for i in range(0, 2):
    f_coeff = sp.diff(cs[i], FESTIM.t)
    for coord in coords:
        f_coeff -= sp.diff(D[i]*sp.diff(cs[i], coord), coord)
    f.append(f_coeff)
parameters["source_term"] = [
    {
        "value": f[0],
        "volumes": [1]
    },
    {
        "value": f[1],
        "volumes": 2
    }
    ]

for t_final in [0.01, 0.06]:
    plt.figure()
    parameters["solving_parameters"]["final_time"] = t_final
    # Solve
    out = FESTIM.generic_simulation.run(parameters)

    # Post processing
    levels = 100
    t_final = parameters["solving_parameters"]["final_time"]
    fontsize = 15
    if len(coords) == 1:
        def u_exact(x):
            if x <= interface:
                val = c1
            elif x > interface:
                val = c2
            val = val.subs(FESTIM.x, x).subs(FESTIM.t, t_final)
            return val

        X1 = np.linspace(0, interface, num=50, endpoint=True)
        X2 = np.linspace(interface + 1e-5, size)

        C1 = [u_exact(val) for val in X1]
        C2 = [u_exact(val) for val in X2]

        plt.plot(X1, C1, color="tab:blue", label="Analytical")
        plt.plot(X2, C2, color="tab:blue")
        fenics.plot(out["u"], linestyle="None", marker="+")
        u_exact_expr = c1*(FESTIM.x < interface) + c2*(FESTIM.x > interface)
        u_exact_expr = fenics.Expression(sp.printing.ccode(u_exact_expr), degree=4, t=t_final)
        V = out["u"].function_space()
        u_exact_expr = fenics.interpolate(u_exact_expr, V)
        print(fenics.errornorm(u_exact_expr, out["u"], "L2"))
        plt.savefig("u_computed.pdf")
    if len(coords) == 2:
        # plot u_computed
        CS = fenics.plot(out["u"], levels=levels)
        fenics.plot(out["u"], mode="contour", colors="white")
        cb = plt.colorbar(CS)
        cb.ax.tick_params(labelsize=fontsize)
        cb.ax.set_title(r"c")
        plt.xlabel(r"$x$", fontsize=fontsize)
        plt.ylabel(r"$y$", fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        for c in CS.collections:  # for avoiding white lines in pdf
            c.set_edgecolor("face")
        plt.savefig("u_computed_t={:.2f}.pdf".format(t_final))
        V = out["u"].function_space()
        # plot u_exact
        u_exact_expr = c1*(FESTIM.x < interface) + c2*(FESTIM.x >= interface)
        u_exact_expr = fenics.Expression(sp.printing.ccode(u_exact_expr), degree=1, t=t_final)
        u_exact_expr = fenics.interpolate(u_exact_expr, V)
        plt.figure()
        CS = fenics.plot(u_exact_expr, levels=levels)
        fenics.plot(u_exact_expr, mode="contour", colors="white")
        cb = plt.colorbar(CS)
        cb.ax.set_title(r"c")
        cb.ax.tick_params(labelsize=fontsize)
        plt.xlabel(r"$x$", fontsize=fontsize)
        plt.ylabel(r"$y$", fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        for c in CS.collections:  # for avoiding white lines in pdf
            c.set_edgecolor("face")
        plt.savefig("u_exact_t={:.2f}.pdf".format(t_final))

        # l2 norm

        print(fenics.errornorm(u_exact_expr, out["u"], "L2"))

        # plot difference
        plt.figure()
        diff = fenics.project(abs(out["u"] - u_exact_expr), V)
        levels = np.linspace(min(diff.vector()), 0.0324, num=100)
        CS = fenics.plot(diff, levels=levels)
        cb = plt.colorbar(CS, ticks=np.linspace(0, 0.0324, num=7))
        cb.ax.set_title(r"|$c-c_M$|")
        cb.ax.tick_params(labelsize=fontsize)
        plt.xlabel(r"$x$", fontsize=fontsize)
        plt.ylabel(r"$y$", fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        for c in CS.collections:  # for avoiding white lines in pdf
            c.set_edgecolor("face")
        plt.savefig("diff_t={:.2f}.pdf".format(t_final))
