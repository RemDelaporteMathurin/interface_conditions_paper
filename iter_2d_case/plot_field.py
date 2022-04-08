from fenics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


meshfile = "mesh/mesh_domains.xdmf"


def scientificNotation(value, pos=0):
    if value == 0:
        return '0'
    else:
        e = np.log10(np.abs(value))
        m = np.sign(value) * 10 ** (e - int(e))
        return r'${:.1f} \times 10^{{{:d}}}$'.format(m, int(e))


fmt_labels = scientificNotation
fmt_colorbar = ticker.FuncFormatter(scientificNotation)

fontsize = 10
inline_spacing = 25
linewidths = 0.5

transparent = True

def load_mesh(meshfile):
    """Loads mesh

    Args:
        meshfile (str): path of the meshfile

    Returns:
        fenics.mesh: the mesh of the solution
    """
    mesh = Mesh()
    XDMFFile(meshfile).read(mesh)
    return mesh


def load_field(mesh, fieldfile, field, counter=-1):
    V = FunctionSpace(mesh, "DG", 1)
    u = Function(V)

    XDMFFile(fieldfile).read_checkpoint(u, field, counter)
    return u


def save(filename, ext, transparent=True):
    """Saves active figure

    Args:
        filename (str): path of the file (without extension)
        ext (str, list): extension(s) of the output file(s)
        (ex: "pdf", "svg", ["png", "pdf"])
        transparent (bool, optional): Sets the background transparent.
        Defaults to True.
    """
    if not isinstance(ext, list):
        ext = [ext]
    for e in ext:
        plt.savefig(filename + "." + e, transparent=transparent)
    return


def create_figure(u):

    fig = plt.figure(figsize=(4.8, 4.8))

    levels = np.linspace(
        -1e22, 1.1e26,
        num=400, endpoint=True)
    CS = plot(u, levels=levels, extend="min", vmin=0)

    # levels_contours = np.linspace(
    #     1e23, u.vector().max(),
    #     num=5)

    # CS2 = plot(
    #     u, mode="contour", levels=levels_contours,
    #     colors="white", linewidths=linewidths)

    # CL = plt.clabel(
    #     CS2, inline=True, fontsize=fontsize,
    #     inline_spacing=inline_spacing,
    #     fmt=fmt_labels)
    # CS2, CL = check_inline_labels(
    #     CS2, CL, thresh=thresh,
    #     inline_spacing=inline_spacing, fontsize=fontsize,
    #     fmt=fmt_labels)

    levels = np.linspace(
        0, 1.1e26,
        num=10, endpoint=True)

    cb = fig.colorbar(CS, ticks=levels, format=fmt_colorbar, extendfrac=0)
    cb.ax.set_title(r'm$^{-3}$')

    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    fig.patch.set_visible(False)
    plt.axis('off')
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    mesh = load_mesh(meshfile)

    # plot retention_mu
    field = "retention"

    folder = "chemical_pot_conservation"

    fieldfile = folder + "/retention.xdmf"
    u = load_field(mesh, fieldfile, "retention", -1)
    max_ret_mu = u.vector().max()
    fig = plt.figure(figsize=(4.8, 4.8))

    levels = np.linspace(
        -1e22, max_ret_mu,
        num=400, endpoint=True)
    CS = plot(u, levels=levels, extend="min", vmin=0)

    levels = np.linspace(
        0, max_ret_mu,
        num=10, endpoint=True)

    cb = fig.colorbar(CS, ticks=levels, format=fmt_colorbar, extendfrac=0)
    cb.ax.set_title(r'm$^{-3}$')

    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    fig.patch.set_visible(False)
    plt.axis('off')
    plt.tight_layout()

    save("retention_mu", "pdf", transparent)

    # plot retention_c
    field = "retention"

    folder = "concentration_continuity"

    fieldfile = folder + "/retention.xdmf"
    u = load_field(mesh, fieldfile, "retention", -1)
    fig = plt.figure(figsize=(4.8, 4.8))

    levels = np.linspace(
        -1e22, max_ret_mu,
        num=400, endpoint=True)
    CS = plot(u, levels=levels, extend="min", vmin=0)

    levels = np.linspace(
        0, max_ret_mu,
        num=10, endpoint=True)

    cb = fig.colorbar(CS, ticks=levels, format=fmt_colorbar, extendfrac=0)
    cb.ax.set_title(r'm$^{-3}$')

    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    fig.patch.set_visible(False)
    plt.axis('off')
    plt.tight_layout()

    save("retention_c", "pdf", transparent)


    # plot solute_mu
    field = "0"

    folder = "chemical_pot_conservation"

    fieldfile = folder + "/0.xdmf"
    u = load_field(mesh, fieldfile, "0", -1)
    max_ret_mu = u.vector().max()
    fig = plt.figure(figsize=(4.8, 4.8))

    levels = np.linspace(
        -1e22, max_ret_mu,
        num=400, endpoint=True)
    CS = plot(u, levels=levels, extend="min", vmin=0)

    levels = np.linspace(
        0, max_ret_mu,
        num=10, endpoint=True)

    cb = fig.colorbar(CS, ticks=levels, format=fmt_colorbar, extendfrac=0)
    cb.ax.set_title(r'm$^{-3}$')

    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    fig.patch.set_visible(False)
    plt.axis('off')
    plt.tight_layout()

    save("solute_mu", "pdf", transparent)
    # plot solute_c
    field = "0"

    folder = "concentration_continuity"

    fieldfile = folder + "/0.xdmf"
    u = load_field(mesh, fieldfile, "0", -1)
    fig = plt.figure(figsize=(4.8, 4.8))

    levels = np.linspace(
        -1e22, max_ret_mu,
        num=400, endpoint=True)
    CS = plot(u, levels=levels, extend="min", vmin=0)

    levels = np.linspace(
        0, max_ret_mu,
        num=10, endpoint=True)

    cb = fig.colorbar(CS, ticks=levels, format=fmt_colorbar, extendfrac=0)
    cb.ax.set_title(r'm$^{-3}$')

    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    fig.patch.set_visible(False)
    plt.axis('off')
    plt.tight_layout()

    save("solute_c", "pdf", transparent)

