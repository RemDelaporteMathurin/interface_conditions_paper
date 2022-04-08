from fenics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


meshfile = "2d_meshes/mesh_domains.xdmf"


def scientificNotation(value, pos=0):
    if value == 0:
        return '0'
    else:
        e = np.log10(np.abs(value))
        m = np.sign(value) * 10 ** (e - int(e))
        return r'${:.1f} \times 10^{{{:d}}}$'.format(m, int(e))


fmt_labels = scientificNotation
fmt_colorbar = ticker.FuncFormatter(scientificNotation)

fontsize = 9
inline_spacing = 50
linewidths = 0.5

ext = "pdf"
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


if __name__ == "__main__":
    mesh = load_mesh(meshfile)
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 4.8), gridspec_kw={'width_ratios': [3, 1]})
    plt.axes(axs[0])
    plot(mesh, linewidth=0.01)
    plt.axes(axs[1])
    plot(mesh, linewidth=0.1)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.tick_right()
    width = 0.002
    translation_factor = 0.7
    origin = (-0.014 + translation_factor*width, 0.0145 - translation_factor*width)
    axs[1].set_xlim(origin[0] - width, origin[0] + width)
    axs[1].set_ylim(origin[1] - width, origin[1] + width)
    for ax in axs:
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
    plt.tight_layout()
    plt.savefig('mesh.pdf')
