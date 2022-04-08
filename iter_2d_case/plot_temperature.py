from fenics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


meshfile = "mesh/mesh_domains.xdmf"
field = "T"

chemical_pot_conservation = False
if chemical_pot_conservation:
    folder = "chemical_pot_conservation"
else:
    folder = "concentration_continuity"

fieldfile = folder + "/T.xdmf"

fmt_labels = '%.0f'
fmt_colorbar = '%.0f'

fontsize = 10
inline_spacing = 25
linewidths = 0.5

ext = ["pdf", "svg"]
ext = "png"
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
    V = FunctionSpace(mesh, "P", 1)
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
        u.vector().min(), u.vector().max(),
        num=400, endpoint=True)
    CS = plot(u, levels=levels)

    levels_contours = np.linspace(
        u.vector().min(), u.vector().max(),
        num=8, endpoint=False)

    CS2 = plot(
        u, mode="contour", levels=levels_contours,
        colors="white", linewidths=linewidths)

    manual_locs = [
        (-0.00736115355427, 0.00520610012083),
        (-0.0068084810055, 0.00734090248789),
        (-0.00891319166176, 0.010180299781),
        (-0.00136339992743, 0.009849632786),
        (-0.00642217935178, 0.0112668146506),
        (-0.00765311779037, 0.012256455834),
        (-0.00138311386894, 0.0132052166893)]
    CL = plt.clabel(
        CS2, inline=True, fontsize=fontsize,
        inline_spacing=inline_spacing,
        fmt=fmt_labels, manual=manual_locs)


    # CS2, CL = check_inline_labels(
    #     CS2, CL, thresh=0.01,
    #     inline_spacing=inline_spacing, fontsize=fontsize,
    #     fmt=fmt_labels)


    levels = np.linspace(
        u.vector().min(), u.vector().max(),
        num=10, endpoint=True)

    cb = fig.colorbar(CS, ticks=levels, format=fmt_colorbar, extendfrac=0)
    cb.ax.set_title(r'K')

    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    fig.patch.set_visible(False)
    plt.axis('off')
    plt.tight_layout()
    return fig


def check_inline_labels(CS, CL, thresh=0.10,
                        inline_spacing=10, fontsize=10, fmt='%.2f'):
    """Checks that the contour labels don't overlapp with axes and remove them
    if so.

    Args:
        CS ([type]): [description]
        CL ([type]): [description]
        thresh (float, optional): [description]. Defaults to 0.10.
        inline_spacing (int, optional): [description]. Defaults to 10.
        fontsize (int, optional): [description]. Defaults to 10.
        fmt (str, optional): [description]. Defaults to '%.2f'.

    Returns:
        [type]: [description]
    """
    # check that the labels don't overlap axes
    # get limits if they're automatic
    xmin, xmax, ymin, ymax = plt.axis()
    Dx = xmax-xmin
    Dy = ymax-ymin

    # check which labels are near a border
    keep_labels = []
    for label in CL:
        lx, ly = label.get_position()
        if xmin+thresh*Dx < lx < xmax-thresh*Dx and \
           ymin+thresh*Dy < ly < ymax-thresh*Dy:
            # inlier, redraw it later
            pass
        elif lx < xmin+thresh*Dx:
            lx = xmin+2*thresh*Dx
        keep_labels.append((lx, ly))

    colors = []
    linewidths = []
    for cline in CS.collections:
        cline.remove()
        # colors.append(cline.get_colors())
        linewidths.append(cline.get_linewidths())
    for label in CL:
        label.remove()
    CS = plot(u, mode="contour", levels=CS.levels, colors="white", linewidths=linewidths)
    CL = plt.clabel(CS, inline=True, fontsize=fontsize,
                    inline_spacing=inline_spacing, manual=keep_labels, fmt=fmt)
    return CS, CL


if __name__ == "__main__":
    mesh = load_mesh(meshfile)
    u = load_field(mesh, fieldfile, field, -1)
    fig = create_figure(u)
    save("temperature_field_2d", "pdf", transparent)
