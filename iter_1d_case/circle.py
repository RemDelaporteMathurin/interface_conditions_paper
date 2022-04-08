# Topics: line, color, LineCollection, cmap, colorline, codex
'''
Defines a function colorline that draws a (multi-)colored 2D line with coordinates x and y.
The color is taken from optional data in z, and creates a LineCollection.

z can be:
- empty, in which case a default coloring will be used based on the position along the input arrays
- a single number, for a uniform color [this can also be accomplished with the usual plt.plot]
- an array of the length of at least the same length as x, to color according to this data
- an array of a smaller length, in which case the colors are repeated along the curve

The function colorline returns the LineCollection created, which can be modified afterwards.

See also: plt.streamplot
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.patches as patches

try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass


def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.
    The angle should be given in radians.
    Args:
        origin (float, float): coordinates of origin point
        point (float, float): coordinates of point to be rotated
        angle (float): rotaton angle in radians (counterclockwise)
    Returns:
        float, float: rotated point coordinates
    """
    ox, oy = origin
    px, py = point

    qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
    qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
    return qx, qy


# Data manipulation:

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments


# Interface to LineCollection:

def colorline(x, y, z=None, cmap=plt.get_cmap('viridis'), norm=plt.Normalize(0.0, 1.0), linewidth=12, alpha=1):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
           
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, norm=norm, linewidth=linewidth, alpha=alpha)
    
    ax = plt.gca()
    ax.add_collection(lc)
    
    return lc


def read_profile(i, N, folder, field='0'):
    data = np.genfromtxt(folder + '/{:.0f}.csv'.format(i), delimiter=',', names=True)
    z = data[field]
    return np.log10(z)


# plot cm

x = np.genfromtxt('profiles_solute_mu_continuity/{:.0f}.csv'.format(0), delimiter=',', names=True)['arc_length']
x = x/x.max() + 0.2

data = np.genfromtxt('chemical_pot_conservation/derived_quantities.csv', delimiter=',', names=True)

fig, ax = plt.subplots()
ax.set_aspect('equal')
z_min, z_max = (np.log10(1e20), np.log10(8.13e23))

thetas = np.linspace(np.pi/20, np.pi, num=len(data["ts"])-1)
folder = "profiles_solute_mu_continuity"
for i, theta in enumerate(thetas):
    points = []
    for point_x in x:
        points.append(rotate(origin=(0, 0), point=(point_x, 0), angle=theta))
    points = np.array(points)
    z = read_profile(i, len(x), folder=folder)
    lc = colorline(points[:, 0], points[:, 1], z, norm=plt.Normalize(z_min, z_max))
cb = plt.colorbar(lc)

cb.set_label("log($c_\mathrm{m}$)")

data = np.genfromtxt('concentration_continuity/derived_quantities.csv', delimiter=',', names=True)
thetas = np.linspace(np.pi/20, np.pi, num=len(data["ts"])-1)
folder = "profiles_solute_c_continuity"
for i, theta in enumerate(-thetas[:-1]):
    points = []
    for point_x in x:
        points.append(rotate(origin=(0, 0), point=(point_x, 0), angle=theta))
    points = np.array(points)
    z = read_profile(i, len(x), folder=folder)
    lc = colorline(points[:, 0], points[:, 1], z, norm=plt.Normalize(z_min, z_max))

plt.text(-0.15, 0, 'Time')
plt.text(0, 1.21, r'$\mu$ continuity')
plt.text(0, -1.3, r'$c_\mathrm{m}$ continuity')

# draw arrows
style = "Simple, tail_width=0.5, head_width=4, head_length=8"
kw = dict(arrowstyle=style, color="k")
a1 = patches.FancyArrowPatch((0.15, 0), (0, 0.15),
                             connectionstyle="arc3,rad=.5", **kw)
a2 = patches.FancyArrowPatch((0.15, 0), (0, -0.15),
                             connectionstyle="arc3,rad=-.5", **kw)
plt.gca().add_patch(a1)
plt.gca().add_patch(a2)

xmax = 1.2
plt.xlim(-xmax, xmax)
plt.ylim(-xmax, xmax)
fig.patch.set_visible(False)
plt.axis('off')
plt.show()

# plot retention

x = np.genfromtxt('profiles_retention_mu_continuity/{:.0f}.csv'.format(0), delimiter=',', names=True)['arc_length']
x = x/x.max() + 0.2
h = 0*0.2
data = np.genfromtxt('chemical_pot_conservation/derived_quantities.csv', delimiter=',', names=True)

fig, ax = plt.subplots()
ax.set_aspect('equal')
z_min, z_max = (np.log10(1e22), np.log10(6.5e25))

thetas = np.linspace(np.pi/20, np.pi, num=len(data["ts"])-1)
folder = "profiles_retention_mu_continuity"
for i, theta in enumerate(thetas):
    points = []
    for point_x in x:
        points.append(rotate(origin=(0, 0), point=(point_x, 0), angle=theta))
    points = np.array(points)
    z = read_profile(i, len(x), folder=folder, field="retention")
    lc = colorline(points[:, 0], points[:, 1] + h, z, norm=plt.Normalize(z_min, z_max))
cb = plt.colorbar(lc)

cb.set_label("log(retention)")

data = np.genfromtxt('concentration_continuity/derived_quantities.csv', delimiter=',', names=True)
thetas = np.linspace(np.pi/20, np.pi, num=len(data["ts"])-2)
folder = "profiles_retention_c_continuity"

for i, theta in enumerate(-thetas):
    points = []
    for point_x in x:
        points.append(rotate(origin=(0, 0), point=(point_x, 0), angle=theta))
    points = np.array(points)
    z = read_profile(i, len(x), folder=folder, field="retention")
    lc = colorline(points[:, 0], points[:, 1] - h, z, norm=plt.Normalize(z_min, z_max))

plt.text(-0.15, 0, 'Time')
plt.text(0, 1.21, r'$\mu$ continuity')
plt.text(0, -1.3, r'$c_\mathrm{m}$ continuity')

# draw arrows
style = "Simple, tail_width=0.5, head_width=4, head_length=8"
kw = dict(arrowstyle=style, color="k")
a1 = patches.FancyArrowPatch((0.15, 0), (0, 0.15),
                             connectionstyle="arc3,rad=.5", **kw)
a2 = patches.FancyArrowPatch((0.15, 0), (0, -0.15),
                             connectionstyle="arc3,rad=-.5", **kw)
plt.gca().add_patch(a1)
plt.gca().add_patch(a2)

xmax = 1.2
plt.xlim(-xmax, xmax)
plt.ylim(-xmax, xmax)
fig.patch.set_visible(False)
plt.axis('off')
plt.show()
