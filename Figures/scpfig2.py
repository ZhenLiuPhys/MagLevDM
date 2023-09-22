from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import art3d, proj3d
import numpy as np

fig = plt.figure(figsize = (6., 6.))
ax = fig.add_subplot(projection = '3d')
ax.view_init(10., 0.)

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 14, 'text.latex.preamble': r'\usepackage{bm}'})

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

#wallx = Rectangle((0, 0), 1, 1, fc = '0.9', ec = 'black')
#ax.add_patch(wallx)
#art3d.pathpatch_2d_to_3d(wallx, z = 0, zdir = 'x')
wally = Rectangle((0, 0), 1, 1, fc = '0.9', ec = 'black')
ax.add_patch(wally)
art3d.pathpatch_2d_to_3d(wally, z = 0, zdir = 'y')
wally2 = Rectangle((0, 0), 1, 1, fc = '0.9', ec = 'black')
ax.add_patch(wally2)
art3d.pathpatch_2d_to_3d(wally2, z = 1, zdir = 'y')
wallz = Rectangle((0, 0), 1, 1, fc = '0.9', ec = 'black')
ax.add_patch(wallz)
art3d.pathpatch_2d_to_3d(wallz, z = 0, zdir = 'z')
wallz2 = Rectangle((0, 0), 1, 1, fc = '0.5', ec = 'black')
ax.add_patch(wallz2)
art3d.pathpatch_2d_to_3d(wallz2, z = 1, zdir = 'z')

r0 = np.array([0.5, 0.2, 0.7])
R = 0.1
h = 0.1
ts = np.linspace(0., 1., 1000)
ax.plot(r0[0] + R * np.cos(2 * np.pi * ts), r0[1] + R * np.sin(2 * np.pi * ts), r0[2] + h + 0 * ts, color = 'black', zorder = 3)
ax.add_artist(Arrow3D([r0[0] + R, r0[0] + R], [r0[1], r0[1] + 0.01], [r0[2] + h, r0[2] + h], mutation_scale = 15., lw = 1.5, arrowstyle='->', color = 'black', zorder = 3))
ax.plot(r0[0] + R * np.cos(2 * np.pi * ts), r0[1] + R * np.sin(2 * np.pi * ts), r0[2] - h + 0 * ts, color = 'black', zorder = 3)
ax.add_artist(Arrow3D([r0[0] + R, r0[0] + R], [r0[1], r0[1] - 0.01], [r0[2] - h, r0[2] - h], mutation_scale = 15., lw = 1.5, arrowstyle='->', color = 'black', zorder = 3))

rad = 0.025
r = r0 + np.array([0., 0.03, 0.])
thetas = np.linspace(0., np.pi, 100)
phis = np.linspace(0., 2 * np.pi, 100)
sphtheta, sphphi = np.meshgrid(thetas, phis)
ax.plot_surface(rad * np.sin(sphtheta) * np.cos(sphphi) + r[0], rad * np.sin(sphtheta) * np.sin(sphphi) + r[1], rad * np.cos(sphtheta) + r[2], color = '0.8', zorder = 3)

ax.add_artist(Arrow3D([0.5, 0.5], [0.5, 0.5], [0.1, 0.9], mutation_scale = 25., lw = 1.5, arrowstyle='->', color = 'royalblue', zorder = 3))
ax.text(0.5, 0.57, 0.7, r'$\bm J_\mathrm{eff}$', color = 'royalblue', ha = 'center', va = 'center', size = 'x-large')
ax.plot(0.5 + 0.4 * np.cos(2 * np.pi * ts), 0.5 + 0.4 * np.sin(2 * np.pi * ts), 0.5 + 0 * ts, color = 'limegreen', zorder = 3)
ax.add_artist(Arrow3D([0.9, 0.9], [0.5, 0.55], [0.5, 0.5], mutation_scale = 25., lw = 1.5, arrowstyle='->', color = 'limegreen', zorder = 3))
ax.text(0.5, 0.6, 0.37, r'$\bm B_\mathrm{DM}$', color = 'limegreen', ha = 'center', va = 'center', size = 'x-large')

ax.set_xlim(0.15, 0.85)
ax.set_ylim(0.15, 0.85)
ax.set_zlim(0.15, 0.85)
ax.set_box_aspect((1., 1., 1.))
ax.axis('off')

fig.tight_layout()
#fig.show()
fig.savefig('scpfig2.png')
