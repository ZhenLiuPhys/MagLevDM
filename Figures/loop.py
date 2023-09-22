from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
import numpy as np
from scipy import integrate

fig = plt.figure(figsize = (6., 6.))
ax = fig.add_subplot(projection = '3d')
ax.view_init(20., 0.)

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'text.latex.preamble': r'\usepackage{bm}'})

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

Bdenom = lambda x, y, z, theta: (x ** 2 + y ** 2 + z ** 2 + 1 - 2 * x * np.cos(theta) - 2 * y * np.sin(theta)) ** 1.5
Bx = lambda x, y, z: integrate.quad(lambda theta: z * np.cos(theta) / Bdenom(x, y, z, theta), 0, 2 * np.pi)[0] / (4 * np.pi)
By = lambda x, y, z: integrate.quad(lambda theta: z * np.sin(theta) / Bdenom(x, y, z, theta), 0, 2 * np.pi)[0] / (4 * np.pi)
Bz = lambda x, y, z: integrate.quad(lambda theta: (1 - x * np.cos(theta) - y * np.sin(theta)) / Bdenom(x, y, z, theta), 0, 2 * np.pi)[0] / (4 * np.pi)

ts = np.linspace(0., 1., 1000)
for i in range(-3, 4):
    xs = integrate.odeint(lambda X, t: 30 * np.array([Bx(*X), By(*X), Bz(*X)]), [0, 0.25 * i, 0], ts)
    ax.plot(xs[:,0], xs[:,1], xs[:,2], color = 'darkviolet', alpha = 0.3, lw = 0.8)
    if i != -3 and i != 3:
        ax.plot(xs[:,0], xs[:,1], -xs[:,2], color = 'darkviolet', alpha = 0.3, lw = 0.8)
        ax.add_artist(Arrow3D([xs[200, 0], xs[210, 0]], [xs[200, 1], xs[210, 1]], [xs[200, 2], xs[210, 2]], mutation_scale = 25., lw = 1.2, arrowstyle='->', color = 'darkviolet', alpha = 0.3))
    else:
        ax.add_artist(Arrow3D([xs[50, 0], xs[60, 0]], [xs[50, 1], xs[60, 1]], [xs[50, 2], xs[60, 2]], mutation_scale = 25., lw = 1.2, arrowstyle='->', color = 'darkviolet', alpha = 0.3))

ax.plot(np.cos(2 * np.pi * ts), np.sin(2 * np.pi * ts), 0 * ts, color = 'black')
ax.add_artist(Arrow3D([1., 1.], [0., 0.1], [0., 0.], mutation_scale = 25., lw = 1.5, arrowstyle='->', color = 'black'))
ax.add_artist(Arrow3D([0., 0.], [0., -1.], [0.075, 0.075], mutation_scale = 5., lw = 2., arrowstyle='|-|', color = 'grey'))
ax.text(1.3, 0.05, 0., r'$I$', color = 'black', ha = 'center', va = 'center', size = 'x-large')
ax.text(0., -0.45, 0.15, r'$R$', color = 'grey', ha = 'center', va = 'center', size = 'x-large')

l = np.array([0.5, np.sqrt(3) / 2, 0.])
r = np.array([-0.5, 1.1, 1.])
ax.add_artist(Arrow3D([0., l[0]], [0., l[1]], [0., l[2]], mutation_scale = 10., lw = 2., arrowstyle='-|>', color = 'crimson'))
ax.add_artist(Arrow3D([0., r[0]], [0., r[1]], [0., r[2]], mutation_scale = 10., lw = 2., arrowstyle='-|>', color = 'limegreen'))
ax.add_artist(Arrow3D([l[0], r[0]], [l[1], r[1]], [l[2], r[2]], mutation_scale = 10., lw = 2., arrowstyle='-|>', color = 'royalblue'))
ax.scatter(r[0], r[1], r[2], color = 'black')

ax.text(1.08 * r[0], 1.08 * r[1], r[2], r'$\bm r$', color = 'black', ha = 'center', va = 'center', size = 'x-large')
ax.text(0.52 * l[0], 0.52 * l[1], l[2] + 0.12, r'$\bm l$', color = 'crimson', ha = 'center', va = 'center', size = 'x-large')
ax.text(0.58 * (l[0] + r[0]), 0.58 * (l[1] + r[1]), 0.5 * (l[2] + r[2]), r'$\bm r - \bm l$', color = 'royalblue', ha = 'center', va = 'center', size = 'x-large')
##ax.text(1.08 * r[0], 1.08 * r[1], r[2], r'$r$', color = 'black', ha = 'center', va = 'center', size = 'x-large')
##ax.text(0.55 * l[0], 0.55 * l[1], l[2] + 0.12, r'$l$', color = 'crimson', ha = 'center', va = 'center', size = 'x-large')
##ax.text(0.58 * (l[0] + r[0]), 0.58 * (l[1] + r[1]), 0.5 * (l[2] + r[2]), r'$r - l$', color = 'royalblue', ha = 'center', va = 'center', size = 'x-large')

ax.set_xlim(-1., 1.)
ax.set_ylim(-1., 1.)
ax.set_zlim(-1., 1.)
ax.set_box_aspect((1., 1., 1.))
ax.axis('off')

fig.tight_layout()
#fig.show()
fig.savefig('loop.pdf')
