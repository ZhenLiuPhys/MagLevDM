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
for i in range(4):
    integrand = lambda X, t: 20 * np.array([Bx(*(X - [0., 0., 1.])) - Bx(*(X - [0., 0., -1.])),
                                            By(*(X - [0., 0., 1.])) - By(*(X - [0., 0., -1.])),
                                            Bz(*(X - [0., 0., 1.])) - Bz(*(X - [0., 0., -1.]))])
    xs = integrate.odeint(integrand, [0., 0.4 * i - 0.6, 1.], ts)
    ax.plot(xs[:,0], xs[:,1], xs[:,2], color = 'darkviolet', alpha = 0.3, lw = 0.8)
    ax.add_artist(Arrow3D([xs[200, 0], xs[210, 0]], [xs[200, 1], xs[210, 1]], [xs[200, 2], xs[210, 2]], mutation_scale = 25., lw = 1.2, arrowstyle='->', color = 'darkviolet', alpha = 0.3))
    ax.plot(xs[:,0], xs[:,1], -xs[:,2], color = 'darkviolet', alpha = 0.3, lw = 0.8)
    ax.add_artist(Arrow3D([xs[200, 0], xs[210, 0]], [xs[200, 1], xs[210, 1]], [-xs[200, 2], -xs[210, 2]], mutation_scale = 25., lw = 1.2, arrowstyle='->', color = 'darkviolet', alpha = 0.3))
    xs = integrate.odeint(integrand, [0., 0.4 * i - 0.6, 1.], -ts)
    ax.plot(xs[:,0], xs[:,1], xs[:,2], color = 'darkviolet', alpha = 0.3, lw = 0.8)
    ax.plot(xs[:,0], xs[:,1], -xs[:,2], color = 'darkviolet', alpha = 0.3, lw = 0.8)

ax.plot(np.cos(2 * np.pi * ts), np.sin(2 * np.pi * ts), 0 * ts + 1, color = 'black')
ax.add_artist(Arrow3D([1., 1.], [0., 0.1], [1., 1.], mutation_scale = 25., lw = 1.5, arrowstyle='->', color = 'black'))
ax.plot(np.cos(2 * np.pi * ts), np.sin(2 * np.pi * ts), 0 * ts - 1, color = 'black')
ax.add_artist(Arrow3D([1., 1.], [0., -0.1], [-1., -1.], mutation_scale = 25., lw = 1.5, arrowstyle='->', color = 'black'))

rs = np.linspace(0., 1., 100)
thetas = np.linspace(0., np.pi, 100)
phis = np.linspace(0., 2 * np.pi, 100)
polarr, polarphi = np.meshgrid(rs, phis)
S1 = np.array([polarr * np.cos(polarphi), polarr * np.sin(polarphi), np.ones_like(polarr)]).T
ax.plot_surface(*S1.T, color = 'crimson', alpha = 0.6)
S2 = np.array([polarr * np.cos(polarphi), polarr * np.sin(polarphi), -np.ones_like(polarr)]).T
ax.plot_surface(*S2.T, color = 'limegreen', alpha = 0.6)
ax.text(0., 0., 1.4, r'$S_1$', color = 'crimson', ha = 'center', va = 'center', size = 'x-large')
ax.text(0., 0.5, -0.55, r'$S_2$', color = 'limegreen', ha = 'center', va = 'center', size = 'x-large')

r = np.array([0.4, 1.2, 0.4])
rad = 0.3
ax.scatter(r[0], r[1], r[2], color = 'black')
ax.text(r[0], r[1] + 0.1, r[2], r'$\bm r$', color = 'black', ha = 'center', va = 'center', size = 'x-large')

project = lambda X: rad * (X - r) / np.linalg.norm(X - r, axis = -1)[...,None] + r
sphtheta, sphphi = np.meshgrid(thetas, phis)
ax.plot_surface(rad * np.sin(sphtheta) * np.cos(sphphi) + r[0], rad * np.sin(sphtheta) * np.sin(sphphi) + r[1], rad * np.cos(sphtheta) + r[2], color = 'black', alpha = 0.05)
ax.plot_surface(*project(S1).T, color = 'crimson', alpha = 0.4)
ax.plot_surface(*project(S2).T, color = 'limegreen', alpha = 0.4)
prpoints = [[0., -1., 1.], [0., 1., 1.], [0., -1., -1.], [0., 1., -1.]]
for prpoint in prpoints:
    ax.plot([prpoint[0], r[0]], [prpoint[1], r[1]], [prpoint[2], r[2]], color = '0.2', lw = 0.5)
ax.text(r[0], r[1] - 0.25, r[2] + 0.35, r'$\Omega_1$', color = 'crimson', ha = 'center', va = 'center', size = 'x-large')
ax.text(r[0], r[1] - 0.25, r[2] - 0.35, r'$\Omega_2$', color = 'limegreen', ha = 'center', va = 'center', size = 'x-large')

ax.set_xlim(-1.25, 1.25)
ax.set_ylim(-1.25, 1.25)
ax.set_zlim(-1.25, 1.25)
ax.set_box_aspect((1., 1., 1.))
ax.axis('off')

fig.tight_layout()
#fig.show()
fig.savefig('omegas.pdf')
