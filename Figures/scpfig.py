from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
import numpy as np
from scipy import integrate, special

fig = plt.figure(figsize = (6., 6.))
ax = fig.add_subplot(projection = '3d')
ax.view_init(20., 0.)

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

Bdenom = lambda x, y, z, theta: (x ** 2 + y ** 2 + z ** 2 + 1 - 2 * x * np.cos(theta) - 2 * y * np.sin(theta)) ** 1.5
Bx = lambda x, y, z: integrate.quad(lambda theta: z * np.cos(theta) / Bdenom(x, y, z, theta), 0, 2 * np.pi)[0] / (4 * np.pi)
By = lambda x, y, z: integrate.quad(lambda theta: z * np.sin(theta) / Bdenom(x, y, z, theta), 0, 2 * np.pi)[0] / (4 * np.pi)
Bz = lambda x, y, z: integrate.quad(lambda theta: (1 - x * np.cos(theta) - y * np.sin(theta)) / Bdenom(x, y, z, theta), 0, 2 * np.pi)[0] / (4 * np.pi)
Bapp = lambda X: np.array([Bx(*(X - [0., 0., 1.])) - Bx(*(X - [0., 0., -1.])), By(*(X - [0., 0., 1.])) - By(*(X - [0., 0., -1.])), Bz(*(X - [0., 0., 1.])) - Bz(*(X - [0., 0., -1.]))])

r = np.array([0., 0.3, 0.])
rad = 0.25
scp = lambda theta, phi: np.array([rad * np.sin(theta) * np.cos(phi) + r[0], rad * np.sin(theta) * np.sin(phi) + r[1], rad * np.cos(theta) + r[2]])
rhat = lambda theta, phi: np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
thetahat = lambda theta, phi: np.array([np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)])
phihat = lambda theta, phi: np.array([-np.sin(phi), np.cos(phi), 0])

alms = []
lnum = 3
for l in range(1, lnum):
    coeffs = []
    coeffs.append(integrate.nquad(lambda theta, phi: np.sin(theta) * np.real(special.sph_harm(0, l, phi, theta))
                                  * np.dot(Bapp(scp(theta, phi)), rhat(theta, phi)), [[0, np.pi], [0, 2 * np.pi]])[0])
    for m in range(1, l + 1):
        re = integrate.nquad(lambda theta, phi: np.sin(theta) * np.real(special.sph_harm(m, l, phi, theta))
                             * np.dot(Bapp(scp(theta, phi)), rhat(theta, phi)), [[0, np.pi], [0, 2 * np.pi]])[0]
        im = integrate.nquad(lambda theta, phi: np.sin(theta) * -np.imag(special.sph_harm(m, l, phi, theta))
                             * np.dot(Bapp(scp(theta, phi)), rhat(theta, phi)), [[0, np.pi], [0, 2 * np.pi]])[0]
        coeffs.append(re + 1j * im)
        coeffs.insert(0, (-1) ** m * (re - 1j * im))
    alms.append(coeffs)

thetacart = lambda X: np.angle(X[2] + 1j * np.sqrt(X[0] ** 2 + X[1] ** 2))
phicart = lambda X: np.angle(X[0] + 1j * X[1])
Ylm = lambda X, l, m: special.sph_harm(m, l, phicart(X), thetacart(X)) * rhat(thetacart(X), phicart(X))
Philm_theta = lambda theta, phi, l, m: -np.sqrt((2 * l + 1) / (4 * np.pi) * special.factorial(l - m) / special.factorial(l + m)) * special.lpmn(m, l, np.cos(theta))[1][m, l] * np.sin(theta) * np.exp(1j * m * phi)
Philm_phi = lambda theta, phi, l, m: 1j * m * special.sph_harm(m, l, phi, theta) / np.sin(theta)
def Psilm(X, l, m):
    if np.abs(thetacart(X)) <= 1e-5:
        return (np.abs(m) == 1) * -np.sqrt((2 * l + 1) * (l + 1) * l / (16 * np.pi)) * np.array([m, 1j, 0.])
    if np.abs(thetacart(X) - np.pi) <= 1e-5:
        return (np.abs(m) == 1) * (-1) ** l * np.sqrt((2 * l + 1) * (l + 1) * l / (16 * np.pi)) * np.array([m, 1j, 0.])
    else:
        return Philm_theta(thetacart(X), phicart(X), l, m) * thetahat(thetacart(X), phicart(X)) + Philm_phi(thetacart(X), phicart(X), l, m) * phihat(thetacart(X), phicart(X))
Bresp = lambda X: np.real(np.sum([alms[l - 1][m + l] * np.linalg.norm((X - r) / rad) ** (-l - 2) * (-Ylm(X - r, l, m) + Psilm(X - r, l, m) / (l + 1)) for l in range(1, lnum) for m in range(-l, l + 1)], axis = 0))
Btot = lambda X: Bapp(X) + Bresp(X)

ts = np.linspace(0., 1., 1000)
for i in range(4):
    xs = integrate.odeint(lambda X, t: 20 * Btot(X), [0., 0.3 * i - 0.45, 1.], ts)
    ax.plot(xs[:,0], xs[:,1], xs[:,2], color = 'darkviolet', alpha = 0.5)
    ax.add_artist(Arrow3D([xs[100, 0], xs[110, 0]], [xs[100, 1], xs[110, 1]], [xs[100, 2], xs[110, 2]], mutation_scale = 25., lw = 1.4, arrowstyle='->', color = 'darkviolet', alpha = 0.5))
    ax.plot(xs[:,0], xs[:,1], -xs[:,2], color = 'darkviolet', alpha = 0.5)
    ax.add_artist(Arrow3D([xs[100, 0], xs[110, 0]], [xs[100, 1], xs[110, 1]], [-xs[100, 2], -xs[110, 2]], mutation_scale = 25., lw = 1.4, arrowstyle='->', color = 'darkviolet', alpha = 0.5))
    xs = integrate.odeint(lambda X, t: 20 * Btot(X), [0., 0.3 * i - 0.45, 1.], -ts)
    ax.plot(xs[:,0], xs[:,1], xs[:,2], color = 'darkviolet', alpha = 0.5)
    ax.plot(xs[:,0], xs[:,1], -xs[:,2], color = 'darkviolet', alpha = 0.5)

ax.plot(np.cos(2 * np.pi * ts), np.sin(2 * np.pi * ts), 0 * ts + 1, color = 'black')
ax.add_artist(Arrow3D([1., 1.], [0., 0.1], [1., 1.], mutation_scale = 25., lw = 1.5, arrowstyle='->', color = 'black'))
ax.plot(np.cos(2 * np.pi * ts), np.sin(2 * np.pi * ts), 0 * ts - 1, color = 'black')
ax.add_artist(Arrow3D([1., 1.], [0., -0.1], [-1., -1.], mutation_scale = 25., lw = 1.5, arrowstyle='->', color = 'black'))

thetas = np.linspace(0., np.pi, 100)
phis = np.linspace(0., 2 * np.pi, 100)
sphtheta, sphphi = np.meshgrid(thetas, phis)
ax.plot_surface(rad * np.sin(sphtheta) * np.cos(sphphi) + r[0], rad * np.sin(sphtheta) * np.sin(sphphi) + r[1], rad * np.cos(sphtheta) + r[2], color = '0.8')
ax.add_artist(Arrow3D([r[0], 0.], [r[1], -0.07], [r[2], 0.], mutation_scale = 25., lw = 2., arrowstyle='-|>', color = 'crimson', zorder = 3))

current_thetas = [np.pi / 15, np.pi / 4, np.pi / 3, 14 * np.pi / 15]
for i in range(4):
    if i != 3:
        current = integrate.odeint(lambda X, t: 20 * np.cross((X - r) / np.linalg.norm(X - r), Btot(X)), r + rad * rhat(current_thetas[i], np.pi / 2), ts)
    else:
        current = integrate.odeint(lambda X, t: 20 * np.cross((X - r) / np.linalg.norm(X - r), Btot(X)), r + rad * rhat(current_thetas[i], np.pi / 2), -ts)
    current = (current - r) * 1.05 + r
    endind = np.where(current[:,0] < 0)[0][0]
    current = current[:endind]
    ax.plot(current[:,0], current[:,1], current[:,2], color = 'deepskyblue', lw = 1.5, zorder = 3)
    arrowstart = [50, 70, 70, 70][i]
    if i != 3:
        ax.add_artist(Arrow3D([current[arrowstart, 0], current[arrowstart + 10, 0]], [current[arrowstart, 1], current[arrowstart + 10, 1]], [current[arrowstart, 2], current[arrowstart + 10, 2]],
                              mutation_scale = 15., lw = 1.9, arrowstyle='->', color = 'deepskyblue', zorder = 3))
    else:
        ax.add_artist(Arrow3D([current[arrowstart, 0], current[arrowstart - 10, 0]], [current[arrowstart, 1], current[arrowstart - 10, 1]], [current[arrowstart, 2], current[arrowstart - 10, 2]],
                              mutation_scale = 15., lw = 1.9, arrowstyle='->', color = 'deepskyblue', zorder = 3))

ax.text(0., 1.2, 0.38, r'$\bm B$', color = 'darkviolet', ha = 'center', va = 'center', size = 'x-large')
ax.text(r[0], r[1] + 1.1 * rad, r[2] + 1.1 * rad, r'$\bm K$', color = 'deepskyblue', ha = 'center', va = 'center', size = 'x-large')
ax.text(0., -0.15, 0., r'$\bm F$', color = 'crimson', ha = 'center', va = 'center', size = 'x-large')

ax.set_xlim(-1., 1.)
ax.set_ylim(-1., 1.)
ax.set_zlim(-0.9, 1.1)
ax.set_box_aspect((1., 1., 1.))
ax.axis('off')

fig.tight_layout()
#fig.show()
fig.savefig('scpfig.png')
