import colorsys
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import FancyArrowPatch, Wedge
import numpy as np

fig, ax = plt.subplots(figsize = (6., 3.))

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'text.latex.preamble': r'\usepackage{bm}'})

def lighten(c, scale_l):
    rgb = colors.to_rgba(c)[:3]
    h, l, s = colorsys.rgb_to_hls(*rgb)
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

width = 0.5
delta = 0.2
ax.add_artist(Wedge((-1., 0.), width, 90, 270, ec = None, fc = lighten('crimson', 0.75), alpha = 0.3))
ax.add_artist(Wedge((1., 0.), width, 270, 90, ec = None, fc = lighten('crimson', 0.75), alpha = 0.3))
ax.add_artist(plt.Rectangle((-1., -width), delta, 2 * width, ec = None, fc = 'crimson', alpha = 0.3))
ax.add_artist(plt.Rectangle((1 - delta, -width), delta, 2 * width, ec = None, fc = 'crimson', alpha = 0.3))
ax.add_artist(plt.Rectangle((-1 + delta, -width), 2 - 2 * delta, 2 * width, ec = None, fc = lighten('crimson', 1.5), alpha = 0.3))
ax.plot([-1., -1.], [-width, width], color = lighten('crimson', 1.5), ls = ':')
ax.plot([-1 + delta, -1 + delta], [-width, width], color = lighten('crimson', 1.5), ls = ':')
ax.plot([1 - delta, 1 - delta], [-width, width], color = lighten('crimson', 1.5), ls = ':')
ax.plot([1., 1.], [-width, width], color = lighten('crimson', 1.5), ls = ':')

r = [0.6, 0.35]
ax.plot([-1., 1.], [0., 0.], color = 'black')
ax.add_artist(FancyArrowPatch((0. - 0.025, 0.), (0., 0.), mutation_scale = 25., arrowstyle='->', color = 'black'))
ax.add_artist(FancyArrowPatch((0., 0.), (r[0], r[1]), mutation_scale = 10., lw = 2., arrowstyle='-|>', color = 'black'))
ax.plot([r[0], r[0]], [0, r[1]], color = 'grey', ls = '--')
ax.text(r[0] + 0.075, r[1] / 2, r'$z$', color = 'grey', ha = 'center', va = 'center', size = 'x-large')
ax.plot([0.015 * r[0] / r[1], r[0]], [0.015, 0.015], color = 'grey', ls = '--')
ax.text(r[0] / 2, -0.075, r'$\rho$', color = 'grey', ha = 'center', va = 'center', size = 'x-large')
ax.add_artist(plt.Circle(r, radius = 0.015, color = 'black', zorder = 3))
ax.text(r[0] + 0.06, r[1] + 0.06, r'$\bm r$', color = 'black', ha = 'center', va = 'center', size = 'x-large')

ax.add_artist(FancyArrowPatch((-0.5, 0.), (-0.5, width), mutation_scale = 3., arrowstyle='|-|', color = 'grey'))
ax.text(-0.57, width / 2, r'$\epsilon$', color = 'grey', ha = 'center', va = 'center', size = 'x-large')
ax.add_artist(FancyArrowPatch((1., 0.), (1 + width * np.cos(0.6), width * np.sin(0.6)), mutation_scale = 3., arrowstyle='|-|', color = '0.25'))
ax.text(1.245, 0.065, r'$\epsilon$', color = '0.25', ha = 'center', va = 'center', size = 'x-large')
ax.add_artist(FancyArrowPatch((1 - delta, -width / 2), (1., -width / 2), mutation_scale = 3., arrowstyle='|-|', color = '0.4'))
ax.text(1 - delta / 2, -width / 2 - 0.1, r'$\delta$', color = '0.4', ha = 'center', va = 'center', size = 'x-large')

ax.text(0, -width - 0.1, r'$\rho<R-\delta$', color = lighten('crimson', 1.5), ha = 'center', va = 'center', size = 'large')
ax.text(-1 + delta / 2, width + 0.1, r'$R-\delta<\rho<R$', color = 'crimson', ha = 'center', va = 'center', size = 'large')
ax.text(1 - delta / 2, width + 0.1, r'$R-\delta<\rho<R$', color = 'crimson', ha = 'center', va = 'center', size = 'large')
ax.text(-1 - width / 2, -width - 0.1, r'$\rho>R$', color = lighten('crimson', 0.75), ha = 'center', va = 'center', size = 'large')
ax.text(1 + width / 2, -width - 0.1, r'$\rho>R$', color = lighten('crimson', 0.75), ha = 'center', va = 'center', size = 'large')

ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-0.8, 0.8)
ax.axis('off')

fig.tight_layout()
#fig.show()
fig.savefig('v1int.pdf')
