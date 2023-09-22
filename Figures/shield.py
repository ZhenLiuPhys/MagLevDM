import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Wedge
import numpy as np

fig, ax = plt.subplots(figsize = (6., 6.))

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'text.latex.preamble': r'\usepackage{bm}'})

ax.add_artist(plt.Rectangle((-0.1, -0.1), 1.2, 1.2, fc = '0.9', ec = 'black'))
ax.add_artist(plt.Rectangle((0., 0.), 1., 1., fc = 'white', ec = 'black'))

r0 = [0.8, 0.7]
R = 0.1
h = 0.1
ax.add_artist(plt.Circle(r0, radius = 0.005, color = 'black'))
ax.text(r0[0] + 0.04, r0[1], r'$\bm r_0$', color = 'black', ha = 'center', va = 'center', size = 'x-large')
ax.plot([r0[0] - R, r0[0] + R], [r0[1] + h, r0[1] + h], color = 'black')
ax.add_artist(FancyArrowPatch((r0[0], r0[1] + h), (r0[0] + 0.025, r0[1] + h), mutation_scale = 15., arrowstyle='->', color = 'black'))
ax.plot([r0[0] - R, r0[0] + R], [r0[1] - h, r0[1] - h], color = 'black')
ax.add_artist(FancyArrowPatch((r0[0], r0[1] - h), (r0[0] - 0.025, r0[1] - h), mutation_scale = 15., arrowstyle='->', color = 'black'))

ax.add_artist(FancyArrowPatch((r0[0] - 0.01, r0[1] + h + 0.035), (r0[0] + R + 0.01, r0[1] + h + 0.035), mutation_scale = 3., arrowstyle='|-|', color = 'grey'))
ax.text(r0[0] + R / 2, r0[1] + h + 0.065, r'$R$', color = 'grey', ha = 'center', va = 'center', size = 'x-large')
ax.add_artist(FancyArrowPatch((r0[0] - 0.03, r0[1] + 0.01), (r0[0] - 0.03, r0[1] - h + 0.003), mutation_scale = 3., arrowstyle='|-|', color = 'grey'))
ax.text(r0[0] - 0.06, r0[1] - h / 2, r'$h$', color = 'grey', ha = 'center', va = 'center', size = 'x-large')
ax.add_artist(FancyArrowPatch((0.03, 0.02), (0.99, 0.02), mutation_scale = 3., arrowstyle='|-|', color = 'grey'))
ax.text(0.5, 0.06, r'$L_x$', color = 'grey', ha = 'center', va = 'center', size = 'x-large')
ax.add_artist(FancyArrowPatch((0.02, 0.03), (0.02, 0.99), mutation_scale = 3., arrowstyle='|-|', color = 'grey'))
ax.text(0.065, 0.5, r'$L_z$', color = 'grey', ha = 'center', va = 'center', size = 'x-large')

width = 0.03
ax.add_artist(Wedge((r0[0] - R, r0[1] + h), width, 90, 270, ec = None, fc = 'crimson', alpha = 0.3))
ax.add_artist(plt.Rectangle((r0[0] - R, r0[1] + h - width), 2 * R, 2 * width, ec = None, fc = 'crimson', alpha = 0.3))
ax.add_artist(Wedge((r0[0] + R, r0[1] + h), width, 270, 90, ec = None, fc = 'crimson', alpha = 0.3))
ax.text(r0[0] - R - width - 0.04, r0[1] + h, r'$V_1$', color = 'crimson', ha = 'center', va = 'center', size = 'x-large')
ax.add_artist(Wedge((r0[0] - R, r0[1] - h), width, 90, 270, ec = None, fc = 'limegreen', alpha = 0.3))
ax.add_artist(plt.Rectangle((r0[0] - R, r0[1] - h - width), 2 * R, 2 * width, ec = None, fc = 'limegreen', alpha = 0.3))
ax.add_artist(Wedge((r0[0] + R, r0[1] - h), width, 270, 90, ec = None, fc = 'limegreen', alpha = 0.3))
ax.text(r0[0] - R - width - 0.04, r0[1] - h, r'$V_2$', color = 'limegreen', ha = 'center', va = 'center', size = 'x-large')
ax.text(0.4, 0.4, r'$V_3$', color = 'black', ha = 'center', va = 'center', size = 'x-large')

ax.set_xlim(-0.2, 1.2)
ax.set_ylim(-0.2, 1.2)
ax.axis('off')

fig.tight_layout()
#fig.show()
fig.savefig('shield.pdf')
