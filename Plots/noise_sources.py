from curlyBrace import curlyBrace
from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 14})

fig, ax = plt.subplots(figsize = (6., 6.))

hbar = 6.62607015e-34 / (2 * np.pi) # in kg*m^2/s
kB = 1.380649e-23 # in kg*m^2/s^2/K
mu0 = 1.256637e-6 # in kg*m/s^2/A^2
c = 2.9979e8 # in m/s

T = 0.01 # in K
m = 0.001 # in kg
gamma = 2 * np.pi * 1e-8 # in Hz
rho = 1e2 # in kg/m^3
kappa = 5
f0 = 10.
eta_res = np.sqrt(4 * kB * T * m * gamma / kappa / hbar)
eta_broad = lambda omega0: np.sqrt(m) * omega0

chi = lambda omega, omega0: 1 / (m * (omega0 ** 2 - omega ** 2 - 1j * gamma * omega))
Bfactor = lambda omega0: 2 * mu0 * rho / (3 * m ** 2 * omega0 ** 2)
SBB_therm = lambda omega, omega0: Bfactor(omega0) * np.full_like(omega, 4 * kB * T * m * gamma)
SBB_imp_res = lambda omega, omega0: Bfactor(omega0) * kappa * hbar / eta_res ** 2 / np.abs(chi(omega, omega0)) ** 2
SBB_imp_broad = lambda omega, omega0: Bfactor(omega0) * kappa * hbar / eta_broad(omega0) ** 2 / np.abs(chi(omega, omega0)) ** 2
SBB_back_res = lambda omega, omega0: Bfactor(omega0) * np.full_like(omega, kappa * hbar * eta_res ** 2)
SBB_back_broad = lambda omega, omega0: Bfactor(omega0) * np.full_like(omega, kappa * hbar * eta_broad(omega0) ** 2)
SBB_res = lambda omega, omega0: SBB_therm(omega, omega0) + SBB_imp_res(omega, omega0) + SBB_back_res(omega, omega0)
SBB_broad = lambda omega, omega0: SBB_therm(omega, omega0) + SBB_imp_broad(omega, omega0) + SBB_back_broad(omega, omega0)

class CustomTicker(ticker.LogFormatterSciNotation): 
    def __call__(self, x, pos = None): 
        if x not in np.concatenate((0.1 * np.arange(1, 10), np.arange(1, 10), 10 * np.arange(1, 10))): 
            return ticker.LogFormatterSciNotation.__call__(self, x, pos = None) 
        else:
            return "{x:g}".format(x = x)

xlim1 = 3.
xlim2 = 30.
ylim1 = 1e-39
ylim2 = 1e-32
ax.set_xlim(xlim1, xlim2)
ax.set_ylim(ylim1, ylim2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$f\,[\mathrm{Hz}]$')
ax.set_ylabel(r'$S_{BB}\,[\mathrm{T}^2/\mathrm{Hz}]$')
ax.xaxis.set_major_locator(ticker.FixedLocator([3, 10, 30]))
ax.xaxis.set_major_formatter(CustomTicker())
ax.xaxis.set_minor_formatter(ticker.NullFormatter())

ax.tick_params(which = 'both', direction = 'in')
secxax = ax.secondary_xaxis('top')
secxax.tick_params(which = 'both', direction = 'in')
plt.setp(secxax.get_xticklabels(), visible = False)
secxax.xaxis.set_minor_formatter(ticker.NullFormatter())
secyax = ax.secondary_yaxis('right')
secyax.tick_params(which = 'both', direction = 'in')
plt.setp(secyax.get_yticklabels(), visible = False)

colors = [(0.921569, 0.494118, 0.431373), (1., 0.721569, 0.219608), (0.317647, 0.654902, 0.752941)]
freqs = np.logspace(np.log10(xlim1), np.log10(xlim2), 1000)
##ax.plot(freqs, SBB_imp_broad(2 * np.pi * freqs, 2 * np.pi * f0), color = 'limegreen', linestyle = '--', label = r'$S_{BB}^\mathrm{imp}$ (broad)')
##ax.plot(freqs, SBB_imp_res(2 * np.pi * freqs, 2 * np.pi * f0), color = 'limegreen', label = r'$S_{BB}^\mathrm{imp}$ (res)')
##ax.plot(freqs, SBB_back_broad(2 * np.pi * freqs, 2 * np.pi * f0), color = 'royalblue', linestyle = '--', label = r'$S_{BB}^\mathrm{back}$ (broad)')
##ax.plot(freqs, SBB_therm(2 * np.pi * freqs, 2 * np.pi * f0), color = 'crimson', label = r'$S_{BB}^\mathrm{th}=S_{BB}^\mathrm{back}$ (res)')
##ax.plot(freqs, SBB_broad(2 * np.pi * freqs, 2 * np.pi * f0), color = 'k', linestyle = '--', label = r'$S_{BB}^\mathrm{tot}$ (broad)')
##ax.plot(freqs, SBB_res(2 * np.pi * freqs, 2 * np.pi * f0), color = 'k', label = r'$S_{BB}^\mathrm{tot}$ (res)')
ax.plot(freqs, SBB_imp_broad(2 * np.pi * freqs, 2 * np.pi * f0), color = colors[2], linestyle = '--', label = r'Imprecision (broad.)')
ax.plot(freqs, SBB_imp_res(2 * np.pi * freqs, 2 * np.pi * f0), color = colors[2], label = r'Imprecision (res.)')
ax.plot(freqs, SBB_back_broad(2 * np.pi * freqs, 2 * np.pi * f0), color = colors[1], linestyle = '--', label = r'Back-action (broad.)')
ax.plot(0, 0, color = 'white', label = 'Back-action (res.)')
ax.plot(freqs, SBB_therm(2 * np.pi * freqs, 2 * np.pi * f0), color = colors[0], label = r'Thermal')
ax.plot(freqs, SBB_broad(2 * np.pi * freqs, 2 * np.pi * f0), color = 'k', linestyle = '--', label = r'Total (broad.)')
ax.plot(freqs, SBB_res(2 * np.pi * freqs, 2 * np.pi * f0), color = 'k', label = r'Total (res.)')

handles, labels = ax.get_legend_handles_labels()
handles[4] = handles[3]
handles.reverse()
labels.reverse()
ax.legend(handles, labels, fontsize = 13, loc = 'lower left').set_zorder(3)
ax.plot([3.25, 3.85], [4.7e-38, 4.7e-38], color = colors[0], zorder = 4)
results = curlyBrace(fig, ax, [4.09, 2.5e-38], [4.09, 8.8e-38], k_r = 0.08, color = 'k', lw = 0.8, zorder = 4)

fig.tight_layout()
#fig.show()
fig.savefig('noise_sources.pdf')
