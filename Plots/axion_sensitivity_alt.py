from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 14})

fig, ax = plt.subplots(figsize = (6., 6.))

hbar = 6.62607015e-34 / (2 * np.pi) # in kg*m^2/s
kB = 1.380649e-23 # in kg*m^2/s^2/K
mu0 = 1.256637e-6 # in kg*m/s^2/A^2
c = 2.9979e8 # in m/s

T = 0.01 # in K
ms = [1e-8, 1e-3] # in kg
gammas = [2 * np.pi * 1e-5, 2 * np.pi * 1e-8] # in Hz
rhos = [1e4, 1e2] # in kg/m^3
kappa = 5
rhoDM = 6e-11 # in T^2
hs = [0.01, 0.1] # in m
Rs = [0.01, 0.1] # in m
Ls = [0.1, 1.] # in m
const = 0.09
Tint = 365 * 86400 # in s
alpha = 1.
f0 = 100.
fres = 10.
SNR = 3

chi = lambda omega, omega0: 1 / (m * (omega0 ** 2 - omega ** 2 - 1j * gamma * omega))
Bfactor = lambda omega0: 2 * mu0 * rho / (3 * m ** 2 * omega0 ** 2)
SBB_therm = lambda omega, omega0: Bfactor(omega0) * np.full_like(omega, 4 * kB * T * m * gamma)
SBB_meas_res = lambda omega, omega0: Bfactor(omega0) * kappa * hbar / eta_res ** 2 / np.abs(chi(omega, omega0)) ** 2
SBB_meas_broad = lambda omega, omega0: Bfactor(omega0) * kappa * hbar / eta_broad(omega0) ** 2 / np.abs(chi(omega, omega0)) ** 2
SBB_back_res = lambda omega, omega0: Bfactor(omega0) * np.full_like(omega, kappa * hbar * eta_res ** 2)
SBB_back_broad = lambda omega, omega0: Bfactor(omega0) * np.full_like(omega, kappa * hbar * eta_broad(omega0) ** 2)
SBB_res = lambda omega, omega0: SBB_therm(omega, omega0) + SBB_meas_res(omega, omega0) + SBB_back_res(omega, omega0)
SBB_broad = lambda omega, omega0: SBB_therm(omega, omega0) + SBB_meas_broad(omega, omega0) + SBB_back_broad(omega, omega0)

cohT = lambda omega: 2 * np.pi * 1e6 / omega
g_res = lambda omega, omega0: np.sqrt(SNR * 3 * SBB_res(omega, omega0) / (rho * rhoDM)) * L ** 3 / (const * 2 * omega0 * (R ** 2 + h ** 2) ** 2.5) / np.sqrt(hbar * c) * 1.602e-10 / (cohT(omega0) * alpha * np.minimum(cohT(omega0) * alpha, cohT(omega))) ** 0.25
g_res_combined = lambda omega, omega0s: (SNR ** 2 / np.sum((rho * rhoDM * const ** 2 * 4 * omega0s[:, None] ** 2 * (R ** 2 + h ** 2) ** 5 * hbar * c * np.sqrt(cohT(omega0s[:, None]) * alpha * np.minimum(cohT(omega0s[:, None]) * alpha, cohT(omega))) / (3 * SBB_res(omega, omega0s[:, None]) * L ** 6 * 1.602e-10 ** 2)) ** 2, axis = 0)) ** 0.25
g_broad = lambda omega, omega0: np.sqrt(SNR * 3 * SBB_broad(omega, omega0) / (rho * rhoDM)) * L ** 3 / (const * 2 * omega0 * (R ** 2 + h ** 2) ** 2.5) / np.sqrt(hbar * c) * 1.602e-10 / (Tint * np.minimum(Tint, cohT(omega))) ** 0.25

xlim1 = 2 * np.pi * 6.582e-16 * 1
xlim2 = 2 * np.pi * 6.582e-16 * 1e3
ylim1 = 1e-14
ylim2 = 1e-3
ax.set_xlim(xlim1, xlim2)
ax.set_ylim(ylim1, ylim2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"$m_a$\,[eV]")
ax.set_ylabel(r'$g_{a\gamma}\,[\mathrm{GeV}^{-1}]$')

masses = np.logspace(np.log10(xlim1), np.log10(xlim2), 10000)
chandra = np.loadtxt('chandra.txt').T
ax.fill_between(chandra[0], chandra[1], 1, color = '0.98', zorder = 0.7)
ax.plot(chandra[0], chandra[1], color = '0.83', zorder = 0.7)
ax.text(1.3e-13, 1.5e-12, r'Chandra', ha = 'center', va = 'center')
ax.fill_between(masses, 5.2304e-12, 1, color = '0.94', zorder = 0.8)
ax.axhline(5.2304e-12, color = '0.79', zorder = 0.8)
ax.text(8e-13, 1.7e-11, r'SN1987A', ha = 'center', va = 'center')
ax.fill_between(masses, 6.6316e-11, 1, color = '0.86', zorder = 0.9)
ax.axhline(6.6316e-11, color = '0.71', zorder = 0.9)
ax.text(2e-12, 1.8e-10, r'CAST', ha = 'center', va = 'center')
snipehunt = np.loadtxt('snipehunt_axion.txt').T
ax.fill_between(2 * np.pi * 6.582e-16 * snipehunt[0], snipehunt[1], 1, color = '0.81', zorder = 1)
ax.plot(2 * np.pi * 6.582e-16 * snipehunt[0], snipehunt[1], color = '0.66', zorder = 1)
ax.text(9e-15, 5e-6, 'SNIPE\nHunt', ha = 'center', va = 'center')

colors = [[(0.317647, 0.654902, 0.752941)], [(1., 0.721569, 0.219608), (0.921569, 0.494118, 0.431373)]]
for i in range(2):
    m = ms[i]
    gamma = gammas[i]
    rho = rhos[i]
    h = hs[i]
    R = Rs[i]
    L = Ls[i]
    eta_res = np.sqrt(4 * kB * T * m * gamma / kappa / hbar)
    eta_broad = lambda omega0: np.sqrt(m) * omega0
    if i == 0:
        ax.plot(masses, g_broad(masses / 6.582e-16, 2 * np.pi * f0), color = colors[i][0], label = 'Broadband (existing)', zorder = 1)
    elif i == 1:
        ax.plot(masses, g_broad(masses / 6.582e-16, 2 * np.pi * f0), color = colors[i][0], label = 'Broadband (improved)', zorder = 1)
##        minmass = 2 * np.pi * 6.582e-16 * 3
##        maxmass = 2 * np.sqrt(2) * gamma * kB * T * 1e-6 / (np.pi * kappa * hbar) * Tint * 6.582e-16 / alpha + minmass
####        opt_res = lambda omega: np.exp(optimize.minimize(lambda omega0: np.log(g_res(omega, omega0)), omega).fun)
##        def scan(mass):
##            if mass <= minmass: return(g_res(mass / 6.582e-16, minmass / 6.582e-16))
##            elif mass <= maxmass: return(g_res(mass / 6.582e-16, mass / 6.582e-16))
##            else: return(g_res(mass / 6.582e-16, maxmass / 6.582e-16))
##        ax.plot(masses, [scan(mass) for mass in masses], color = colors[i][1], label = 'Scanning (improved)', zorder = 1)
        delomega = lambda omega: 4 * np.sqrt(2) * gamma * kB * T / (kappa * hbar * omega)
        omega0s = [2 * np.pi * 3]
        totalT = cohT(omega0s[0]) * alpha
        while totalT <= Tint:
            omega0s.append(omega0s[-1] + delomega(omega0s[-1]))
            totalT += cohT(omega0s[-1]) * alpha
        omega0s = np.array(omega0s)
        ax.plot(masses, g_res_combined(masses / 6.582e-16, omega0s), color = colors[i][1], label = 'Scanning (improved)', zorder = 1)
        ax.plot(masses, g_res(masses / 6.582e-16, 2 * np.pi * fres), color = colors[i][1], linestyle = '--', label = 'Single experiment', zorder = 1)
ax.legend(loc = 'lower left', ncol = 2, fontsize = 12)

class CustomTicker(ticker.LogFormatterSciNotation): 
    def __call__(self, x, pos = None): 
        if x not in [0.1, 1, 10]: 
            return ticker.LogFormatterSciNotation.__call__(self, x, pos = None) 
        else: 
            return "{x:g}".format(x = x)

ax.tick_params(which = 'both', direction = 'in')
ax.yaxis.set_major_formatter(CustomTicker())
#plt.setp(ax.get_yticklabels()[-2], visible = False)
ax.yaxis.set_minor_locator(ticker.FixedLocator(np.logspace(-15, 0, 16)))
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
secxax = ax.secondary_xaxis('top', functions = (lambda x: x / (2 * np.pi * 6.582e-16), lambda x: 2 * np.pi * 6.582e-16 * x), zorder = 1)
secxax.tick_params(which = 'both', direction = 'in')
secxax.xaxis.set_major_formatter(CustomTicker())
secxax.set_xlabel(r"$f_a$\,[Hz]")
secyax = ax.secondary_yaxis('right', zorder = 1)
secyax.tick_params(which = 'both', direction = 'in')
secyax.yaxis.set_minor_locator(ticker.FixedLocator(np.logspace(-15, 0, 16)))
secyax.yaxis.set_minor_formatter(ticker.NullFormatter())
plt.setp(secyax.get_yticklabels(), visible = False)

fig.tight_layout()
#fig.show()
fig.savefig('axion_sensitivity_alt.pdf')
