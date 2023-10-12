import numpy as np
from scipy import integrate
import time

np.seterr(divide = 'ignore')
t0 = time.time()
lnum = 10

R = 1.
h = 1.
I = 1.
Lx = 10.
Ly = 10.
Lz = 10.
x0 = 7.
y0 = 8.
z0 = 5.

TEnorm = lambda m, n, p: (m ** 2 / Lx ** 2 + n ** 2 / Ly ** 2) * Lx * Ly * Lz / 8 * 2 ** ((m == 0) + (n == 0))
TE_Efield = lambda m, n, p, x, y, z: np.nan_to_num(1 / np.sqrt(TEnorm(m, n, p))) * np.array([
    n / Ly * np.cos(m * np.pi * x / Lx) * np.sin(n * np.pi * y / Ly) * np.sin(p * np.pi * z / Lz),
    -m / Lx * np.sin(m * np.pi * x / Lx) * np.cos(n * np.pi * y / Ly) * np.sin(p * np.pi * z / Lz),
    0.])

TMnorm = lambda m, n, p: (m ** 2 * p ** 2 / Lx ** 2 / Lz ** 2 + n ** 2 * p ** 2 / Ly ** 2 / Lz ** 2 + (m ** 2 / Lx ** 2 + n ** 2 / Ly ** 2) ** 2) * Lx * Ly * Lz / 8 * 2 ** (p == 0)
TM_Efield = lambda m, n, p, x, y, z: np.nan_to_num(1 / np.sqrt(TMnorm(m, n, p))) * np.array([
    m * p / Lx / Lz * np.cos(m * np.pi * x / Lx) * np.sin(n * np.pi * y / Ly) * np.sin(p * np.pi * z / Lz),
    n * p / Ly / Lz * np.sin(m * np.pi * x / Lx) * np.cos(n * np.pi * y / Ly) * np.sin(p * np.pi * z / Lz),
    -(m ** 2 / Lx ** 2 + n ** 2 / Ly ** 2) * np.sin(m * np.pi * x / Lx) * np.sin(n * np.pi * y / Ly) * np.cos(p * np.pi * z / Lz)])

Vquad = lambda x, y, z: I * R ** 2 * h / (2 * (x ** 2 + y ** 2 + z ** 2) ** 1.5) * (1 - 3 * z ** 2 /  (x ** 2 + y ** 2 + z ** 2))
Vexact = lambda x, y, z: I / (4 * np.pi) * integrate.nquad(lambda r, theta: r * (z + h) / ((x - r * np.cos(theta)) ** 2 + (y - r * np.sin(theta)) ** 2 + (z + h) ** 2) ** 1.5
                                                           - r * (z - h) / ((x - r * np.cos(theta)) ** 2 + (y - r * np.sin(theta)) ** 2 + (z - h) ** 2) ** 1.5, [(0, R), (0, 2 * np.pi)])[0]

TE_boundint = lambda m, n, p: (integrate.nquad(lambda y, z: TE_Efield(m, n, p, Lx, y, z)[0] * Vquad(Lx - x0, y - y0, z - z0), [(0, Ly), (0, Lz)])[0]
                               - integrate.nquad(lambda y, z: TE_Efield(m, n, p, 0, y, z)[0] * Vquad(-x0, y - y0, z - z0), [(0, Ly), (0, Lz)])[0]
                               + integrate.nquad(lambda x, z: TE_Efield(m, n, p, x, Ly, z)[1] * Vquad(x - x0, Ly - y0, z - z0), [(0, Lx), (0, Lz)])[0]
                               - integrate.nquad(lambda x, z: TE_Efield(m, n, p, x, 0, z)[1] * Vquad(x - x0, -y0, z - z0), [(0, Lx), (0, Lz)])[0])
TM_boundint = lambda m, n, p: (integrate.nquad(lambda y, z: TM_Efield(m, n, p, Lx, y, z)[0] * Vquad(Lx - x0, y - y0, z - z0), [(0, Ly), (0, Lz)])[0]
                               - integrate.nquad(lambda y, z: TM_Efield(m, n, p, 0, y, z)[0] * Vquad(-x0, y - y0, z - z0), [(0, Ly), (0, Lz)])[0]
                               + integrate.nquad(lambda x, z: TM_Efield(m, n, p, x, Ly, z)[1] * Vquad(x - x0, Ly - y0, z - z0), [(0, Lx), (0, Lz)])[0]
                               - integrate.nquad(lambda x, z: TM_Efield(m, n, p, x, 0, z)[1] * Vquad(x - x0, -y0, z - z0), [(0, Lx), (0, Lz)])[0]
                               + integrate.nquad(lambda x, y: TM_Efield(m, n, p, x, y, Lz)[2] * Vquad(x - x0, y - y0, Lz - z0), [(0, Lx), (0, Ly)])[0]
                               - integrate.nquad(lambda x, y: TM_Efield(m, n, p, x, y, 0)[2] * Vquad(x - x0, y - y0, -z0), [(0, Lx), (0, Ly)])[0])
TM_fluxint = lambda m, n, p: I * (integrate.nquad(lambda r, theta: r * TM_Efield(m, n, p, x0 + r * np.cos(theta), y0 + r * np.sin(theta), z0 + h)[2], [(0, R), (0, 2 * np.pi)])[0]
                                  - integrate.nquad(lambda r, theta: r * TM_Efield(m, n, p, x0 + r * np.cos(theta), y0 + r * np.sin(theta), z0 - h)[2], [(0, R), (0, 2 * np.pi)])[0])

TE_boundint_exact = lambda m, n, p: (integrate.nquad(lambda y, z: TE_Efield(m, n, p, Lx, y, z)[0] * Vexact(Lx - x0, y - y0, z - z0), [(0, Ly), (0, Lz)])[0]
                                     - integrate.nquad(lambda y, z: TE_Efield(m, n, p, 0, y, z)[0] * Vexact(-x0, y - y0, z - z0), [(0, Ly), (0, Lz)])[0]
                                     + integrate.nquad(lambda x, z: TE_Efield(m, n, p, x, Ly, z)[1] * Vexact(x - x0, Ly - y0, z - z0), [(0, Lx), (0, Lz)])[0]
                                     - integrate.nquad(lambda x, z: TE_Efield(m, n, p, x, 0, z)[1] * Vexact(x - x0, -y0, z - z0), [(0, Lx), (0, Lz)])[0])
TM_boundint_exact = lambda m, n, p: (integrate.nquad(lambda y, z: TM_Efield(m, n, p, Lx, y, z)[0] * Vexact(Lx - x0, y - y0, z - z0), [(0, Ly), (0, Lz)])[0]
                                     - integrate.nquad(lambda y, z: TM_Efield(m, n, p, 0, y, z)[0] * Vexact(-x0, y - y0, z - z0), [(0, Ly), (0, Lz)])[0]
                                     + integrate.nquad(lambda x, z: TM_Efield(m, n, p, x, Ly, z)[1] * Vexact(x - x0, Ly - y0, z - z0), [(0, Lx), (0, Lz)])[0]
                                     - integrate.nquad(lambda x, z: TM_Efield(m, n, p, x, 0, z)[1] * Vexact(x - x0, -y0, z - z0), [(0, Lx), (0, Lz)])[0]
                                     + integrate.nquad(lambda x, y: TM_Efield(m, n, p, x, y, Lz)[2] * Vexact(x - x0, y - y0, Lz - z0), [(0, Lx), (0, Ly)])[0]
                                     - integrate.nquad(lambda x, y: TM_Efield(m, n, p, x, y, 0)[2] * Vexact(x - x0, y - y0, -z0), [(0, Lx), (0, Ly)])[0])

def overlap(params, exact = False):
    mode, m, n, p = params
    if not exact:
        if mode == 'TE':
            overlap = TE_boundint(m, n, p)
        if mode == 'TM':
            overlap = TM_boundint(m, n, p) + TM_fluxint(m, n, p)
        #print(params, time.time() - t0)
        return overlap
    else:
        if mode == 'TE':
            overlap = TE_boundint_exact(m, n, p)
        if mode == 'TM':
            overlap = TM_boundint_exact(m, n, p) + TM_fluxint(m, n, p)
        #print(params, time.time() - t0)
        return overlap

TEcoeff = np.zeros((lnum, lnum, lnum))
TMcoeff = np.zeros((lnum, lnum, lnum))
for m in range(1, lnum):
    for n in range(1, lnum):
        TEcoeff[0, m, n] = overlap(['TE', 0, m, n])
        TEcoeff[m, 0, n] = overlap(['TE', m, 0, n])
        TMcoeff[m, n, 0] = overlap(['TM', m, n, 0])
for m in range(1, lnum):
    for n in range(1, lnum):
        for p in range(1, lnum):
            TEcoeff[m, n, p] = overlap(['TE', m, n, p])
            TMcoeff[m, n, p] = overlap(['TM', m, n, p])

TE_Bfield = lambda m, n, p, x, y, z: -1j * np.nan_to_num(1 / np.sqrt(TEnorm(m, n, p)) / np.sqrt(m ** 2 / Lx ** 2 + n ** 2 / Ly ** 2 + p ** 2 / Lz ** 2)) * np.array([
    m * p / Lx / Lz * np.sin(m * np.pi * x / Lx) * np.cos(n * np.pi * y / Ly) * np.cos(p * np.pi * z / Lz),
    n * p / Ly / Lz * np.cos(m * np.pi * x / Lx) * np.sin(n * np.pi * y / Ly) * np.cos(p * np.pi * z / Lz),
    -(m ** 2 / Lx ** 2 + n ** 2 / Ly ** 2) * np.cos(m * np.pi * x / Lx) * np.cos(n * np.pi * y / Ly) * np.sin(p * np.pi * z / Lz)])

TM_Bfield = lambda m, n, p, x, y, z: 1j * np.nan_to_num(1 / np.sqrt(TMnorm(m, n, p)) * np.sqrt(m ** 2 / Lx ** 2 + n ** 2 / Ly ** 2 + p ** 2 / Lz ** 2)) * np.array([
    n / Ly * np.sin(m * np.pi * x / Lx) * np.cos(n * np.pi * y / Ly) * np.cos(p * np.pi * z / Lz),
    -m / Lx * np.cos(m * np.pi * x / Lx) * np.sin(n * np.pi * y / Ly) * np.cos(p * np.pi * z / Lz),
    0.])

Bresp = lambda x, y, z, num: -np.sum(TEcoeff[:num,:num,:num,None] * [[[TE_Bfield(m, n, p, x, y, z) * np.nan_to_num(1 / np.pi / np.sqrt(m ** 2 / Lx ** 2 + n ** 2 / Ly ** 2 + p ** 2 / Lz ** 2)) for p in range(num)] for n in range(num)] for m in range(num)]
       + TMcoeff[:num,:num,:num,None] * [[[TM_Bfield(m, n, p, x, y, z) * np.nan_to_num(1 / np.pi / np.sqrt(m ** 2 / Lx ** 2 + n ** 2 / Ly ** 2 + p ** 2 / Lz ** 2)) for p in range(num)] for n in range(num)] for m in range(num)], axis = (0, 1, 2))

print(Bresp(x0, y0, z0, lnum))
print(lnum, time.time() - t0)
