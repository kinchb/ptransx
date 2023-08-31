import numpy as np

import h5py

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/pdflatex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

M = 10.0

dir = './h3d_10_a0_01/'

def get_spec(flnm):
    infile = h5py.File(dir + flnm, 'r')

    spec = np.array(infile['/data']).T

    infile.close()

    lum = np.zeros(len(spec[0]))

    for i in range(0, len(spec[0])):
        for j in range(0, len(spec)):
            lum[i] += spec[j][i]

    j = 0
    lum = np.zeros(len(spec[0]))
    for i in range(0, len(spec[0])):
        lum[i] = spec[j][i]

    e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)

    return e_grid, lum * 2.417989e14

def count_photons(e, f):
    e_min = 2.0e3
    for i in range(0, len(e)):
        if (e[i] > e_min):
            break
    return np.trapz(f[i:]/e[i:], e[i:])

# e_grid, spec = get_spec('no_refl_scat_spec.1000.0.h5')

# plt.plot(e_grid/1000.0, e_grid * spec, 'g-')

# print(count_photons(e_grid, spec))

# e_grid, spec = get_spec('resp_scat_spec.1000.0.h5')

# print(count_photons(e_grid, spec))

# plt.plot(e_grid/1000.0, e_grid * spec, 'k-')

# e_grid, spec = get_spec('no_scat_spec.1000.0.h5')

# print(count_photons(e_grid, spec))

# plt.plot(e_grid/1000.0, e_grid * spec, 'k-')

# e_grid, spec = get_spec('scat_spec.1000.0.h5')

# plt.plot(e_grid/1000.0, e_grid * spec, 'k-')

# e_grid, spec = get_spec('scat_spec.1000.1.h5')

# plt.plot(e_grid/1000.0, e_grid * spec, 'k-')

e_grid, spec = get_spec('scat_spec.1000.5.h5')

plt.plot(e_grid/1000.0, e_grid * spec, 'k-')

plt.loglog()

plt.xlim([0.1, 300.0])
plt.ylim([1.0e33, 1.0e36])

plt.xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\varepsilon L_\varepsilon \left[\mathrm{erg}\ \mathrm{s}^{-1}\ \mathrm{cm}^{-2}\right]$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('flux_3.pdf', bbox_inches = 'tight')

plt.show()
