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

def get_spec(flnm):
    infile = h5py.File(flnm, 'r')

    spec = np.array(infile['/data']).T

    infile.close()

    lum = np.zeros(len(spec[0]))

    for i in range(0, len(spec[0])):
        for j in range(0, len(spec)):
            lum[i] += spec[j][i]

    e_grid = np.logspace(0, 8, len(spec[0]), endpoint = True)

    return e_grid, lum * 2.417989e14

specs = []

for i in range(50, 201, 10):
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_c1.5/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_1 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 201, 10):
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_redo/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_2 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 201, 10):
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_c0.5/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_3 = mean_spec.copy()

# --- --- ---

plt.plot(e_grid/1000.0, e_grid * mean_spec_1, 'r-', label = r'$B^2/\rho > 1.5$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_2, 'k-', label = r'$B^2/\rho > 1$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_3, 'b-', label = r'$B^2/\rho > 0.5$')

plt.loglog()

plt.xlim([0.1, 1.0e5])
plt.ylim([1.0e32, 2.0e37])

plt.xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\varepsilon L_\varepsilon \left[\mathrm{erg}\ \mathrm{s}^{-1}\ \mathrm{cm}^{-2}\right]$')

plt.legend(frameon = False, loc = 'lower left')

plt.tight_layout()

f = plt.gcf()
f.savefig('cutoff_obs_spec.pdf', bbox_inches = 'tight')

plt.show()

