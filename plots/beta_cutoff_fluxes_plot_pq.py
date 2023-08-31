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
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_higher_beta/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_0 = mean_spec.copy()

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

mean_spec_1 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 201, 10):
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_lower_beta/scat_specs/scat_spec.' + repr(i) + '.0.h5')
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
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_lowest_beta/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_3 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 201, 10):
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_even_lower_beta/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_4 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 201, 10):
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_super_low_beta/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_5 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 201, 10):
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_1o2000_beta/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_6 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 201, 10):
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_1o5000_beta/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_7 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 201, 10):
    tmp = get_spec('/mnt/archive/h3d_10_a0_01_2cuts/scat_specs/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_8 = mean_spec.copy()

# --- --- ---

plt.plot(e_grid/1000.0, e_grid * mean_spec_0, 'r-', label = r'$1/\beta_\mathrm{min} = 50$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_1, 'g-', label = r'$1/\beta_\mathrm{min} = 100$')
# plt.plot(e_grid/1000.0, e_grid * mean_spec_2, 'b-', label = r'$1/\beta_\mathrm{min} = 150$')
# plt.plot(e_grid/1000.0, e_grid * mean_spec_3, 'c-', label = r'$1/\beta_\mathrm{min} = 300$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_4, 'b-', label = r'$1/\beta_\mathrm{min} = 500$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_5, 'c-', label = r'$1/\beta_\mathrm{min} = 1000$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_8, 'c--', label = r'$1/\beta_\mathrm{min} = 1000$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_6, 'y-', label = r'$1/\beta_\mathrm{min} = 2000$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_7, 'm-', label = r'$1/\beta_\mathrm{min} = 5000$')


plt.loglog()

plt.xlim([0.1, 1.0e3])
plt.ylim([1.0e36, 1.0e37])

plt.xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\varepsilon L_\varepsilon \left[\mathrm{erg}\ \mathrm{s}^{-1}\ \mathrm{cm}^{-2}\right]$')

plt.legend(frameon = False, loc = 'best', fontsize = 14)

plt.tight_layout()

f = plt.gcf()
f.savefig('beta_cutoff_obs_spec.pdf', bbox_inches = 'tight')

plt.show()

