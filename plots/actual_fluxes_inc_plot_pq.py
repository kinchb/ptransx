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

def print_angle(ndx):
    angle = np.arccos(np.linspace(1, -1, 41, endpoint = True)) * (180.0/np.pi)
    print(angle[ndx])

def get_spec(flnm, inc_ndx):


    infile = h5py.File(flnm, 'r')

    spec = np.array(infile['/data']).T

    infile.close()

    lum = np.zeros(len(spec[0]))

    for i in range(0, len(spec[0])):
        for j in range(0, len(spec)):
            lum[i] += spec[j][i]

    j = inc_ndx
    lum = np.zeros(len(spec[0]))
    for i in range(0, len(spec[0])):
        lum[i] = spec[j][i]

    e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)

    return e_grid, lum * 2.417989e14

"""
specs = []

for i in range(0, 51, 1):
    tmp = get_spec('../panptx_data/h3d_10_a0_01_old/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_1 = mean_spec.copy()
"""

specs = []

for i in range(50, 2001, 50):
    tmp = get_spec('../panptx_data/h3d_10_a0_01/scat_spec.' + repr(i) + '.0.h5', 0)
    e_grid = tmp[0]
    specs.append(tmp[1])
print_angle(0)

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_1 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 2001, 50):
    tmp = get_spec('../panptx_data/h3d_10_a0_01/scat_spec.' + repr(i) + '.0.h5', 3)
    e_grid = tmp[0]
    specs.append(tmp[1])
print_angle(3)

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_2 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 2001, 50):
    tmp = get_spec('../panptx_data/h3d_10_a0_01/scat_spec.' + repr(i) + '.0.h5', 10)
    e_grid = tmp[0]
    specs.append(tmp[1])
print_angle(10)

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_3 = mean_spec.copy()

# --- --- ---

specs = []

for i in range(50, 2001, 50):
    tmp = get_spec('../panptx_data/h3d_10_a0_01/scat_spec.' + repr(i) + '.0.h5', 20)
    e_grid = tmp[0]
    specs.append(tmp[1])
print_angle(20)

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_4 = mean_spec.copy()

# --- --- ---

plt.plot(e_grid/1000.0, e_grid * mean_spec_1, 'g-', label = r'$i = 0^\circ$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_2, 'b-', label = r'$i = 30^\circ$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_3, 'c-', label = r'$i = 60^\circ$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_4, 'm-', label = r'$i = 90^\circ$')

plt.loglog()

plt.xlim([0.1, 300.0])
plt.ylim([3.0e32, 5.0e35])

plt.xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\varepsilon F_\varepsilon \left[\mathrm{arbitrary\ units}\right]$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dist_obs_spec_inc.pdf', bbox_inches = 'tight')

plt.show()
