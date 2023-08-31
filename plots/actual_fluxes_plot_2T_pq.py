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

inc_ndx = -1

def pl_fit(e, f, e_min, e_max):
    for i in range(0, len(e)):
        if (e[i] > e_min):
            break
    min_ndx = i - 1

    for i in range(0, len(e)):
        if (e[i] > e_max):
            break
    max_ndx = i + 1

    e = e[min_ndx:max_ndx]
    f = f[min_ndx:max_ndx]

#   print 'gamma = ' + repr(-np.polyfit(np.log10(e), np.log10(f/e), 1)[0])

    return -np.polyfit(np.log10(e), np.log10(f/e), 1)[0]

#   p = np.poly1d(np.polyfit(np.log10(e), np.log10(f), 1))

#   return e, 10.0**p(np.log10(e)), f

"""
def get_spec(flnm):
    infile = h5py.File(flnm, 'r')

    spec = np.array(infile['/data']).T

    infile.close()

    lum = np.zeros(len(spec[0]))

    for i in range(0, len(spec[0])):
        for j in range(0, len(spec)):
            lum[i] += spec[j][i]

    e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)

    return e_grid, lum * 2.417989e14
"""

def print_angle(ndx):
    angle = np.arccos(np.linspace(1, -1, 41, endpoint = True)) * (180.0/np.pi)
    print(angle[ndx])

print_angle(inc_ndx)

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

specs = []

for i in range(50, 1001, 50):
#   tmp = get_spec('../panptx_data/h3d_10_a0_01/scat_spec.' + repr(i) + '.0.h5')
    tmp = get_spec('../panptx_data/h3d_10_a0_01/scat_spec.' + repr(i) + '.0.h5', inc_ndx)
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_2 = mean_spec.copy()

specs = []

for i in range(50, 1001, 50):
#   tmp = get_spec('../panptx_data/h3d2T_10_a0_01/scat_spec.' + repr(i) + '.0.h5')
    tmp = get_spec('../panptx_data/h3d2T_10_a0_01/scat_spec.' + repr(i) + '.0.h5', inc_ndx)
    e_grid = tmp[0]
    specs.append(tmp[1])

mean_spec = np.zeros(len(e_grid))

for i in range(0, len(specs)):
    mean_spec += specs[i]

mean_spec /= len(specs)

mean_spec_3 = mean_spec.copy()

plt.plot(e_grid/1000.0, e_grid * mean_spec_2, 'r-', label = r'$\mathrm{1T}$')
plt.plot(e_grid/1000.0, e_grid * mean_spec_3, 'k-', label = r'$\mathrm{2T}$')

# plt.plot(e_grid/1000.0, mean_spec_3/e_grid, 'k-', label = r'$\mathrm{2T}$')
# plt.plot(e_grid/1000.0, mean_spec_2/e_grid, 'r-', label = r'$\mathrm{1T}$')

gamma_3 = pl_fit(e_grid, mean_spec_3, 3.0e3, 70.0e3)
gamma_2 = pl_fit(e_grid, mean_spec_2, 3.0e3, 70.0e3)

print(gamma_3)
print(gamma_2)

plt.loglog()

plt.xlim([0.1, 300.0])
plt.ylim([1.0e34, 1.0e36])

plt.xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\varepsilon F_\varepsilon \left[\mathrm{arbitrary\ units}\right]$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dist_obs_spec_inc_2T.pdf', bbox_inches = 'tight')

plt.show()
