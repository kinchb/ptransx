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

    e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)

    return e_grid, lum * 2.417989e14

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

e_min = 2.0
e_max = 30.0

specs = []

for i in range(0, 51, 1):
    tmp = get_spec('../panptx_data/h3d_10_a0_01_old/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

gammas = np.zeros(len(specs))

for i in range(0, len(gammas)):
    gammas[i] = pl_fit(e_grid, specs[i], e_min * 1.0e3, e_max * 1.0e3)

t_old = np.load('../panptx_data/h3d_10_a0_01_old/old_a0_data_times.npz')['times']
gammas_old = gammas.copy()

# print gammas_old.std()/gammas_old.mean()
print(gammas_old.mean())

# --- --- ---

t = []
specs = []

times = np.load('../panptx_data/h3d_10_a0_01/h3d_10_a0_01_times.npz')['times']

for i in range(50, 2001, 50):
    t.append(times[i])
    tmp = get_spec('../panptx_data/h3d_10_a0_01/scat_spec.' + repr(i) + '.0.h5')
    e_grid = tmp[0]
    specs.append(tmp[1])

t = np.array(t)

gammas = np.zeros(len(specs))

for i in range(0, len(gammas)):
    gammas[i] = pl_fit(e_grid, specs[i], e_min * 1.0e3, e_max * 1.0e3)

# print gammas.std()/gammas.mean()
print(gammas.mean())

plt.plot([0.0, 0.0], [-1.0, 10.0], 'k--')

plt.plot(t_old - 10000.0, gammas_old, 'k-')

plt.plot(t - 10000.0, gammas, 'k-')

plt.xlim([-1000, 2000])
plt.ylim([2.0, 3.0])

plt.xlabel(r'$t/M$')
plt.ylabel(r'$\Gamma$')

plt.tight_layout()

f = plt.gcf()
f.savefig('gammas.pdf', bbox_inches = 'tight')

plt.show()
