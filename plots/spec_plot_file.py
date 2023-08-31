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

e_min = 2.0
e_max = 30.0

def get_spec(spec_flnm):
    infile    = h5py.File(spec_flnm, 'r')
    spec      = np.array(infile['/data']).T
    infile.close()
    lum    = np.zeros(len(spec[0]))
    e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)
    for i in range(0, len(spec[0])):
        for j in range(0, len(spec)):
            lum[i] += spec[j][i]
    return (e_grid, lum * 2.417989e14)

spec_flnm = '/Users/kinch/spin_analysis/data/h3d_10_a09_10_gcc/scat_specs/v2/scat_spec.1000.0.h5'

nu, Lnu = get_spec(spec_flnm)

plt.plot(nu/1.0e3, nu * Lnu, 'k-')

plt.xlim([1.0e-1, 1.0e4])
plt.ylim([1.0e35, 1.0e38])

plt.xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\varepsilon L_\varepsilon \left[\mathrm{erg}\ \mathrm{s}^{-1}\ \mathrm{cm}^{-2}\right]$')

plt.loglog()

plt.tight_layout()

plt.savefig('spec_lum_no_hi_Te.pdf')

plt.show()

"""
cdf = np.zeros(len(nu))
for i in range(1, len(nu)):
    cdf[i] = cdf[i-1] + np.trapz(Lnu[:i], nu[:i])
cdf /= cdf[-1]

plt.plot(nu/1.0e3, cdf)
plt.semilogx()

plt.show()
"""