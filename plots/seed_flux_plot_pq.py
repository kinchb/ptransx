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

def get_seed(flnm):
    infile = h5py.File(dir + flnm, 'r')

    seed_spec_top = np.array(infile['/seed_spec_top'])/(4.135668e-15/np.pi)
    seed_spec_bot = np.array(infile['/seed_spec_bot'])/(4.135668e-15/np.pi)

    Ne   = len(seed_spec_top)
    Nphi = len(seed_spec_top[0])
    Nr   = len(seed_spec_top[0][0])

    new_seed_spec_top = np.zeros((Nr, Nphi, Ne))
    new_seed_spec_bot = np.zeros((Nr, Nphi, Ne))

    for i in range(0, Nr):
        for j in range(0, Nphi):
            for k in range(0, Ne):
                new_seed_spec_top[i][j][k] = seed_spec_top[k][j][i]
                new_seed_spec_bot[i][j][k] = seed_spec_bot[k][j][i]

    e_grid = np.logspace(1, 5, Ne, endpoint = True)

    return e_grid, new_seed_spec_top

e_grid, c_seed = get_seed('spec_out_c_hr.h5')
e_grid, l_seed = get_seed('spec_out_l_hr.h5')

r_ndx = 70

seed = c_seed[r_ndx][0] + l_seed[r_ndx][0]

plt.plot(e_grid/1000.0, seed, 'k-')

plt.loglog()

plt.xlim([0.1, 300.0])
# plt.ylim([1.0e33, 1.0e36])

plt.xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\varepsilon L_\varepsilon \left[\mathrm{erg}\ \mathrm{s}^{-1}\ \mathrm{cm}^{-2}\right]$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

# f = plt.gcf()
# f.savefig('flux_2.pdf', bbox_inches = 'tight')

plt.show()
