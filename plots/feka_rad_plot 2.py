import numpy as np

import h5py

import matplotlib.pyplot as plt
custom_preamble = {
    "text.usetex": True,
    "text.latex.preamble": [
        r"\usepackage{amsmath}", # for the align enivironment
        ],
    }
plt.rcParams.update(custom_preamble)

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/pdflatex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

def isco_radius(a):
    z1 = 1.0 + ((1.0 - a*a)**(1.0/3.0))*((1.0 + a)**(1.0/3.0) + (1.0 - a)**(1.0/3.0))
    z2 = np.sqrt(3.0*a*a + z1*z1)
    r  = 3.0 + z2 - np.sqrt((3.0 - z1)*(3.0 + z1 + 2.0*z2))
    return r

def get_feka(flnm):
    infile = h5py.File(flnm, 'r')

    seed_spec_top = np.array(infile['/seed_spec_top'])/(4.135668e-15/np.pi)
    seed_spec_bot = np.array(infile['/seed_spec_top'])/(4.135668e-15/np.pi)

    Ne   = len(seed_spec_top)
    Nphi = len(seed_spec_top[0])
    Nr   = len(seed_spec_top[0][0])

    new_seed_spec_top = np.zeros((Nr, Nphi, Ne))
    new_seed_spec_bot = np.zeros((Nr, Nphi, Ne))

    for i in range(0, Nr):
        for j in range(0, Nphi):
            for k in range(0, Ne):
                new_seed_spec_top[i][j][k] = seed_spec_top[k][j][i]

    e_grid = np.logspace(1, 5, Ne, endpoint = True)

    feka = np.zeros(Nr)

    for i in range(0, Nr):
        for j in range(0, Nphi):
            feka[i] += (1.0/(2*Nphi))*(np.trapz(6.2415091e11*new_seed_spec_top[i][j]/e_grid, e_grid) + np.trapz(6.2415091e11*new_seed_spec_bot[i][j]/e_grid, e_grid))

    return feka

a = 0.9

harm_file = '../data/h3d_10_a09_10_gcc/postprocess/run_1000/RADFLUX_panptx_1000.h5'
line_file = '../data/h3d_10_a09_10_gcc/postprocess/run_1000/spec_out_l_hr.h5'

harm_data = h5py.File(harm_file, 'r')
r = np.array(harm_data['/r'])
harm_data.close()

feka = get_feka(line_file)

trunc_r    = np.zeros(int(len(r)/3))
trunc_feka = np.zeros(int(len(feka)/3))

j = 0
for i in range(0, len(r), 3):
    trunc_r[j]    = r[i]
    trunc_feka[j] = feka[i]
    j += 1

plt.plot([isco_radius(a), isco_radius(a)], [1.0e24, 1.0e29], 'k--')

plt.plot(trunc_r, trunc_feka, 'k-', label = r'$a = 0.9$')

# --- --- ---

plt.loglog()

plt.xlim([1, 50])
plt.ylim([1.0e25, 1.0e28])

# plt.legend(loc = 'upper right', bbox_to_anchor = (1.12, 1.06), framealpha = 1.0, handletextpad = 0.2)

plt.xlabel(r'$r/M$')
plt.ylabel(r'$\mathrm{Fe\ K}\alpha\ \mathrm{flux}\ \left[\mathrm{cm}^{-2}\ \mathrm{s}^{-1}\right]$')

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((1, 2, 3, 4, 5, 6, 8, 10, 20, 30, 40, 50))
ax.set_xticklabels(('1', '2', '3', '4', '5', '6', '8', '10', '20', '30', '40', '50'))

# f = plt.gcf()
# f.savefig('feka_rad.pdf', bbox_inches = 'tight', dpi = 400)

plt.tight_layout()
plt.show()
