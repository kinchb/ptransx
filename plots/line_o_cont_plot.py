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

def get_spec(flnm, ndx):
    infile = h5py.File(flnm, 'r')

    spec = np.array(infile['/data']).T

    angle = np.arccos(np.linspace(1, -1, len(spec), endpoint = True)) * (180.0/np.pi)

    print(angle[ndx])

    e_grid = np.logspace(1, 5, len(spec[0]), endpoint = True)

    return e_grid, spec[ndx] * 2.417989e14 * 1.0e3 * (2.0/len(spec))

# --- --- --- line/continuum plot --- --- ---

flnm_l = '../data/h3d_10_a09_10_gcc/postprocess/run_1000/scat_spec_l_hr.h5'
flnm_c = '../data/h3d_10_a09_10_gcc/postprocess/run_1000/scat_spec_c_hr.h5'

e_grid, spec_l = get_spec(flnm_l, 0)
e_grid, spec_c = get_spec(flnm_c, 0)

plt.plot(e_grid/1.0e3, 1.0 + spec_l/spec_c, 'r-', label = r'$i = 0^\circ$')

e_grid, spec_l = get_spec(flnm_l, 3)
e_grid, spec_c = get_spec(flnm_c, 3)

plt.plot(e_grid/1.0e3, 1.0 + spec_l/spec_c, 'g-', label = r'$i = 30^\circ$')

e_grid, spec_l = get_spec(flnm_l, 6)
e_grid, spec_c = get_spec(flnm_c, 6)

plt.plot(e_grid/1.0e3, 1.0 + spec_l/spec_c, 'b-', label = r'$i = 45^\circ$')

e_grid, spec_l = get_spec(flnm_l, 17)
e_grid, spec_c = get_spec(flnm_c, 17)

plt.plot(e_grid/1.0e3, 1.0 + spec_l/spec_c, 'c-', label = r'$i = 80^\circ$')

plt.xlim([0, 14])
plt.ylim([1.0, 1.5])

ax = plt.gca()
ax.xaxis.set_ticks(np.arange(0, 15, 2))

plt.legend(loc = 'best', frameon = False)

plt.xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')

plt.tight_layout()
plt.show()
