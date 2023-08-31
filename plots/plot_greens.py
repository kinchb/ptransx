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

from h5py_wrapper import GetFile

greens_data_filename = './h3d_10_a0_01/greens_0.h5'
with GetFile(greens_data_filename, open_type = 'r') as f:
    pan_d_Nphi              = int(f.read('Nphi'))
    pan_d_Nr                = int(f.read('Nr'))
    pan_phi_list            = f.read('phi_list')
    pan_r_list              = f.read('r_list')
    pan_slab_exists         = f.read('slab_exists')
    pan_e_coarse            = f.read('e_coarse')
    pan_r_frac_top_disk     = f.read('r_frac_top_disk')
    pan_r_frac_bot_disk     = f.read('r_frac_bot_disk')
    pan_refl_profs_top_disk = f.read('refl_profs_top_disk')
    pan_refl_profs_bot_disk = f.read('refl_profs_bot_disk')
    e_fine = f.read('e_fine')

prof = pan_refl_profs_top_disk[0, 20]

Ne = len(e_fine)-1
A  = 10.0**((1.0/Ne)*np.log10(e_fine[-1]/e_fine[0]))
new_eta_grid = np.zeros(2*Ne+1)
for i in range(0, 2*Ne+1):
    new_eta_grid[i] = A**(i-Ne)
bnd_eta_grid = np.zeros(2*Ne+2)
for i in range(1, 2*Ne+1):
    bnd_eta_grid[i] = 0.5*(new_eta_grid[i] + new_eta_grid[i-1])
bnd_eta_grid[0]  = new_eta_grid[0] - 0.5*(new_eta_grid[1] - new_eta_grid[0])
bnd_eta_grid[-1] = new_eta_grid[-1] + 0.5*(new_eta_grid[-1] - new_eta_grid[-2])

cdf = prof[22]

pdf = np.zeros(len(cdf)-1)
for i in range(0, len(pdf)):
    pdf[i] = (cdf[i+1] - cdf[0])/(bnd_eta_grid[i+1] - bnd_eta_grid[i])

plt.plot(new_eta_grid, pdf/np.trapz(pdf, new_eta_grid), 'k-')
plt.xlabel(r'$\nu/\nu^\prime$', fontsize = 22)
plt.ylabel(r'$G\left(\nu^\prime \to \nu\right)$', fontsize = 22)
plt.ylim([0, 0.07])
plt.xlim([1.0e-3, 1.0e3])
plt.text(1.0e1, 0.05, r'$\nu^\prime = 65\ \mathrm{keV}$', fontsize = 22)
plt.semilogx()
plt.tight_layout()
plt.savefig('greens.pdf', bbox_inches = 'tight')
plt.show()

r_frac = pan_r_frac_top_disk[0, 20]
plt.plot(e_fine/1000.0, r_frac, 'k-')
plt.xlim([1.0e-2, 1.0e3])
plt.ylim([0.0, 1.03])
plt.xlabel(r'$\nu^\prime$', fontsize = 22)
plt.ylabel(r'$f\left(\nu^\prime\right)$', fontsize = 22)
plt.tight_layout()
plt.semilogx()
plt.savefig('albedo.pdf', bbox_inches = 'tight')
plt.show()
