import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

M     = 10.0
L_edd = (1.26e38) * M

def get_data(flnm):
    T_grid    = np.load(flnm)['T_grids'][0]/511.0e3
    rho_pdfs  = np.load(flnm)['rho_pdfs']
    cool_pdfs = np.load(flnm)['cool_pdfs']

    # rho_pdf  = np.mean(rho_pdfs,  axis = 0)
    # cool_pdf = np.mean(cool_pdfs, axis = 0)

    rho_pdfs  = np.nan_to_num(rho_pdfs)
    cool_pdfs = np.nan_to_num(cool_pdfs)

    rho_pdf  = rho_pdfs[-1]#/np.trapz(rho_pdfs[-1], np.log10(T_grid))
    cool_pdf = cool_pdfs[-1]#/np.trapz(cool_pdfs[-1], np.log10(T_grid))

    print(np.trapz(rho_pdf, np.log10(T_grid)))
    print(np.trapz(cool_pdf, np.log10(T_grid)))

    return (T_grid, rho_pdf, cool_pdf)

flnm = '../data/h3d_10_a09_10_gcc/dLdT_radflux_data_a0.9.npz'

T_grid, rho_pdf, cool_pdf = get_data(flnm)

plt.plot(T_grid, rho_pdf, 'k-', label = r'\textsc{harm3d}')

flnm = '../data/h3d_10_a09_10_gcc/dLdT_pan_data_a0.9.npz'

T_grid, rho_pdf, cool_pdf = get_data(flnm)

plt.plot(T_grid, rho_pdf, 'r-', label = r'\textsc{pandurata}')

flnm = '../data/h3d_10_a09_10_gcc/dLdT_pan_data_v2_a0.9.npz'

T_grid, rho_pdf, cool_pdf = get_data(flnm)

plt.plot(T_grid, rho_pdf, 'k--', label = r'\textsc{pandurata}')

plt.semilogx()

plt.legend(loc = 'best', frameon = False)

plt.xlabel(r'$\Theta_e$')
plt.ylabel(r'$\partial \rho_\mathrm{cor} / \partial \log T_e$')

plt.tight_layout()

plt.savefig('dMdT_cutphoto.pdf')

plt.show()

# --- --- ---

flnm = '../data/h3d_10_a09_10_gcc/dLdT_radflux_data_a0.9.npz'

T_grid, rho_pdf, cool_pdf = get_data(flnm)

plt.plot(T_grid, cool_pdf, 'k-', label = r'\textsc{harm3d}')

flnm = '../data/h3d_10_a09_10_gcc/dLdT_pan_data_a0.9.npz'

T_grid, rho_pdf, cool_pdf = get_data(flnm)

plt.plot(T_grid, cool_pdf, 'r-', label = r'\textsc{pandurata}')

plt.semilogx()

plt.legend(loc = 'best', frameon = False)

plt.xlabel(r'$\Theta_e$')
plt.ylabel(r'$\partial L_\mathrm{IC} / \partial \log T_e$')

plt.tight_layout()

plt.savefig('dLdT_cutphoto.pdf')

plt.show()