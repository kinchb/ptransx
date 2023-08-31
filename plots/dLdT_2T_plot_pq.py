import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

M = 10.0

L_edd = (1.26e38) * M

T_grids = np.load('../harm_data/dLdT_2T_pan_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_data/dLdT_2T_pan_data_a0.0.npz')['cool_pdfs']
times   = np.load('../harm_data/dLdT_2T_pan_data_a0.0.npz')['times']

plt.plot(T_grids[0], pdfs[0]/L_edd, 'k-', label = r'$\mathrm{2T}$')

a = np.trapz(pdfs[0], np.log(T_grids))

T_grids = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['cool_pdfs']
times   = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['times']

plt.plot(T_grids[0], pdfs[0]/L_edd, 'r-', label = r'$\mathrm{1T}$')

b = np.trapz(pdfs[0], np.log(T_grids))

plt.semilogx()

plt.xlim([1.0e3, 1.0e7])
plt.ylim([0, 0.0035])

plt.xlabel(r'$T_e\ \left[\mathrm{eV}\right]$')
plt.ylabel(r'$\partial \left(L_\mathrm{IC}/L_\mathrm{Edd}\right) / \partial \log T_e$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dLdT_2T.pdf', bbox_inches = 'tight')

plt.show()

# --- --- --- --- ---

M = 10.0

L_edd = (1.26e38) * M

T_grids = np.load('../harm_data/dLdT_2T_pan_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_data/dLdT_2T_pan_data_a0.0.npz')['pdfs']
times   = np.load('../harm_data/dLdT_2T_pan_data_a0.0.npz')['times']

plt.plot(T_grids[0], pdfs[0]/np.trapz(pdfs[0][1:], np.log10(T_grids[0][1:])), 'k-', label = r'$\mathrm{2T}$')

T_grids = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['pdfs']
times   = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['times']

plt.plot(T_grids[0], pdfs[0]/np.trapz(pdfs[0][1:], np.log10(T_grids[0][1:])), 'r-', label = r'$\mathrm{1T}$')

plt.semilogx()

plt.xlim([1.0e3, 1.0e6])
plt.ylim([0, 1.75])

plt.xlabel(r'$T_e\ \left[\mathrm{eV}\right]$')
plt.ylabel(r'$\partial L_\mathrm{IC} / \partial \log T_e\ (\mathrm{normalized})$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dLdT_2T_norm.pdf', bbox_inches = 'tight')

plt.show()
