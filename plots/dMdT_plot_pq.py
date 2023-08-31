import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

T_grid = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['T_grids'][0]
pdf    = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['harm_rho_pdfs'][0]

pdf[0] = 0.0

T_harm_1000   = T_grid.copy()
pdf_harm_1000 = pdf.copy()

plt.plot(T_grid, pdf, 'k-', label = '$1000M$')

T_grid = np.load('../harm_data/dLdT_old_pan_data_a0.0.npz')['T_grids'][0]
pdf    = np.load('../harm_data/dLdT_old_pan_data_a0.0.npz')['harm_rho_pdfs'][0]

pdf[0] = 0.0

T_harm_m1   = T_grid.copy()
pdf_harm_m1 = pdf.copy()

plt.plot(T_grid, pdf, 'r-', label = '$-1M$')

plt.semilogx()

plt.xlim([1.0e3, 1.0e6])
plt.ylim([0, 5.0e16])

plt.xlabel(r'$T_e\ \left[\mathrm{eV}\right]$')
plt.ylabel(r'$\partial \rho / \partial \log T_e$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dMdT.pdf', bbox_inches = 'tight')

plt.show()

# --- --- ---

T_grids = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['rho_pdfs']
times   = np.load('../harm_data/dLdT_pan_data_a0.0.npz')['times']

T_pan_1000   = T_grids[0].copy()
pdf_pan_1000 = pdfs[0].copy()
pdf_pan_1000[0] = 0.0

T_grids = np.load('../harm_data/dLdT_old_pan_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_data/dLdT_old_pan_data_a0.0.npz')['rho_pdfs']
times   = np.load('../harm_data/dLdT_old_pan_data_a0.0.npz')['times']

T_pan_m1   = T_grids[0].copy()
pdf_pan_m1 = pdfs[0].copy()
pdf_pan_m1[0] = 0.0

"""
pdf_harm_1000 = pdf_harm_1000/pdf_harm_1000.max() #np.trapz(T_harm_1000 * pdf_harm_1000, T_harm_1000)
pdf_harm_m1   = pdf_harm_m1/pdf_harm_m1.max() #np.trapz(T_harm_m1 * pdf_harm_m1, T_harm_m1)
pdf_pan_1000  = pdf_pan_1000/pdf_pan_1000.max() #np.trapz(T_pan_1000 * pdf_pan_1000, T_pan_1000)
pdf_pan_m1    = pdf_pan_m1/pdf_pan_m1.max() #np.trapz(T_pan_m1 * pdf_pan_m1, T_pan_m1)
"""

pdf_harm_1000 = np.nan_to_num(pdf_harm_1000)
pdf_harm_m1   = np.nan_to_num(pdf_harm_m1)
pdf_pan_1000  = np.nan_to_num(pdf_pan_1000)
pdf_pan_m1    = np.nan_to_num(pdf_pan_m1)

pdf_harm_1000 = pdf_harm_1000/np.trapz(pdf_harm_1000, np.log10(T_harm_1000))
pdf_harm_m1   = pdf_harm_m1/np.trapz(pdf_harm_m1, np.log10(T_harm_m1))
pdf_pan_1000  = pdf_pan_1000/np.trapz(pdf_pan_1000, np.log10(T_pan_1000))
pdf_pan_m1    = pdf_pan_m1/np.trapz(pdf_pan_m1, np.log10(T_pan_m1))

print(np.trapz(pdf_harm_1000, np.log10(T_harm_1000)))
print(np.trapz(pdf_pan_1000, np.log10(T_pan_1000)))

plt.plot(T_pan_1000,  pdf_pan_1000,  'k-',  label = r'$1000M\ \textsc{pandurata}$')
plt.plot(T_harm_1000, pdf_harm_1000, 'k--', label = r'$1000M\ \textsc{harm3d}$')
plt.plot(T_pan_m1,    pdf_pan_m1,    'r-',  label = r'$-1M\ \textsc{pandurata}$')
plt.plot(T_harm_m1,   pdf_harm_m1,   'r--', label = r'$-1M\ \textsc{harm3d}$')

plt.semilogx()

plt.xlim([1.0e3, 1.0e6])
plt.ylim([0.0, 1.75])

plt.xlabel(r'$T_e\ \left[\mathrm{eV}\right]$')
plt.ylabel(r'$\partial \rho / \partial \log T_e$')

ax = plt.gca()

ax.legend(bbox_to_anchor = (0.65, 0.6), frameon = True, fontsize = 14, framealpha = 1.0)

plt.tight_layout()

f = plt.gcf()
f.savefig('dMdT_comp.pdf', bbox_inches = 'tight')

plt.show()
