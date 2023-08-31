import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

radii  = np.load('../harm_data/Te_errors_radial_rho_wght.npz')['radii']
errors = np.load('../harm_data/Te_errors_radial_rho_wght.npz')['errors']

plt.plot([1.0, 100.0], [1.0, 1.0], 'k--')
plt.plot([1.0, 100.0], [2.0, 2.0], 'r--')
plt.plot([1.0, 100.0], [0.5, 0.5], 'r--')

plt.plot(radii, errors, 'k-')

plt.xlim([2, 65])
plt.ylim([0.4, 40.0])

plt.loglog()

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 65))
ax.set_xticklabels(('2', '3', '4', '5', '6', '8', '10', '15', '20', '30', '40', '50', '65'))

plt.xlabel(r'$r/M$')
plt.ylabel(r'$\langle T_e^{\textsc{harm}} / T_e^{\textsc{pandurata}} \rangle$')

# plt.legend(frameon = False, loc = 'upper left')

plt.tight_layout()

f = plt.gcf()
f.savefig('Te_errors_radial.pdf', bbox_inches = 'tight')

plt.show()
