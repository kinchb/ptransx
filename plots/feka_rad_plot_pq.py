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
    print(a, r)
    return r

def get_k_edge(flnm):
    infile = h5py.File(flnm, 'r')
    spec_top = np.array(infile['data'])
    infile.close()

    flnm = flnm.replace('diskt', 'diskb')
    infile = h5py.File(flnm, 'r')
    spec_bot = np.array(infile['data'])
    infile.close()

    Ne   = len(spec_top)
    Nphi = len(spec_top[0])
    Nr   = len(spec_top[0][0])

    new_spec_top = np.zeros((Nr, Nphi, Ne))
    new_spec_bot = np.zeros((Nr, Nphi, Ne))

    for i in range(0, Nr):
        for j in range(0, Nphi):
            for k in range(0, Ne):
                new_spec_top[i][j][k] = spec_top[k][j][i]
                new_spec_bot[i][j][k] = spec_bot[k][j][i]

    e_grid = np.logspace(0.0, 7.0, Ne, endpoint = True)

    e_min = 7.0e3
    for i in range(0, len(e_grid)):
        if (e_grid[i] > e_min):
            break
    e_min_ndx = i

#   plt.plot(e_grid, new_seed_spec_top[50][0])
#   plt.semilogx()
#   plt.xlabel('energy [eV]')
#   plt.ylabel('disk surface flux')
#   plt.show()

    e_grid = e_grid[i:]
    new_spec_top = new_spec_top[::, ::, i:]
    new_spec_bot = new_spec_bot[::, ::, i:]

    edge = np.zeros(Nr)

    for i in range(0, Nr):
        for j in range(0, Nphi):
            edge[i] += (1.0/(2*Nphi))*(np.trapz(6.2415091e11*new_spec_top[i][j] * (e_grid**-3), e_grid) + np.trapz(6.2415091e11*new_spec_bot[i][j] * (e_grid**-3), e_grid))

    return edge

def get_feka(flnm, do_cont):
    infile = h5py.File(flnm, 'r')

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

    e_min = 3.0e3
    for i in range(0, len(e_grid)):
        if (e_grid[i] > e_min):
            break
    e_min_ndx = i

#   plt.plot(e_grid, new_seed_spec_top[50][0])
#   plt.semilogx()
#   plt.xlabel('energy [eV]')
#   plt.ylabel('disk surface flux')
#   plt.show()

#   e_grid = e_grid[i:]
#   new_seed_spec_top = new_seed_spec_top[::, ::, i:]
#   new_seed_spec_bot = new_seed_spec_bot[::, ::, i:]

    feka = np.zeros(Nr)

    for i in range(0, Nr):
        for j in range(0, Nphi):
            if do_cont:
                feka[i] += (1.0/(2*Nphi))*(np.trapz(6.2415091e11*new_seed_spec_top[i][j], e_grid) + np.trapz(6.2415091e11*new_seed_spec_bot[i][j], e_grid))
            else:
                feka[i] += (1.0/(2*Nphi))*(np.trapz(6.2415091e11*new_seed_spec_top[i][j], e_grid) + np.trapz(6.2415091e11*new_seed_spec_bot[i][j], e_grid))

    return feka

root = '/Users/kinch/panptx/panptx_data/survey/'

mdot = '10'

do_cont = False

# --- --- ---
a = 0.0

dir = root + 'h3d_10_a0_' + mdot + '/'

harm_file = dir + 'RADFLUX_panptx_1000.h5'
line_file = dir + 'spec_out_l_hr.h5'
cont_file = dir + 'scat_diskt.h5'

harm_data = h5py.File(harm_file, 'r')
r = np.array(harm_data['/r'])
harm_data.close()

if do_cont:
    feka = get_k_edge(cont_file)
else:
    feka = get_feka(line_file, do_cont)

trunc_r    = []
trunc_feka = []
j = 0
for i in range(0, len(r), 3):
    if (r[i] > 50.0):
        break
    trunc_r.append(r[i])
    trunc_feka.append(feka[i])
trunc_r    = np.array(trunc_r)
trunc_feka = np.array(trunc_feka)

# plt.plot([isco_radius(a), isco_radius(a)], [1.0e24, 1.0e29], 'r--')

plt.plot(trunc_r/isco_radius(a), trunc_feka, 'r-', label = r'$a = 0$')

# --- --- ---

a = 0.5

dir = root + 'h3d_10_a05_' + mdot + '/'

harm_file = dir + 'RADFLUX_panptx_1000.h5'
line_file = dir + 'spec_out_l_hr.h5'
cont_file = dir + 'scat_diskt.h5'

harm_data = h5py.File(harm_file, 'r')
r = np.array(harm_data['/r'])
harm_data.close()

if do_cont:
    feka = get_k_edge(cont_file)
else:
    feka = get_feka(line_file, do_cont)

trunc_r    = []
trunc_feka = []
j = 0
for i in range(0, len(r), 3):
    if (r[i] > 50.0):
        break
    trunc_r.append(r[i])
    trunc_feka.append(feka[i])
trunc_r    = np.array(trunc_r)
trunc_feka = np.array(trunc_feka)

# plt.plot([isco_radius(a), isco_radius(a)], [1.0e24, 1.0e29], 'g--')

plt.plot(trunc_r/isco_radius(a), trunc_feka, 'g-', label = r'$a = 0.5$')

# --- --- ---

a = 0.9

dir = root + 'h3d_10_a09_' + mdot + '/'

harm_file = dir + 'RADFLUX_panptx_1000.h5'
line_file = dir + 'spec_out_l_hr.h5'
cont_file = dir + 'scat_diskt.h5'

harm_data = h5py.File(harm_file, 'r')
r = np.array(harm_data['/r'])
harm_data.close()

if do_cont:
    feka = get_k_edge(cont_file)
else:
    feka = get_feka(line_file, do_cont)

trunc_r    = []
trunc_feka = []
j = 0
for i in range(0, len(r), 3):
    if (r[i] > 50.0):
        break
    trunc_r.append(r[i])
    trunc_feka.append(feka[i])
trunc_r    = np.array(trunc_r)
trunc_feka = np.array(trunc_feka)

# plt.plot([isco_radius(a), isco_radius(a)], [1.0e24, 1.0e29], 'b--')

# plt.plot(trunc_r, trunc_cont, 'b--',)
plt.plot(trunc_r/isco_radius(a), trunc_feka, 'b-', label = r'$a = 0.9$')

# --- --- ---

plt.loglog()

plt.xlim([0.5, 30])
# plt.ylim([3.0e5, 3.0e8])
plt.ylim([1.0e28, 1.0e31])

if mdot == '10':
#   plt.legend(loc = 'upper right', bbox_to_anchor = (1.12, 1.06), framealpha = 1.0, handletextpad = 0.2)
    plt.legend(loc = 'lower center', frameon = False)

plt.xlabel(r'$r/r_\mathrm{ISCO}$')
if do_cont:
    plt.ylabel(r'$\mathrm{flux\ above\ K\ edge}$')
else:
    plt.ylabel(r'$\mathrm{Fe\ K}\alpha\ \mathrm{flux}\ \left[\mathrm{erg}\ \mathrm{cm}^{-2}\ \mathrm{s}^{-1}\right]$')

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((0.5, 1, 2, 3, 4, 5, 6, 8, 10, 20, 30))
ax.set_xticklabels(('0.5', '1', '2', '3', '4', '5', '', '', '10', '20', '30'))

mass = 10.0
if mdot == '01':
    mdot = 0.01
else:
    mdot = 0.1
mass_label  = r'M/M_\odot &= ' + r'{:g}'.format(mass)    + r'\\[-1ex]'
mdot_label  = r'\dot{m} &= '   + r'{:g}'.format(mdot)    + r'\\'
text        = r'\begin{align*} ' + mass_label + mdot_label + r'\end{align*}'
text = r'$\dot{m} = ' + r'{:g}'.format(mdot) + r'$'
ax.text(0.8, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

f = plt.gcf()
# f.savefig('feka_rad_mdot_0.01.pdf', bbox_inches = 'tight', dpi = 400)

plt.tight_layout()

if (mdot == 0.1):
    plt.savefig('feka_rad_10.pdf', bbox_inches = 'tight')
else:
    plt.savefig('feka_rad_01.pdf', bbox_inches = 'tight')

plt.show()
