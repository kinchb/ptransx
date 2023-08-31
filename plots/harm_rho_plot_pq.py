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

from scipy.interpolate import interp1d

# simulation parameters
# --- --- --- --- --- --- --- --- --- ---
M    = 10.0
a    = 0.0
mdot = 0.01
# --- --- --- --- --- --- --- --- --- ---

dir = '../panptx_data/survey/h3d_10_a0_01'

# useful constants (cgs)
# --- --- --- --- --- --- --- --- --- ---
m_bh   = M * 2.0e33
G      = 6.6726e-8
c      = 3.0e10
kappa  = 0.4
# --- --- --- --- --- --- --- --- --- ---

if (a == 0.0):
    eta = 0.0572
if (a == 0.5):
    eta = 0.0821
if (a == 0.9):
    eta = 0.1558
if (a == 0.99):
    eta = 0.2640

harm_data_1 = h5py.File(dir + '/KDHARM0.001500.h5', 'r')

x1 = np.array(harm_data_1['x1'])
x2 = np.array(harm_data_1['x2'])
x3 = np.array(harm_data_1['x3'])

r   = np.zeros(len(x1))
th  = np.zeros(len(x2[0]))
phi = np.zeros(len(x3[0][0]))

for i in range(0, len(x1)):
    r[i]   = x1[i][0][0]
for i in range(0, len(x2[0])):
    th[i]  = x2[0][i][0]
for i in range(0, len(x3[0][0])):
    phi[i] = x3[0][0][i]

max_r = 20.0
for i in range(0, len(r)):
    if (r[i] > max_r):
        break
max_r_index = i
r = r[:max_r_index]

def get_mean_photo_thetas(rho):
    # convert rho and coolfunc to cgs units
    rho *= (4 * np.pi * c**2 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)

    # determine tau = 1 surface location and related quantities
    th_b = np.zeros(len(th)+1)
    for i in range(1, len(th_b)-1):
        th_b[i] = 0.5*(th[i-1] + th[i])
    th_b[0]  = 0.0
    th_b[-1] = np.pi

    tau_from_top = np.zeros((len(r), len(th_b), len(phi)))
    tau_from_bot = np.zeros((len(r), len(th_b), len(phi)))

    for i in range(0, len(r)):
        for j in range(1, len(th_b)):
            for k in range(0, len(phi)):
                tau_from_top[i][j][k] = tau_from_top[i][j-1][k] + kappa * rho[i][j-1][k] * ((G*m_bh)/(c*c)) * np.sqrt(r[i]**2 + a**2 * np.cos(th[j-1])**2) * (th_b[j] - th_b[j-1])

    for i in range(0, len(r)):
        for j in range(len(th_b)-2, -1, -1):
            for k in range(0, len(phi)):
                tau_from_bot[i][j][k] = tau_from_bot[i][j+1][k] + kappa * rho[i][j][k] * ((G*m_bh)/(c*c)) * np.sqrt(r[i]**2 + a**2 * np.cos(th[j])**2) * (th_b[j+1] - th[j])

    emtop_ik = (np.pi/2) * np.ones((len(r), len(phi)))
    embot_ik = (np.pi/2) * np.ones((len(r), len(phi)))

    for i in range(0, len(r)-1):
        for k in range(0, len(phi)):
            if ((tau_from_top[i,::,k].max() > 1.0) and (tau_from_bot[i,::,k].max() > 1.0)):
                emtop_ik[i][k] = interp1d(tau_from_top[i,::,k], th_b)(1.0)
                embot_ik[i][k] = interp1d(tau_from_bot[i,::,k][::-1], th_b[::-1])(1.0)
                if (emtop_ik[i][k] > embot_ik[i][k]):
                    emtop_ik[i][k] = np.pi/2
                    embot_ik[i][k] = np.pi/2
    for k in range(0, len(phi)):
        emtop_ik[len(r)-1][k] = emtop_ik[len(r)-2][k]
        embot_ik[len(r)-1][k] = embot_ik[len(r)-2][k]

    emtop = np.mean(emtop_ik, axis = 1)
    embot = np.mean(embot_ik, axis = 1)

    r_abbr = r.copy()

    start_ndx = 0
    while ((emtop[start_ndx] == np.pi/2) or (embot[start_ndx] == np.pi/2)):
        start_ndx += 1

    return (r_abbr[start_ndx:], emtop[start_ndx:], embot[start_ndx:])

rho = np.array(harm_data_1['rho'])

rho = rho[:max_r_index,::,::]

r_abbr, emtop, embot = get_mean_photo_thetas(rho.copy())

rho = np.mean(rho, axis = 2)
# rho = rho[:max_r_index,::,32]

slice_data = rho

rad, theta = np.meshgrid(r, th)

fig = plt.figure(figsize = (4.5, 8))
ax  = Axes3D(fig)

rad, theta = np.meshgrid(r, th)

ax = plt.subplot(projection = 'polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_rlabel_position(80)
ax.set_xticks(np.pi/180.0 * np.array([0, 30, 60, 90, 120, 150, 180]))
ax.set_yticks(np.array([0, 10, 20, 30, 40, 50]))

cmap = plt.get_cmap('jet')

rho_conv = (4 * np.pi * c**2 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)

slice_data *= rho_conv

plt.plot(emtop, r_abbr, 'w-')
plt.plot(embot, r_abbr, 'w-')

np.savez(dir + '/photo_locs.npz', r_abbr = r_abbr, emtop = emtop, embot = embot)

plt.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = -10.0, vmax = -5.0)
plt.colorbar(label = r'$\log \rho \left[\mathrm{g}\ \mathrm{cm}^{-3} \right]$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()
f.savefig('harm_rho_' + dir.split('/')[-1] + '.png', bbox_inches = 'tight', dpi = 600)

plt.show()
