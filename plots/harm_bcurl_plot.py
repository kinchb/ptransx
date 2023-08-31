import numpy as np

from scipy.interpolate import interp1d

import h5py

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc

import imageio

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 15)

# simulation parameters
# --- --- --- --- --- --- --- --- --- ---
M    = 10.0
a    = 0.99
mdot = 0.01
# --- --- --- --- --- --- --- --- --- ---

# useful constants (cgs)
# --- --- --- --- --- --- --- --- --- ---
m_bh   = M * 2.0e33
G      = 6.6726e-8
c      = 3.0e10
kappa  = 0.4
# --- --- --- --- --- --- --- --- --- ---

m_e = 9.1093826e-28
m_i = 1.67262171e-24

if (a == 0.0):
    eta = 0.0572
if (a == 0.5):
    eta = 0.0821
if (a == 0.9):
    eta = 0.1558
if (a == 0.99):
    eta = 0.2640

def calc_b_curl(r, th, phi, dx1, dx2, dx3, B1, B2, B3):
    A = np.zeros((len(r), len(th), len(phi), 4))

    for i in range(0, len(r)):
        for j in range(0, len(th)):
            for k in range(0, len(phi)):
                A[i][j][k][1] = B1[i][j][k]
                A[i][j][k][2] = B2[i][j][k]
                A[i][j][k][3] = B3[i][j][k]

    f1 = 0.25/dx1
    f2 = 0.25/dx2
    f3 = 0.25/dx3

    curl_B = np.zeros((len(r), len(th), len(phi)))

    B1 = 1
    B2 = 2
    B3 = 3

    for i in range(1, len(r)-1):
        print(i)
        for j in range(1, len(th)-1):
            for k in range(1, len(phi)-1):
                b1 = (
                f2*(
                    A[i  ][j+1][k  ][B3] - A[i  ][j  ][k  ][B3] +
                    A[i+1][j+1][k  ][B3] - A[i+1][j  ][k  ][B3] +
                    A[i  ][j+1][k+1][B3] - A[i  ][j  ][k+1][B3] +
                    A[i+1][j+1][k+1][B3] - A[i+1][j  ][k+1][B3]
                    ) -
                f3*(
                    A[i  ][j  ][k+1][B2] - A[i  ][j  ][k  ][B2] +
                    A[i+1][j  ][k+1][B2] - A[i+1][j  ][k  ][B2] +
                    A[i  ][j+1][k+1][B2] - A[i  ][j+1][k  ][B2] +
                    A[i+1][j+1][k+1][B2] - A[i+1][j+1][k  ][B2]
                    )
                )
                b2 = (
                f3*(
                    A[i  ][j  ][k+1][B1] - A[i  ][j  ][k  ][B1] +
                    A[i+1][j  ][k+1][B1] - A[i+1][j  ][k  ][B1] +
                    A[i  ][j+1][k+1][B1] - A[i  ][j+1][k  ][B1] +
                    A[i+1][j+1][k+1][B1] - A[i+1][j+1][k  ][B1]
                    ) -
                f1*(
                    A[i+1][j  ][k  ][B3] - A[i  ][j  ][k  ][B3] +
                    A[i+1][j+1][k  ][B3] - A[i  ][j+1][k  ][B3] +
                    A[i+1][j  ][k+1][B3] - A[i  ][j  ][k+1][B3] +
                    A[i+1][j+1][k+1][B3] - A[i  ][j+1][k+1][B3]
                    )
                )
                b3 =  (
        	    f1*(
                    A[i+1][j  ][k  ][B2] - A[i  ][j  ][k  ][B2] +
                    A[i+1][j  ][k+1][B2] - A[i  ][j  ][k+1][B2] +
                    A[i+1][j+1][k  ][B2] - A[i  ][j+1][k  ][B2] +
                    A[i+1][j+1][k+1][B2] - A[i  ][j+1][k+1][B2]
                ) -
                f2*(
                    A[i  ][j+1][k  ][B1] - A[i  ][j  ][k  ][B1] +
                    A[i+1][j+1][k  ][B1] - A[i+1][j  ][k  ][B1] +
                    A[i  ][j+1][k+1][B1] - A[i  ][j  ][k+1][B1] +
                    A[i+1][j+1][k+1][B1] - A[i+1][j  ][k+1][B1]
                    )
                )
                curl_B[i][j][k] = np.sqrt(b1*b1 + b2*b2 + b3*b3)

    return curl_B

def pad(r, th, data):
    new_r    = np.zeros(len(r)+1)
    new_data = np.zeros((len(r)+1, len(th)))

    for i in range(1, len(r)+1):
        new_r[i] = r[i-1]
        for j in range(0, len(th)):
            new_data[i][j] = data[i-1][j]

    return (new_r, new_data)

dir = '../data/h3d_10_a099_01_gcc_later/dumps'

harm_data_1 = h5py.File(dir + '/KDHARM0.000520.h5', 'r')
harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.010020.h5', 'r')

cool = np.array(harm_data_2['coolfunc'])
rho  = np.array(harm_data_1['rho'])
bsq  = np.array(harm_data_1['bsq'])
B1   = np.array(harm_data_1['B1'])
B2   = np.array(harm_data_1['B2'])
B3   = np.array(harm_data_1['B3'])
dx1  = np.array(harm_data_1['/Header/Grid/dx1'])[0]
dx2  = np.array(harm_data_1['/Header/Grid/dx2'])[0]
dx3  = np.array(harm_data_1['/Header/Grid/dx3'])[0]

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

max_r = 51.0
for i in range(0, len(r)):
    if (r[i] > max_r):
        break
max_r_index = i

bcurl = calc_b_curl(r, th, phi, dx1, dx2, dx3, B1, B2, B3)

slice_data = np.mean(bcurl, axis = 2)

r = r[:max_r_index]

th = th[1:-1]

r, slice_data = pad(r, th, slice_data)

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

plt.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap)#, vmin = -1.0, vmax = 1.0)
plt.colorbar(label = r'$| \nabla \times \mathbf{B} |$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()
f.savefig('harm_curl_B.png', bbox_inches = 'tight', dpi = 600)

plt.show()
