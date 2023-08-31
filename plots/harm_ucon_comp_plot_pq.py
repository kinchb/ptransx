import numpy as np

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
a    = 0.0
mdot = 0.01
# --- --- --- --- --- --- --- --- --- ---

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

gdump_fname = '../harm_data/KDHARM0.gdump.a0.h5'

num_phi = 64

grid_data = h5py.File('../harm_data/KDHARM0.000500.h5', 'r')

dx_dxp10 = np.array(grid_data['/dx_dxp10'])
dx_dxp11 = np.array(grid_data['/dx_dxp11'])
dx_dxp12 = np.array(grid_data['/dx_dxp12'])
dx_dxp13 = np.array(grid_data['/dx_dxp13'])
dx_dxp20 = np.array(grid_data['/dx_dxp20'])
dx_dxp21 = np.array(grid_data['/dx_dxp21'])
dx_dxp22 = np.array(grid_data['/dx_dxp22'])
dx_dxp23 = np.array(grid_data['/dx_dxp23'])
dx_dxp30 = np.array(grid_data['/dx_dxp30'])
dx_dxp31 = np.array(grid_data['/dx_dxp31'])
dx_dxp32 = np.array(grid_data['/dx_dxp32'])
dx_dxp33 = np.array(grid_data['/dx_dxp33'])

x1 = np.array(grid_data['/x1'])
x2 = np.array(grid_data['/x2'])
x3 = np.array(grid_data['/x3'])

grid_data.close()

r   = np.zeros(len(x1))
th  = np.zeros(len(x2[0]))
phi = np.zeros(len(x3[0][0]))

rr = x1

for i in range(0, len(x1)):
    r[i]   = x1[i][0][0]
for i in range(0, len(x2[0])):
    th[i]  = x2[0][i][0]
for i in range(0, len(x3[0][0])):
    phi[i] = x3[0][0][i]

def get_metric(gdump_fname):
    gdump_file = h5py.File(gdump_fname, 'r')

    gcov00 = np.dstack([np.squeeze(np.array(gdump_file['/gcov300']))] * num_phi)
    gcov01 = np.dstack([np.squeeze(np.array(gdump_file['/gcov301']))] * num_phi)
    gcov02 = np.dstack([np.squeeze(np.array(gdump_file['/gcov302']))] * num_phi)
    gcov03 = np.dstack([np.squeeze(np.array(gdump_file['/gcov303']))] * num_phi)
    gcov11 = np.dstack([np.squeeze(np.array(gdump_file['/gcov311']))] * num_phi)
    gcov12 = np.dstack([np.squeeze(np.array(gdump_file['/gcov312']))] * num_phi)
    gcov13 = np.dstack([np.squeeze(np.array(gdump_file['/gcov313']))] * num_phi)
    gcov22 = np.dstack([np.squeeze(np.array(gdump_file['/gcov322']))] * num_phi)
    gcov23 = np.dstack([np.squeeze(np.array(gdump_file['/gcov323']))] * num_phi)
    gcov33 = np.dstack([np.squeeze(np.array(gdump_file['/gcov333']))] * num_phi)
    gdet   = np.dstack([np.squeeze(np.array(gdump_file['/gdet3']))]   * num_phi)

    gdump_file.close()

    gcon00 =  gcov11*(gcov22*gcov33 - gcov23*gcov23) - gcov12*(gcov12*gcov33 - gcov13*gcov23) + gcov13*(gcov12*gcov23 - gcov13*gcov22)
    gcon01 = -gcov01*(gcov22*gcov33 - gcov23*gcov23) + gcov02*(gcov12*gcov33 - gcov13*gcov23) - gcov03*(gcov12*gcov23 - gcov13*gcov22)
    gcon02 =  gcov01*(gcov12*gcov33 - gcov23*gcov13) - gcov02*(gcov11*gcov33 - gcov13*gcov13) + gcov03*(gcov11*gcov23 - gcov13*gcov12)
    gcon03 = -gcov01*(gcov12*gcov23 - gcov22*gcov13) + gcov02*(gcov11*gcov23 - gcov12*gcov13) - gcov03*(gcov11*gcov22 - gcov12*gcov12)
    gcon11 =  gcov00*(gcov22*gcov33 - gcov23*gcov23) - gcov02*(gcov02*gcov33 - gcov03*gcov23) + gcov03*(gcov02*gcov23 - gcov03*gcov22)
    gcon12 = -gcov00*(gcov12*gcov33 - gcov23*gcov13) + gcov02*(gcov01*gcov33 - gcov03*gcov13) - gcov03*(gcov01*gcov23 - gcov03*gcov12)
    gcon13 =  gcov00*(gcov12*gcov23 - gcov22*gcov13) - gcov02*(gcov01*gcov23 - gcov02*gcov13) + gcov03*(gcov01*gcov22 - gcov02*gcov12)
    gcon22 =  gcov00*(gcov11*gcov33 - gcov13*gcov13) - gcov01*(gcov01*gcov33 - gcov03*gcov13) + gcov03*(gcov01*gcov13 - gcov03*gcov11)
    gcon23 = -gcov00*(gcov11*gcov23 - gcov12*gcov13) + gcov01*(gcov01*gcov23 - gcov02*gcov13) - gcov03*(gcov01*gcov12 - gcov02*gcov11)
    gcon33 =  gcov00*(gcov11*gcov22 - gcov12*gcov12) - gcov01*(gcov01*gcov22 - gcov02*gcov12) + gcov02*(gcov01*gcov12 - gcov02*gcov11)

    det = gcov00*gcon00 + gcov01*gcon01 + gcov02*gcon02 + gcov03*gcon03

    det[det==0]       +=  1.e-10
    gcon00[gcon00==0] += -1.e-10

    inv_det = 1.0/det

    gcon00 *= inv_det
    gcon01 *= inv_det
    gcon02 *= inv_det
    gcon03 *= inv_det
    gcon11 *= inv_det
    gcon12 *= inv_det
    gcon13 *= inv_det
    gcon22 *= inv_det
    gcon23 *= inv_det
    gcon33 *= inv_det

    alpha = np.sqrt(-1.0/gcon00)
    beta1 = -gcon01/gcon00
    beta2 = -gcon02/gcon00
    beta3 = -gcon03/gcon00

    metric = {'gdet':gdet,
              'gcov00':gcov00, 'gcov01':gcov01, 'gcov02':gcov02, 'gcov03':gcov03,
              'gcov11':gcov11, 'gcov12':gcov12, 'gcov13':gcov13,
              'gcov22':gcov22, 'gcov23':gcov23,
              'gcov33':gcov33,
              'gcon00':gcon00, 'gcon01':gcon01, 'gcon02':gcon02, 'gcon03':gcon03,
              'gcon11':gcon11, 'gcon12':gcon12, 'gcon13':gcon13,
              'gcon22':gcon22, 'gcon23':gcon23,
              'gcon33':gcon33,
              'alpha':alpha, 'beta1':beta1, 'beta2':beta2, 'beta3':beta3}

    return metric

def uconp2ucon_dynamic(data_file, ucon0, ucon1, ucon2, ucon3):
    ucontt = ucon0
    uconrr = dx_dxp10 * ucon0 + dx_dxp11 * ucon1 + dx_dxp12 * ucon2 + dx_dxp13 * ucon3
    uconth = dx_dxp20 * ucon0 + dx_dxp21 * ucon1 + dx_dxp22 * ucon2 + dx_dxp23 * ucon3
    uconph = dx_dxp30 * ucon0 + dx_dxp31 * ucon1 + dx_dxp32 * ucon2 + dx_dxp33 * ucon3

    return ucontt, uconrr, uconth, uconph

def ks2bl_con2(data_file, uconttKS, uconrrKS, uconthKS, uconphKS):
    spin = a

    delta = rr*rr - 2.0*rr + spin*spin

    dtBL_drKS   = -2.0*rr/delta
    dphiBL_drKS = -spin/delta

    ucon0BL = uconttKS + dtBL_drKS*uconrrKS
    ucon1BL = uconrrKS
    ucon2BL = uconthKS
    ucon3BL = uconphKS + dphiBL_drKS*uconrrKS

    return ucon0BL, ucon1BL, ucon2BL, ucon3BL

def calc_ucon_bl(data_file, metric):
    v1 = np.array(data_file['/v1'])
    v2 = np.array(data_file['/v2'])
    v3 = np.array(data_file['/v3'])

    vsq = metric['gcov11']*v1*v1 + metric['gcov22']*v2*v2 + metric['gcov33']*v3*v3 + 2.0*(v1*(metric['gcov12']*v2 + metric['gcov13']*v3) + metric['gcov23']*v2*v3)

    gamma = np.sqrt(1.0 + vsq)
    alpha = metric['alpha']
    beta1 = metric['beta1']
    beta2 = metric['beta2']
    beta3 = metric['beta3']

    ucon0 = gamma/alpha
    ucon1 = v1 - ucon0 * beta1
    ucon2 = v2 - ucon0 * beta2
    ucon3 = v3 - ucon0 * beta3

    ucon0_ks, ucon1_ks, ucon2_ks, ucon3_ks = uconp2ucon_dynamic(data_file, ucon0, ucon1, ucon2, ucon3)

    ucon0_bl, ucon1_bl, ucon2_bl, ucon3_bl = ks2bl_con2(data_file, ucon0_ks, ucon1_ks, ucon2_ks, ucon3_ks)

    return vsq, ucon1_bl, ucon2_bl, ucon3_bl

metric = get_metric(gdump_fname)

dir = '../harm_data/'

harm_data_1 = h5py.File(dir + '/KDHARM0.001500.h5', 'r')

x1 = np.array(harm_data_1['x1'])
x2 = np.array(harm_data_1['x2'])
x3 = np.array(harm_data_1['x3'])

harm_data_1.close()

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
r = r[:max_r_index]

harm_data = h5py.File(dir + '/KDHARM0.RADFLUX.011000.h5', 'r')

u0_ic, u1_ic, u2_ic, u3_ic = calc_ucon_bl(harm_data, metric)

harm_data.close()

dir = '../harm_data/old_data'

harm_data = h5py.File(dir + '/KDHARM0.RADFLUX.011000.h5', 'r')

u0_tt, u1_tt, u2_tt, u3_tt = calc_ucon_bl(harm_data, metric)

harm_data.close()

slice_data = u0_ic[:max_r_index,::,32]/u0_tt[:max_r_index,::,32]

# slice_data = np.mean(slice_data, axis = 2)

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

cmap = plt.get_cmap('RdBu')

plt.pcolormesh(theta, rad, slice_data.T, cmap = cmap, vmin = 0.9, vmax = 1.1)
plt.colorbar(label = r'$[(\mathcal{L}/\rho)_\mathrm{IC} - (\mathcal{L}/\rho)_\mathrm{TT}]/(\mathcal{L}/\rho)_\mathrm{TT}$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

# plt.title('IC in corona')

# f = plt.gcf()
# f.savefig('harm_cool_o_rho_comp_ic.png', bbox_inches = 'tight', dpi = 600)

plt.show()
