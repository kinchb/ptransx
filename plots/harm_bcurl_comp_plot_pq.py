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

num_phi = 64

gdump_fname = '../harm_data/KDHARM0.gdump.a0.h5'

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

# metric = get_metric(gdump_fname)

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

def calc_incor(rho, r, th, phi):
    rho *= (4 * np.pi * c**2 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)

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

    emtop_ik    = (np.pi/2) * np.ones((len(r), len(phi)))
    embot_ik    = (np.pi/2) * np.ones((len(r), len(phi)))
    diskbody_ik = np.zeros((len(r), len(phi)), dtype = int)

    for i in range(0, len(r)-1):
        for k in range(0, len(phi)):
            if ((tau_from_top[i,::,k].max() > 1.0) and (tau_from_bot[i,::,k].max() > 1.0)):
                emtop_ik[i][k] = interp1d(tau_from_top[i,::,k], th_b)(1.0)
                embot_ik[i][k] = interp1d(tau_from_bot[i,::,k][::-1], th_b[::-1])(1.0)
                if (emtop_ik[i][k] > embot_ik[i][k]):
                    emtop_ik[i][k] = np.pi/2
                    embot_ik[i][k] = np.pi/2
                else:
                    diskbody_ik[i][k] = 2
    for k in range(0, len(phi)):
        emtop_ik[len(r)-1][k]    = emtop_ik[len(r)-2][k]
        embot_ik[len(r)-1][k]    = embot_ik[len(r)-2][k]
        diskbody_ik[len(r)-1][k] = 2

    incor = np.ones((len(r), len(th), len(phi)), dtype = int)

    for i in range(0, len(r)):
        for j in range(0, len(th)):
            for k in range(0, len(phi)):
                if (diskbody_ik[i][k] != 0):
                    if ((th[j] > emtop_ik[i][k]) and (th[j] < embot_ik[i][k])):
                        incor[i][j][k] = 0

    return incor

dir = '../harm_data'
a   = 0.0

harm_data_1 = h5py.File(dir + '/KDHARM0.002500.h5', 'r')

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

dV_cgs = np.load('../harm_data/dV_cgs.npz')['dV_cgs']

harm_data_1 = h5py.File(dir + '/KDHARM0.002500.h5', 'r')
harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.012000.h5', 'r')

cool = np.array(harm_data_2['coolfunc'])
rho  = np.array(harm_data_1['rho'])
bsq  = np.array(harm_data_1['bsq'])
B1   = np.array(harm_data_1['B1'])
B2   = np.array(harm_data_1['B2'])
B3   = np.array(harm_data_1['B3'])
dx1  = np.array(harm_data_1['/Header/Grid/dx1'])[0]
dx2  = np.array(harm_data_1['/Header/Grid/dx2'])[0]
dx3  = np.array(harm_data_1['/Header/Grid/dx3'])[0]

cool_new = cool.copy()
rho_new  = rho.copy()

incor_new = calc_incor(rho.copy(), r, th, phi)

bcurl_new = calc_b_curl(r, th, phi, dx1, dx2, dx3, B1, B2, B3)

dir = '../harm_data/old_data'

harm_data_1 = h5py.File(dir + '/KDHARM0.000600.h5', 'r')
harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.012000.h5', 'r')

cool = np.array(harm_data_2['coolfunc'])
rho  = np.array(harm_data_1['rho'])
bsq  = np.array(harm_data_1['bsq'])
B1   = np.array(harm_data_1['B1'])
B2   = np.array(harm_data_1['B2'])
B3   = np.array(harm_data_1['B3'])
dx1  = np.array(harm_data_1['/Header/Grid/dx1'])[0]
dx2  = np.array(harm_data_1['/Header/Grid/dx2'])[0]
dx3  = np.array(harm_data_1['/Header/Grid/dx3'])[0]

cool_old = cool.copy()
rho_old  = rho.copy()

incor_old = calc_incor(rho.copy(), r, th, phi)

bcurl_old = calc_b_curl(r, th, phi, dx1, dx2, dx3, B1, B2, B3)

int_tt = (bcurl_old * incor_old * rho_old * dV_cgs).sum()/(incor_old * rho_old * dV_cgs).sum()
int_ic = (bcurl_new * incor_new * rho_new * dV_cgs).sum()/(incor_new * rho_new * dV_cgs).sum()

print(int_ic/int_tt)

int_tt = (bcurl_old * (1.0 - incor_old) * rho_old * dV_cgs).sum()/((1.0 - incor_old) * rho_old * dV_cgs).sum()
int_ic = (bcurl_new * (1.0 - incor_new) * rho_new * dV_cgs).sum()/((1.0 - incor_new) * rho_new * dV_cgs).sum()

print(int_ic/int_tt)

int_tt = (bcurl_old * incor_old * dV_cgs).sum()
int_ic = (bcurl_new * incor_new * dV_cgs).sum()

print(int_ic/int_tt)

int_tt = (bcurl_old * rho_old * dV_cgs).sum()/(rho_old * dV_cgs).sum()
int_ic = (bcurl_new * rho_new * dV_cgs).sum()/(rho_new * dV_cgs).sum()

print(int_ic/int_tt)

cor_select  = np.zeros((len(r), len(th), len(phi)))

frac = 35.0
for i in range(0, len(r)):
    for j in range(0, len(th)):
        for k in range(0, len(phi)):
            if ((th[j] > (np.pi - (frac/180.0)*np.pi)) and (th[j] < (np.pi + (frac/180.0)*np.pi)) and (r[i] < 25.0)):
                cor_select[i][j][k]  = 1.0
            else:
                cor_select[i][j][k]  = 0.0

int_tt = (bcurl_old * incor_old * cor_select * dV_cgs).sum()
int_ic = (bcurl_new * incor_new * cor_select * dV_cgs).sum()

print(int_ic/int_tt)

int_tt = (bcurl_old * incor_old * cool_new * dV_cgs).sum()/(incor_old * cool_new * dV_cgs).sum()
int_ic = (bcurl_new * incor_new * cool_new * dV_cgs).sum()/(incor_new * cool_new * dV_cgs).sum()

print(int_ic/int_tt)

slice_data_new = bcurl_new
slice_data_old = bcurl_old

slice_data = np.mean((slice_data_new[1:max_r_index,1:-1,1:-1] - slice_data_old[1:max_r_index,1:-1,1:-1])/slice_data_old[1:max_r_index,1:-1,1:-1], axis = 2)

# slice_data = (np.mean(slice_data_new[1:max_r_index,1:-1,1:-1], axis = 2)  - np.mean(slice_data_old[1:max_r_index,1:-1,1:-1], axis = 2))/np.mean(slice_data_old[1:max_r_index,1:-1,1:-1], axis = 2)

r = r[:max_r_index]

th = th[1:-1]

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

plt.pcolormesh(theta, rad, slice_data.T, cmap = cmap, vmin = -1.0, vmax = 1.0)
plt.colorbar(label = r'$\mathrm{fractional\ change\ in\ } | \nabla \times \mathbf{B} |$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()
f.savefig('harm_curl_B_comp.png', bbox_inches = 'tight', dpi = 600)

plt.show()
