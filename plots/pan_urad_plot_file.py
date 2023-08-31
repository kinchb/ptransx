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
a    = 0.9
mdot = 0.1
# --- --- --- --- --- --- --- --- --- ---

# useful constants (cgs)
# --- --- --- --- --- --- --- --- --- ---
m_bh   = M * 2.0e33
G      = 6.6726e-8
c      = 3.0e10
kappa  = 0.4
# --- --- --- --- --- --- --- --- --- ---

def deriv(x):
    dx = np.zeros(len(x))

    for i in range(0, len(x)):
        if (i == 0):
            dx[i] = (-3*x[0] + 4*x[1] - x[2])/2
        elif (i == len(x)-1):
            dx[i] = (3*x[len(x)-1] - 4*x[len(x)-2] + x[len(x)-3])/2
        else:
            dx[i] = (x[i+1] - x[i-1])/2

    return dx

def calc_dV_cgs(radflux_file):
    data = h5py.File(radflux_file, 'r')

    r   = np.array(data['r'])
    th  = np.array(data['th'])
    phi = np.array(data['phi'])

    data.close()

    Nr   = len(r)
    Nth  = len(th)
    Nphi = len(phi)

    rr = np.zeros((Nr, Nth, Nphi))
    tt = np.zeros((Nr, Nth, Nphi))
    pp = np.zeros((Nr, Nth, Nphi))

    for i in range(0, Nr):
        for j in range(0, Nth):
            for k in range(0, Nphi):
                rr[i][j][k] = r[i]
                tt[i][j][k] = th[j]
                pp[i][j][k] = phi[k]

    R_hor = 1 + np.sqrt(1 - a**2)
    Sigma = rr*rr + a*a*np.cos(tt)*np.cos(tt)
    Delta = rr*rr - 2*rr + a*a
    alpha = np.sqrt(Sigma*Delta/(Sigma*Delta + 2*rr*(a*a + rr*rr)))
    omega = 2*rr*a/(Sigma*Delta + 2*rr*(a*a + rr*rr))
    varph = np.sqrt((Sigma*Delta + 2*rr*(a*a + rr*rr))/Sigma*np.sin(tt)*np.sin(tt))

    g = np.zeros((Nr, Nth, 4, 4))
    for i in range(0, Nr):
        for j in range(0, Nth):
            g[i][j][0][0] = -alpha[i][j][0]**2 + (omega[i][j][0]**2)*(varph[i][j][0]**2)
            g[i][j][0][3] = -omega[i][j][0]*(varph[i][j][0]**2)
            g[i][j][1][1] = Sigma[i][j][0]/Delta[i][j][0]
            g[i][j][2][2] = Sigma[i][j][0]
            g[i][j][3][3] = varph[i][j][0]**2
            g[i][j][3][0] = g[i][j][0][3]

    dr   = deriv(r)
    dth  = deriv(th)
    dphi = deriv(phi)
    dV   = np.zeros((Nr, Nth, Nphi))
    tmp  = np.outer(dr, dth)
    for i in range(0, Nr):
        print(i)
        for j in range(0, Nth):
            for k in range(0, Nphi):
                dV[i][j][k] = np.sqrt(np.abs(g[i][j][1][1] * g[i][j][2][2] * g[i][j][3][3]))*tmp[i][j]*dphi[k]
#               dV[i][j][k] = np.sqrt(-np.linalg.det(g[i][j]))*tmp[i][j]*dphi[k]
                if np.isnan(dV[i][j][k]):
                    dV[i][j][k] = 0.0

    dV_cgs = dV*(((G*m_bh)/(c*c))**3)

    return (r, th, phi, dV_cgs)

dir = '../data/h3d_10_a09_10_gcc'

flnm = dir + '/radfluxes/RADFLUX_panptx_1000.h5'

r, th, phi, dV_cgs = calc_dV_cgs(flnm)

panptx_data = h5py.File(flnm, 'r')

harm_urad = np.array(panptx_data['urad'])
rho       = np.array(panptx_data['rho'])

diskbody_ijk = np.array(panptx_data['diskbody_ijk'])

panptx_data.close()

max_r = 51.0
for i in range(0, len(r)):
    if (r[i] > max_r):
        break
max_r_index = i
r = r[:max_r_index]

flnm = dir + '/tes/te.1000.0.h5'

pan_data = h5py.File(flnm, 'r')

T_e = np.array(pan_data['T_e'])

pan_data.close()

flnm = dir + '/scat_cpows/scat_cpow.1000.0.h5'

pan_data = h5py.File(flnm, 'r')

pan_cool = 2.42e17 * np.array(pan_data['data'])

pan_data.close()

incor = 0.0 * pan_cool.copy()

for ii in range(0, len(r)):
    for jj in range(0, len(th)):
        for kk in range(0, len(phi)):
            if ((r[ii] > 1.0 + np.sqrt(1.0 - a**2)) and (diskbody_ijk[ii][jj][kk] == 0)):
                incor[ii][jj][kk] = 1.0

pan_cool *= incor

dV_cgs *= incor

rho *= incor

# --- --- ---

Theta_e = (8.617e-5 * T_e)/511.0e3

Theta_e *= incor

pan_cool = pan_cool/dV_cgs

coef = (4 * 6.652e-25 * 3.0e10 * 1.21)/1.673e-24

pan_urad = pan_cool/(coef * rho * Theta_e * (1 + 4 * Theta_e))

for ii in range(0, len(r)):
    for jj in range(0, len(th)):
        for kk in range(0, len(phi)):
            if incor[ii][jj][kk] == 0.0:
                pan_urad[ii][jj][kk] = 0.0

# --- --- ---

slice_data = pan_urad[:max_r_index,::,::]

slice_data = np.mean(slice_data, axis = 2)

slice_data = np.ma.masked_where(slice_data == 0.0 , slice_data)

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

plt.title(r'\textsc{pandurata}')

plt.pcolormesh(theta, rad, np.log10(slice_data.T/4.0), cmap = cmap, vmin = 9.0, vmax = 13.0)
plt.colorbar(label = r'$\log u_\mathrm{rad} \left[\mathrm{erg}\ \mathrm{cm}^{-3} \right]$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()
f.savefig('pan_urad.png', bbox_inches = 'tight', dpi = 600)

plt.show()
