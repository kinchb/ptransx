import sys

import numpy as np

from scipy.interpolate import interp1d

import h5py

time_step = int(sys.argv[2])

# simulation parameters
# --- --- --- --- --- --- --- --- --- ---
M        = 10.0
a        = 0.0
mdot     = 0.01
Fe_abund = 1.0
# --- --- --- --- --- --- --- --- --- ---

input_filename  = sys.argv[1]
output_filename = './RADFLUX_panptx_' + repr(time_step) + '.h5'

if (a == 0.0):
    eta = 0.0572
if (a == 0.5):
    eta = 0.0821
if (a == 0.9):
    eta = 0.1558
if (a == 0.99):
    eta = 0.2640

# useful constants (cgs)
# --- --- --- --- --- --- --- --- --- ---
m_bh   = M * 2.0e33
G      = 6.6726e-8
c      = 3.0e10
kappa  = 0.4
sig_sb = 5.6704e-5
m_e    = 9.1093826e-28
m_i    = 1.67262171e-24
# --- --- --- --- --- --- --- --- --- ---

tau_photo = 1.0  # optical depth of upper/lower disk surface
dtau_tau  = 0.1  # logarithmic step size
dtau_max  = 0.25 # maximum step size
tau_max   = 10.0 # 1/2 maximum thickness of a 'whole' slab

mass_frac_list     = np.array([0.704, 0.28, 6.1E-11, 1.58E-10, 3.78E-09, 0.00278, 0.000814, 0.00756, 0.000000418, 0.00169, 0.0000343, 0.000645, 0.0000556, 0.000696, 0.0000061, 0.000479, 0.00000783, 0.0000701, 0.0000036, 0.0000641, 4.64E-08, 0.0000035, 0.000000356, 0.000017, 0.00000942, 0.00123, 0.00000342, 0.0000729, 0.00000072, 0.00000182])
mass_frac_list[25] = Fe_abund * mass_frac_list[25]
mass_frac_list     = mass_frac_list/mass_frac_list.sum()

amu_list = np.array([1.0079, 4.0026, 6.941, 9.0122, 10.811, 12.0107, 14.0067, 15.9994, 18.9984, 20.1797, 22.9897, 24.305, 26.9815, 28.0855, 30.9738, 32.065, 35.453, 39.948, 39.0983, 40.078, 44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845, 58.9332, 58.6934, 63.546, 65.39])

Z_list = np.linspace(1, 30, num = 30, endpoint = True)

num_frac_list = mass_frac_list/amu_list
num_frac_list = num_frac_list/num_frac_list.sum()

xee_init = (num_frac_list * Z_list).sum()

def get_slab_type(r, th, rho, phi_index, r_index, emtop_ik, embot_ik):
    dens = np.zeros(len(th))
    for i in range(0, len(th)):
        dens[i] = rho[r_index][i][phi_index]

    dens *= (mass_frac_list/amu_list).sum()/1.661e-24

    th_top = emtop_ik[phi_index][r_index]
    th_bot = embot_ik[phi_index][r_index]

    tmp_th   = np.linspace(th_top, th_bot, 1e5, endpoint = True)
    tmp_dens = interp1d(th, dens, kind = 'linear')(tmp_th)

    th   = tmp_th
    dens = tmp_dens

    s = np.zeros(len(th))

    for i in range(1, len(s)):
        s[i] = s[i-1] + ((G*m_bh)/(c*c)) * np.sqrt(r**2 + a**2 * np.cos(0.5*(th[i] + th[i-1]))**2) * (th[i] - th[i-1])

    tau_top    = np.zeros(len(s))
    tau_top[0] = tau_photo

    for i in range(1, len(s)):
        tau_top[i] = tau_top[i-1] + xee_init * 6.6524e-25 * 0.5 * (dens[i] + dens[i-1]) * (s[i] - s[i-1])

    tau_bot     = np.zeros(len(s))
    tau_bot[-1] = tau_photo

    for i in range(len(s)-2, -1, -1):
        tau_bot[i] = tau_bot[i+1] + xee_init * 6.6524e-25 * 0.5 * (dens[i] + dens[i+1]) * (s[i+1] - s[i])

    eq_index = 0
    curr_min = abs(tau_top[0] - tau_bot[0])
    for i in range(0, len(s)):
        if (abs(tau_top[i] - tau_bot[i]) < curr_min):
            curr_min = abs(tau_top[i] - tau_bot[i])
            eq_index = i

    mid_tau = min(tau_top[eq_index], tau_bot[eq_index])
    if (mid_tau < tau_max):
        to_cleave = False
    else:
        to_cleave = True

    tau_steps = []

    tau_steps.append(tau_photo)

    if not to_cleave:
        next_step = tau_steps[len(tau_steps)-1] * dtau_tau
        while ((next_step < dtau_max) and (tau_steps[len(tau_steps)-1] + next_step < mid_tau)):
            tau_steps.append(tau_steps[len(tau_steps)-1] + next_step)
            next_step = tau_steps[len(tau_steps)-1] * dtau_tau
        next_step = dtau_max
        while (tau_steps[len(tau_steps)-1] + next_step < mid_tau):
            tau_steps.append(tau_steps[len(tau_steps)-1] + next_step)
        if (len(tau_steps) < 2):
            print('too few steps')
            return -1
        tau_steps.insert(1, tau_photo + 0.0001)
    else:
        next_step = tau_steps[len(tau_steps)-1] * dtau_tau
        while ((next_step < dtau_max) and (tau_steps[len(tau_steps)-1] + next_step < tau_max)):
            tau_steps.append(tau_steps[len(tau_steps)-1] + next_step)
            next_step = tau_steps[len(tau_steps)-1] * dtau_tau
        next_step = dtau_max
        while (tau_steps[len(tau_steps)-1] + next_step < tau_max):
            tau_steps.append(tau_steps[len(tau_steps)-1] + next_step)
        if (tau_steps[len(tau_steps)-1] < tau_max - 0.0001):
            tau_steps.append(tau_max - 0.0001)
        tau_steps.append(tau_max)
        tau_steps.insert(1, tau_photo + 0.0001)

    tau_steps = np.array(tau_steps)

    tau_top_indices = []
    tau_top_indices.append(0)
    for i in range(1, len(tau_steps)):
        index = tau_top_indices[i-1] + 1
        curr_min = abs(tau_top[index] - tau_steps[i])
        for j in range(index+1, len(tau_top)):
            if (abs(tau_top[j] - tau_steps[i]) < curr_min):
                curr_min = abs(tau_top[j] - tau_steps[i])
                index = j
        tau_top_indices.append(index)

    tau_bot_indices = []
    tau_bot_indices.append(len(tau_bot)-1)
    for i in range(1, len(tau_steps)):
        index = tau_bot_indices[i-1] - 1
        curr_min = abs(tau_bot[index] - tau_steps[i])
        for j in range(index-1, -1, -1):
            if (abs(tau_bot[j] - tau_steps[i]) < curr_min):
                curr_min = abs(tau_bot[j] - tau_steps[i])
                index = j
        tau_bot_indices.append(index)

    tau_bot_indices = tau_bot_indices[::-1]

    if (tau_top_indices[-1] >= tau_bot_indices[0]):
        tau_top_indices = tau_top_indices[:-1]
        tau_bot_indices = tau_bot_indices[1:]

    if ((not to_cleave) and (tau_top_indices[-1] < eq_index) and (tau_bot_indices[0] > eq_index)):
        tau_top_indices.append(eq_index)

    # check that the combined indices are sorted and without duplicates
    if (not np.all(np.array(tau_top_indices + tau_bot_indices) == np.sort(np.array(tau_top_indices + tau_bot_indices)))):
        print('array poorly sorted')
        return -1
    if (len(np.unique(np.array(tau_top_indices + tau_bot_indices))) != len(np.array(tau_top_indices + tau_bot_indices))):
        print('array contains duplicates')
        return -1

    if not to_cleave:
        combined = tau_top_indices + tau_bot_indices

        z = s[combined]

        zc = np.zeros(len(z)-1)

        for i in range(0, len(zc)):
            zc[i] = 0.5 * (z[i] + z[i+1])

        dens = np.zeros(len(zc))

        for i in range(0, len(dens)):
            dens[i] = (tau_top[combined[i+1]] - tau_top[combined[i]])/(xee_init * 6.6524e-25 * (z[i+1] - z[i]))

        dens[0]  = dens[1]
        dens[-1] = dens[-2]

        return 'whole'
    else:
        z_u = s[tau_top_indices]
        z_l = s[tau_bot_indices]

        zc_u = np.zeros(len(z_u)-1)
        zc_l = np.zeros(len(z_l)-1)

        for i in range(0, len(zc_u)):
            zc_u[i] = 0.5 * (z_u[i] + z_u[i+1])
        for i in range(0, len(zc_l)):
            zc_l[i] = 0.5 * (z_l[i] + z_l[i+1])

        dens_u = np.zeros(len(zc_u))
        dens_l = np.zeros(len(zc_l))

        for i in range(0, len(dens_u)):
            dens_u[i] = (tau_top[tau_top_indices[i+1]] - tau_top[tau_top_indices[i]])/(xee_init * 6.6524e-25 * (z_u[i+1] - z_u[i]))
        for i in range(0, len(dens_l)):
            dens_l[i] = (tau_bot[tau_bot_indices[i]] - tau_bot[tau_bot_indices[i+1]])/(xee_init * 6.6524e-25 * (z_l[i+1] - z_l[i]))

        dens_u[0]  = dens_u[1]
        dens_u[-1] = dens_u[-2]
        dens_l[0]  = dens_l[1]
        dens_l[-1] = dens_l[-2]

        return 'cleave'

# read in relevant data
harm_file = h5py.File(input_filename, 'r')

r   = np.array(harm_file['r'])
th  = np.array(harm_file['th'])
phi = np.array(harm_file['phi'])

rho      = np.array(harm_file['rho_'      + repr(time_step)])
uu       = np.array(harm_file['uu_'       + repr(time_step)])
urad     = np.array(harm_file['urad_'     + repr(time_step)])
T_C      = np.array(harm_file['T_C_'      + repr(time_step)])
coolfunc = np.array(harm_file['coolfunc_' + repr(time_step)])

T_e = 511.0e3 * (1.0/(1.0 + 1.21)) * ((5.0/3.0) - 1) * (m_i/m_e) * (uu/rho)

ucon0 = np.array(harm_file['ucon0_' + repr(time_step)])
ucon1 = np.array(harm_file['ucon1_' + repr(time_step)])
ucon2 = np.array(harm_file['ucon2_' + repr(time_step)])
ucon3 = np.array(harm_file['ucon3_' + repr(time_step)])

harm_file.close()

# convert rho and coolfunc to cgs units
rho      *= (4 * np.pi * c**2 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)
urad     *= (4 * np.pi * c**4 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)
coolfunc *= (4 * np.pi * c**7 * mdot)/(kappa * G**2 * m_bh**2 * 3.0e-4 * eta)

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

diskbody_ijk = np.zeros((len(r), len(th), len(phi)), dtype = int)

for i in range(0, len(r)):
    for j in range(0, len(th)):
        for k in range(0, len(phi)):
            if (diskbody_ik[i][k] != 0):
                if ((th[j] > emtop_ik[i][k]) and (th[j] < embot_ik[i][k])):
                    diskbody_ijk[i][j][k] = 1

# determine disk effective temperature
Tdisk_ik = np.zeros((len(r), len(phi)))
int_flux = np.zeros((len(r), len(phi)))

for i in range(0, len(r)):
    for j in range(1, len(th_b)):
        for k in range(0, len(phi)):
            int_flux[i][k] += diskbody_ijk[i][j-1][k] * coolfunc[i][j-1][k] * ((G*m_bh)/(c*c)) * np.sqrt(r[i]**2 + a**2 * np.cos(th[j-1])**2) * (th_b[j] - th_b[j-1])
Tdisk_ik = (0.5*int_flux/sig_sb)**0.25

slab_type_ik = np.zeros((len(phi), len(r)), dtype = int)
k = 0
for i in range(0, len(phi), 8):
    for j in range(0, len(r), 3):
        print(k)
        k += 1
        if ((diskbody_ik[j][i] == 0) or (r[j] < 1.0 + np.sqrt(1.0 - a**2)) or (r[j] > 50.0)):
            continue
        slab_type = get_slab_type(r[j], th, rho, i, j, emtop_ik.T, embot_ik.T)
        if (slab_type == -1):
            print('should not get here')
        if (slab_type == 'whole'):
            slab_type_ik[i][j] = int(1)
        if (slab_type == 'cleave'):
            slab_type_ik[i][j] = int(2)

# output to file
f = h5py.File(output_filename, 'w')

f['/M']        = np.array([M])
f['/a']        = np.array([a])
f['/mdot']     = np.array([mdot])
f['/eta']      = np.array([eta])
f['/Fe_abund'] = np.array([Fe_abund])

f['/r']            = r
f['/th']           = th
f['/phi']          = phi
f['/rho']          = rho
f['/urad']         = urad
f['/T_C']          = T_C
f['/T_e']          = T_e
f['/coolfunc']     = coolfunc
f['/ucon0']        = ucon0
f['/ucon1']        = ucon1
f['/ucon2']        = ucon2
f['/ucon3']        = ucon3
f['/diskbody_ik']  = diskbody_ik.T
f['/diskbody_ijk'] = diskbody_ijk
f['/Tdisk_ik']     = Tdisk_ik.T
f['/emtop_ik']     = emtop_ik.T
f['/embot_ik']     = embot_ik.T
f['/slab_type_ik'] = slab_type_ik

f.close()
