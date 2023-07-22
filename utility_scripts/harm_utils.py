import os

import numpy as np

import h5py

def gather_files():
    dirs = []
    i = 0
    while True:
        dir = 'run' + '{:03}'.format(i)
        if dir in os.listdir('./'):
            print dir
            dirs.append(dir)
        else:
            break
        i += 1
    dump_list    = []
    radflux_list = []
    curr_dir_index = 0
    file_index     = 0
    while True:
        dump_flnm    = '/KDHARM0.' + '{:06}'.format(500 + file_index) + '.h5'
        radflux_flnm = '/KDHARM0.RADFLUX.' + '{:06}'.format(10000 + file_index) + '.h5'
        for dir_index in range(len(dirs)-1, curr_dir_index-1, -1):
            if (os.path.isfile(dirs[dir_index] + dump_flnm) and os.path.isfile(dirs[dir_index] + radflux_flnm)):
                dump_list.append(dirs[dir_index] + dump_flnm)
                radflux_list.append(dirs[dir_index] + radflux_flnm)
                print radflux_list[-1]
                curr_dir_index = dir_index
                found = True
                break
            else:
                found = False
        file_index += 1
        if (not found):
            break
    return dump_list, radflux_list

def get_metric(gdump_fname):
    num_phi = 64

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

def calc_ucon0(data_file, metric):
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

    return ucon0

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

def calc_dV_cgs(harm_data):
    # useful constants (cgs)
    # --- --- --- --- --- --- --- --- --- ---
    G    = 6.6726e-8
    c    = 3.0e10
    M    = 10.0
    print('WARNING: assuming 10 solar mass black hole!')
    m_bh = M * 2.0e33
    # --- --- --- --- --- --- --- --- --- ---

    if os.path.isfile('./dV_cgs.npz'):
        r      = np.load('dV_cgs.npz')['r']
        th     = np.load('dV_cgs.npz')['th']
        phi    = np.load('dV_cgs.npz')['phi']
        dV_cgs = np.load('dV_cgs.npz')['dV_cgs']
    else:
        harm_data = h5py.File(harm_data, 'r')

        a = harm_data['/Header/Grid/a'][0]

        x1 = np.array(harm_data['x1'])
        x2 = np.array(harm_data['x2'])
        x3 = np.array(harm_data['x3'])

        r   = np.zeros(len(x1))
        th  = np.zeros(len(x2[0]))
        phi = np.zeros(len(x3[0][0]))

        for i in range(0, len(x1)):
            r[i]   = x1[i][0][0]
        for i in range(0, len(x2[0])):
            th[i]  = x2[0][i][0]
        for i in range(0, len(x3[0][0])):
            phi[i] = x3[0][0][i]

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
#                   dV[i][j][k] = np.sqrt(np.abs(g[i][j][1][1] * g[i][j][2][2] * g[i][j][3][3]))*tmp[i][j]*dphi[k]
                    dV[i][j][k] = np.sqrt(-np.linalg.det(g[i][j]))*tmp[i][j]*dphi[k]
                    if np.isnan(dV[i][j][k]):
                        dV[i][j][k] = 0.0

        dV_cgs = dV*(((G*m_bh)/(c*c))**3)

        np.savez('dV_cgs.npz', r = r, th = th, phi = phi, dV_cgs = dV_cgs)

    return (r, th, phi, dV_cgs)
