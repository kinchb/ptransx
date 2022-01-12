# command line execution:
# mpirun -np <numper of processors> python nr_iterate_corona.py <'panptx' input file> <run ID> <grand iteration> <disk spectrum file> <'normal' or 'hires'> >pan.out 2>pan.err

import sys

if (sys.argv[5] == 'hires'):
    run_hr = True
else:
    run_hr = False

import os

# machine-specific task-allocation options
# nppc must be a multiple of skip!

machine    = ''
nppc       = 0
skip       = 0
if ('bw.txt' in os.listdir('.')):
    machine = 'bw'
    nppc    = 16
    if run_hr:
        skip = 4
    else:
        skip = 2
if ('marcc.txt' in os.listdir('.')):
    machine = 'marcc'
    nppc    = 24
    if run_hr:
        skip = 4
    else:
        skip = 1
    sys.path.insert(0, '/home-1/bkinch1@jhu.edu/.local/lib/python2.7/site-packages') # our custom-compiled mpi4py which only works with intelmpi/2018
if ('lanl.txt' in os.listdir('.')):
    machine = 'lanl'
    nppc    = 36
    if run_hr:
        skip = 3
    else:
        skip = 1
if ('protagoras.txt' in os.listdir('.')):
    machine = 'protagoras'
    nppc    = 6
    if run_hr:
        skip = 2
    else:
        skip = 1

num_procs_per_node = 36
if run_hr:
    num_procs_to_use = 13
else:
    num_procs_to_use = 26

cwd = os.path.dirname(os.path.realpath(__file__))

sys.path.insert(0, cwd + '/ptxlib/pandurata')

import time

import numpy as np

from scipy.interpolate import interp1d

from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD

from h5py_wrapper import GetFile

rank = COMM_WORLD.Get_rank()
size = COMM_WORLD.Get_size()

num_nodes = int(size/num_procs_per_node)

head_proc_ranks = []
for i in range(0, size, num_procs_per_node):
    head_proc_ranks.append(i)

full_proc_list = []
proc_lists = []
for i in range(0, num_nodes):
    tmp = []
    full_proc_list.append(head_proc_ranks[i])
    for j in range(head_proc_ranks[i] + 1, head_proc_ranks[i] + num_procs_to_use):
        tmp.append(j)
        full_proc_list.append(j)
    proc_lists.append(tmp)

if (rank == 0):
    print('size = ' + repr(size))
    print('num_nodes = ' + repr(num_nodes))
    print('len(full_proc_list) = ' + repr(len(full_proc_list)))
    print('head_proc_ranks = ' + repr(head_proc_ranks))
    print('full_proc_list  = ' + repr(full_proc_list))

# useful constants (cgs)
# --- --- --- --- --- --- --- --- --- ---
G = 6.6726e-8
c = 3.0e10
# --- --- --- --- --- --- --- --- --- ---

theta_thresh = 1.0

def get_r(flnm):
    with GetFile(flnm, open_type = 'r') as f:
        r = float(f.read('r'))
    return r

def collect_data(flnm):
    with GetFile(flnm, open_type = 'r') as f:
        e_coarse = f.read('e_coarse')
        e_fine   = f.read('e_fine')
        eta_grid = f.read('eta')

        slab_type = flnm.split('/')[-1].split('_')[1]

        if (slab_type == 'whole'):
            r_frac_top     = f.read('r_frac_top')
            r_frac_bot     = f.read('r_frac_bot')
            refl_profs_top = f.read('refl_profs_top')
            refl_profs_bot = f.read('refl_profs_bot')
        else:
            r_frac_top     = f.read('r_frac_top')
            refl_profs_top = f.read('refl_profs_top')

    if (slab_type != 'whole'):
        with GetFile(flnm.replace('upper', 'lower'), open_type = 'r') as f:
            r_frac_bot     = f.read('r_frac_bot')
            refl_profs_bot = f.read('refl_profs_bot')
    
    Ne = len(e_fine)-1
    A  = 10.0**((1.0/Ne)*np.log10(e_fine[-1]/e_fine[0]))
    new_eta_grid = np.zeros(2*Ne+1)
    for i in range(0, 2*Ne+1):
        new_eta_grid[i] = A**(i-Ne)
    bnd_eta_grid = np.zeros(2*Ne+2)
    for i in range(1, 2*Ne+1):
        bnd_eta_grid[i] = 0.5*(new_eta_grid[i] + new_eta_grid[i-1])
    bnd_eta_grid[0]  = new_eta_grid[0] - 0.5*(new_eta_grid[1] - new_eta_grid[0])
    bnd_eta_grid[-1] = new_eta_grid[-1] + 0.5*(new_eta_grid[-1] - new_eta_grid[-2])

    new_refl_profs_top = np.zeros((len(e_coarse), 2*Ne+2))
    new_refl_profs_bot = np.zeros((len(e_coarse), 2*Ne+2))

    for i in range(0, len(e_coarse)):
        f = interp1d(eta_grid, refl_profs_top[i], bounds_error = False, fill_value = (0.0, 1.0))
        new_refl_profs_top[i] = f(bnd_eta_grid)
        f = interp1d(eta_grid, refl_profs_bot[i], bounds_error = False, fill_value = (0.0, 1.0))
        new_refl_profs_bot[i] = f(bnd_eta_grid)

    return r_frac_top, r_frac_bot, new_refl_profs_top, new_refl_profs_bot

def rsp2pan(dir_name):
    dir = dir_name

    phi_index_list = []
    r_index_list   = []

    phi_list = []
    r_list   = []

    for flnm in os.listdir(dir):
        if (flnm[:2] == 'mc'):
            phi_index = int(flnm.split('_')[2])
            if phi_index not in phi_index_list:
                phi_index_list.append(phi_index)
                phi_list.append((phi_index/64.0) * 0.5 * np.pi)
            r_index = int(flnm.split('_')[3].split('.')[0])
            if r_index not in r_index_list:
                r_index_list.append(r_index)
                r_list.append(get_r('./' + dir + '/' + flnm))

    phi_index_list.sort()
    r_index_list.sort()

    phi_list.sort()
    r_list.sort()

    Nphi = len(phi_index_list)
    Nr   = len(r_index_list)

    # get energy grids, presumably the same across all files

    for flnm in os.listdir(dir):
        if (flnm[:2] == 'mc'):
            with GetFile('./' + dir + '/' + flnm) as f:
                e_coarse = f.read('e_coarse')
                e_fine   = f.read('e_fine')
                break

    # not all pairs of (phi_index, r_index) will have slabs;
    # record existence in slab_exists: 1 indicates existence

    slab_exists = np.zeros((Nphi, Nr), dtype = 'i')

    for i in range(0, Nphi):
        for j in range(0, Nr):
            # an extant slab will either have a "whole" file or an "upper" *and* "lower"
            flnm = 'mc_whole_' + repr(phi_index_list[i]) + '_' + repr(r_index_list[j]) + '.h5'
            if flnm in os.listdir(dir):
                slab_exists[i][j] = 1
                r_frac_top, r_frac_bot, refl_profs_top, refl_profs_bot = collect_data('./' + dir + '/' + flnm)
            flnm1 = 'mc_lower_' + repr(phi_index_list[i]) + '_' + repr(r_index_list[j]) + '.h5'
            flnm2 = 'mc_upper_' + repr(phi_index_list[i]) + '_' + repr(r_index_list[j]) + '.h5'
            if (flnm1 in os.listdir(dir)) and (flnm2 in os.listdir(dir)):
                slab_exists[i][j] = 1
                r_frac_top, r_frac_bot, refl_profs_top, refl_profs_bot = collect_data('./' + dir + '/' + flnm2)
            # only make disk-covering data structures the one time
            if (slab_exists.sum() == 1):
                r_frac_top_disk     = np.zeros((Nphi, Nr, len(e_fine)))
                r_frac_bot_disk     = np.zeros((Nphi, Nr, len(e_fine)))
                refl_profs_top_disk = np.zeros((Nphi, Nr, len(e_coarse), 2*(len(e_fine)-1)+2))
                refl_profs_bot_disk = np.zeros((Nphi, Nr, len(e_coarse), 2*(len(e_fine)-1)+2))
            if (slab_exists[i][j] > 0):
                r_frac_top_disk[i][j]     = r_frac_top
                r_frac_bot_disk[i][j]     = r_frac_bot
                refl_profs_top_disk[i][j] = refl_profs_top
                refl_profs_bot_disk[i][j] = refl_profs_bot

    print(Nphi, Nr, slab_exists.sum())

    out_flnm = dir.split('/')[-1].replace('run', 'greens') + '.h5'

    with GetFile(out_flnm, open_type = 'w') as f:        
        f.write('Nphi',                Nphi)
        f.write('Nr',                  Nr)
        f.write('slab_exists',         slab_exists)
        f.write('phi_list',            np.array(phi_list))
        f.write('r_list',              np.array(r_list))
        f.write('e_coarse',            e_coarse)
        f.write('e_fine',              e_fine)
        f.write('r_frac_top_disk',     r_frac_top_disk)
        f.write('r_frac_bot_disk',     r_frac_bot_disk)
        f.write('refl_profs_top_disk', refl_profs_top_disk)
        f.write('refl_profs_bot_disk', refl_profs_bot_disk)

grand_iter = int(sys.argv[3])
if (grand_iter > 0):
    if (rank == 0):
        if run_hr:
            rsp2pan('run_hr')
        else:
            rsp2pan('run_' + repr(grand_iter-1))
    else:
        if run_hr:
            time.sleep(600)
        else:
            time.sleep(180)

# load immutable data needed by pandurata:
if (rank in full_proc_list):
    input_filename = sys.argv[1]

    with GetFile('./' + input_filename, open_type = 'r') as f:
        pan_rr          = f.read('r')
        pan_tt          = f.read('th')
        pan_rr          = f.read('r')
        pan_tt          = f.read('th')
        pan_pp          = f.read('phi')
        pan_rho_ijk     = f.read('rho')
        pan_ut_ijk      = f.read('ucon0')
        pan_ur_ijk      = f.read('ucon1')
        pan_uz_ijk      = f.read('ucon2')
        pan_up_ijk      = f.read('ucon3')
        pan_diskbody_ik = f.read('diskbody_ik')
        pan_Tdisk_ik    = f.read('Tdisk_ik')
        pan_emtop_ik    = f.read('emtop_ik') # same as reftop_ik
        pan_embot_ik    = f.read('embot_ik') # same as refbot_ik
        # zero out density (so no scattering) where coolfunc is zero
        pan_rho_ijk = np.where(f.read('coolfunc') == 0.0, 0.0, pan_rho_ijk)

    if run_hr:
        compton_data_filename = '/lustre/scratch3/turquoise/kinch/compy/compton_data_hr.h5'
    else:
        compton_data_filename = '/lustre/scratch3/turquoise/kinch/compy/compton_data_pan.h5'

    with GetFile(compton_data_filename, open_type = 'r') as f:
        pan_T_list = f.read('T_list')
        pan_E_list = f.read('E_list')
        pan_ratio  = f.read('ratio')
        pan_cdf    = f.read('cdf')

    grand_iter = int(sys.argv[3])

    if (grand_iter > 0):
        if run_hr:
            greens_data_filename = './greens_hr.h5'
        else:
            greens_data_filename = './greens_' + repr(grand_iter-1) + '.h5'
        with GetFile(greens_data_filename, open_type = 'r') as f:
            pan_d_Nphi              = int(f.read('Nphi'))
            pan_d_Nr                = int(f.read('Nr'))
            pan_phi_list            = f.read('phi_list')
            pan_r_list              = f.read('r_list')
            pan_slab_exists         = f.read('slab_exists')
            pan_e_coarse            = f.read('e_coarse')
            pan_r_frac_top_disk     = f.read('r_frac_top_disk')
            pan_r_frac_bot_disk     = f.read('r_frac_bot_disk')
            pan_refl_profs_top_disk = f.read('refl_profs_top_disk')
            pan_refl_profs_bot_disk = f.read('refl_profs_bot_disk')
    else:
        # value are not used; objects need to exist for passing
        pan_d_Nphi              = 0
        pan_d_Nr                = 0
        pan_phi_list            = np.zeros(1)
        pan_r_list              = np.zeros(1)
        pan_slab_exists         = np.zeros(1, dtype = int)
        pan_e_coarse            = np.zeros(1)
        pan_r_frac_top_disk     = np.zeros(1)
        pan_r_frac_bot_disk     = np.zeros(1)
        pan_refl_profs_top_disk = np.zeros(1)
        pan_refl_profs_bot_disk = np.zeros(1)

    if (grand_iter > 0):
        spec_flnm = sys.argv[4]
        with GetFile(spec_flnm, open_type = 'r') as f:
            pan_xstar_fluxtot_top = f.read('seed_spec_top')
            pan_xstar_fluxtot_bot = f.read('seed_spec_bot')
    else:
        pan_xstar_fluxtot_top = np.zeros(1)
        pan_xstar_fluxtot_bot = np.zeros(1)

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

def calc_dV_cgs(r, th, phi, M, a):
    m_bh = M * 2.0e33

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
        for j in range(0, Nth):
            for k in range(0, Nphi):
                dV[i][j][k] = np.sqrt(np.abs(g[i][j][1][1] * g[i][j][2][2] * g[i][j][3][3]))*tmp[i][j]*dphi[k]
#               dV[i][j][k] = np.sqrt(-np.linalg.det(g[i][j]))*tmp[i][j]*dphi[k]
                if ((dV[i][j][k] < 0.0) or np.isnan(dV[i][j][k]) or (r[i] < R_hor)):
                    dV[i][j][k] = 0.0

    dV_cgs = dV*(((G*m_bh)/(c*c))**3)

    return dV_cgs

def replace_line(flnm, line_num, repl_line):
#   orig_file = open(flnm, 'r')
    with open(flnm, 'r') as orig_file:
        lines = []
        for line in orig_file:
            lines.append(line)
#   orig_file.close()

#   new_file = open(flnm, 'w')
    with open(flnm, 'w') as new_file:
        for i in range(0, len(lines)):
            if (i+1 != line_num):
                new_file.write(lines[i])
            else:
                new_file.write(repl_line + '\n')
#   new_file.close()

def Tsec2Tcor(cor_sectors, T_sec, Nr, Nth, Nphi):
    T_cor = np.zeros((Nr, Nth, Nphi))

    for i in range(0, len(cor_sectors)):
        for cell in cor_sectors[i]:
            T_cor[cell[0]][cell[1]][cell[2]] = T_sec[i]

    return T_cor

def run_pandurata_full(run_id, spec_flnm, cor_sectors, T_sec, Nr, Nth, Nphi, Nsec, run_type = 'full'):
    T_cor = Tsec2Tcor(cor_sectors, T_sec, Nr, Nth, Nphi)

    if (run_type == 'threshold'):
        thresh = theta_thresh
        """
        for i in range(0, len(cor_sectors)):
            for cell in cor_sectors[i]:
                if (T_cor[cell[0]][cell[1]][cell[2]] * (8.617e-8/511.0) > theta_thresh):
                    pan_rho_ijk[cell[0]][cell[1]][cell[2]] = 0.0
        """
    else:
        thresh = 1.0e20

    if (num_nodes == 1):
        for i in full_proc_list[1:]:
            COMM_WORLD.send([run_type, run_id, spec_flnm, T_cor, cor_sectors, rank], dest = i, tag = 0)
        start = time.time()
        print('about to run pandurata')
        scat_spec, scat_diskt, scat_diskb, scat_cpow = _pandurata.pandurata(rank,
                                                       pan_rr, pan_tt, pan_pp, pan_rho_ijk, pan_ut_ijk, pan_ur_ijk, pan_uz_ijk, pan_up_ijk, pan_diskbody_ik, pan_Tdisk_ik, pan_emtop_ik, pan_embot_ik,
                                                       pan_T_list, pan_E_list, pan_ratio, pan_cdf,
                                                       pan_d_Nphi, pan_d_Nr,
                                                       pan_phi_list, pan_r_list, pan_slab_exists, pan_e_coarse, pan_r_frac_top_disk, pan_r_frac_bot_disk, pan_refl_profs_top_disk, pan_refl_profs_bot_disk,
                                                       T_cor,
                                                       pan_xstar_fluxtot_top, pan_xstar_fluxtot_bot,
                                                       thresh)
        end = time.time()
        print('at rank = ' + repr(rank) + ', subprocess took ' + '{0:.2f}'.format(end - start) + ' s')
        finished = 1
        while (finished < len(full_proc_list)):
            results    = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)
            scat_spec  += results[0]
            scat_diskt += results[1]
            scat_diskb += results[2]
            scat_cpow  += results[3]
            finished   += 1
    else:
        for i in proc_lists[0]:
            COMM_WORLD.send([run_type, run_id, spec_flnm, T_cor, cor_sectors, rank], dest = i, tag = 0)
        for i in head_proc_ranks[1:]:
            COMM_WORLD.send([run_type, run_id, spec_flnm, T_cor, cor_sectors, rank], dest = i, tag = 0)
        start = time.time()
        scat_spec, scat_diskt, scat_diskb, scat_cpow = _pandurata.pandurata(rank,
                                                       pan_rr, pan_tt, pan_pp, pan_rho_ijk, pan_ut_ijk, pan_ur_ijk, pan_uz_ijk, pan_up_ijk, pan_diskbody_ik, pan_Tdisk_ik, pan_emtop_ik, pan_embot_ik,
                                                       pan_T_list, pan_E_list, pan_ratio, pan_cdf,
                                                       pan_d_Nphi, pan_d_Nr,
                                                       pan_phi_list, pan_r_list, pan_slab_exists, pan_e_coarse, pan_r_frac_top_disk, pan_r_frac_bot_disk, pan_refl_profs_top_disk, pan_refl_profs_bot_disk,
                                                       T_cor,
                                                       pan_xstar_fluxtot_top, pan_xstar_fluxtot_bot,
                                                       thresh)
        end = time.time()
        print('at rank = ' + repr(rank) + ', subprocess took ' + '{0:.2f}'.format(end - start) + ' s')
        finished = 1
        while (finished < len(proc_lists[0]) + len(head_proc_ranks[1:]) + 1):
            results     = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)
            scat_spec  += results[0]
            scat_diskt += results[1]
            scat_diskb += results[2]
            scat_cpow  += results[3]
            finished   += 1

    scat_spec  /= len(full_proc_list)
    scat_diskt /= len(full_proc_list)
    scat_diskb /= len(full_proc_list)
    scat_cpow  /= len(full_proc_list)

    return scat_spec, scat_diskt, scat_diskb, scat_cpow

# ----------+---------- root process ----------+----------
if (rank == 0):
    input_filename = sys.argv[1]

    # simulation parameters
    with GetFile('./' + input_filename, open_type = 'r') as f:
        M = float(f.read('M'))
        a = float(f.read('a'))

    # command-line parameters
    run_id     = int(sys.argv[2])
    grand_iter = int(sys.argv[3])
    spec_flnm  = sys.argv[4]

    # other file(s)
    panhead = './ptxlib/pandurata/panhead.h'

    # read in grid data and construct physical volume elements:
    with GetFile('./' + input_filename, open_type = 'r') as f:
        r   = f.read('r')
        th  = f.read('th')
        phi = f.read('phi')

    Nr   = len(r)
    Nth  = len(th)
    Nphi = len(phi)

#   if os.path.isfile('./dV_cgs.npz'):
#       dV_cgs = np.load('dV_cgs.npz')['dV_cgs']
#   else:
#       dV_cgs = calc_dV_cgs(r, th, phi, M, a)
#       np.savez('dV_cgs.npz', dV_cgs = dV_cgs)
    dV_cgs = calc_dV_cgs(r, th, phi, M, a)

    # construct incor: True if (i, j, k) is in corona
    with GetFile('./' + input_filename, open_type = 'r') as f:
        diskbody_ijk = f.read('diskbody_ijk')
    incor = np.full((Nr, Nth, Nphi), False, dtype = bool)
    for i in range(0, Nr):
        for j in range(0, Nth):
            for k in range(0, Nphi):
                if ((diskbody_ijk[i][j][k] == 0) and (r[i] > 1.0 + np.sqrt(1.0 - a**2))):
                    incor[i][j][k] = True

    # construct harmcor: cooling data per cell (erg/s), in corona cells only
    # (factor of 4 because harm3d data is from a single quadrant simulation)
    with GetFile('./' + input_filename, open_type = 'r') as f:
        coolfunc = f.read('coolfunc')
    harm_cor = np.zeros((Nr, Nth, Nphi))
    for i in range(0, Nr):
        for j in range(0, Nth):
            for k in range(0, Nphi):
                if incor[i][j][k]:
                    harm_cor[i][j][k] = 4 * coolfunc[i][j][k] * dV_cgs[i][j][k]
                    if (np.isnan(harm_cor[i][k][k]) or np.isinf(harm_cor[i][j][k])):
                        harm_cor[i][j][k] = 0.0

    with GetFile('./' + input_filename, open_type = 'r') as f:
        T_C  = f.read('T_C')
        urad = f.read('urad')

    # sector divisions
    del_r   = 3
    del_th  = 4
    del_phi = 4

    # construct cor_sectors:
    # cor_map[i][j][k] --> sector number; cor_sectors[sector number] --> list of cells (i, j, k) in sector
    cor_sectors = []
    cor_map     = -np.ones((Nr, Nth, Nphi), dtype = int)
    for i in range(0, int(Nr/del_r)):
        for j in range(0, int(Nth/del_th)):
            for k in range(0, int(Nphi/del_phi)):
                tmp = []
                for ii in range(i*del_r, (i+1)*del_r):
                    for jj in range(j*del_th, (j+1)*del_th):
                        for kk in range(k*del_phi, (k+1)*del_phi):
                            if incor[ii][jj][kk]:
                                tmp.append((ii,jj,kk))
                                cor_map[ii][jj][kk] = len(cor_sectors)
                if (len(tmp) > 0):
                    cor_sectors.append(tmp)
    Nsec = len(cor_sectors)

    # if not first "grand" iterations, read in last T_e map
    if (grand_iter > 0):
        te_filename = './te.' + repr(run_id) + '.' + repr(grand_iter-1) + '.h5'
        with GetFile(te_filename, open_type = 'r') as f:
            T_cor = f.read('T_e')
    # otherwise, initialize T_cor (at, say, ~ 50 keV everywhere in corona)
    else:
        T_cor = 6.0e8 * np.ones((Nr, Nth, Nphi)) * incor.astype(float)
        """
        T_sec = np.zeros(Nsec)
        for i in range(0, Nsec):
            mean_T_C  = 0.0
            den1      = 0.0
            mean_rho  = 0.0
            mean_cool = 0.0
            mean_urad = 0.0
            den2      = 0.0
            for cell in cor_sectors[i]:
                ii, jj, kk  = cell
                mean_T_C   += T_C[ii][jj][kk] * pan_rho_ijk[ii][jj][kk] * dV_cgs[ii][jj][kk]
                den1       += pan_rho_ijk[ii][jj][kk] * dV_cgs[ii][jj][kk]
                mean_rho   += pan_rho_ijk[ii][jj][kk] * dV_cgs[ii][jj][kk]
                mean_cool  += 4 * coolfunc[ii][jj][kk] * dV_cgs[ii][jj][kk]
                mean_urad  += urad[ii][jj][kk] * dV_cgs[ii][jj][kk]
                den2       += dV_cgs[ii][jj][kk]
            mean_T_C  /= den1
            mean_rho  /= den2
            mean_cool /= den2
            mean_urad /= den2
            if (mean_cool == 0.0):
                tmp_T = 11605.0 * 511.0e3 * mean_T_C
            else:
                tmp_A = 4.0 * 6.652e-25 * 2.998e10 * 1.2 * (mean_rho / 1.673e-24) * mean_urad
                tmp_T = 11605.0 * 511.0e3 * (-1.0/8.0 + np.sqrt((1.0/64.0) + (1.0/4.0)*mean_T_C + (1.0/4.0)*(mean_cool/tmp_A)))
                if (tmp_A <= 0.0):
                    print tmp_A
                    print mean_rho
                    print mean_urad
                    print mean_T_C
                    print '------'
            if ((tmp_T > 0.0) and (not np.isnan(tmp_T)) and (not np.isinf(tmp_T))):
                T_sec[i] = tmp_T
            elif ((mean_T_C > 0.0) and (not np.isnan(mean_T_C)) and (not np.isinf(mean_T_C))):
                T_sec[i] = 11605.0 * 511.0e3 * mean_T_C
            else:
                T_sec[i] = 1.0e8
#           print 'T_sec[' + repr(i) + '] = ' + repr(T_sec[i])
        T_cor = Tsec2Tcor(cor_sectors, T_sec, Nr, Nth, Nphi)
        """

    # sectorize initial T_cor
    T_sec = np.zeros(Nsec)
    for i in range(0, Nsec):
        for cell in cor_sectors[i]:
            T_sec[i] += T_cor[cell[0]][cell[1]][cell[2]]
        T_sec[i] /= len(cor_sectors[i])
    mean_T = T_sec.mean()
    for i in range(0, Nsec):
        if (T_sec[i] == 0.0):
            T_sec[i] = mean_T

    # sectorize coronal power for comparison
    harm_cor_sec = np.zeros(Nsec)
    for i in range(0, Nr):
        for j in range(0, Nth):
            for k in range(0, Nphi):
                l = cor_map[i][j][k]
                if (l != -1):
                    if (not np.isnan(harm_cor[i][j][k])):
                        harm_cor_sec[l] += harm_cor[i][j][k]
                    else:
                        harm_cor[i][j][k] = 0.0
    for i in range(0, Nsec):
        if (harm_cor_sec[i] == 0.0):
            T_sec[i] = 1.0e7

    # perform initial pandurata run at this temperature:
    # set all parameters as needed in panhead
    if run_hr:
        # high spectral resolution run
        replace_line(panhead, 20, '#define GFILE \"./greens_hr.h5\"')
        replace_line(panhead, 44, '#define N 6')
        replace_line(panhead, 46, '#define Ne 400')
        replace_line(panhead, 74, '#define CFILE \"./compton_data_hr.h5\"')
        replace_line(panhead, 77, '#define hires 1')
    else:
        # regular, iterating run
        replace_line(panhead, 20, '#define GFILE \"./greens_' + repr(grand_iter-1) + '.h5\"')
        replace_line(panhead, 44, '#define N 2')
        replace_line(panhead, 46, '#define Ne 140')
        replace_line(panhead, 74, '#define CFILE \"./compton_data_pan.h5\"')
        replace_line(panhead, 77, '#define hires 0')
    replace_line(panhead, 19, '#define HFILE \"./' + input_filename + '\"')
    replace_line(panhead, 23, '#define Mstar ' + '{0:.3}'.format(M/3.0))
    replace_line(panhead, 30, '#define aa ' + '{0:.3}'.format(a))
    replace_line(panhead, 41, '#define RUN_ID ' + repr(run_id))
    replace_line(panhead, 43, '#define grand_iter ' + repr(grand_iter))
    replace_line(panhead, 47, '#define Nr ' + repr(Nr-1))
    replace_line(panhead, 48, '#define Nth ' + repr(Nth-1))
    replace_line(panhead, 49, '#define Nph ' + repr(Nphi-1))
    replace_line(panhead, 75, '#define TFILE \"./data/te.' + repr(run_id) + '.h5\"')
    replace_line(panhead, 79, '#define SFILE \"./' + spec_flnm + '\"')

    time.sleep(15)

    # rebuild pandurata
    os.chdir('./ptxlib/pandurata')
    os.system('./build.sh')
    os.chdir('..')
    os.chdir('..')

    print('about to import pandurata')

    time.sleep(15)

    import _pandurata

    # actually run pandurata
    if run_hr:
        scat_spec, scat_diskt, scat_diskb, scat_cpow = run_pandurata_full(run_id, spec_flnm, cor_sectors, T_sec, Nr, Nth, Nphi, Nsec, run_type = 'threshold')
    else:
        scat_spec, scat_diskt, scat_diskb, scat_cpow = run_pandurata_full(run_id, spec_flnm, cor_sectors, T_sec, Nr, Nth, Nphi, Nsec)

    if run_hr:
        scat_spec_filename = 'scat_spec.' + repr(run_id) + '.' + repr(grand_iter) + '.h5'
        with GetFile(scat_spec_filename, open_type = 'w') as f:
            f.write('data', scat_spec)
        # all we needed was the one run, so exit early
        for i in range(1, size):
            COMM_WORLD.send(-1, dest = i, tag = 0)
        sys.exit()

    # ray-tracing power; 2.42e17 accounts for a Hz^-1 to keV^-1 factor
    rt_cor = 2.42e17 * scat_cpow

    # zero out any non-coronal cells
    rt_cor *= incor.astype(float)

    # sectorize coronal power for comparison
    rt_cor_sec   = np.zeros(Nsec)
    harm_cor_sec = np.zeros(Nsec)
    for i in range(0, Nr):
        for j in range(0, Nth):
            for k in range(0, Nphi):
                l = cor_map[i][j][k]
                if (l != -1):
                    if (not np.isnan(rt_cor[i][j][k])):
                        rt_cor_sec[l] += rt_cor[i][j][k]
                    else:
                        rt_cor[i][j][k] = 0.0
                    if (not np.isnan(harm_cor[i][j][k])):
                        harm_cor_sec[l] += harm_cor[i][j][k]
                    else:
                        harm_cor[i][j][k] = 0.0

    T_hist  = []
    r_hist  = []
    rt_hist = []

    drdT = np.zeros(Nsec)

    last_ratio = 0.0
    ratio      = 0.0

    # pandurata temperature balance iterations:
    # --- --- --- --- --- --- --- --- --- ---
    itr = 0
    while True:
        # construct and print all sorts of various ratios, scores, etc., for measuring how well the convergence process is going
        print('ratio 1 = ' + repr(rt_cor_sec.sum()/harm_cor_sec.sum()))
        print('ratio 2 = ' + repr(rt_cor.sum()/harm_cor.sum()))
        harm_in_zero_secs = 0.0
        for i in range(0, Nsec):
            if (rt_cor_sec[i] == 0.0):
                harm_in_zero_secs += harm_cor_sec[i]
        print('ratio 3 = ' + repr(rt_cor_sec.sum()/(harm_cor_sec.sum() - harm_in_zero_secs)))
        if ((itr > 0) and (itr % 2 == 1)):
            last_ratio = ratio
            ratio      = rt_cor_sec.sum()/(harm_cor_sec.sum() - harm_in_zero_secs)
        score_1 = 0.0
        for i in range(0, Nsec):
            if (harm_cor_sec[i] != 0.0 or rt_cor_sec[i] != 0.0):
                score_1 += 2*(abs(harm_cor_sec[i] - rt_cor_sec[i])/(abs(harm_cor_sec[i]) + abs(rt_cor_sec[i])))*abs(harm_cor_sec[i])
        score_1 /= abs(harm_cor_sec).sum()
        score_2 = 0.0
        for i in range(0, Nr):
            for j in range(0, Nth):
                for k in range(0, Nphi):
                    if (harm_cor[i][j][k] != 0.0 or rt_cor[i][j][k] != 0.0):
                        score_2 += 2*(abs(harm_cor[i][j][k] - rt_cor[i][j][k])/(abs(harm_cor[i][j][k]) + abs(rt_cor[i][j][k])))*abs(harm_cor[i][j][k])
        score_2 /= abs(harm_cor).sum()
        print('score 1 = ' + repr(score_1))
        print('score 2 = ' + repr(score_2))

        T_hist.append(T_sec.copy())
        r_hist.append(rt_cor_sec.copy() - harm_cor_sec.copy())
        rt_hist.append(rt_cor_sec.copy())

#       for i in range(0, Nsec):
#           print 'T_sec[' + repr(i) + '] = ' + "{:.6E}".format(T_sec[i]) + ',    r  = ' + "{:.6E}".format(rt_cor_sec[i] - harm_cor_sec[i]) + ', ' + repr(rt_cor_sec[i]/harm_cor_sec[i])
#           print 'T_sec[' + repr(i) + '] = ' + "{:.6E}".format(T_sec[i]) + ',    rt = ' + "{:.6E}".format(rt_cor_sec[i])

        cnv_sectors = 0
        nz_sectors  = 0
        for i in range(0, Nsec):
            if ((rt_cor_sec[i] != 0.0) and (harm_cor_sec[i] != 0.0)):
                nz_sectors += 1
                if ((rt_cor_sec[i]/harm_cor_sec[i] > 0.99) and (rt_cor_sec[i]/harm_cor_sec[i] < 1.01)):
                    cnv_sectors += 1
        print('sectors converged: ' + repr(cnv_sectors) + '/' + repr(nz_sectors) + ', ' + repr(float(cnv_sectors)/float(nz_sectors)))

        # the actual convergence criterion: does the total ray-tracing power match the total harm power---excluding sectors which see no scattering events---to within 1%?
        if ((float(cnv_sectors)/float(nz_sectors) > 0.95) and (abs((rt_cor_sec.sum()/(harm_cor_sec.sum() - harm_in_zero_secs)) - 1.0) < 0.01)):
            break
        if (float(cnv_sectors)/float(nz_sectors) > 0.99):
            break
        if (itr > 100):
            break

        print('itr = ' + repr(itr))

        """
        if (itr == 4):
            for i in range(0, Nsec):
                if ((rt_cor_sec[i] != 0.0) and ((harm_cor_sec[i] == 0.0) or (not ((rt_cor_sec[i]/harm_cor_sec[i] > 0.999) and (rt_cor_sec[i]/harm_cor_sec[i] < 1.001))))):
                    Theta_matrix = np.array([[(8.617e-5/511.0e3)*T_hist[0][i], ((8.617e-5/511.0e3)*T_hist[0][i])**2, -1.0],
                                             [(8.617e-5/511.0e3)*T_hist[2][i], ((8.617e-5/511.0e3)*T_hist[2][i])**2, -1.0],
                                             [(8.617e-5/511.0e3)*T_hist[4][i], ((8.617e-5/511.0e3)*T_hist[4][i])**2, -1.0]])
                    L_vec = np.array([rt_hist[0][i], rt_hist[2][i], rt_hist[4][i]])
                    try:
                        a, b, c   = np.linalg.solve(Theta_matrix, L_vec)
                        new_Theta = (np.sqrt(a**2 + 4*b*(c + harm_cor_sec[i])) - a)/(2*b)
                        if ((not np.isnan(new_Theta)) and (not np.isinf(new_Theta)) and (new_Theta > 0.0)):
                            T_sec[i] = (511.0e3/8.617e-5)*new_Theta
                    except:
                        print 'singular matrix: ' + repr(Theta_matrix)
        # on even iterations, randomly adjust non-converged sectors up or down in temperature by a little bit
        """
        if (itr % 2 == 0):
            for i in range(0, Nsec):
                if ((harm_cor_sec[i] == 0.0) or (not ((rt_cor_sec[i]/harm_cor_sec[i] > 0.999) and (rt_cor_sec[i]/harm_cor_sec[i] < 1.001)))):
                    if (np.random.rand() < 0.5):
                        T_sec[i] += 0.00001*T_sec[i]
                    else:
                        T_sec[i] -= 0.00001*T_sec[i]
        # on odd iterations, Newton-Raphson your way to a new guess at the temperature in each sector (without letting T change too much);
        # drdT should always be positive---raising the temperature should cause the electrons to cool more---so if it's negative,
        # it's probably some Monte Carlo issue and not a real effect
        else:
            for i in range(0, Nsec):
                if ((rt_cor_sec[i] != 0.0) and (T_hist[-1][i] != T_hist[-2][i])):
                    drdT[i] = (r_hist[-1][i] - r_hist[-2][i])/(T_hist[-1][i] - T_hist[-2][i])
#                   print 'drdT[' + repr(i) + '] = ' + "{:.6E}".format(drdT[i])
                    if (r_hist[-1][i] * r_hist[-2][i] < 0.0):
                        T_sec[i] = 0.5*(T_hist[-1][i] + T_hist[-2][i])
                    elif ((itr >= 5) and (T_hist[-2][i] != T_hist[-4][i]) and (r_hist[-2][i] * r_hist[-4][i] < 0.0)):
                        m = (r_hist[-2][i] - r_hist[-4][i])/(T_hist[-2][i] - T_hist[-4][i])
                        b = r_hist[-2][i] - m * T_hist[-2][i]
                        T_sec[i] = -b/m
                    elif (drdT[i] > 0.0):
                        target_T = T_hist[-2][i] - (r_hist[-2][i]/drdT[i])
                        if (np.isnan(target_T) or np.isinf(target_T)):
                            T_sec[i] = T_hist[-2][i]
                        elif (target_T < 0.1*T_hist[-2][i]):
                            T_sec[i] = 0.1*T_hist[-2][i]
                        elif (target_T > 10*T_hist[-2][i]):
                            T_sec[i] = 10*T_hist[-2][i]
                        else:
                            T_sec[i] = target_T
                    else:
                        T_sec[i] = T_hist[-2][i]

        # run pandurata again with the revised temperature map, and iterate
        scat_spec, scat_diskt, scat_diskb, scat_cpow = run_pandurata_full(run_id, spec_flnm, cor_sectors, T_sec, Nr, Nth, Nphi, Nsec)

        rt_cor = 2.42e17 * scat_cpow

        rt_cor *= incor.astype(float)

        rt_cor_sec = np.zeros(Nsec)
        for i in range(0, Nr):
            for j in range(0, Nth):
                for k in range(0, Nphi):
                    l = cor_map[i][j][k]
                    if (l != -1):
                        if (not np.isnan(rt_cor[i][j][k])):
                            rt_cor_sec[l] += rt_cor[i][j][k]
                        else:
                            rt_cor[i][j][k]  = 0.0

        itr += 1

    # one final run with sectors with T_e above some threshold zeroed out
    scat_spec, scat_diskt, scat_diskb, scat_cpow = run_pandurata_full(run_id, spec_flnm, cor_sectors, T_sec, Nr, Nth, Nphi, Nsec, run_type = 'threshold')

    T_cor = Tsec2Tcor(cor_sectors, T_sec, Nr, Nth, Nphi)

    te_filename = 'te.' + repr(run_id) + '.' + repr(grand_iter) + '.h5'
    with GetFile(te_filename, open_type = 'w') as f:
        f.write('T_e', T_cor)

#   te_file = h5py.File('te.' + repr(run_id) + '.' + repr(grand_iter) + '.h5', 'w')
#   te_file.create_dataset('T_e', data = T_cor)
#   te_file.close()

    scat_spec_filename = 'scat_spec.' + repr(run_id) + '.' + repr(grand_iter) + '.h5'
    with GetFile(scat_spec_filename, open_type = 'w') as f:
        f.write('data', scat_spec)

#   scat_spec_file = h5py.File('scat_spec.' + repr(run_id) + '.' + repr(grand_iter) + '.h5', 'w')
#   scat_spec_file.create_dataset('data', data = scat_spec)
#   scat_spec_file.close()

    scat_diskt_filename = 'scat_diskt.h5'
    with GetFile(scat_diskt_filename, open_type = 'w') as f:
        f.write('data', scat_diskt)

#   scat_diskt_file = h5py.File('scat_diskt.h5', 'w')
#   scat_diskt_file.create_dataset('data', data = scat_diskt)
#   scat_diskt_file.close()

    scat_diskb_filename = 'scat_diskb.h5'
    with GetFile(scat_diskb_filename, open_type = 'w') as f:
        f.write('data', scat_diskb)

#   scat_diskb_file = h5py.File('scat_diskb.h5', 'w')
#   scat_diskb_file.create_dataset('data', data = scat_diskb)
#   scat_diskb_file.close()

    scat_cpow_filename = 'scat_cpow.' + repr(run_id) + '.' + repr(grand_iter) + '.h5'
    with GetFile(scat_cpow_filename, open_type = 'w') as f:
        f.write('data', scat_cpow)

#   scat_cpow_file = h5py.File('scat_cpow.' + repr(run_id) + '.' + repr(grand_iter) + '.h5', 'w')
#   scat_cpow_file.create_dataset('data', data = scat_cpow)
#   scat_cpow_file.close()

    # kill *all* processors
    for i in range(1, size):
        COMM_WORLD.send(-1, dest = i, tag = 0)

# ----------+---------- worker process ----------+----------
else:
    loaded = 0
    while True:
        # receive job
        tmp = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)

        if (loaded == 0):
            import _pandurata
            loaded = 1

        # break out if sent done signal
        if (tmp == -1):
            break

        run_type = tmp[0]
        if (run_type == 'full'):
            run_id, spec_flnm, T_cor, cor_sectors, orig_rank = tmp[1:]
            thresh = 1.0e20
        else:
            run_id, spec_flnm, T_cor, cor_sectors, orig_rank = tmp[1:]
            thresh = theta_thresh
            """
            for i in range(0, len(cor_sectors)):
                for cell in cor_sectors[i]:
                    if (T_cor[cell[0]][cell[1]][cell[2]] * (8.617e-8/511.0) > theta_thresh):
                        pan_rho_ijk[cell[0]][cell[1]][cell[2]] = 0.0
            """

        # are we a "head" core?
        if (rank in head_proc_ranks):
            for i in proc_lists[int(rank/num_procs_per_node)]:
                COMM_WORLD.send([run_type, run_id, spec_flnm, T_cor, cor_sectors, rank], dest = i, tag = 0)
            start = time.time()
            scat_spec, scat_diskt, scat_diskb, scat_cpow = _pandurata.pandurata(rank,
                                                           pan_rr, pan_tt, pan_pp, pan_rho_ijk, pan_ut_ijk, pan_ur_ijk, pan_uz_ijk, pan_up_ijk, pan_diskbody_ik, pan_Tdisk_ik, pan_emtop_ik, pan_embot_ik,
                                                           pan_T_list, pan_E_list, pan_ratio, pan_cdf,
                                                           pan_d_Nphi, pan_d_Nr,
                                                           pan_phi_list, pan_r_list, pan_slab_exists, pan_e_coarse, pan_r_frac_top_disk, pan_r_frac_bot_disk, pan_refl_profs_top_disk, pan_refl_profs_bot_disk,
                                                           T_cor,
                                                           pan_xstar_fluxtot_top, pan_xstar_fluxtot_bot,
                                                           thresh)
            end = time.time()
            print('at rank = ' + repr(rank) + ', subprocess took ' + '{0:.2f}'.format(end - start) + ' s')
            finished = 1
            while (finished < len(proc_lists[int(rank/num_procs_per_node)]) + 1):
                results     = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)
                scat_spec  += results[0]
                scat_diskt += results[1]
                scat_diskb += results[2]
                scat_cpow  += results[3]
                finished   += 1
            # reply to root
            COMM_WORLD.send((scat_spec, scat_diskt, scat_diskb, scat_cpow), dest = 0, tag = 0)
        else:
            start = time.time()
            scat_spec, scat_diskt, scat_diskb, scat_cpow = _pandurata.pandurata(rank,
                                                           pan_rr, pan_tt, pan_pp, pan_rho_ijk, pan_ut_ijk, pan_ur_ijk, pan_uz_ijk, pan_up_ijk, pan_diskbody_ik, pan_Tdisk_ik, pan_emtop_ik, pan_embot_ik,
                                                           pan_T_list, pan_E_list, pan_ratio, pan_cdf,
                                                           pan_d_Nphi, pan_d_Nr,
                                                           pan_phi_list, pan_r_list, pan_slab_exists, pan_e_coarse, pan_r_frac_top_disk, pan_r_frac_bot_disk, pan_refl_profs_top_disk, pan_refl_profs_bot_disk,
                                                           T_cor,
                                                           pan_xstar_fluxtot_top, pan_xstar_fluxtot_bot,
                                                           thresh)
            end = time.time()
            print('at rank = ' + repr(rank) + ', subprocess took ' + '{0:.2f}'.format(end - start) + ' s')
            # reply to local head
            COMM_WORLD.send((scat_spec, scat_diskt, scat_diskb, scat_cpow), dest = orig_rank, tag = 0)
