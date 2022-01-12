import sys
import time
import numpy as np
import h5py
from scipy.interpolate import interp1d

class GetFile:
    def __init__(self, flnm, open_type = 'r+'):
        self.flnm      = flnm
        self.open_type = open_type
    def __enter__(self):
        while True:
            try:
                self.f = h5py.File(self.flnm, self.open_type)
            except:
                print('failed to open file ' + self.flnm + ', waiting to try again...')
                time.sleep(5)
                continue
            return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        (self.f).close()
    def write(self, dsetname, dset):
        (self.f)[dsetname] = dset
    def read(self, dsetname):
        return np.array((self.f)[dsetname])
    def delete(self, dsetname):
        if (dsetname in self.f):
            del (self.f)[dsetname]
        else:
            print('cannot delete ' + dsetname + '; does not exist')
    def has_dset(self, dsetname):
        return dsetname in self.f

def logtable_ndx(value_list, value):
    if (value <= value_list[0]):
        ndx   = 0
        value = value_list[0]
        return (ndx, value)
    elif (value >= value_list[-1]):
        ndx   = len(value_list)-2
        value = value_list[-1]
        return (ndx, value)
    else:
        ndx = int(np.log(value/value_list[0])/np.log(value_list[1]/value_list[0]))
        if ((value > value_list[ndx]) and (value <= value_list[ndx+1])):
            return (ndx, value)
        else:
            for ndx in range(0, len(value_list)-1):
                if ((value > value_list[ndx]) and (value <= value_list[ndx+1])):
                    return (ndx, value)
            print('in logtable_ndx: something went horribly wrong')
            print('in logtable_ndx: value = ' + repr(value) + ', value_list[0] = ' + repr(value_list[0]) + ', value_list[-1] = ' + repr(value_list[-1]))

def construct_cdf(T_list, E_list, eta_grid, T, E, f):
    i, T = logtable_ndx(T_list, T)
    j, E = logtable_ndx(E_list, E)

    x  = T
    x1 = T_list[i]
    x2 = T_list[i+1]
    y  = E
    y1 = E_list[j]
    y2 = E_list[j+1]

    q11 = f.read('cdf_' + repr(i)   + '_' + repr(j))
    q12 = f.read('cdf_' + repr(i)   + '_' + repr(j+1))
    q21 = f.read('cdf_' + repr(i+1) + '_' + repr(j))
    q22 = f.read('cdf_' + repr(i+1) + '_' + repr(j+1))

    cdf = (1.0/((x2 - x1)*(y2 - y1)))*(q11*(x2 - x)*(y2 - y) + q21*(x - x1)*(y2 - y) + q12*(x2 - x)*(y - y1) + q22*(x - x1)*(y - y1))

    e_grid = E * eta_grid

    cdf = interp1d(e_grid, cdf, bounds_error = False, fill_value = (0.0, 1.0))

    return cdf

def prep_for_pandurata(input_flnm, output_flnm, Ne, A):
    with GetFile(input_flnm, open_type = 'r') as f:
        T_list   = f.read('T_list')
        E_list   = f.read('E_list')
        eta_grid = f.read('eta_grid')
        sa_ratio = f.read('sa_ratio')

        prepped_cdf = np.zeros((len(T_list), len(E_list), 2*Ne+2))

        for i in range(0, len(T_list)):
            print(i)
            for j in range(0, len(E_list)):
                Ec_tmp     = np.zeros(2*Ne+1)
                Ec_tmp[Ne] = E_list[j]
                for k in range(Ne+1, 2*Ne+1):
                    Ec_tmp[k] = A*Ec_tmp[k-1]
                for k in range(Ne-1, -1):
                    Ec_tmp[k] = Ec_tmp[k+1]/A
                Eb_tmp = np.zeros(2*Ne+2)
                for k in range(1, 2*Ne+1):
                    Eb_tmp[k] = 0.5*(Ec_tmp[k-1] + Ec_tmp[k])
                Eb_tmp[-1] = 1.0e50
                prepped_cdf[i][j] = construct_cdf(T_list, E_list, eta_grid, T_list[i], E_list[j], f)(Eb_tmp)

    with GetFile(output_flnm) as f:
        f.write('T_list', T_list)
        f.write('E_list', E_list)
        f.write('ratio',  sa_ratio)
        f.write('cdf',    prepped_cdf)

compton_data_flnm = sys.argv[1]

Ne = 140
tmp_e_grid = np.logspace(0, 7, endpoint = True, num = Ne+1)
A = tmp_e_grid[1]/tmp_e_grid[0]

prep_for_pandurata(compton_data_flnm, 'compton_data_pan.h5', Ne, A)

Ne = 400
tmp_e_grid = np.logspace(1, 5, endpoint = True, num = Ne+1)
A = tmp_e_grid[1]/tmp_e_grid[0]

prep_for_pandurata(compton_data_flnm, 'compton_data_hr.h5', Ne, A)
