import time
import numpy as np
import h5py

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
