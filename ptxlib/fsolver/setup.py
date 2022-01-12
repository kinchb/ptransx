from distutils.core import setup, Extension
import numpy.distutils.misc_util

c_ext = Extension("_fsolver",
                  ["_fsolver.c", "fsolver.c"]
)

setup(
    ext_modules  = [c_ext],
    include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
)

