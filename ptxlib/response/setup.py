from distutils.core import setup, Extension
import numpy.distutils.misc_util

c_ext = Extension("_response",
                  ["_response.c", "response.c"]
)

setup(
    ext_modules  = [c_ext],
    include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
)

