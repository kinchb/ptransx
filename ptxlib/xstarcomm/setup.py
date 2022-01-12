from distutils.core import setup, Extension
import numpy.distutils.misc_util
import json

with open('../../config.json', 'r') as f:
    config = json.load(f)

c_ext = Extension("_xstarcomm",
                  ["_xstarcomm.c", "xstarcomm.c"],
                  library_dirs = [value for key, value in config['xstar lib dirs'].items()],
                  runtime_library_dirs = [value for key, value in config['xstar lib dirs'].items()],
                  libraries = [value for key, value in config['xstar libs'].items()],
                  extra_objects = ["xstarsub.o"],
                  extra_link_args = ['-lgfortran']
)

setup(
    ext_modules  = [c_ext],
    include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
)
