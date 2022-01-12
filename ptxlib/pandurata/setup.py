from distutils.core import setup, Extension
import numpy.distutils.misc_util

c_ext = Extension("_pandurata",
                  ["_pandurata.c", "pandurata.c"],
                  extra_objects = ["accel.o", "calc_g.o", "cashkarp.o", "time_keeper.o", "vector_math.o", "tensor_math.o", "calc_scat_angles.o", "calc_tetrad.o", "nt_spectrum.o", "chandra.o", "get_harm3d_data.o", "lookup_data.o"]
)

setup(
    ext_modules  = [c_ext],
    include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs(),
)

