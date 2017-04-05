import os
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy

os.environ["CC"] = "c++"

sourcefiles = ['nm_conn_molID.pyx']#, 'ctest.c']

extensions = [Extension("nm_conn_molID", sourcefiles, language = 'c++', extra_compile_args=["-fopenmp"],
                        extra_link_args=["-fopenmp"], include_dirs=[numpy.get_include()])]



setup(
    name = "nm_conn_molID",
    ext_modules = cythonize(extensions),
    )




