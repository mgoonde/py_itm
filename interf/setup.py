from setuptools import setup, Extension
from Cython.Build import cythonize
from os.path import dirname,abspath,join
from inspect import getsourcefile
import numpy


# path to this file
mypath=dirname(abspath(getsourcefile(lambda:0)))
# one dir up
mypath = dirname(mypath)
libpath = join( mypath, "lib/libitm.so" )

extensions = [
    Extension(
              "py_itm", ["py_itm.pyx"],
              extra_objects=[libpath],
              define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
              )
]


setup(
    name="py_itm",
    ext_modules=cythonize(extensions, include_path = [numpy.get_include()]),
    )

