"""
Setup file to Cythonize Column Cruncher
"""
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name="Odin Column Cruncher",
    ext_modules=cythonize("Column_Cruncher.pyx", annotate=True)
)
