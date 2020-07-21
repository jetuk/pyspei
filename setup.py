#!/usr/bin/env python
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension('pyspei._spei', ['pyspei/_spei.pyx'],
              include_dirs=[np.get_include()],),
]

setup(
    name='pyspei',
    version='0.1',
    description='Python (Cython) Wrapper of the SPEI C library.',
    author='James Tomlinson',
    author_email='tomo.bbe@gmail.com',
    packages=['pyspei'],
    ext_modules=cythonize(extensions),
    entry_points={'console_scripts': ['pyspei=pyspei.cli:cli']},
)
