#!/usr/bin/env python

# --------------------------------------------------------------- #
# Phydelity                                                       #
# Author: Alvin X. Han & Edyth Parker                             #
# --------------------------------------------------------------- #

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np


ext_modules=[
    Extension("phyilpx",
              ["phyilpx/*.pyx"],
              include_dirs=[np.get_include()],
              extra_link_args=['-L/usr/lib/x86_64-linux-gnu/']
              )
]

setup(
    name = 'Phydelity',
    version = '2.0',
    author = 'Alvin X. Han',
    author_email = 'hanxc@bii.a-star.edu.sg',
    description = ('Phydelity.'),
    keywords = 'Phylogenetics, Clustering',
    url = 'https://github.com/alvinxhan/Phydelity',
    ext_modules=cythonize("phyilpx/*.pyx"),
    include_dirs=[np.get_include()],
    packages = ['phyilpd'],
    scripts = ['bin/phydelity.py'],
)
