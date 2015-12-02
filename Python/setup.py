import os
import re
import sys
import numpy

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

from Cython.Build import cythonize

libraries = []
if os.name == "posix":
    libraries.append("m")
include_dirs = [
    "../C",
    numpy.get_include(),
]

ext = cythonize([
    Extension("ttvfaster._ttvfaster",
              sources=["../C/ttvfaster.c", "ttvfaster/_ttvfaster.pyx"],
              libraries=libraries, include_dirs=include_dirs)
])

setup(
    name="ttvfaster",
    version="0.0.1",
    author="Eric Agol, Kat Deck, Daniel Foreman-Mackey",
    url="https://github.com/ericagol/TTVFaster",
    license="MIT",
    packages=["ttvfaster", ],
    ext_modules=ext,
    # description="Blazingly fast Gaussian Processes for regression.",
    # long_description=open("README.rst").read(),
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
