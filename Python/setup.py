import os
import numpy

from setuptools import setup, Extension

from Cython.Build import cythonize

libraries = []
if os.name == "posix":
    libraries.append("m")
include_dirs = [
    "ttvfaster/C",
    numpy.get_include(),
]

ext = cythonize(
    [
        Extension(
            "ttvfaster._ttvfaster",
            sources=["ttvfaster/C/ttvfaster.c", "ttvfaster/_ttvfaster.pyx"],
            libraries=libraries,
            include_dirs=include_dirs,
        )
    ],
    language_level=3,
)

setup(
    name="ttvfaster",
    version="0.1.0",
    author="Eric Agol, Kat Deck, Daniel Foreman-Mackey",
    author_email="foreman.mackey@gmail.com",
    url="https://github.com/ericagol/TTVFaster",
    license="MIT",
    packages=["ttvfaster"],
    ext_modules=ext,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    install_requires=["numpy", "Cython"],
)
