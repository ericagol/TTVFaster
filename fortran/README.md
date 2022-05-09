
5/9/2022

This directory contains a translation of TTVFaster from Python
into FORTRAN by David Kipping.

To run the example code, from the shell command prompt type:

> gfortran -O3 -c ttvfaster.f90 && gfortran -O3 -o example example.f90 ttvfaster.o

> ./example

This reads the input parameters from kepler62ef_planets.txt, creates
the transit times, and writes these to two files: inner_ttv.txt and
outer_ttv.txt.  David has verified
