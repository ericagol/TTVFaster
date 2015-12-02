This is an initial Python version of TTVFaster
(Agol & Deck 2015, http://arxiv.org/abs/1509.01623)

This calls the C implementation of the code using Cython.

Building
--------

To build this code, run:

```
python setup.py install
```

and then navigate to a different directory before you try to use it. If you
want to run code from within this same directory, you should instead run

```
python setup.py build_ext --inplace
```

Usage
-----

There's only one function ``run_ttvfaster``. To look at the docs, run:

```
from ttvfaster import run_ttvfaster
print(run_ttvfaster.__doc__)
```

to see something that looks like:

```
Run TTVFaster by calling the C implementation via Cython.

:param n_planets:
    The number of planets.

:param param_vec:
    The parameter vector in the same format as used by the C
    implementation::

        mstar m1 p1 e1*cos(arg peri1) i1 Omega1 e1*sin(arg peri1) TT1
        + repeat for remaining planets

:param tmin:
    The initial time.

:param tmax:
    The final time.

:param j_max:
    The maximum j to evaluate.

:returns:
    A list of numpy arrays (one for each planet) giving the list of
    transit times for each planet.
```

``demo.py`` gives an example for how you might want to run this.
