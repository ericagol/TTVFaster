cimport cython
from libc.stdlib cimport free

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


cdef extern from "ttvfaster.h":

    ctypedef struct CalcTransit:
        int planet
        int epoch
        double time

    cdef CalcTransit* ttvfaster (
        # Inputs:
        int n_planets,
        double* params,
        double t0,
        double tf,
        int m_max,

        # Outputs:
        int* n_events_out
    )


def run_ttvfaster(int n_planets, param_vec, double tmin, double tmax, int j_max):
    """
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

    """
    cdef np.ndarray[DTYPE_t, mode='c'] params = np.atleast_1d(param_vec)

    cdef int npar = 1 + 7 * n_planets
    if (npar != params.shape[0]):
        raise ValueError("incorrect number of parameters")

    # Run TTVFaster.
    cdef int n_events
    cdef CalcTransit* model = ttvfaster(n_planets, <double*>params.data, tmin, tmax, j_max, &n_events);

    # Convert the output into something Python can use.
    times = [[] for _ in range(n_planets)]
    cdef int i
    for i in range(n_events):
        times[model[i].planet].append(model[i].time)
    for i in range(n_planets):
        times[i] = np.array(times[i], dtype=float)

    # Clean up the memory allocated by TTVFaster.
    free(model)
    return times
