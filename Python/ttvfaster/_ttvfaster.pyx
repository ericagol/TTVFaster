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


def run_ttvfaster(int n_planets, np.ndarray[DTYPE_t, mode='c'] params,
                  double tmin, double tmax, int m_max):
    cdef int npar = 1 + 7 * n_planets
    if (npar != params.shape[0]):
        raise ValueError("incorrect number of parameters")

    cdef int n_events
    cdef CalcTransit* model = ttvfaster(n_planets, <double*>params.data, tmin, tmax, m_max, &n_events);

    times = [[] for _ in range(n_planets)]
    cdef int i
    for i in range(n_events):
        times[model[i].planet].append(model[i].time)

    for i in range(n_planets):
        times[i] = np.array(times[i], dtype=float)

    free(model)
    return times
