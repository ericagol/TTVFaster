#ifndef _TTVFASTER_H_
#define _TTVFASTER_H_

typedef struct {
    int planet;
    int epoch;
    double time;
} CalcTransit;

CalcTransit* ttvfaster (
    // Inputs:
    int n_planets,
    double* params,
    double t0,
    double tf,
    int m_max,

    // Outputs:
    int* n_events_out
);

#endif
