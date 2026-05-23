/**
 * domain.c — Data structure implementations for SAFB
 *
 * Translated from: Sg_init/domain.py
 */

#include "domain.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ========================================================================
 * LatticeInfo
 * ======================================================================== */

LatticeInfo lattice_info_new(double a, double b, double c,
                             double alpha, double beta, double gamma, int dim)
{
    LatticeInfo l;
    l.a = a; l.b = b; l.c = c;
    l.alpha = alpha; l.beta = beta; l.gamma = gamma;
    l.dim = dim;
    return l;
}

LatticeInfo lattice_info_new_2d(double a, double b, double gamma)
{
    return lattice_info_new(a, b, 1.0, 90.0, 90.0, gamma, 2);
}

/* ========================================================================
 * ScatteringProfile
 * ======================================================================== */

ScatteringProfile scattering_profile_new(int N)
{
    ScatteringProfile sp;
    sp.hkl = (int32_t*)calloc((size_t)N * 3, sizeof(int32_t));
    sp.q   = (double*)  calloc((size_t)N,    sizeof(double));
    sp.intensity = (double*) calloc((size_t)N, sizeof(double));
    sp.num_peaks = N;
    return sp;
}

void scattering_profile_free(ScatteringProfile* sp)
{
    if (!sp) return;
    free(sp->hkl);
    free(sp->q);
    free(sp->intensity);
    sp->hkl = NULL;
    sp->q = NULL;
    sp->intensity = NULL;
    sp->num_peaks = 0;
}
