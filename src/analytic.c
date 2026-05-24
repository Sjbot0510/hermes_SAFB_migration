/**
 * analytic.c — Analytical calculations for SAFB
 *
 * Translated from: Sg_init/Analytic.py
 *
 * Contains:
 *   - calculate_square_norm — numerical integration of |f|² over the unit cell
 */

#define _POSIX_C_SOURCE 200809L
#include "analytic.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_TAU
#define M_TAU (2.0 * M_PI)
#endif

/* ========================================================================
 * evaluate_analytic_real — compute f(X,Y,Z) for a real-valued Fourier series
 *
 * The AnalyticField represents:
 *   f(X,Y,Z) = Σ_j c_j · cos(h_j·X + k_j·Y + l_j·Z)
 *
 * This matches the Python code's sp.re(phi) for closed stars, which gives:
 *   Σ c_j · cos(φ_j)
 *
 * Parameters:
 *   field — the Fourier series to evaluate
 *   X, Y, Z — evaluation point in [0, 2π]
 *   out_val — returns f(X,Y,Z)
 * ======================================================================== */

static void evaluate_analytic_real(const AnalyticField *field,
                                    double X, double Y, double Z,
                                    double *out_val)
{
    *out_val = 0.0;

    for (int j = 0; j < field->count; j++) {
        const AnalyticTerm *t = &field->terms[j];
        double phi = (double)t->hkl[0] * X +
                     (double)t->hkl[1] * Y +
                     (double)t->hkl[2] * Z;
        *out_val += t->real_part * cos(phi);
    }
}

/* ========================================================================
 * calculate_square_norm — numerical integration of f² over [0,2π]³
 *
 * Computes the mean squared value:
 *   <|f|²> = (1 / (2π)³) ∫₀^{2π} ∫₀^{2π} ∫₀^{2π} f(X,Y,Z)² dX dY dZ
 *
 * where f(X,Y,Z) = Σ_j c_j · cos(h_j·X + k_j·Y + l_j·Z)
 *
 * This matches the Python code's:
 *   sq_expr = sp.expand(basis_expr * sp.conjugate(basis_expr))
 *   integral_val = sp.integrate(sq_expr, (X, 0, 2π), ...)
 *   mean_sq_val = integral_val / (2π)³
 *
 * For real-valued f, conjugate(f) = f, so |f|² = f².
 *
 * Numerical integration uses a grid of quadrature points.
 * For orthogonal cosine terms, this converges to:
 *   <|f|²> = ½ Σ c_j²  (for non-DC terms)
 *   <|f|²> = c₀²       (for DC term only)
 *
 * Parameters:
 *   field     — the Fourier series to integrate
 *   grid_size — number of quadrature points per axis (default 64)
 *
 * Returns:
 *   The mean squared value (e.g., 0.5 for cos(X) over [0,2π]).
 *   The normalization coefficient is 1/sqrt(this_value).
 * ======================================================================== */

double calculate_square_norm(const AnalyticField *field, int grid_size)
{
    if (!field || grid_size <= 0) return 0.0;
    if (field->count == 0) return 0.0;

    double total_sq = 0.0;

    for (int ix = 0; ix < grid_size; ix++) {
        double X = M_TAU * (double)ix / (double)grid_size;
        for (int iy = 0; iy < grid_size; iy++) {
            double Y = M_TAU * (double)iy / (double)grid_size;
            for (int iz = 0; iz < grid_size; iz++) {
                double Z = M_TAU * (double)iz / (double)grid_size;

                double f_val;
                evaluate_analytic_real(field, X, Y, Z, &f_val);

                total_sq += f_val * f_val;
            }
        }
    }

    /* Average over all grid points — equals integral / (2π)³ */
    double N = (double)(grid_size * grid_size * grid_size);
    return total_sq / N;
}
