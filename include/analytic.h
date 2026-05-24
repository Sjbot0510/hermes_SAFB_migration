/**
 * analytic.h — Analytical calculations for SAFB
 *
 * Translated from: Sg_init/Analytic.py
 *
 * Contains:
 *   - calculate_square_norm — numerical integration of |f|² over the unit cell
 */

#ifndef ANALYTIC_H
#define ANALYTIC_H

#include <stdint.h>
#include <stddef.h>
#include "initializers.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * AnalyticField — a Fourier series defined by wave vectors + coefficients
 * ======================================================================== */

#define MAX_ANALYTIC_TERMS 256

typedef struct {
    int hkl[3];             /* wave vector (h, k, l) */
    double real_part;       /* real coefficient */
    double imag_part;       /* imaginary coefficient */
} AnalyticTerm;

typedef struct {
    AnalyticTerm terms[MAX_ANALYTIC_TERMS];
    int count;               /* number of terms */
} AnalyticField;

/* ========================================================================
 * calculate_square_norm — numerical integration of |f(X,Y,Z)|² over [0,2π]³
 *
 * Computes the mean squared value:
 *   <|f|²> = (1 / (2π)³) ∫₀^{2π} ∫₀^{2π} ∫₀^{2π} |f(X,Y,Z)|² dX dY dZ
 *
 * Numerical integration uses an Nx×Ny×Nz grid of quadrature points.
 * The Fourier series is evaluated as:
 *   f(X,Y,Z) = Σ_j c_j · exp(i · (h_j·X + k_j·Y + l_j·Z))
 * where c_j = real_part + i·imag_part.
 *
 * For open stars (non-closed under inversion), the Python code computes
 * |f|² = f · conj(f) = f · f* to handle the conjugate pair structure.
 *
 * Parameters:
 *   field     — the Fourier series to integrate
 *   grid_size — number of quadrature points per axis (default 64)
 *
 * Returns:
 *   The mean squared value (e.g., 0.375 for cos(X) over [0,2π]).
 *   The normalization coefficient is 1/sqrt(this_value).
 * ======================================================================== */

double calculate_square_norm(const AnalyticField *field, int grid_size);

#ifdef __cplusplus
}
#endif

#endif /* ANALYTIC_H */
