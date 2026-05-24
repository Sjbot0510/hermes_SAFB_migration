/**
 * analytic.h — Analytical calculations for SAFB
 *
 * Translated from: Sg_init/Analytic.py
 *
 * Contains:
 *   - calculate_square_norm — numerical integration of |f|² over the unit cell
 *   - derive_analytical_star_function — build AnalyticField from star coeffs + compute norms
 *   - extract_basis — numerical normalizer (divide by leading coefficient)
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
 * AnalyticStarResult — output of derive_analytical_star_function
 * ======================================================================== */

#define MAX_FAMILY_KEYS 64
#define MAX_STARS 256

typedef struct {
    char family_key[64];     /* string key like "{222}" */
    AnalyticField field;     /* the AnalyticField for this star */
    double sq_norm;          /* |field|² mean value */
    double norm_factor;      /* 1.0 / sqrt(sq_norm) for normalization */
    int star_closed;         /* whether the star is closed under inversion */
} AnalyticStarResult;

typedef struct {
    AnalyticStarResult stars[MAX_STARS];
    int count;
} AnalyticBasis;

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

/* ========================================================================
 * evaluate_analytic_real — compute f(X,Y,Z) for a real-valued Fourier series
 *
 * Evaluates f(X,Y,Z) = Σ_j c_j · cos(h_j·X + k_j·Y + l_j·Z)
 *
 * Parameters:
 *   field — the Fourier series to evaluate
 *   X, Y, Z — evaluation point
 *   out_val — returns f(X,Y,Z)
 * ======================================================================== */

void evaluate_analytic_real(const AnalyticField *field,
                            double X, double Y, double Z,
                            double *out_val);

/* ========================================================================
 * derive_analytical_star_function — build AnalyticField from star coeffs
 *
 * Translated from: Sg_init/Analytic.py derive_analytical_star_function()
 *
 * For each family_key in the result, constructs an AnalyticField from the
 * complex coefficients, computes the square norm, and returns normalized
 * results. This is the numerical equivalent of the SymPy symbolic pipeline.
 *
 * Closed stars: f = Σ c_j · cos(h_j·X + k_j·Y + l_j·Z)  (real part only)
 * Open stars:   handles even/odd recombination via conjugate pairs
 *
 * Parameters:
 *   result     — FullInitializationResult with coefficients and modes
 *   grid_size  — quadrature points per axis for square norm (default 64)
 *   out        — output AnalyticBasis (caller-allocated)
 * ======================================================================== */

int derive_analytical_star_function(const FullInitializationResult *result,
                                     int grid_size,
                                     AnalyticBasis *out);

/* ========================================================================
 * extract_basis — numerical normalizer
 *
 * Translated from: Sg_init/Analytic.py extract_basis()
 *
 * Divides the field by its leading coefficient magnitude to get
 * "integer-like" normalized terms. For C, this returns a copy of the
 * field with all coefficients divided by the first term's magnitude.
 *
 * In the Python code this uses SymPy's nsimplify for rational snapping.
 * In C we use numerical tolerance-based snapping to clean noisy coefficients.
 *
 * Parameters:
 *   src   — input AnalyticField
 *   out   — output AnalyticField (caller-allocated)
 * ======================================================================== */

void extract_basis(const AnalyticField *src, AnalyticField *out);

#ifdef __cplusplus
}
#endif

#endif /* ANALYTIC_H */
