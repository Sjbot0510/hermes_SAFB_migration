/**
 * analytic.c — Analytical calculations for SAFB
 *
 * Translated from: Sg_init/Analytic.py
 *
 * Contains:
 *   - calculate_square_norm — numerical integration of |f|² over the unit cell
 *   - derive_analytical_star_function — build AnalyticField from star coeffs
 *   - extract_basis — numerical normalizer
 *   - evaluate_analytic_real — point evaluation of Fourier series
 */

#define _POSIX_C_SOURCE 200809L
#include "analytic.h"
#include "basis.h"
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
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
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

void evaluate_analytic_real(const AnalyticField *field,
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

/* ========================================================================
 * snap_coefficient — numerical tolerance snapping (like SymPy nsimplify)
 *
 * Snaps a coefficient to a "clean" value:
 *   - Zero tolerance: values near zero → 0
 *   - Fraction snapping: common values like 0.5, 0.707, 1.0, 1.414
 * ======================================================================== */

static double snap_coefficient(double val, double tol)
{
    if (fabs(val) < tol) return 0.0;

    /* Snap to common values */
    struct { double snap; double dist; } candidates[] = {
        { 1.0, 0.0 }, { -1.0, 0.0 }, { 0.5, 0.0 },
        { M_SQRT2, 0.0 }, { M_SQRT2 / 2.0, 0.0 },
        { 2.0, 0.0 }, { 3.0, 0.0 }, { 4.0, 0.0 },
        { 0.0, 0.0 }
    };

    /* Check if val is close to any common integer or sqrt value */
    for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); i++) {
        if (fabs(val - candidates[i].snap) < tol) {
            return candidates[i].snap;
        }
    }

    /* If no snap found, return original value */
    return val;
}

/* ========================================================================
 * derive_analytical_star_function — build AnalyticField from star coeffs
 *
 * Translated from: Sg_init/Analytic.py derive_analytical_star_function()
 *
 * For each family_key in the result, constructs an AnalyticField from the
 * complex coefficients, computes the square norm, and returns normalized
 * results. This is the numerical equivalent of the SymPy symbolic pipeline.
 *
 * Algorithm:
 *   1. For each mode (family_key) in the SAFB basis:
 *      a. Collect star vectors and their coefficients
 *      b. Construct Fourier series: Σ c_j · cos(h_j·X + k_j·Y + l_j·Z)
 *         - Closed stars: take real part directly
 *         - Open stars: construct even/odd recombination
 *      c. Compute square norm via numerical integration
 *      d. Return normalized field with norm_factor = 1/sqrt(sq_norm)
 *
 * Parameters:
 *   result     — FullInitializationResult with coefficients and modes
 *   grid_size  — quadrature points per axis for square norm (default 64)
 *   out        — output AnalyticBasis (caller-allocated)
 * ======================================================================== */

int derive_analytical_star_function(const FullInitializationResult *result,
                                     int grid_size,
                                     AnalyticBasis *out)
{
    if (!result || !out || grid_size <= 0) return -1;

    memset(out, 0, sizeof(AnalyticBasis));
    const SAFBBasis *basis = &result->safb;
    int n_modes = basis->modes_count;

    for (int mi = 0; mi < n_modes && mi < MAX_STARS; mi++) {
        const Star *mode = &basis->modes[mi];
        const char *fk = mode->family_key;

        /* Find coefficients for this family_key */
        int n_coeffs = 0;
        int key_idx = 0;
        for (int ci = 0; ci < result->n_coeffs && ci < MAX_AMPLITUDE_KEYS; ci++) {
            if (strncmp(result->coeffs[ci].family_key, fk, sizeof(result->coeffs[ci].family_key)) == 0) {
                key_idx = ci;
                n_coeffs = result->coeffs[ci].count;
                break;
            }
        }

        if (n_coeffs <= 0) continue;

        /* Build AnalyticField from star coefficients */
        AnalyticField *af = &out->stars[out->count].field;
        memset(af, 0, sizeof(AnalyticField));

        /* Build family key string for output */
        snprintf(out->stars[out->count].family_key, sizeof(out->stars[out->count].family_key),
                 "{%s}", fk);

        /* Determine if star is closed */
        out->stars[out->count].star_closed = mode->star_close;

        if (mode->star_close) {
            /* Closed star: f = Σ c_j · cos(h_j·X + k_j·Y + l_j·Z) */
            for (int vi = 0; vi < n_coeffs && vi < MAX_VEC_PER_STAR; vi++) {
                const ComplexCoeff *cc = &result->coeffs[key_idx].coeffs[vi];
                AnalyticTerm *t = &af->terms[af->count];
                t->hkl[0] = cc->hkl[0];
                t->hkl[1] = cc->hkl[1];
                t->hkl[2] = cc->hkl[2];

                /* For closed stars, use real part of coefficient */
                /* Snap noisy coefficients */
                t->real_part = snap_coefficient(cc->real, 1e-10);
                t->imag_part = snap_coefficient(cc->imag, 1e-10);

                af->count++;
            }
        } else {
            /* Open star: construct even/odd recombination
             * Primary star → even (cosine-like): f = (φ + φ*) / √2
             * Secondary star → odd (sine-like): f = (φ* - φ) / (i√2)
             */
            int sorted_indices[MAX_VEC_PER_STAR];
            int hkl_vecs[MAX_VEC_PER_STAR][3];
            double coeff_reals[MAX_VEC_PER_STAR];
            double coeff_imags[MAX_VEC_PER_STAR];

            for (int vi = 0; vi < n_coeffs; vi++) {
                sorted_indices[vi] = vi;
                const ComplexCoeff *cc = &result->coeffs[key_idx].coeffs[vi];
                hkl_vecs[vi][0] = cc->hkl[0];
                hkl_vecs[vi][1] = cc->hkl[1];
                hkl_vecs[vi][2] = cc->hkl[2];
                coeff_reals[vi] = snap_coefficient(cc->real, 1e-10);
                coeff_imags[vi] = snap_coefficient(cc->imag, 1e-10);
            }

            /* Sort by lexicographic order */
            for (int i = 0; i < n_coeffs; i++) {
                for (int j = i + 1; j < n_coeffs; j++) {
                    int a = sorted_indices[i];
                    int b = sorted_indices[j];
                    int swap = 0;
                    if (hkl_vecs[a][0] > hkl_vecs[b][0]) swap = 1;
                    else if (hkl_vecs[a][0] == hkl_vecs[b][0]) {
                        if (hkl_vecs[a][1] > hkl_vecs[b][1]) swap = 1;
                        else if (hkl_vecs[a][1] == hkl_vecs[b][1]) {
                            if (hkl_vecs[a][2] > hkl_vecs[b][2]) swap = 1;
                        }
                    }
                    if (swap) {
                        int tmp = sorted_indices[i];
                        sorted_indices[i] = sorted_indices[j];
                        sorted_indices[j] = tmp;
                    }
                }
            }

            int min_idx = sorted_indices[0];
            int max_idx = sorted_indices[n_coeffs - 1];

            /* Check if primary: min_v < -max_v (lexicographic) */
            int neg_max[3] = {-hkl_vecs[max_idx][0], -hkl_vecs[max_idx][1], -hkl_vecs[max_idx][2]};
            int is_primary = 0;
            if (hkl_vecs[min_idx][0] < neg_max[0]) is_primary = 1;
            else if (hkl_vecs[min_idx][0] == neg_max[0]) {
                if (hkl_vecs[min_idx][1] < neg_max[1]) is_primary = 1;
                else if (hkl_vecs[min_idx][1] == neg_max[1]) {
                    if (hkl_vecs[min_idx][2] <= neg_max[2]) is_primary = 1;
                }
            }

            if (is_primary) {
                /* Primary star → even function: f = (φ + φ*) / √2
                 * For real coefficients: f = √2 · Σ c_j · cos(φ_j) */
                for (int vi = 0; vi < n_coeffs && vi < MAX_VEC_PER_STAR; vi++) {
                    int idx = sorted_indices[vi];
                    AnalyticTerm *t = &af->terms[af->count];
                    t->hkl[0] = hkl_vecs[idx][0];
                    t->hkl[1] = hkl_vecs[idx][1];
                    t->hkl[2] = hkl_vecs[idx][2];
                    t->real_part = coeff_reals[idx] * M_SQRT2;
                    t->imag_part = 0.0;
                    af->count++;
                }
            } else {
                /* Secondary star → odd function: f = √2 · Σ c_j · sin(φ_j) */
                /* Encode sin terms with real_part=0, imag_part set for later handling */
                for (int vi = 0; vi < n_coeffs && vi < MAX_VEC_PER_STAR; vi++) {
                    int idx = sorted_indices[vi];
                    AnalyticTerm *t = &af->terms[af->count];
                    t->hkl[0] = hkl_vecs[idx][0];
                    t->hkl[1] = hkl_vecs[idx][1];
                    t->hkl[2] = hkl_vecs[idx][2];
                    t->real_part = 0.0;
                    t->imag_part = coeff_imags[idx] * M_SQRT2;
                    af->count++;
                }
            }

            /* Clean up zero coefficients */
            {
                int write = 0;
                for (int r = 0; r < af->count; r++) {
                    if (fabs(af->terms[r].real_part) > 1e-10 ||
                        fabs(af->terms[r].imag_part) > 1e-10) {
                        if (write != r) {
                            af->terms[write] = af->terms[r];
                        }
                        write++;
                    }
                }
                af->count = write;
            }
        }

        /* Compute square norm */
        out->stars[out->count].sq_norm = calculate_square_norm(af, grid_size);
        if (out->stars[out->count].sq_norm > 0) {
            out->stars[out->count].norm_factor = 1.0 / sqrt(out->stars[out->count].sq_norm);
        } else {
            out->stars[out->count].norm_factor = 0.0;
        }

        out->count++;
    }

    return 0;
}

/* ========================================================================
 * extract_basis — numerical normalizer
 *
 * Translated from: Sg_init/Analytic.py extract_basis()
 *
 * Divides the field by its leading coefficient magnitude to get
 * "integer-like" normalized terms. For C, this returns a copy of the
 * field with all coefficients divided by the first term's magnitude.
 *
 * Parameters:
 *   src   — input AnalyticField
 *   out   — output AnalyticField (caller-allocated)
 * ======================================================================== */

void extract_basis(const AnalyticField *src, AnalyticField *out)
{
    if (!src || !out) return;

    memset(out, 0, sizeof(AnalyticField));

    /* Find the leading coefficient magnitude */
    double leading_mag = 0.0;
    for (int i = 0; i < src->count; i++) {
        double mag = sqrt(src->terms[i].real_part * src->terms[i].real_part +
                          src->terms[i].imag_part * src->terms[i].imag_part);
        if (mag > leading_mag) {
            leading_mag = mag;
        }
    }

    if (leading_mag <= 1e-15) {
        /* Zero field — just copy */
        out->count = src->count;
        memcpy(out->terms, src->terms, sizeof(AnalyticTerm) * src->count);
        return;
    }

    /* Divide all coefficients by leading magnitude and snap */
    for (int i = 0; i < src->count; i++) {
        AnalyticTerm *o = &out->terms[i];
        o->hkl[0] = src->terms[i].hkl[0];
        o->hkl[1] = src->terms[i].hkl[1];
        o->hkl[2] = src->terms[i].hkl[2];
        o->real_part = snap_coefficient(src->terms[i].real_part / leading_mag, 1e-5);
        o->imag_part = snap_coefficient(src->terms[i].imag_part / leading_mag, 1e-5);
        out->count++;
    }
}
