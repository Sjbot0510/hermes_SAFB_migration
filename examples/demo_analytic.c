/**
 * demo_analytic.c — Analytical field utilities
 *
 * Demonstrates:
 *   1. Building an AnalyticField from scratch
 *   2. Calculating square norm (|f|² mean value)
 *   3. Point evaluation of the Fourier series
 *   4. Deriving analytical fields from a FullInitializationResult
 *
 * These utilities let you analyze field properties without generating
 * the full 3D grid — useful for energy calculations and validation.
 *
 * Build:
 *   source /sandbox/setup_build.sh && gcc -o demo_analytic demo_analytic.c \
 *     ../src/domain.c ../src/symmetry_ops.c ../src/basis.c ../src/initializers.c \
 *     ../src/field.c ../src/engine.c ../src/analytic.c \
 *     -std=c11 -Wall -O2 -I../include -lm -lfftw3 \
 *     -L/sandbox/miniforge3/envs/build/lib
 *
 * Run:
 *   ./demo_analytic
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/engine.h"
#include "../include/analytic.h"
#include "../include/domain.h"

#define IA3D_PATH "../examples/space_groups_3d/Cubic/I_a_-3_d.txt"

/* ================================================================
 * Part 1: Square norm of a hand-crafted AnalyticField
 * ================================================================
 *
 * For a single cosine term cos(n·x) with amplitude A, the
 * mean squared value is: <A²·cos²> = A²/2
 *
 * For N orthogonal cosines: <Σ Aᵢcosᵢ>² = Σ Aᵢ²/2
 */
static void demo_square_norm(void)
{
    printf("=== Part 1: Square Norm Calculation ===\n\n");

    /* {222} star: 8 vectors (±2,±2,±2) with equal amplitude */
    AnalyticField field;
    memset(&field, 0, sizeof(field));

    int vectors[][3] = {
        { 2, 2, 2}, {-2,-2,-2}, { 2,-2,-2}, {-2, 2,-2},
        {-2,-2, 2}, {-2, 2, 2}, { 2, 2,-2}, { 2,-2, 2}
    };

    double amp = sqrt(1.0 / 8.0);  /* normalization so total = 1.0 */
    field.count = 8;
    for (int i = 0; i < 8; i++) {
        field.terms[i].hkl[0] = vectors[i][0];
        field.terms[i].hkl[1] = vectors[i][1];
        field.terms[i].hkl[2] = vectors[i][2];
        field.terms[i].real_part = amp;
        field.terms[i].imag_part = 0.0;
    }

    double sq_norm = calculate_square_norm(&field, 64);  /* 64³ grid */

    printf("{222} star: 8 vectors, amp = sqrt(1/8) = %.4f\n", amp);
    printf("Square norm <|f|²> = %.6f (expected ≈ 1.0)\n", sq_norm);
    printf("Norm factor 1/√<|f|²> = %.6f\n\n", 1.0 / sqrt(sq_norm));

    /* --- A single cosine: <cos²(x)> = 0.5 --- */
    AnalyticField simple;
    memset(&simple, 0, sizeof(simple));
    simple.count = 1;
    simple.terms[0].hkl[0] = 1;
    simple.terms[0].hkl[1] = 0;
    simple.terms[0].hkl[2] = 0;
    simple.terms[0].real_part = 1.0;
    simple.terms[0].imag_part = 0.0;

    double simple_norm = calculate_square_norm(&simple, 64);
    printf("Simple cos(x):  square norm = %.6f (expected 0.500000)\n", simple_norm);
    printf("Norm factor     = %.6f\n\n", 1.0 / sqrt(simple_norm));

    /* --- Two cosines: <cos²(x) + cos²(y)> = 0.5 + 0.5 = 1.0 --- */
    AnalyticField two;
    memset(&two, 0, sizeof(two));
    two.count = 2;
    two.terms[0].hkl[0] = 1; two.terms[0].hkl[1] = 0; two.terms[0].hkl[2] = 0;
    two.terms[0].real_part = 1.0; two.terms[0].imag_part = 0.0;
    two.terms[1].hkl[0] = 0; two.terms[1].hkl[1] = 1; two.terms[1].hkl[2] = 0;
    two.terms[1].real_part = 1.0; two.terms[1].imag_part = 0.0;

    double two_norm = calculate_square_norm(&two, 64);
    printf("cos(x) + cos(y): square norm = %.6f (expected 1.000000)\n\n", two_norm);
}

/* ================================================================
 * Part 2: Point evaluation of Fourier series
 * ================================================================ */
static void demo_point_evaluation(void)
{
    printf("=== Part 2: Point Evaluation ===\n\n");

    /* f(x,y,z) = cos(x) + cos(y) */
    AnalyticField field;
    memset(&field, 0, sizeof(field));
    field.count = 2;
    field.terms[0].hkl[0] = 1; field.terms[0].hkl[1] = 0; field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 1.0; field.terms[0].imag_part = 0.0;
    field.terms[1].hkl[0] = 0; field.terms[1].hkl[1] = 1; field.terms[1].hkl[2] = 0;
    field.terms[1].real_part = 1.0; field.terms[1].imag_part = 0.0;

    /* Evaluate at several points */
    struct { double x, y, z; } pts[] = {
        {0.0, 0.0, 0.0},    /* cos(0)+cos(0) = 2.0 */
        {3.14159265, 0.0, 0.0},  /* cos(π)+cos(0) = 0.0 */
        {1.57079633, 1.57079633, 0.0},  /* cos(π/2)+cos(π/2) = 0.0 */
        {3.14159265, 3.14159265, 0.0},  /* cos(π)+cos(π) = -2.0 */
    };
    const char *labels[] = {
        "(0, 0, 0)         → 2cos(0)      = 2.000000",
        "(π, 0, 0)         → cos(π)+1     = 0.000000",
        "(π/2, π/2, 0)     → 0+0          = 0.000000",
        "(π, π, 0)         → -1+-1        = -2.000000"
    };

    printf("f(X,Y,Z) = cos(X) + cos(Y)\n\n");
    for (int i = 0; i < 4; i++) {
        double val;
        evaluate_analytic_real(&field, pts[i].x, pts[i].y, pts[i].z, &val);
        printf("  %-35s → %.6f\n", labels[i], val);
    }
    printf("\n");
}

/* ================================================================
 * Part 3: Derive analytical fields from initialization result
 * ================================================================ */
static void demo_derive_from_init(void)
{
    printf("=== Part 3: Derive Analytical Fields from Initialization ===\n\n");

    /* Build a proper Ia-3d initialization */
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d", IA3D_PATH,
                            &(LatticeInfo){4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3},
                            5);
    if (ret != 0) {
        fprintf(stderr, "ERROR: engine_create failed\n");
        return;
    }

    const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
    double amps[] = {1.0, 1.0, 1.0, 1.0, 1.0};

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    ret = engine_manual_init(&ctx, (const char *const *)keys, amps, 5, &result);
    if (ret != 0) {
        fprintf(stderr, "ERROR: manual_init failed\n");
        engine_free(&ctx);
        return;
    }

    printf("Built Ia-3d basis with 5 modes, equal amplitudes = 1.0\n\n");

    /* Derive analytical fields */
    AnalyticBasis ab;
    memset(&ab, 0, sizeof(ab));

    ret = derive_analytical_star_function(&result, 64, &ab);
    if (ret != 0) {
        fprintf(stderr, "ERROR: derive_analytical_star_function failed\n");
        engine_free(&ctx);
        return;
    }

    printf("Derived %d analytical star fields:\n\n", ab.count);
    printf("  %-8s  %-10s  %-14s  %-10s\n",
           "Family", "Closed", "Square Norm", "Norm Factor");
    printf("  %-8s  %-10s  %-14s  %-10s\n",
           "------", "------", "-----------", "-----------");

    for (int i = 0; i < ab.count; i++) {
        printf("  %-8s  %-10s  %-14.6f  %-10.6f\n",
               ab.stars[i].family_key,
               ab.stars[i].star_closed ? "yes" : "no",
               ab.stars[i].sq_norm,
               ab.stars[i].norm_factor);
    }
    printf("\n");

    printf("Each star field's norm_factor = 1/sqrt(<|f|²>)\n");
    printf("Multiply coefficients by norm_factor for unit-normalized modes\n\n");

    /* Cleanup */
    engine_free(&ctx);
}

/* ================================================================
 * Part 4: extract_basis — numerical coefficient cleanup
 * ================================================================ */
static void demo_extract_basis(void)
{
    printf("=== Part 4: extract_basis (Coefficient Cleanup) ===\n\n");

    /* Simulate a field with noisy coefficients (like from numerical iFFT) */
    AnalyticField noisy;
    memset(&noisy, 0, sizeof(noisy));
    noisy.count = 2;

    /* "Ideal" coefficient is 1.0, but numerical noise gives 1.000003 */
    noisy.terms[0].hkl[0] = 1; noisy.terms[0].hkl[1] = 0; noisy.terms[0].hkl[2] = 0;
    noisy.terms[0].real_part = 1.000003;
    noisy.terms[0].imag_part = 0.000001;

    noisy.terms[1].hkl[0] = 0; noisy.terms[1].hkl[1] = 1; noisy.terms[1].hkl[2] = 0;
    noisy.terms[1].real_part = 0.999997;
    noisy.terms[1].imag_part = -0.000002;

    printf("Before extract_basis:\n");
    for (int i = 0; i < noisy.count; i++) {
        printf("  term[%d]: hkl=(%d,%d,%d)  c = %.7f + %.7fi\n",
               i, noisy.terms[i].hkl[0], noisy.terms[i].hkl[1], noisy.terms[i].hkl[2],
               noisy.terms[i].real_part, noisy.terms[i].imag_part);
    }

    AnalyticField cleaned;
    memset(&cleaned, 0, sizeof(cleaned));
    extract_basis(&noisy, &cleaned);

    printf("\nAfter extract_basis:\n");
    for (int i = 0; i < cleaned.count; i++) {
        printf("  term[%d]: hkl=(%d,%d,%d)  c = %.7f + %.7fi\n",
               i, cleaned.terms[i].hkl[0], cleaned.terms[i].hkl[1], cleaned.terms[i].hkl[2],
               cleaned.terms[i].real_part, cleaned.terms[i].imag_part);
    }
    printf("\nSmall noise (< 1e-6) is snapped to zero; all coeffs scaled by first term's magnitude\n\n");
}

/* ================================================================
 * Main
 * ================================================================ */
int main(void)
{
    printf("=== SAFB Demo: Analytical Field Utilities ===\n\n");

    demo_square_norm();
    demo_point_evaluation();
    demo_derive_from_init();
    demo_extract_basis();

    printf("Done. ✅\n");
    return 0;
}
