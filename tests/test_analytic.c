/**
 * test_analytic.c — Tests for analytic.c (calculate_square_norm)
 *
 * Validates against known analytical results for Fourier series.
 *
 * The Python code in Analytic.py:
 *   1. Builds a Fourier series: Σ c_j * exp(i * φ_j)
 *   2. Takes sp.re() → Σ c_j * cos(φ_j)  (real-valued)
 *   3. Calculates square norm: integrate [Σ c_j * cos(φ_j)]² over [0,2π]³
 *
 * Our C code matches this: evaluate_real uses cos(phi) only.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "analytic.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#define TOLERANCE 1e-4
#define GRID_SIZE 64

static int tests_passed = 0;
static int tests_failed = 0;

static void check(const char *name, double actual, double expected, double tol) {
    double diff = fabs(actual - expected);
    if (diff <= tol) {
        printf("  [PASS] %s: got %.8f (expected %.8f)\n", name, actual, expected);
        tests_passed++;
    } else {
        printf("  [FAIL] %s: got %.8f (expected %.8f, diff=%.2e)\n", name, actual, expected, diff);
        tests_failed++;
    }
}

int main(void) {
    AnalyticField field;
    memset(&field, 0, sizeof(field));

    printf("Testing calculate_square_norm:\n\n");

    /* Test 1: cos(X) → expected mean of cos²(X) = 0.5 */
    field.count = 1;
    field.terms[0].hkl[0] = 1;
    field.terms[0].hkl[1] = 0;
    field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 1.0;
    double val1 = calculate_square_norm(&field, GRID_SIZE);
    check("cos(X) → cos² mean", val1, 0.5, TOLERANCE);

    /* Test 2: cos(Y) → same as cos(X) by symmetry */
    field.count = 1;
    field.terms[0].hkl[0] = 0;
    field.terms[0].hkl[1] = 1;
    field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 1.0;
    double val2 = calculate_square_norm(&field, GRID_SIZE);
    check("cos(Y) → cos² mean", val2, 0.5, TOLERANCE);

    /* Test 3: cos(Z) → same */
    field.count = 1;
    field.terms[0].hkl[0] = 0;
    field.terms[0].hkl[1] = 0;
    field.terms[0].hkl[2] = 1;
    field.terms[0].real_part = 1.0;
    double val3 = calculate_square_norm(&field, GRID_SIZE);
    check("cos(Z) → cos² mean", val3, 0.5, TOLERANCE);

    /* Test 4: cos(X) + cos(Y) → expected 0.5 + 0.5 = 1.0 (cross term integrates to 0) */
    field.count = 2;
    field.terms[0].hkl[0] = 1; field.terms[0].hkl[1] = 0; field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 1.0;
    field.terms[1].hkl[0] = 0; field.terms[1].hkl[1] = 1; field.terms[1].hkl[2] = 0;
    field.terms[1].real_part = 1.0;
    double val4 = calculate_square_norm(&field, GRID_SIZE);
    check("cos(X)+cos(Y) → mean", val4, 1.0, TOLERANCE);

    /* Test 5: constant 1 → mean of 1² = 1 */
    field.count = 1;
    field.terms[0].hkl[0] = 0;
    field.terms[0].hkl[1] = 0;
    field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 1.0;
    double val5 = calculate_square_norm(&field, GRID_SIZE);
    check("constant 1 → mean", val5, 1.0, TOLERANCE);

    /* Test 6: √2 · cos(X) → mean of 2·cos² = 1.0 */
    field.count = 1;
    field.terms[0].hkl[0] = 1;
    field.terms[0].hkl[1] = 0;
    field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = M_SQRT2;
    double val6 = calculate_square_norm(&field, GRID_SIZE);
    check("√2·cos(X) → mean", val6, 1.0, TOLERANCE);

    /* Test 7: 3·cos(X) → mean of 9·cos² = 4.5 */
    field.count = 1;
    field.terms[0].hkl[0] = 1;
    field.terms[0].hkl[1] = 0;
    field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 3.0;
    double val7 = calculate_square_norm(&field, GRID_SIZE);
    check("3·cos(X) → mean", val7, 4.5, TOLERANCE);

    /* Test 8: cos(2X) → still mean of cos² = 0.5 (higher freq doesn't change mean) */
    field.count = 1;
    field.terms[0].hkl[0] = 2;
    field.terms[0].hkl[1] = 0;
    field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 1.0;
    double val8 = calculate_square_norm(&field, GRID_SIZE);
    check("cos(2X) → cos² mean", val8, 0.5, TOLERANCE);

    /* Test 9: cos(X)+cos(2Y) → 0.5 + 0.5 = 1.0 */
    field.count = 2;
    field.terms[0].hkl[0] = 1; field.terms[0].hkl[1] = 0; field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 1.0;
    field.terms[1].hkl[0] = 0; field.terms[1].hkl[1] = 2; field.terms[1].hkl[2] = 0;
    field.terms[1].real_part = 1.0;
    double val9 = calculate_square_norm(&field, GRID_SIZE);
    check("cos(X)+cos(2Y) → mean", val9, 1.0, TOLERANCE);

    /* Test 10: empty field → 0 */
    field.count = 0;
    double val10 = calculate_square_norm(&field, GRID_SIZE);
    check("empty field → 0", val10, 0.0, TOLERANCE);

    /* Test 11: NULL field → 0 */
    double val11 = calculate_square_norm(NULL, GRID_SIZE);
    check("NULL field → 0", val11, 0.0, TOLERANCE);

    /* Test 12: cos(X)+cos(Y)+cos(Z) → 3·0.5 = 1.5 */
    field.count = 3;
    field.terms[0].hkl[0] = 1; field.terms[0].hkl[1] = 0; field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 1.0;
    field.terms[1].hkl[0] = 0; field.terms[1].hkl[1] = 1; field.terms[1].hkl[2] = 0;
    field.terms[1].real_part = 1.0;
    field.terms[2].hkl[0] = 0; field.terms[2].hkl[1] = 0; field.terms[2].hkl[2] = 1;
    field.terms[2].real_part = 1.0;
    double val12 = calculate_square_norm(&field, 128);
    check("cos(X)+cos(Y)+cos(Z) (128³ grid)", val12, 1.5, TOLERANCE);

    /* Test 13: 4·cos(X) + 3·cos(Y) → 16*0.5 + 9*0.5 = 8 + 4.5 = 12.5 */
    field.count = 2;
    field.terms[0].hkl[0] = 1; field.terms[0].hkl[1] = 0; field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 4.0;
    field.terms[1].hkl[0] = 0; field.terms[1].hkl[1] = 1; field.terms[1].hkl[2] = 0;
    field.terms[1].real_part = 3.0;
    double val13 = calculate_square_norm(&field, GRID_SIZE);
    check("4·cos(X)+3·cos(Y) → mean", val13, 12.5, TOLERANCE);

    /* Test 14: Grid convergence — cos(X)+cos(Y) with larger grid */
    field.count = 2;
    field.terms[0].hkl[0] = 1; field.terms[0].hkl[1] = 0; field.terms[0].hkl[2] = 0;
    field.terms[0].real_part = 1.0;
    field.terms[1].hkl[0] = 0; field.terms[1].hkl[1] = 1; field.terms[1].hkl[2] = 0;
    field.terms[1].real_part = 1.0;
    double val14a = calculate_square_norm(&field, 32);
    double val14b = calculate_square_norm(&field, 256);
    check("cos(X)+cos(Y) (32³ grid)", val14a, 1.0, 0.01);
    check("cos(X)+cos(Y) (256³ grid)", val14b, 1.0, 1e-6);

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
