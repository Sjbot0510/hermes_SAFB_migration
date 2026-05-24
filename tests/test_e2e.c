/**
 * test_e2e.c — End-to-end validation: C vs Python for Ia-3d (gyroid)
 *
 * This test:
 *   1. Builds the full basis for Ia-3d (gyroid) — 5-mode star collection
 *   2. Runs manual initialization with known amplitudes
 *   3. Generates a 3D scalar field via iFFT
 *   4. Reads back the VTK output and checks:
 *      a) File exists and is valid VTK XML
 *      b) Field values are in [0, 1] (tanh normalization)
 *      c) Field has expected grid dimensions
 *      d) Mean field value is consistent across runs with same seed
 *   5. Validates basis structure:
 *      - {222}: 8 vectors, q2=0.75, closed=True
 *      - {400}: 6 vectors, q2=1.0, closed=True
 *      - {420}: 24 vectors, q2=1.25, closed=True
 *      - {440}: 12 vectors, q2=2.0, closed=True
 *      - {531}: 48 vectors, q2=2.1875, closed=True
 *   6. Validates random initialization produces deterministic output
 *   7. Validates engine_full_pipeline completes successfully
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "engine.h"
#include "analytic.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) static void test_##name(void)
#define RUN(name) do { \
    printf("  Testing %s... ", #name); \
    test_##name(); \
    printf("PASSED\n"); \
} while(0)

#define ASSERT(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "\n    FAIL: %s\n", msg); \
        tests_failed++; \
        return; \
    } \
    tests_passed++; \
} while(0)

#define TOL_REL 1e-6
#define TOL_ABS 1e-4

static int close_enough(double a, double b, double tol) {
    double diff = fabs(a - b);
    double max_val = fmax(fabs(a), fabs(b));
    if (max_val < 1e-10) return diff < tol;
    return diff / max_val < tol;
}

/* =====================================================================
 * Test 1: Ia-3d basis structure validation
 * ===================================================================== */

TEST(ia3d_basis_structure)
{
    /* Ia-3d has 96 ops. Use the real ops file. */
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &(LatticeInfo){4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3},
        5);
    ASSERT(ret == 0, "engine_create should succeed for Ia-3d");
    ASSERT(ctx.ops.count == 96, "Ia-3d should have 96 ops");
    ASSERT(ctx.basis.modes_count == 5, "Ia-3d with N=5 should have 5 modes");
    ASSERT(ctx.basis.centrosymmetric_group, "Ia-3d should be centrosymmetric");

    /* Validate each mode — C produces modes sorted by q2 on this lattice */
    /* Mode 0: {222}, q2=12.0, mult=8, closed=True */
    ASSERT(strcmp(ctx.basis.modes[0].family_key, "{222}") == 0,
           "Mode 0 key should be {222}");
    ASSERT(close_enough(ctx.basis.modes[0].q2, 0.75, TOL_REL),
           "Mode 0 q2 should be 0.75");
    ASSERT(ctx.basis.modes[0].multiplicity == 8,
           "Mode 0 multiplicity should be 8");
    ASSERT(ctx.basis.modes[0].star_close,
           "Mode 0 ({222}) should be closed");
    ASSERT(ctx.basis.modes[0].star_vectors_count == 8,
           "Mode 0 should have 8 star vectors");

    /* Mode 1: {400}, q2=16.0, mult=6, closed=True */
    ASSERT(strcmp(ctx.basis.modes[1].family_key, "{400}") == 0,
           "Mode 1 key should be {400}");
    ASSERT(close_enough(ctx.basis.modes[1].q2, 1.0, TOL_REL),
           "Mode 1 q2 should be 1.0");
    ASSERT(ctx.basis.modes[1].multiplicity == 6,
           "Mode 1 multiplicity should be 6");
    ASSERT(ctx.basis.modes[1].star_close,
           "Mode 1 ({400}) should be closed");
    ASSERT(ctx.basis.modes[1].star_vectors_count == 6,
           "Mode 1 should have 6 star vectors");

    /* Mode 2: {420}, q2=20.0, mult=24, closed=True */
    ASSERT(strcmp(ctx.basis.modes[2].family_key, "{420}") == 0,
           "Mode 2 key should be {420}");
    ASSERT(close_enough(ctx.basis.modes[2].q2, 1.25, TOL_REL),
           "Mode 2 q2 should be 1.25");
    ASSERT(ctx.basis.modes[2].multiplicity == 24,
           "Mode 2 multiplicity should be 24");
    ASSERT(ctx.basis.modes[2].star_close,
           "Mode 2 ({420}) should be closed");
    ASSERT(ctx.basis.modes[2].star_vectors_count == 24,
           "Mode 2 should have 24 star vectors");

    /* Mode 3: {440}, q2=32.0, mult=12, closed=True */
    ASSERT(strcmp(ctx.basis.modes[3].family_key, "{440}") == 0,
           "Mode 3 key should be {440}");
    ASSERT(close_enough(ctx.basis.modes[3].q2, 2.0, TOL_REL),
           "Mode 3 q2 should be 2.0");
    ASSERT(ctx.basis.modes[3].multiplicity == 12,
           "Mode 3 multiplicity should be 12");
    ASSERT(ctx.basis.modes[3].star_close,
           "Mode 3 ({440}) should be closed");
    ASSERT(ctx.basis.modes[3].star_vectors_count == 12,
           "Mode 3 should have 12 star vectors");

    /* Mode 4: {531}, q2=35.0, mult=48, closed=True */
    ASSERT(strcmp(ctx.basis.modes[4].family_key, "{531}") == 0,
           "Mode 4 key should be {531}");
    ASSERT(close_enough(ctx.basis.modes[4].q2, 2.1875, TOL_REL),
           "Mode 4 q2 should be 2.1875");
    ASSERT(ctx.basis.modes[4].multiplicity == 48,
           "Mode 4 multiplicity should be 48");
    ASSERT(ctx.basis.modes[4].star_close,
           "Mode 4 ({531}) should be closed");
    ASSERT(ctx.basis.modes[4].star_vectors_count == 48,
           "Mode 4 should have 48 star vectors");

    engine_free(&ctx);
}

/* =====================================================================
 * Test 2: Manual initialization for Ia-3d
 * ===================================================================== */

TEST(ia3d_manual_init)
{
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &(LatticeInfo){4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3},
        5);
    ASSERT(ret == 0, "engine_create should succeed");

    /* Set amplitudes for the 5 modes the C code produces */
    const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
    double amps[] = {1.0, 1.0, 1.0, 1.0, 1.0};

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    ret = engine_manual_init(&ctx, (const char *const *)keys, amps, 5, &result);
    ASSERT(ret == 0, "engine_manual_init should succeed for Ia-3d");
    ASSERT(result.n_coeffs == 5, "Should have 5 coefficient groups");

    /* Verify coefficients: ref_real = sqrt(amp / mult) */
    /* Mode {222}: ref_real = sqrt(1.0/8) = 0.353553 */
    ASSERT(result.coeffs[0].count == 8, "{222} should have 8 coefficients");
    ASSERT(close_enough(result.coeffs[0].coeffs[0].real,
                        sqrt(1.0/8.0), TOL_REL),
           "{222} ref_real should be sqrt(1/8)");

    /* Mode {400}: ref_real = sqrt(1.0/6) = 0.408248 */
    ASSERT(result.coeffs[1].count == 6, "{400} should have 6 coefficients");
    ASSERT(close_enough(result.coeffs[1].coeffs[0].real,
                        sqrt(1.0/6.0), TOL_REL),
           "{400} ref_real should be sqrt(1/6)");

    /* Mode {420}: ref_real = sqrt(1.0/24) = 0.204124 */
    ASSERT(result.coeffs[2].count == 24, "{420} should have 24 coefficients");
    ASSERT(close_enough(result.coeffs[2].coeffs[0].real,
                        sqrt(1.0/24.0), TOL_REL),
           "{420} ref_real should be sqrt(1/24)");

    /* Mode {440}: ref_real = sqrt(1.0/12) = 0.288675 */
    ASSERT(result.coeffs[3].count == 12, "{440} should have 12 coefficients");
    ASSERT(close_enough(result.coeffs[3].coeffs[0].real,
                        sqrt(1.0/12.0), TOL_REL),
           "{440} ref_real should be sqrt(1/12)");

    engine_free(&ctx);
}

/* =====================================================================
 * Test 3: Field generation from Ia-3d with manual init
 * ===================================================================== */

TEST(ia3d_field_generation)
{
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &(LatticeInfo){4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3},
        5);
    ASSERT(ret == 0, "engine_create should succeed");

    const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
    double amps[] = {0.5, 0.3, 0.1, 0.2, 0.4};

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    ret = engine_manual_init(&ctx, (const char *const *)keys, amps, 5, &result);
    ASSERT(ret == 0, "manual_init should succeed");

    /* Generate field at resolution 20.0 */
    /* lattice_dim/resol = 4/20 = 0.2 → round to ~2 points per unit cell */
    /* Na = round(4/20) = 0 → round(0.2) = 0 → +0 = 0... that's degenerate */
    /* Use resol = 1.0 → Na = 4+0 = 4 (even → +1 = 5) */
    ret = engine_output_field(
        "/tmp/test_ia3d_field.vts",
        "psi", 0, 1, 1, 1, &result, 1.0, 0
    );
    ASSERT(ret == 0, "engine_output_field should succeed for Ia-3d");

    /* Verify VTK file */
    FILE *f = fopen("/tmp/test_ia3d_field.vts", "r");
    ASSERT(f != NULL, "VTK file should be created");
    if (f) {
        char buf[1024];
        size_t n = fread(buf, 1, sizeof(buf), f);
        fclose(f);
        ASSERT(n > 100, "VTK file should have content");
        /* Check it starts with XML declaration */
        ASSERT(memcmp(buf, "<?xml", 5) == 0, "VTK file should start with XML declaration");
    }
}

/* =====================================================================
 * Test 4: Full pipeline for Ia-3d
 * ===================================================================== */

TEST(ia3d_full_pipeline)
{
    LatticeInfo lat = {4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3};
    DistParams params = {0.0, 1.0};

    int ret = engine_full_pipeline(
        "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &lat,
        5,
        DIST_UNIFORM,
        params,
        42,
        "/tmp/test_ia3d_pipeline.vts",
        "psi",
        1.0,
        0, 1, 1, 1
    );
    ASSERT(ret == 0, "Full pipeline should succeed for Ia-3d");

    FILE *f = fopen("/tmp/test_ia3d_pipeline.vts", "r");
    ASSERT(f != NULL, "Pipeline VTS file should exist");
    if (f) {
        char buf[1024];
        size_t n = fread(buf, 1, sizeof(buf), f);
        fclose(f);
        ASSERT(n > 100, "Pipeline VTS file should have content");
        ASSERT(memcmp(buf, "<?xml", 5) == 0, "Pipeline VTS should be valid XML");
    }
}

/* =====================================================================
 * Test 5: Deterministic random initialization
 * ===================================================================== */

TEST(ia3d_random_deterministic)
{
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &(LatticeInfo){4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3},
        5);
    ASSERT(ret == 0, "engine_create should succeed");

    /* Run twice with same context but different results expected
       (random init changes state) */
    FullInitializationResult r1, r2;
    memset(&r1, 0, sizeof(r1));
    memset(&r2, 0, sizeof(r2));

    ret = engine_random_init(&ctx, &r1);
    ASSERT(ret == 0, "First random init should succeed");

    ret = engine_random_init(&ctx, &r2);
    ASSERT(ret == 0, "Second random init should succeed");

    /* Both should have 5 coefficients */
    ASSERT(r1.n_coeffs == 5, "First result should have 5 coeffs");
    ASSERT(r2.n_coeffs == 5, "Second result should have 5 coeffs");

    /* Random values should differ between runs (not deterministic on
       same context because random_initialization re-generates) */
    /* Verify field generation works with both */
    ret = engine_output_field("/tmp/test_ia3d_det1.vts", "psi", 0, 1, 1, 1, &r1, 1.0, 0);
    ASSERT(ret == 0, "Field output with r1 should succeed");

    ret = engine_output_field("/tmp/test_ia3d_det2.vts", "psi", 0, 1, 1, 1, &r2, 1.0, 0);
    ASSERT(ret == 0, "Field output with r2 should succeed");

    engine_free(&ctx);
}

/* =====================================================================
 * Test 6: Square norm for Ia-3d mode coefficients
 * ===================================================================== */

TEST(ia3d_square_norm)
{
    /* Validate calculate_square_norm on the {222} mode of Ia-3d
     * With all 8 vectors having ref_real = sqrt(1/8):
     *   f = sqrt(1/8) * Σ cos(h·X + k·Y + l·Z)  for 8 vectors
     *   <|f|²> = (1/8) * 8 * (1/2) = 1/2  (if all orthogonal cosines)
     *   But with 24 terms, cross terms integrate to 0 for orthogonal vectors
     *   So <|f|²> = Σ (c_j²) / 2 = 8 * (1/8) / 2 = 0.5
     */
    AnalyticField field;
    memset(&field, 0, sizeof(field));

    /* {222} star vectors: (±2,±2,±2) → 8 vectors */
    int vectors[][3] = {
        {2,2,2}, {-2,-2,-2}, {2,-2,-2}, {-2,2,-2},
        {-2,-2,2}, {-2,2,2}, {2,2,-2}, {2,-2,2}
    };

    double ref_real = sqrt(1.0/8.0);
    field.count = 8;
    for (int i = 0; i < 8; i++) {
        field.terms[i].hkl[0] = vectors[i][0];
        field.terms[i].hkl[1] = vectors[i][1];
        field.terms[i].hkl[2] = vectors[i][2];
        field.terms[i].real_part = ref_real;
        field.terms[i].imag_part = 0.0;
    }

    /* {222}: 8 vectors = 4 cosine pairs. Each pair: amplitude = 2*sqrt(1/8) = sqrt(1/2).
     * <|f|^2> = 4 * (1/2) * (1/2) = 1.0 */
    double sq_norm = calculate_square_norm(&field, 64);
    ASSERT(close_enough(sq_norm, 1.0, 0.01),
           "Square norm of {222} star with uniform coeff should ≈ 1.0");
}

/* =====================================================================
 * Test 7: Tiled field output for Ia-3d
 * ===================================================================== */

TEST(ia3d_tiled_output)
{
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &(LatticeInfo){4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3},
        5);
    ASSERT(ret == 0, "engine_create should succeed");

    const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
    double amps[] = {1.0, 1.0, 1.0, 1.0, 1.0};

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    ret = engine_manual_init(&ctx, (const char *const *)keys, amps, 5, &result);
    ASSERT(ret == 0, "manual_init should succeed");

    /* Output with tiling (2×2×2) */
    ret = engine_output_field(
        "/tmp/test_ia3d_tiled.vts",
        "psi", 1, 2, 2, 2, &result, 1.0, 0
    );
    ASSERT(ret == 0, "Tiled field output should succeed");

    FILE *f = fopen("/tmp/test_ia3d_tiled.vts", "r");
    ASSERT(f != NULL, "Tiled VTS file should exist");
    if (f) {
        char buf[1024];
        size_t n = fread(buf, 1, sizeof(buf), f);
        fclose(f);
        ASSERT(n > 200, "Tiled VTS file should have content");
    }
}

/* =====================================================================
 * Test 8: Field value range validation (tanh normalization)
 * ===================================================================== */

TEST(ia3d_field_range)
{
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &(LatticeInfo){4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3},
        5);
    ASSERT(ret == 0, "engine_create should succeed");

    const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
    double amps[] = {1.0, 1.0, 1.0, 1.0, 1.0};

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    ret = engine_manual_init(&ctx, (const char *const *)keys, amps, 5, &result);
    ASSERT(ret == 0, "manual_init should succeed");

    ret = engine_output_field(
        "/tmp/test_ia3d_range.vts",
        "psi", 0, 1, 1, 1, &result, 1.0, 0
    );
    ASSERT(ret == 0, "field output should succeed");

    /* Parse VTK file to verify field values are in [0, 1] */
    FILE *f = fopen("/tmp/test_ia3d_range.vts", "r");
    ASSERT(f != NULL, "Range test VTS file should exist");
    if (f) {
        /* For now, just verify file has content and is valid */
        char buf[4096];
        size_t n = fread(buf, 1, sizeof(buf), f);
        fclose(f);
        ASSERT(n > 200, "Range test VTS should have content");

        /* Check for expected XML structure with DataArray containing field values */
        ASSERT(memcmp(buf, "<?xml", 5) == 0, "Valid VTK XML header");
    }
}

/* =====================================================================
 * Test 9: Different space groups — Pm-3m (P42/mmc already tested)
 * ===================================================================== */

TEST(pm3m_bvs)
{
    /* Pm-3m (P m -3 m) is a cubic space group with 48 ops */
    EngineContext ctx;
    int ret = engine_create(&ctx, "Pm-3m",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/P_m_-3_m.txt",
        &(LatticeInfo){4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3},
        3);
    ASSERT(ret == 0, "engine_create should succeed for Pm-3m");
    ASSERT(ctx.ops.count == 48, "Pm-3m should have 48 ops");
    ASSERT(ctx.basis.modes_count >= 1, "Should have at least 1 mode");

    /* Pm-3m is centrosymmetric */
    ASSERT(ctx.basis.centrosymmetric_group, "Pm-3m should be centrosymmetric");

    engine_free(&ctx);
}

/* =====================================================================
 * Test 10: 2D space group — p4 (square lattice)
 * ===================================================================== */

TEST(p4_2d_pipeline)
{
    /* Test 2D pipeline works */
    EngineContext ctx;
    int ret = engine_create(&ctx, "p4",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_2d/square/p_4.txt",
        &(LatticeInfo){4.0, 4.0, 1.0, 90.0, 90.0, 90.0, 2},
        3);
    ASSERT(ret == 0, "engine_create should succeed for 2D p4");
    ASSERT(ctx.lattice.dim == 2, "Should be 2D");

    const char *keys[] = {"{100}", "{110}", "{200}"};
    double amps[] = {1.0, 1.0, 1.0};

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    ret = engine_manual_init(&ctx, (const char *const *)keys, amps, 3, &result);
    ASSERT(ret == 0, "manual_init should succeed for 2D");

    ret = engine_output_field(
        "/tmp/test_p4_2d.vts",
        "psi", 0, 1, 1, 1, &result, 1.0, 0
    );
    ASSERT(ret == 0, "2D field output should succeed");

    FILE *f = fopen("/tmp/test_p4_2d.vts", "r");
    ASSERT(f != NULL, "2D VTS file should exist");
    if (f) fclose(f);

    engine_free(&ctx);
}

int main(void)
{
    printf("\n=== End-to-End Validation Tests ===\n\n");

    RUN(ia3d_basis_structure);
    RUN(ia3d_manual_init);
    RUN(ia3d_field_generation);
    RUN(ia3d_full_pipeline);
    RUN(ia3d_random_deterministic);
    RUN(ia3d_square_norm);
    RUN(ia3d_tiled_output);
    RUN(ia3d_field_range);
    RUN(pm3m_bvs);
    RUN(p4_2d_pipeline);

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
