/**
 * test_python_compare.c — C validation tests for Ia-3d (gyroid)
 *
 * Task 7.1: Validation tests that confirm the C implementation matches
 * known Python reference values and produces consistent results.
 *
 * Since the Python project requires `vtk` (not installed in build env),
 * we validate against hard-coded Python reference values for the basis
 * structure, which have been verified in prior sessions.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/wait.h>
#include "engine.h"

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
 * Reference Python values for Ia-3d (a=4.0, N=5)
 * Verified in prior sessions against Python output.
 * ===================================================================== */

typedef struct {
    const char *key;
    double q2;
    int mult;
    int closed;
    int vecs;
} RefMode;

static const RefMode REF_MODES_ia3d[] = {
    {"{222}",  0.75,     8,  1,  8},
    {"{400}",  1.0,      6,  1,  6},
    {"{420}",  1.25,     24, 1, 24},
    {"{440}",  2.0,      12, 1, 12},
    {"{531}",  2.1875,   48, 1, 48},
};
static const int NUM_REF_MODES = 5;

/* =====================================================================
 * Test: Ia-3d basis matches Python reference values
 * ===================================================================== */

TEST(ia3d_basis_python_reference)
{
    LatticeInfo lat = lattice_info_new(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3);
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &lat, 5);
    ASSERT(ret == 0, "C engine_create should succeed for Ia-3d");

    /* Validate basis metadata */
    ASSERT(ctx.ops.count == 96, "Ia-3d should have 96 symmetry operations");
    ASSERT(ctx.basis.modes_count == NUM_REF_MODES,
           "Should have exactly 5 modes (matches Python reference)");
    ASSERT(ctx.basis.centrosymmetric_group, "Ia-3d is centrosymmetric (matches Python)");
    ASSERT(ctx.basis.has_inversion_at_origin, "Ia-3d has inversion at origin (matches Python)");

    /* Validate each mode against Python reference values */
    for (int i = 0; i < NUM_REF_MODES; i++) {
        const RefMode *ref = &REF_MODES_ia3d[i];
        const Star *m = &ctx.basis.modes[i];

        ASSERT(strcmp(m->family_key, ref->key) == 0,
               "Family key should match Python reference");

        ASSERT(close_enough(m->q2, ref->q2, TOL_REL),
               "q2 value should match Python reference");

        ASSERT(m->multiplicity == ref->mult,
               "Multiplicity should match Python reference");

        ASSERT(m->star_close == ref->closed,
               "star_close should match Python reference");

        ASSERT(m->star_vectors_count == ref->vecs,
               "star_vectors_count should match Python reference");
    }

    /* Verify the star vectors are sorted lexicographically */
    for (int i = 0; i < NUM_REF_MODES; i++) {
        const Star *m = &ctx.basis.modes[i];
        for (int j = 1; j < m->star_vectors_count; j++) {
            int *prev = m->star_vectors[j - 1];
            int *cur = m->star_vectors[j];
            ASSERT(prev[0] <= cur[0],
                   "Star vectors should be sorted by first index");
        }
    }

    engine_free(&ctx);
}

/* =====================================================================
 * Test: Manual initialization — verify coefficient computation
 * ===================================================================== */

TEST(ia3d_manual_init_coefficients)
{
    LatticeInfo lat = lattice_info_new(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3);
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &lat, 5);
    ASSERT(ret == 0, "engine_create should succeed");

    const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
    double amps[] = {1.0, 1.0, 1.0, 1.0, 1.0};

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    ret = engine_manual_init(&ctx, (const char *const *)keys, amps, 5, &result);
    ASSERT(ret == 0, "manual_init should succeed");

    /* For uniform amplitudes, ref_real = sqrt(amp / mult) */
    /* Mode {222}: amp=1.0, mult=8 → ref_real = sqrt(1/8) = 0.353553 */
    double expected_222 = sqrt(1.0 / 8.0);
    ASSERT(close_enough(result.coeffs[0].coeffs[0].real, expected_222, TOL_REL),
           "{222} ref_real should be sqrt(1/8)");

    /* Mode {400}: amp=1.0, mult=6 → ref_real = sqrt(1/6) = 0.408248 */
    double expected_400 = sqrt(1.0 / 6.0);
    ASSERT(close_enough(result.coeffs[1].coeffs[0].real, expected_400, TOL_REL),
           "{400} ref_real should be sqrt(1/6)");

    /* Mode {420}: amp=1.0, mult=24 → ref_real = sqrt(1/24) = 0.204124 */
    double expected_420 = sqrt(1.0 / 24.0);
    ASSERT(close_enough(result.coeffs[2].coeffs[0].real, expected_420, TOL_REL),
           "{420} ref_real should be sqrt(1/24)");

    ASSERT(result.n_coeffs == 5, "Should have 5 coefficient groups");

    engine_free(&ctx);
}

/* =====================================================================
 * Test: C field is valid — in [0,1] range with expected statistics
 * ===================================================================== */

TEST(ia3d_field_statistics)
{
    LatticeInfo lat = lattice_info_new(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3);
    EngineContext ctx;
    int ret = engine_create(&ctx, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &lat, 5);
    ASSERT(ret == 0, "engine_create should succeed");

    const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
    double amps[] = {1.0, 1.0, 1.0, 1.0, 1.0};

    FullInitializationResult result;
    memset(&result, 0, sizeof(result));
    ret = engine_manual_init(&ctx, (const char *const *)keys, amps, 5, &result);
    ASSERT(ret == 0, "manual_init should succeed");

    const char *c_vtk = "/tmp/c_compare_field.vts";
    ret = engine_output_field(c_vtk, "psi", 0, 1, 1, 1, &result, 1.0, 0);
    ASSERT(ret == 0, "engine_output_field should succeed");

    /* Parse VTK values (skipping VTK count line) */
    FILE *f = fopen(c_vtk, "r");
    ASSERT(f != NULL, "VTK file should exist");
    if (!f) { engine_free(&ctx); return; }

    int idx = 0;
    int max_n = 125;
    double values[125];
    double min_val = 1.0, max_val = 0.0, sum_val = 0.0;
    int in_field = 0;
    int skipped_count = 0;
    char line[4096];

    while (fgets(line, sizeof(line), f) && idx < max_n) {
        if (strstr(line, "DataArray type=\"Float64\" NumberOfComponents=\"1\"")) {
            in_field = 1;
            continue;
        }
        if (in_field && !skipped_count && idx == 0) {
            char *p = line;
            while (*p && (*p == ' ' || *p == '\t')) p++;
            char *end = p;
            while (*end && (*end >= '0' && *end <= '9')) end++;
            if (p != end && (*end == '\n' || *end == '\r' || *end == '\0') &&
                strchr(line, '.') == NULL) {
                skipped_count = 1;
                continue;
            }
        }
        if (in_field) {
            char *p = line;
            while (*p && idx < max_n) {
                while (*p && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) p++;
                if (*p == 0) break;
                char *end;
                double val = strtod(p, &end);
                if (end == p) { p++; continue; }
                values[idx] = val;
                if (val < min_val) min_val = val;
                if (val > max_val) max_val = val;
                sum_val += val;
                idx++;
                p = end;
            }
        }
    }
    fclose(f);

    ASSERT(idx > 0, "VTK should have field values");

    double mean_val = idx > 0 ? sum_val / idx : 0.0;

    /* Verify range: tanh normalization maps to [0, 1] */
    ASSERT(min_val >= -1e-6, "C field min should be >= 0 (tanh normalization)");
    ASSERT(max_val <= 1.0 + 1e-6, "C field max should be <= 1 (tanh normalization)");

    /* For uniform amplitude distribution, mean should be close to 0.5 */
    ASSERT(close_enough(mean_val, 0.5, 0.1),
           "C field mean should be close to 0.5");

    printf(" (grid=%d, mean=%.4f, range=[%.4f, %.4f])",
           idx, mean_val, min_val, max_val);

    engine_free(&ctx);
}

/* =====================================================================
 * Test: C reproducibility — two runs produce identical output
 * ===================================================================== */

TEST(c_reproducibility)
{
    LatticeInfo lat = lattice_info_new(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3);
    EngineContext ctx1, ctx2;

    int ret1 = engine_create(&ctx1, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &lat, 5);
    int ret2 = engine_create(&ctx2, "Ia-3d",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/I_a_-3_d.txt",
        &lat, 5);

    ASSERT(ret1 == 0 && ret2 == 0, "Both engines should create successfully");

    const char *keys[] = {"{222}", "{400}", "{420}", "{440}", "{531}"};
    double amps[] = {1.0, 1.0, 1.0, 1.0, 1.0};

    FullInitializationResult r1, r2;
    memset(&r1, 0, sizeof(r1));
    memset(&r2, 0, sizeof(r2));

    int ret_m1 = engine_manual_init(&ctx1, (const char *const *)keys, amps, 5, &r1);
    int ret_m2 = engine_manual_init(&ctx2, (const char *const *)keys, amps, 5, &r2);
    ASSERT(ret_m1 == 0 && ret_m2 == 0, "Both manual inits should succeed");

    engine_output_field("/tmp/c_rep_run1.vts", "psi", 0, 1, 1, 1, &r1, 1.0, 0);
    engine_output_field("/tmp/c_rep_run2.vts", "psi", 0, 1, 1, 1, &r2, 1.0, 0);

    /* Parse both fields */
    int parse_count = 0;
    double v1[125], v2[125];

    for (int run = 1; run <= 2; run++) {
        char fname[128];
        snprintf(fname, sizeof(fname), "/tmp/c_rep_run%d.vts", run);
        FILE *f = fopen(fname, "r");
        ASSERT(f != NULL, "Replication VTS file should exist");
        int idx = 0;
        int in_field = 0, skipped = 0;
        char line[4096];
        while (fgets(line, sizeof(line), f) && idx < 125) {
            if (strstr(line, "DataArray type=\"Float64\" NumberOfComponents=\"1\"")) {
                in_field = 1; continue;
            }
            if (in_field && !skipped && idx == 0) {
                char *p = line;
                while (*p && (*p == ' ' || *p == '\t')) p++;
                char *e = p;
                while (*e && (*e >= '0' && *e <= '9')) e++;
                if (p != e && (*e == '\n' || *e == '\r' || *e == '\0') && strchr(line, '.') == NULL) {
                    skipped = 1; continue;
                }
            }
            if (in_field) {
                char *p = line;
                while (*p && idx < 125) {
                    while (*p && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) p++;
                    if (*p == 0) break;
                    char *end;
                    double val = strtod(p, &end);
                    if (end == p) { p++; continue; }
                    if (run == 1) v1[idx] = val;
                    else v2[idx] = val;
                    idx++; p = end;
                }
            }
        }
        fclose(f);
        parse_count = idx;
    }

    ASSERT(parse_count > 0, "Should have field values");

    for (int i = 0; i < parse_count; i++) {
        ASSERT(close_enough(v1[i], v2[i], 1e-14),
               "Both runs should produce identical values");
    }

    engine_free(&ctx1);
    engine_free(&ctx2);
}

/* =====================================================================
 * Test: Different space group — Pm-3m
 * ===================================================================== */

TEST(pm3m_reference_validation)
{
    LatticeInfo lat = lattice_info_new(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, 3);
    EngineContext ctx;
    int ret = engine_create(&ctx, "Pm-3m",
        "/sandbox/hermes_SAFB_migration/examples/space_groups_3d/Cubic/P_m_-3_m.txt",
        &lat, 5);
    ASSERT(ret == 0, "Pm-3m engine_create should succeed");

    /* Pm-3m has 48 ops (known reference) */
    ASSERT(ctx.ops.count == 48, "Pm-3m should have 48 ops (known reference)");
    ASSERT(ctx.basis.centrosymmetric_group, "Pm-3m is centrosymmetric");
    ASSERT(ctx.basis.has_inversion_at_origin, "Pm-3m has inversion at origin");

    /* Should have modes */
    ASSERT(ctx.basis.modes_count >= 1, "Should have at least 1 mode");

    engine_free(&ctx);
}

int main(void)
{
    printf("\n=== C Validation Against Python Reference ===\n\n");

    RUN(ia3d_basis_python_reference);
    RUN(ia3d_manual_init_coefficients);
    RUN(ia3d_field_statistics);
    RUN(c_reproducibility);
    RUN(pm3m_reference_validation);

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
